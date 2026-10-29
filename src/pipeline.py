"""One-target orchestration and reproducible output artifacts."""
import json
from pathlib import Path
import platform
from datetime import datetime, timezone
from importlib.metadata import version
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .download import download_light_curve
from .preprocess import clean_and_flatten, transit_mask
from .search import search, search_multiple
from .vetting import vet


def save_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False), encoding='utf-8')


def run_pipeline(target_id, mission='Kepler', output='results', curves=None, duration=.15,
                 min_period=.5, max_period=50., max_segments=4, fit=True, steps=1500,
                 burn=500, stellar_radius=None, stellar_radius_error=0., stellar_density=None,
                 seed=42, multiple=True):
    folder = Path(output) / ''.join(c if c.isalnum() or c in '-_' else '_' for c in str(target_id))
    folder.mkdir(parents=True, exist_ok=True)
    def status(message):
        print(f'[{target_id}] {message}', flush=True)
    if curves is None:
        status('Finding and loading archive light curves (cached files are reused)...')
        curves = download_light_curve(target_id, mission, max_segments=max_segments)
    status(f'Loaded {len(curves)} segments. Cleaning and detrending...')
    raw, flat = clean_and_flatten(curves, duration)
    status(f'First BLS search on {len(flat):,} cleaned cadences...')
    candidate, _ = search(flat, min_period, max_period)
    status('Protecting detected transits and repeating detrending and BLS...')
    # Refit the trend while protecting the first-pass transit ephemeris.
    raw, flat = clean_and_flatten(curves, duration,
        (candidate['period'], candidate['t0'], candidate['duration']*1.5))
    candidate, power = search(flat, min_period, max_period)
    status(f'Primary period: {candidate["period"]:.8f} days. Vetting and searching for additional signals...')
    result = dict(target=str(target_id), mission=mission, **candidate)
    result['created_utc'] = datetime.now(timezone.utc).isoformat()
    result['data_segments'] = [dict(label=str(c.meta.get('LABEL', c.meta.get('OBJECT', 'unknown'))),
        filename=str(c.meta.get('FILENAME', 'in-memory')), cadences=len(c),
        start=float(c.time.value.min()), stop=float(c.time.value.max()), time_format=c.time.format)
        for c in curves]
    result['vetting'] = vet(flat, candidate, stellar_density)
    result['signals'] = search_multiple(flat, min_period=min_period, max_period=max_period) if multiple else [candidate]
    result['configuration'] = dict(duration=duration, min_period=min_period, max_period=max_period,
                                   max_segments=max_segments, seed=seed)
    result['versions'] = {name: version(name) for name in ['numpy','astropy','lightkurve','scipy']}
    result['versions']['python'] = platform.python_version()
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    axes[0].plot(raw.time.value, raw.flux.value, '.', ms=1)
    axes[1].plot(flat.time.value, flat.flux.value, '.', ms=1)
    marks = transit_mask(flat.time.value, candidate['period'], candidate['t0'], candidate['duration'])
    axes[1].plot(flat.time.value[marks], flat.flux.value[marks], '.', ms=2, color='tomato')
    axes[0].set_ylabel('Normalized flux'); axes[1].set_ylabel('Flattened flux')
    axes[1].set_xlabel(f'Time ({flat.time.format}, days)')
    fig.tight_layout(); fig.savefig(folder/'before_after_flatten.png', dpi=150); plt.close(fig)
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(power.period, power.power); ax.set(xlabel='Period (days)', ylabel='BLS power')
    fig.tight_layout(); fig.savefig(folder/'bls_periodogram.png', dpi=150); plt.close(fig)
    phase = (flat.time.value-candidate['t0']+candidate['period']/2)%candidate['period']-candidate['period']/2
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(phase, flat.flux.value, '.', ms=2, alpha=.4)
    if fit:
        from .fit import fit_transit, transit_model, PARAMETERS
        import corner
        status('Starting physical parameter fitting...')
        summary, samples = fit_transit(flat, candidate, steps=steps, burn=burn, seed=seed,
            stellar_radius=stellar_radius, stellar_radius_error=stellar_radius_error, progress=True)
        result['fit'] = summary
        np.save(folder/'posterior.npy', samples)
        status('Saving fitted model and corner plot...')
        theta = np.array([summary[k]['median'] for k in PARAMETERS])
        model_phase = np.linspace(-candidate['duration']*2, candidate['duration']*2, 1000)
        ax.plot(model_phase, transit_model(theta[1]+model_phase, theta,
                exposure=np.median(np.diff(flat.time.value))), color='tomato')
        cf = corner.corner(samples, labels=PARAMETERS, show_titles=True)
        cf.savefig(folder/'corner.png', dpi=130); plt.close(cf)
    ax.set(xlim=(-candidate['duration']*2, candidate['duration']*2), xlabel='Time from transit (days)', ylabel='Normalized flux')
    fig.tight_layout(); fig.savefig(folder/'folded_transit.png', dpi=150); plt.close(fig)
    np.savez_compressed(folder/'lightcurve.npz', time=flat.time.value, flux=flat.flux.value,
                        flux_err=flat.flux_err.value, time_format=flat.time.format)
    save_json(folder/'result.json', result)
    status(f'Finished. Results saved to {folder.resolve()}')
    return result
