"""Injection recovery, labeled evaluation, and catalog-blind batch searches."""
from pathlib import Path
import numpy as np
import pandas as pd
import lightkurve as lk
import matplotlib.pyplot as plt
from .fit import transit_model
from .preprocess import clean_and_flatten
from .search import search
from .pipeline import run_pipeline, save_json


def synthetic_curve(seed=42, period=3., radius_ratio=.08, days=27., noise=.0003):
    rng = np.random.default_rng(seed)
    t = np.arange(0, days, 1/48)
    flux = transit_model(t, [period, .7, radius_ratio, 10., .3, 0.])
    flux += .001*np.sin(2*np.pi*t/13) + rng.normal(0, noise, len(t))
    return lk.LightCurve(time=t, flux=flux, flux_err=np.full(len(t), noise))


def injection_recovery(curves, periods, radii, output='results/injection', repeats=5,
                       stellar_radius=1., stellar_mass=1., seed=42, min_period=.5, max_period=50.):
    if repeats < 1 or stellar_radius <= 0 or stellar_mass <= 0:
        raise ValueError('Positive repeats and stellar properties required')
    rng, rows = np.random.default_rng(seed), []
    for p in periods:
        for radius in radii:
            for trial in range(repeats):
                t0 = min(lc.time.value.min() for lc in curves) + rng.uniform(0, p)
                rp = radius/(109.076*stellar_radius)
                a = 4.2083 * stellar_mass**(1/3)*p**(2/3)/stellar_radius
                injected = []
                for curve in curves:
                    lc = curve.copy()
                    lc.flux *= transit_model(lc.time.value, [p, t0, rp, a, .3, 0.])
                    injected.append(lc)
                _, flat = clean_and_flatten(injected, duration=.15)
                c, _ = search(flat, min_period, max_period)
                # Strict period AND epoch requirement; harmonics are recorded separately.
                phase_error = abs((c['t0']-t0+p/2)%p-p/2)
                exact = abs(c['period']/p-1) < .01
                recovered = exact and phase_error < c['duration'] and c['depth_snr'] >= 7 and c['n_transits'] >= 3
                rows.append(dict(period=p, radius_earth=radius, trial=trial, t0=t0,
                    recovered_period=c['period'], snr=c['depth_snr'], recovered=bool(recovered),
                    alias=bool(min(abs(c['period']/p-.5), abs(c['period']/p-2)) < .01)))
    out = Path(output); out.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame(rows); table.to_csv(out/'trials.csv', index=False)
    grouped = table.groupby(['radius_earth','period']).recovered.agg(['mean','sum','count'])
    # Wilson binomial interval: include uncertainty even for all-success bins.
    n, k = grouped['count'], grouped['sum']; z = 1.96
    center = (k/n+z*z/(2*n))/(1+z*z/n)
    half = z*np.sqrt(k/n*(1-k/n)/n+z*z/(4*n*n))/(1+z*z/n)
    grouped['lower_95'], grouped['upper_95'] = center-half, center+half
    grouped.to_csv(out/'completeness.csv')
    grid = grouped['mean'].unstack()
    fig, ax = plt.subplots(); im = ax.imshow(grid.values, origin='lower', vmin=0, vmax=1, aspect='auto')
    ax.set_xticks(range(len(grid.columns)), labels=grid.columns)
    ax.set_yticks(range(len(grid.index)), labels=grid.index)
    ax.set(xlabel='Period (days)', ylabel='Radius (Earth radii)')
    fig.colorbar(im, ax=ax, label='Recovery fraction'); fig.tight_layout()
    fig.savefig(out/'completeness.png', dpi=150); plt.close(fig)
    save_json(out/'configuration.json', dict(seed=seed, repeats=repeats, stellar_radius=stellar_radius,
        stellar_mass=stellar_mass, baseline=[dict(label=str(c.meta.get('LABEL',c.meta.get('OBJECT','unknown'))),
        cadences=len(c), start=float(c.time.value.min()), stop=float(c.time.value.max())) for c in curves],
        caveat='Conditional on selected stars and fixed impact 0.3; injections before detrending.'))
    return table


def batch(manifest, output='results/batch', fit=False, **kwargs):
    """Freeze results before catalog matching; continue after individual failures."""
    entries = pd.read_csv(manifest)
    results, failures = [], []
    for _, row in entries.iterrows():
        print(f'Analyzing {row.target}', flush=True)
        try:
            results.append(run_pipeline(str(row.target), mission=str(row.get('mission','Kepler')),
                                        output=output, fit=fit, **kwargs))
        except Exception as exc:
            failures.append(dict(target=str(row.target), error=f'{type(exc).__name__}: {exc}'))
            print(f'Failed {row.target}: {exc}', flush=True)
        out = Path(output); out.mkdir(parents=True, exist_ok=True)
        save_json(out/'frozen_results.json', results)
        save_json(out/'failures.json', failures)
    out = Path(output); out.mkdir(parents=True, exist_ok=True)
    save_json(out/'frozen_results.json', results); save_json(out/'failures.json', failures)
    ranked = pd.DataFrame([dict(target=r['target'], period=r['period'], depth=r['depth'], duration=r['duration'],
        snr=r['depth_snr'], passed=r['vetting']['passed'], incomplete=r['vetting']['incomplete']) for r in results])
    if len(ranked):
        ranked = ranked.sort_values(['passed','snr'], ascending=False)
    ranked.to_csv(out/'ranked.csv', index=False)
    return results


def evaluate_labeled(results_file, labels_file, output='results/vetting_evaluation.json'):
    import json
    results = json.loads(Path(results_file).read_text())
    labels = pd.read_csv(labels_file)  # target, disposition; no catalog features in vetting
    predicted = pd.DataFrame([dict(target=r['target'], flagged=not r['vetting']['passed']) for r in results], columns=['target','flagged'])
    merged = labels.merge(predicted, on='target', validate='one_to_one')
    merged = merged[merged.disposition.isin(['CONFIRMED','FALSE POSITIVE'])]
    truth = merged.disposition.eq('FALSE POSITIVE')
    tp = int((truth & merged.flagged).sum()); fn = int((truth & ~merged.flagged).sum())
    fp = int((~truth & merged.flagged).sum()); tn = int((~truth & ~merged.flagged).sum())
    metrics = dict(true_positive=tp, false_negative=fn, false_positive=fp, true_negative=tn,
        false_positive_recall=tp/(tp+fn) if tp+fn else None,
        planet_retention=tn/(tn+fp) if tn+fp else None, matched=len(merged), requested=len(labels))
    save_json(output, metrics)
    return metrics
