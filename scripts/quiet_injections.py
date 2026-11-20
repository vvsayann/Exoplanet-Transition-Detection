"""Inspect two low-variability Kepler baselines and run injection recovery."""
import argparse
from pathlib import Path
from urllib.request import urlopen

import lightkurve as lk
import numpy as np

from src.download import quarantine_incomplete
from src.experiments import injection_recovery
from src.pipeline import save_json, search_transits
import matplotlib.pyplot as plt


# Central asteroseismic estimates and their source.
STARS = {
    '6106415': (1.078, 1.213, 'https://academic.oup.com/mnras/article/479/4/4416/5043232'),
    '6603624': (1.09, 1.181, 'https://doi.org/10.1088/0067-0049/214/2/27'),
}


def run(kic, repeats):
    mass, radius, source = STARS[kic]
    kid = kic.zfill(9)
    filename = f'kplr{kid}-2009166043257_llc.fits'
    url = f'https://stpubdata.s3.amazonaws.com/kepler/public/lightcurves/{kid[:4]}/{kid}/{filename}'
    path = Path('data/kepler_validation') / filename
    path.parent.mkdir(parents=True, exist_ok=True)
    quarantine_incomplete(path)
    if not path.exists():
        with urlopen(url, timeout=30) as response:
            contents = response.read()
        path.write_bytes(contents)
        if quarantine_incomplete(path):
            raise ValueError(f'Incomplete Kepler FITS download for KIC {kic}')
    curve = lk.read(path, flux_column='pdcsap_flux', quality_bitmask='default')

    output = Path('results/research/quiet_injection') / f'KIC_{kic}'
    output.mkdir(parents=True, exist_ok=True)
    raw, flat, peak, _ = search_transits([curve], .15, .5, 8.)
    scatter = []
    fig, axes = plt.subplots(2, 1, figsize=(9, 5), sharex=True)
    for lc, ax, label in zip((raw, flat), axes, ('Normalized PDCSAP', 'Detrended')):
        ppm = (lc.flux.value - np.median(lc.flux.value)) * 1e6
        scatter.append(float(1.4826 * np.median(abs(ppm))))
        ax.plot(lc.time.value, ppm, '.', ms=2)
        ax.set_ylabel(f'{label}\n(ppm)')
    axes[-1].set_xlabel('Time (BKJD days)')
    fig.suptitle(f'KIC {kic} · Kepler quarter 1')
    fig.tight_layout()
    fig.savefig(output / 'baseline.png', dpi=150)
    plt.close(fig)

    flux_error = float(np.median(flat.flux_err.value) * 1e6)
    save_json(output / 'baseline.json', dict(
        target=f'KIC {kic}', quarter=int(curve.quarter), archive_product=url, stellar_source=source,
        stellar_mass_solar=mass, stellar_radius_solar=radius, input_cadences=len(curve),
        usable_cadences=len(flat), baseline_days=float(np.ptp(flat.time.value)),
        normalized_mad_ppm=scatter[0], detrended_mad_ppm=scatter[1], median_flux_error_ppm=flux_error,
        uninjected_best_period=peak['period'], uninjected_best_depth_ppm=peak['depth'] * 1e6,
        uninjected_best_formal_snr=peak['depth_snr'],
        uninjected_snr_scaled_by_mad=peak['depth_snr'] * flux_error / scatter[1],
    ))
    print(f'KIC {kic}: {repeats} trials per bin on inspected quarter 1 photometry', flush=True)
    injection_recovery([curve], [1., 3., 7.], [2., 4., 8.], output=output, repeats=repeats,
                       stellar_mass=mass, stellar_radius=radius, max_period=8., seed=int(kic))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('kic', choices=STARS)
    parser.add_argument('--repeats', type=int, default=20)
    args = parser.parse_args()
    run(args.kic, args.repeats)
