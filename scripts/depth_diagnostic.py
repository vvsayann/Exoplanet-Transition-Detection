"""Measure depth-definition and detrending effects on the known-target fits."""
import argparse
import json
from pathlib import Path

import lightkurve as lk
import numpy as np
import pandas as pd

from src.download import download_light_curve
from src.fit import PARAMETERS, transit_model
from src.preprocess import clean_and_flatten, transit_mask
from src.search import search


def wasp18_occultation(raw, theta):
    """Estimate the eclipse from per-orbit local baselines in normalized PDCSAP."""
    t, flux, error = (np.asarray(v, dtype=float) for v in (raw.time.value, raw.flux.value, raw.flux_err.value))
    period, t0 = theta[:2]
    time_since_occultation = t - (t0 + period / 2)
    epochs = np.rint(time_since_occultation / period).astype(int)
    phase_days = time_since_occultation - epochs * period
    depths = []
    for epoch in np.unique(epochs):
        event = epochs == epoch
        core = event & (abs(phase_days) < .033)
        before = event & (phase_days >= -.14) & (phase_days <= -.065)
        after = event & (phase_days >= .065) & (phase_days <= .14)
        if min(core.sum(), before.sum(), after.sum()) < 15:
            continue
        baseline = before | after
        coef = np.polynomial.polynomial.polyfit(phase_days[baseline], flux[baseline], 2, w=1 / error[baseline])
        prediction = np.polynomial.polynomial.polyval(phase_days[core], coef)
        depths.append(np.average(prediction - flux[core], weights=error[core]**-2) * 1e6)
    return dict(occultation_events=len(depths),
                occultation_depth_ppm=float(np.mean(depths)) if depths else None,
                occultation_se_ppm=float(np.std(depths, ddof=1) / np.sqrt(len(depths))) if len(depths) > 1 else None)


def measure(result):
    cfg = result['configuration']
    files = [Path(segment.get('filename', '')) for segment in result['data_segments']]
    if all(path.is_file() for path in files):
        curves = [lk.read(path, flux_column='pdcsap_flux', quality_bitmask='default') for path in files]
    else:
        curves = download_light_curve(result['target'], result['mission'], max_segments=cfg['max_segments'])

    _, first_pass = clean_and_flatten(curves, cfg['duration'])
    first, _ = search(first_pass, cfg['min_period'], cfg['max_period'])
    mask_ephemeris = (first['period'], first['t0'], first['duration'] * 1.5)
    theta = [result['fit'][key]['median'] for key in PARAMETERS]
    theta[-1] = 0.
    raw, flat = clean_and_flatten(curves, cfg['duration'], mask_ephemeris)
    exposure = np.median(np.diff(flat.time.value))
    # Divide out the fitted transit, then compare both branches with identical
    # masking and weights. Their difference measures the response to that model.
    divided = []
    for curve in curves:
        profile = transit_model(curve.time.value, theta, exposure=exposure)
        baseline = curve.copy()
        baseline.flux /= profile
        baseline.flux_err /= profile
        divided.append(baseline)
    raw_base, flat_base = clean_and_flatten(divided, cfg['duration'], mask_ephemeris)
    time, observed, base = np.intersect1d(flat.time.value, flat_base.time.value, return_indices=True)
    inside = transit_mask(time, result['period'], result['t0'], result['duration'])
    weights = flat.flux_err.value[observed]**-2

    def box_depth(flux):
        return np.average(flux[~inside], weights=weights[~inside]) - np.average(flux[inside], weights=weights[inside])

    model = transit_model(time, theta, exposure=exposure)
    dense_time = theta[1] + np.linspace(-result['duration'], result['duration'], 1001)
    model_peak = 1 - transit_model(dense_time, theta, exposure=exposure).min()
    raw_response = box_depth(raw.flux.value[observed] - raw_base.flux.value[base])
    flat_response = box_depth(flat.flux.value[observed] - flat_base.flux.value[base])
    metrics = dict(target=result['target'], cadences_compared=len(time),
                   model_peak_ppm=model_peak * 1e6, model_box_ppm=box_depth(model) * 1e6,
                   raw_response_ppm=raw_response * 1e6, detrended_response_ppm=flat_response * 1e6,
                   detrending_loss_percent=100 * (1 - flat_response / raw_response))
    if result['target'] == 'WASP-18':
        metrics.update(wasp18_occultation(raw, theta))
    return metrics


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    parser.add_argument('--output', default='results/research/depth_diagnostic.csv')
    args = parser.parse_args()
    results = json.loads(args.results.read_text(encoding='utf-8'))
    pd.DataFrame([measure(result) for result in results]).to_csv(args.output, index=False)
