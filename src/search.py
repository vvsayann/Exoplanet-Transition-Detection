"""Weighted Box Least Squares search with alias diagnostics."""
import numpy as np
from astropy.timeseries import BoxLeastSquares
from .preprocess import transit_mask


def search(lc, min_period=0.5, max_period=50., durations=None):
    t, y, dy = lc.time.value, lc.flux.value, lc.flux_err.value
    maximum = min(max_period, np.ptp(t) / 2)
    if maximum <= min_period:
        raise ValueError('Baseline must cover at least two minimum periods')
    durations = np.array([.04, .08, .12, .20, .30] if durations is None else durations)
    durations = durations[(durations > 0) & (durations < min_period)]
    if not len(durations):
        raise ValueError('Trial durations must be shorter than minimum period')
    model = BoxLeastSquares(t, y, dy)
    power = model.autopower(durations, minimum_period=min_period,
                           maximum_period=maximum, frequency_factor=3)
    i = np.argmax(power.power)
    candidate = {key: float(power[key][i]) for key in
                 ('period', 'transit_time', 'duration', 'depth', 'depth_snr', 'power')}
    candidate['t0'] = candidate.pop('transit_time')
    mask = transit_mask(t, candidate['period'], candidate['t0'], candidate['duration'])
    candidate['n_transits'] = int(len(np.unique(np.round((t[mask] - candidate['t0']) / candidate['period']))))
    candidate['aliases'] = []
    for factor in [.5, 1., 2.]:
        p = candidate['period'] * factor
        result = model.power(p, durations[durations < p])
        j = np.argmax(result.power)
        candidate['aliases'].append({'period': p, 'power': float(result.power[j]),
                                     'snr': float(result.depth_snr[j])})
    return candidate, power


def search_multiple(lc, max_signals=3, min_snr=7., **kwargs):
    remaining, signals = lc.copy(), []
    for _ in range(max_signals):
        candidate, _ = search(remaining, **kwargs)
        if candidate['depth_snr'] < min_snr or candidate['n_transits'] < 3:
            break
        signals.append(candidate)
        remaining = remaining[~transit_mask(remaining.time.value, candidate['period'],
                                          candidate['t0'], candidate['duration'] * 1.5)]
    return signals
