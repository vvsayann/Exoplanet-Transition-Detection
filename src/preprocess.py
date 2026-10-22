"""Segment-wise normalization and transit-preserving detrending."""
import numpy as np
import lightkurve as lk


def transit_mask(time, period, t0, duration):
    return abs((np.asarray(time) - t0 + period / 2) % period - period / 2) < duration / 2


def clean_and_flatten(curves, duration=0.15, ephemeris=None):
    if duration <= 0:
        raise ValueError('duration must be positive (days)')
    raw, flat = [], []
    for curve in curves:
        lc = curve.remove_nans().normalize()
        good = np.isfinite(lc.flux_err.value) & (lc.flux_err.value > 0)
        if hasattr(lc, 'quality'):
            good &= lc.quality == 0
        lc = lc[good].remove_outliers(sigma_lower=float('inf'), sigma_upper=5)
        lc = lc[np.argsort(lc.time.value)]
        if len(lc) < 9:
            raise ValueError('Insufficient valid cadences in a segment')
        cadence = np.median(np.diff(lc.time.value))
        if cadence <= 0:
            raise ValueError('Cadences must have distinct times')
        window = max(5, int(np.ceil(3 * duration / cadence)) | 1)
        window = min(window, len(lc) if len(lc) % 2 else len(lc) - 1)
        mask = None if ephemeris is None else transit_mask(lc.time.value, *ephemeris)
        detrended = lc.flatten(window_length=window, mask=mask, break_tolerance=5)
        raw.append(lc)
        flat.append(detrended)
    raw = lk.LightCurveCollection(raw).stitch()
    flat = lk.LightCurveCollection(flat).stitch()
    order = np.argsort(flat.time.value)
    raw, flat = raw[order], flat[order]
    unique = np.r_[True, np.diff(flat.time.value) > 0]
    return raw[unique], flat[unique]
