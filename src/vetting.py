"""Screening diagnostics; a pass does not establish planetary nature."""
import numpy as np
from scipy.stats import chi2
from .preprocess import transit_mask


def vet(lc, c, stellar_density=None):
    t, y, e = lc.time.value, lc.flux.value, lc.flux_err.value
    p, t0, d = (c[k] for k in ['period', 't0', 'duration'])
    phase = (t-t0+p/2) % p-p/2
    inside = abs(phase) < d/2
    outside = abs(phase) > d
    baseline = np.median(y[outside]) if outside.any() else 1.
    def depth(mask):
        if mask.sum() < 3:
            return None
        weights = 1/e[mask]**2
        return float(baseline-np.average(y[mask], weights=weights)), float(1/np.sqrt(weights.sum()))
    epoch = np.round((t-t0)/p).astype(int)
    odd, even = depth(inside & (epoch % 2 == 1)), depth(inside & (epoch % 2 == 0))
    z = None if odd is None or even is None else abs(odd[0]-even[0])/np.hypot(odd[1], even[1])
    secondary = depth(transit_mask(t, p, t0+p/2, d))
    secondary_z = None if secondary is None else secondary[0]/secondary[1]
    depths = [depth(inside & (epoch == k)) for k in np.unique(epoch[inside])]
    depths = [v for v in depths if v is not None]
    consistency_p = None
    if len(depths) >= 3:
        values, errors = np.array(depths).T
        mean = np.average(values, weights=errors**-2)
        consistency_p = float(chi2.sf(np.sum(((values-mean)/errors)**2), len(values)-1))
    core, edge = depth(abs(phase)<d/6), depth((abs(phase)>d/3) & inside)
    shape = None if core is None or edge is None or core[0] <= 0 else edge[0]/core[0]
    # Circular central transit upper-duration approximation, density in g/cm^3.
    expected = None
    if stellar_density is not None:
        if stellar_density <= 0:
            raise ValueError('stellar_density must be positive')
        a = (6.6743e-11 * stellar_density*1000 * (p*86400)**2/(3*np.pi))**(1/3)
        expected = float(p/np.pi*np.arcsin(min(1., (1+np.sqrt(c['depth']))/a)))
    metrics = dict(odd_even_sigma=z, secondary_sigma=secondary_z,
                   depth_consistency_p=consistency_p, edge_core_ratio=shape,
                   central_duration_days=expected, individual_depths=depths)
    flags = dict(odd_even=None if z is None else bool(z < 3),
                 secondary=None if secondary_z is None else bool(secondary_z < 5),
                 depth_consistency=None if consistency_p is None else bool(consistency_p > .01),
                 shape=None if shape is None else bool(shape > .25),
                 duration=None if expected is None else bool(d < expected*1.5),
                 detection=bool(c['depth_snr'] >= 7 and c['n_transits'] >= 3))
    return dict(metrics=metrics, flags=flags, passed=all(v for v in flags.values() if v is not None),
                incomplete=any(v is None for v in flags.values()))
