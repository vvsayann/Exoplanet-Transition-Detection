import unittest
import numpy as np
import lightkurve as lk
from src.experiments import synthetic_curve
from src.preprocess import clean_and_flatten, transit_mask
from src.search import search
from src.vetting import vet
from src.fit import fit_transit


class PipelineTests(unittest.TestCase):
    def test_period_and_uncertainties(self):
        _, flat = clean_and_flatten([synthetic_curve()], duration=.15)
        c, _ = search(flat, max_period=8)
        self.assertLess(abs(c['period']/3-1), .01)
        self.assertGreater(c['depth_snr'], 7)
        summary, samples = fit_transit(flat, c, steps=80, burn=30)
        self.assertEqual(samples.shape, (1600, 6))
        self.assertGreater(summary['radius_ratio']['plus'], 0)
        self.assertFalse(summary['diagnostics']['converged'])

    def test_downward_dips_preserved(self):
        t = np.arange(0,10,.02); flux = np.ones(len(t))
        flux[100:104] = .95; flux[300] = 1.3
        lc = lk.LightCurve(time=t, flux=flux, flux_err=np.full(len(t), .001))
        raw, flat = clean_and_flatten([lc], ephemeris=(5,2.03,.16))
        self.assertLess(raw.flux.value.min(), .96)
        self.assertLess(raw.flux.value.max(), 1.1)

    def test_binary_odd_even(self):
        t = np.arange(0,30,.02); p = 2.; t0 = .5
        mask = transit_mask(t,p,t0,.16)
        flux = np.ones(len(t)); epoch = np.round((t-t0)/p).astype(int)
        flux[mask] -= np.where(epoch[mask]%2==0,.02,.005)
        lc = lk.LightCurve(time=t, flux=flux, flux_err=np.full(len(t), .001))
        c = dict(period=p,t0=t0,duration=.16,depth=.012,depth_snr=100,n_transits=15)
        self.assertFalse(vet(lc,c)['flags']['odd_even'])

    def test_short_baseline_rejected(self):
        lc = synthetic_curve(days=.8)
        with self.assertRaises(ValueError):
            search(lc)


if __name__ == '__main__':
    unittest.main()
