import unittest
import numpy as np
import lightkurve as lk
from astropy.utils.masked import Masked
from src.experiments import synthetic_curve
from src.preprocess import clean_and_flatten, transit_mask
from src.search import search
from src.vetting import vet
from src.fit import circular_log_likelihood, fit_transit, transit_model


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

    def test_secondary_depth_limit_still_rejects_deep_eclipse(self):
        t = np.arange(0, 30, .02)
        primary = transit_mask(t, 2., .5, .16)
        secondary = transit_mask(t, 2., 1.5, .16)
        c = dict(period=2., t0=.5, duration=.16, depth=.02, depth_snr=100, n_transits=15)
        for secondary_depth in (.0003, .003):
            flux = 1 - .02 * primary - secondary_depth * secondary
            lc = lk.LightCurve(time=t, flux=flux, flux_err=np.full(len(t), .0001))
            with self.subTest(secondary_depth=secondary_depth):
                self.assertFalse(vet(lc, c)['flags']['secondary'])
                self.assertEqual(vet(lc, c, secondary_depth_max=.0004)['flags']['secondary'],
                                 secondary_depth < .0004)

    def test_sparse_likelihood_matches_full_exposure_integrated_model(self):
        t = np.arange(0., 13., .02)
        flux = np.random.default_rng(21).normal(1., .002, len(t))
        error = .002 + .001 * t
        likelihood = circular_log_likelihood(t, Masked(flux), Masked(error), .02)
        for theta in ([2., .31, .14, 6., .5, .001], [.51, .013, .48, 1.11, .98, -.005]):
            with self.subTest(theta=theta):
                full = -.5 * np.sum(((flux - transit_model(t, theta, exposure=.02)) / error)**2)
                self.assertAlmostEqual(likelihood(theta), full, places=8)

    def test_short_baseline_rejected(self):
        lc = synthetic_curve(days=.8)
        with self.assertRaises(ValueError):
            search(lc)


if __name__ == '__main__':
    unittest.main()
