"""Circular, quadratic limb-darkened model; conditional Bayesian inference."""
import numpy as np
import batman
import emcee
from time import perf_counter
from scipy.optimize import least_squares

PARAMETERS = ['period_days', 't0_days', 'radius_ratio', 'a_over_rstar', 'impact', 'offset']


def _parameters(theta):
    period, t0, rp, a, b, offset = theta
    params = batman.TransitParams()
    params.t0, params.per, params.rp, params.a = t0, period, rp, a
    params.inc = np.degrees(np.arccos(np.clip(b / a, -1, 1)))
    params.ecc, params.w = 0., 90.
    params.u, params.limb_dark = [.3, .2], 'quadratic'
    return params


def transit_model(time, theta, exposure=None):
    params = _parameters(theta)
    exposure = float(exposure or np.median(np.diff(np.sort(time))))
    model = batman.TransitModel(params, np.asarray(time), supersample_factor=5, exp_time=exposure)
    return model.light_curve(params) + theta[-1]


def summarize_radius_earth(samples, stellar_radius, stellar_radius_error, rng):
    if stellar_radius <= 0 or stellar_radius_error < 0:
        raise ValueError('Invalid stellar radius or uncertainty')
    radii = samples[:, 2] * rng.normal(stellar_radius, stellar_radius_error, len(samples)) * 109.076
    lo, mid, hi = np.percentile(radii, [16, 50, 84])
    return dict(median=float(mid), minus=float(mid-lo), plus=float(hi-mid))


def circular_log_likelihood(time, flux, flux_error, exposure):
    """Evaluate the same Gaussian likelihood using batman only near transits."""
    time = np.asarray(time, dtype=float)
    residual = np.asarray(flux, dtype=float) - 1.
    weights = np.asarray(flux_error, dtype=float)**-2
    baseline_chi2 = np.dot(weights, residual**2)
    baseline_sum = np.dot(weights, residual)
    weight_sum = weights.sum()

    def log_likelihood(theta):
        period, t0, radius, a, impact, offset = theta
        dt = time - t0
        phase = dt - period * np.rint(dt / period)
        # From z² = b² + (a² - b²) sin²(2π phase / P) at first contact.
        ratio = ((1 + radius)**2 - impact**2) / (a**2 - impact**2)
        half_duration = period / (2 * np.pi) * np.arcsin(np.sqrt(min(1., ratio)))
        affected = np.abs(phase) <= half_duration + exposure / 2 + 1e-10
        chi2 = baseline_chi2 - 2 * offset * baseline_sum + offset**2 * weight_sum
        if affected.any():
            params = _parameters(theta)
            model = batman.TransitModel(params, time[affected], supersample_factor=5, exp_time=exposure)
            dip = 1 - model.light_curve(params)
            chi2 += np.dot(weights[affected], dip * (dip + 2 * (residual[affected] - offset)))
        return -.5 * chi2

    return log_likelihood


def fit_transit(lc, candidate, steps=1500, burn=500, walkers=32, seed=42,
                stellar_radius=None, stellar_radius_error=0., progress=False,
                differential_evolution=False):
    if steps <= burn or walkers < 12:
        raise ValueError('Require steps > burn and at least 12 walkers')
    rng = np.random.default_rng(seed)
    p, t0, dur = (candidate[k] for k in ['period', 't0', 'duration'])
    t, y, dy = lc.time.value.copy(), lc.flux.value, lc.flux_err.value
    # Shift the origin for numerical conditioning; restore it in exported samples.
    origin = float(t.min())
    t -= origin
    t0 -= origin
    guess = [p, t0, np.sqrt(max(candidate['depth'], 1e-6)), max(2., p / (np.pi * dur)), .3, 0.]
    lower = np.array([p*.95, t0-dur, .0001, 1.1, 0., -.05])
    upper = np.array([p*1.05, t0+dur, .5, 500., .99, .05])
    guess = np.clip(guess, lower + 1e-7, upper - 1e-7)
    exposure = np.median(np.diff(t))
    model = batman.TransitModel(_parameters(guess), t, supersample_factor=5, exp_time=exposure)
    def residual(theta):
        return (y - model.light_curve(_parameters(theta)) - theta[-1]) / dy
    if progress:
        print(f'Physical fit: optimizing 6 parameters using {len(t):,} cadences...', flush=True)
    optimization_start = perf_counter()
    opt = least_squares(residual, guess, bounds=(lower, upper), x_scale='jac')
    if not opt.success:
        raise RuntimeError(f'Initial transit optimization failed: {opt.message}')
    if progress:
        print(f'Initial optimization finished in {perf_counter()-optimization_start:.1f}s.', flush=True)
    log_likelihood = circular_log_likelihood(t, y, dy, exposure)

    def log_prob(theta):
        if np.any(theta <= lower) or np.any(theta >= upper):
            return -np.inf
        return log_likelihood(theta)

    moves = None
    if differential_evolution:
        moves = [(emcee.moves.DEMove(), .8), (emcee.moves.DESnookerMove(), .2)]
    # emcee uses NumPy's legacy RNG; explicitly seed its own random state.
    sampler = emcee.EnsembleSampler(walkers, 6, log_prob, moves=moves)
    sampler.random_state = np.random.RandomState(seed).get_state()
    start = opt.x + rng.normal(size=(walkers, 6)) * (upper-lower) * 1e-5
    start = np.clip(start, lower+1e-9, upper-1e-9)
    if progress:
        print(f'MCMC: {steps:,} steps, {walkers} walkers; progress and estimated remaining time follow.', flush=True)
    sampler.run_mcmc(start, steps, progress=progress)
    if progress:
        print('MCMC finished. Computing posterior summaries and convergence diagnostics...', flush=True)
    samples = sampler.get_chain(discard=burn, flat=True)
    samples[:, 1] += origin
    summary = {}
    for i, name in enumerate(PARAMETERS):
        lo, mid, hi = np.percentile(samples[:, i], [16, 50, 84])
        summary[name] = {'median': float(mid), 'minus': float(mid-lo), 'plus': float(hi-mid)}
    if stellar_radius is not None:
        summary['radius_earth'] = summarize_radius_earth(samples, stellar_radius, stellar_radius_error, rng)
    try:
        tau = sampler.get_autocorr_time(discard=burn)
    except emcee.autocorr.AutocorrError as exc:
        tau = exc.tau
    finite_tau = np.isfinite(tau) & (tau > 0)
    converged = bool(finite_tau.all() and np.all((steps-burn) > 50*tau))
    summary['diagnostics'] = dict(converged=converged,
                                  autocorrelation_time=[float(x) if valid else None for x, valid in zip(tau, finite_tau)],
                                  autocorrelation_basis='Retained chain after burn-in',
                                  acceptance_fraction=float(np.mean(sampler.acceptance_fraction)),
                                  steps=steps, burn=burn, walkers=walkers,
                                  moves='Differential evolution and snooker' if differential_evolution else 'Stretch',
                                  assumptions='Circular orbit; fixed limb darkening [0.3,0.2]; independent Gaussian noise; uniform bounded priors')
    return summary, samples
