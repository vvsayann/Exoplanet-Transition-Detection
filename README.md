# Exoplanet Transit Detection

A Python research pipeline for Kepler and TESS PDCSAP photometry: data cleaning, BLS detection, physical transit fitting, Bayesian uncertainties, candidate screening, injection recovery, and catalog comparison.

## Setup

Use **Python 3.12** (Windows binary wheels are available for batman). The original Python 3.14 environment requires Microsoft C++ build tools for batman.

```powershell
py -3.12 -m venv .venv-etd
.\.venv-etd\Scripts\python.exe -m pip install -r requirements.txt
.\.venv-etd\Scripts\python.exe -m unittest discover -s tests -v
```

Run commands from the repository root. The locally prepared `.venv-etd` can be used directly.

## Run

```powershell
# Offline end-to-end demonstration (explicitly synthetic)
.\.venv-etd\Scripts\python.exe -m src demo
# Real data, default four quarters/sectors, 0.5–50 day BLS search
.\.venv-etd\Scripts\python.exe -m src run Kepler-10 --stellar-radius 1.06 --stellar-radius-error 0.02
.\.venv-etd\Scripts\python.exe -m src run "TIC 100100827" --mission TESS
# Help and options
.\.venv-etd\Scripts\python.exe -m src --help
```

The radius arguments above illustrate usage; retrieve and cite appropriate stellar measurements for a scientific analysis. `run_pipeline(target_id)` in `src.pipeline` returns period, depth, duration, depth SNR, vetting, and fit summaries. Each output folder contains JSON, cleaned arrays, posterior samples, and before/after, periodogram, folded-model, and corner figures. Times retain the input mission convention (BKJD/BTJD for typical archive data); never compare epochs without converting conventions.

```powershell
# Known-target validation (three Kepler hosts and a TESS host)
.\.venv-etd\Scripts\python.exe -m src batch config/known_targets.csv --output results/known
.\.venv-etd\Scripts\python.exe -m src catalog pscomppars --output data/planets.csv
.\.venv-etd\Scripts\python.exe -m src compare results/known/frozen_results.json data/planets.csv
# Freeze blind detections before fetching labels
.\.venv-etd\Scripts\python.exe -m src batch config/blind_targets.csv --no-fit --output results/blind
# Full archive workflow: known, blind, catalog crossmatch, labeled screening evaluation
.\.venv-etd\Scripts\python.exe scripts/research.py
# Offline completeness demonstration
.\.venv-etd\Scripts\python.exe -m src inject synthetic --repeats 5
# Injection recovery on two inspected Kepler baselines, 20 trials per bin for each star
.\.venv-etd\Scripts\python.exe -m scripts.quiet_injections 6106415
.\.venv-etd\Scripts\python.exe -m scripts.quiet_injections 6603624
```

For real completeness, call `src.experiments.injection_recovery` with original downloaded segments from inspected quiet stars. The `inject path/to/lightcurve.npz` convenience command uses an already detrended baseline, so it measures conditional search completeness and does not measure the full detrending loss. Set stellar mass/radius to measured values in the Python API. Archive runs need working MAST and NASA archive access and can be lengthy. All failed batch targets are recorded in `failures.json`; inspect coverage before interpreting tables. The blind manifest is an identifier sample, not a guaranteed observed or unbiased survey sample.

## Method and limitations

Each segment is quality-filtered, normalized, clipped upward only, and flattened with a window about three times an assumed duration. A preliminary search supplies a protected transit mask for a second detrending pass. Gaps are handled by Lightkurve's segment splitting. The period search is limited to half the observed baseline, with at least three observed events required for detection. Half/double-period BLS scores are exported for review; no automatic alias is silently selected.

Fitting uses exposure-integrated batman models, bounded least squares, then emcee. Parameters are period, epoch, radius ratio, scaled semimajor axis, impact parameter, and flux offset. Reported 16/50/84 percentiles are conditional on circular orbits, fixed quadratic limb darkening (0.3, 0.2), bounded uniform priors, and independent Gaussian errors. Radius conversion propagates supplied stellar radius uncertainty. Check convergence diagnostics and extend chains; short chains do not establish reliable credible intervals. Strong correlated noise, grazing geometry, and uncertain stellar properties require more sophisticated modeling.

Vetting returns odd/even and per-event depth consistency with an empirical error floor, phase-0.5 secondary significance, a heuristic shape flag, and a circular duration check with stellar density. The known-target manifest supplies host densities and a published occultation bound for WASP-18. Supply these for individual targets with `--stellar-density` in g/cm³ and `--secondary-depth-max` as a fractional depth. Missing diagnostics remain `null` and set the `incomplete` flag. The flags are screening heuristics. Unmatched signals require independent vetting.

Injection trials record exact period/epoch recovery and harmonics, with a recovery heatmap and Wilson 95% binomial intervals. Sensitivity depends on selected stars, noise, radius/mass assumptions, phase, impact parameter and search settings. When the label file supplies `catalog_period`, KOI evaluation requires period agreement within 1% and reports requested, analyzed and matched counts by class alongside recall and retention.

## Reproduction and deliverables

Six staged notebooks accompany the modules; the original exploratory notebook and figures are preserved. See `docs/report.md` for the research write-up and `docs/completion.md` for verified and outstanding deliverables. Generated measurements live in `results/`; cached FITS files and environments are ignored. Freeze the verified environment with `pip freeze` before a publication run and retain catalog query metadata, manifests, seeds, and configuration JSON.

References: [Lightkurve detrending](https://lightkurve.github.io/lightkurve/reference/api/lightkurve.LightCurve.flatten.html), [batman](https://lkreidberg.github.io/batman/), [emcee](https://emcee.readthedocs.io/), [NASA archive TAP](https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html). BLS: Kovács, Zucker & Mazeh (2002); physical model: Kreidberg (2015); ensemble sampling: Foreman-Mackey et al. (2013).

## Interrupted archive downloads

The downloader selects 1800-second Kepler and 120-second TESS SPOC products. Before using cached files it checks FITS structure and expected data length. Incomplete files are preserved with an `.incomplete-...` suffix so a subsequent download can replace them. An archive or network error can still interrupt a fresh transfer; rerunning the command triggers the same cache check. Wait for an existing download to finish before starting another run for the same target and cache.
