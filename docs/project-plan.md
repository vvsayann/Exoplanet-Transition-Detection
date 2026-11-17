# Exoplanet Transit Detection Pipeline: Complete Project Plan

Detect, characterize, and vet exoplanet transit signals in Kepler and TESS light curves, and measure how well the pipeline performs.

**Duration:** 6 weeks (4-week core + 2-week Phase 2)
**Level:** Advanced undergraduate / early masters
**Primary language:** Python

---

## 1. Project Summary

When a planet passes in front of its star, the star's brightness dips slightly and periodically. This project builds an end-to-end pipeline that:

1. Downloads and cleans light curves from Kepler and TESS
2. Searches for periodic transit signals using Box Least Squares (BLS)
3. Fits a physical transit model to get planet parameters **with uncertainties**
4. Vets candidates to reject false positives such as eclipsing binaries
5. Quantifies its own sensitivity with injection-recovery tests
6. Runs a blind search on stars where the answer is not known in advance

What makes it a detection project rather than a tutorial: steps 3 to 6.

---

## 2. Learning Goals

- Time-series processing and detrending of real astronomical data
- Periodogram-based signal detection
- Bayesian parameter estimation (MCMC)
- Distinguishing real signals from astrophysical false positives
- Measuring pipeline completeness and reliability
- Reproducible, documented scientific code

---

## 3. Tech Stack

| Purpose | Tool |
|---|---|
| Data access and light curve handling | `lightkurve` |
| BLS periodogram | `astropy.timeseries.BoxLeastSquares` |
| Alternative detrending (optional) | `wotan` |
| Transit model | `batman-package` |
| MCMC fitting | `emcee`, `corner` |
| Catalog access | `astroquery` (NASA Exoplanet Archive, MAST) |
| Analysis and plots | `numpy`, `pandas`, `matplotlib` |
| Optional search upgrade | `transitleastsquares` (TLS) |

```bash
pip install lightkurve astropy astroquery batman-package emcee corner wotan pandas matplotlib
```

---

## 4. Suggested Repository Structure

```
exoplanet-transit-detection/
├── README.md
├── requirements.txt
├── notebooks/
│   ├── 01_data_and_preprocessing.ipynb
│   ├── 02_bls_search.ipynb
│   ├── 03_parameter_fitting.ipynb
│   ├── 04_vetting.ipynb
│   ├── 05_injection_recovery.ipynb
│   └── 06_blind_search.ipynb
├── src/
│   ├── download.py
│   ├── preprocess.py
│   ├── search.py
│   ├── fit.py
│   ├── vetting.py
│   └── pipeline.py
├── data/            # cached light curves (gitignored)
├── results/
│   ├── figures/
│   └── tables/
└── docs/
```

---

## 5. Timeline

### Phase 1: Core Pipeline (Weeks 1 to 4)

#### Week 1: Data and Preprocessing

**Tasks**
- Use `lightkurve` to search MAST and download light curves for a known planet host (start with Kepler-10 or a bright TESS Object of Interest).
- Use PDCSAP flux (systematics-corrected) rather than raw SAP flux.
- Remove NaNs and bad-quality cadences, and sigma-clip upward outliers only (clipping downward can remove real transits).
- Normalize flux, stitch quarters (Kepler) or sectors (TESS).
- Detrend stellar variability with `flatten()` or `wotan` (biweight filter). Choose a window length of about 3 times the expected transit duration.
- Plot raw vs. flattened for sanity.

**Pitfall:** An over-aggressive window flattens the transit itself. Always overlay the known transit times on the flattened curve.

**Checkpoint:** Clean, flattened light curve for one known planet, with a before/after plot.

#### Week 2: Transit Search

**Tasks**
- Run BLS across a period grid (e.g. 0.5 to 50 days) with a sensible set of trial durations.
- Plot the periodogram and identify the highest peak.
- Phase-fold at the best period and inspect the folded transit shape.
- Compare recovered period, depth, and duration against published values (NASA Exoplanet Archive).
- Check for aliases: peaks at half or double the true period are common.

**Checkpoint:** Periodogram with a clear peak, folded light curve matching the published period within a small tolerance.

#### Week 3: Parameter Fitting

**Tasks**
- Build a transit model with `batman`: period, mid-transit time, planet-to-star radius ratio, scaled semi-major axis, inclination, limb darkening.
- Use BLS output as the initial guess; fit with `scipy.optimize` first, then run `emcee` for posteriors.
- Produce a corner plot and report median and 68% credible intervals.
- Convert radius ratio to planet radius using the stellar radius from the TESS Input Catalog or Kepler stellar table.
- Compare against published values and note any discrepancies.

**Checkpoint:** Fitted parameters with uncertainties for at least one planet, including a corner plot and model-over-data plot.

#### Week 4: Generalization and Packaging (Core)

**Tasks**
- Wrap the steps into a function: `run_pipeline(target_id) -> period, depth, duration, SNR, plots`.
- Test on 3 to 5 known planets spanning different sizes and periods (e.g. a hot Jupiter, a super-Earth, a TESS TOI).
- Write the first README with method and example plots.

**Checkpoint:** One-command pipeline that works on multiple targets, with a results table comparing recovered vs. published parameters.

---

### Phase 2: Making It a Real Detection Project (Weeks 5 to 6)

#### Week 5: Vetting and False-Positive Rejection

**Tasks**
- **Odd/even test:** Fold odd and even transits separately. A significant depth difference suggests an eclipsing binary at twice the period.
- **Secondary eclipse search:** Look for a dip at phase 0.5, which indicates a self-luminous companion.
- **Transit shape and duration check:** Compare duration to what the stellar density predicts; implausibly long or V-shaped dips are suspicious.
- **Depth consistency:** Compare depth across individual transits.
- **Multi-planet search:** Mask the detected transits and rerun BLS to look for additional planets in the same system.
- Build a labeled test set using Kepler KOI dispositions from the NASA Exoplanet Archive (confirmed planets vs. false positives) and measure how well your vetting tests separate them.

**Checkpoint:** Vetting module that outputs pass/fail flags per test, evaluated on a labeled sample of confirmed planets and known false positives.

#### Week 6: Injection-Recovery and Blind Search

**Tasks**
- **Injection-recovery:** Inject synthetic `batman` transits into real light curves from quiet stars across a grid of planet radii and periods. Run the pipeline and record whether each is recovered.
- Plot recovery rate as a heatmap (radius vs. period). This is your pipeline's completeness map.
- **Blind search:** Run the pipeline on a set of stars (e.g. a few dozen TESS targets or Kepler stars) without looking at known results first. Rank signals by SNR and vetting flags.
- Afterward, cross-match against known TOIs/KOIs: how many did you recover, and are any flagged signals not in the catalogs?
- Finalize README, figures, and a short written report.

**Checkpoint:** Completeness heatmap, blind-search results table cross-matched with catalogs, and a finished write-up.

---

## 6. Deliverables

- Clean GitHub repository with the structure above
- README covering method, results, limitations, and how to reproduce
- Results table: recovered vs. published parameters for test planets
- Figures: raw vs. flattened curve, periodogram, folded transit, corner plot, completeness heatmap
- Short report (3 to 5 pages) written like a mini research paper: Introduction, Method, Results, Discussion
- Optional: Streamlit app where a user enters a TIC/KIC ID and gets the periodogram and folded curve

---

## 7. Success Metrics

- Recovers published periods to within about 1% for known test planets
- Reported parameters include uncertainties and agree with literature within errors
- Vetting tests correctly flag a majority of known false positives in the labeled set
- Injection-recovery gives a quantified completeness (e.g. "recovers more than 90% of 2 Earth-radius planets with P < 10 d around bright stars")
- Anyone can clone the repo and reproduce the main figures

---

## 8. Common Pitfalls

- Flattening away the transit with too short a detrending window
- Using SAP instead of PDCSAP flux without correcting systematics
- Trusting a BLS peak without checking aliases (P/2, 2P)
- Ignoring TESS sector gaps and momentum-dump artifacts
- Reporting parameters without uncertainties
- Calling something a "new planet" without vetting; at most it is a candidate

---

## 9. Optional Extensions

- Swap BLS for Transit Least Squares (TLS) and compare sensitivity
- Train a small classifier (random forest or 1D CNN) on folded light curves to separate planets from false positives
- Analyze transit timing variations (TTVs) for multi-planet systems
- Add centroid analysis using target pixel files to test for background contamination

---

## 10. Key References and Resources

- Lightkurve documentation: https://lightkurve.github.io/lightkurve/
- NASA Exoplanet Archive (KOI and TOI tables): https://exoplanetarchive.ipac.caltech.edu/
- MAST archive: https://archive.stsci.edu/
- Kovács, Zucker & Mazeh (2002), the original BLS paper
- Hippke & Heller (2019), Transit Least Squares
- Kreidberg (2015), batman: Bad-Ass Transit Model cAlculatioN
- Foreman-Mackey et al. (2013), emcee

---

## 11. How to Begin Today

```python
import lightkurve as lk

search = lk.search_lightcurve("Kepler-10", mission="Kepler", author="Kepler")
print(search)

lc = search[0].download().remove_nans().normalize()
lc.plot()
```

Download one quarter, plot it, and you have started Week 1.
