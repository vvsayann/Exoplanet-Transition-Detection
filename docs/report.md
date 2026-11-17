# A reproducible pipeline for detecting and screening exoplanet transits

**Implementation and verification report — 22 October 2026**

## Abstract

This project implements a modular transit detection workflow for Kepler and TESS photometry. It combines systematics-corrected archive light curves, segment-wise cleaning, Box Least Squares (BLS), exposure-integrated physical transit models, ensemble Markov chain Monte Carlo, and candidate screening. Experiment drivers support injection recovery, labeled false-positive evaluation, and searches whose detections are frozen before catalog matching. An offline numerical demonstration recovered an injected 3-day period at 2.99943475 days, an error of approximately 0.019%. A real Kepler-10 smoke test recovered a 0.83771248-day signal. Twenty-seven synthetic injection trials demonstrate how the completeness products are generated. These checks establish useful software behavior, but do not establish survey completeness, catalog-wide reliability, converged parameter uncertainties, or discovery of a new planet. The distinction between implemented capability and completed observational evidence is maintained throughout.

## 1. Introduction

Transit photometry searches for repeated decreases in stellar brightness when an orbiting object crosses the stellar disk. Periodicity alone cannot establish planetary nature: stellar companions, background eclipsing binaries, instrumental artifacts, and activity can create similar signatures. A useful detection project therefore requires more than a periodogram. It must preserve shallow dips during preprocessing, quantify parameter uncertainty, examine plausible false-positive explanations, and measure the probability of detecting signals under specified observing conditions.

The initial repository contained a single-target demonstration, a notebook, and three figures. It searched a narrow period range and used symmetric outlier rejection, which could remove downward transit events. The expanded implementation separates archive access, detrending, search, fitting, screening, and experiments. This makes individual assumptions inspectable and allows the same analysis to be applied to multiple targets without manually editing a script.

The main intended use is an undergraduate or early postgraduate research workflow. It produces candidate measurements and diagnostic evidence. Confirmation, occurrence-rate inference, and population-level completeness require additional observations, stronger noise models, explicit sample selection, and substantially larger validation experiments.

## 2. Method

### 2.1 Data and preprocessing

The downloader requests Kepler or SPOC TESS light curves through Lightkurve, selects PDCSAP flux, and caches downloaded FITS products. Each segment is processed separately before stitching. Invalid fluxes and uncertainties are rejected; quality flags are checked; normalization removes segment-dependent scaling. Sigma clipping removes upward outliers while retaining downward dips. Segments are sorted, duplicate times removed, and long gaps passed to Lightkurve's gap-aware flattening.

The detrending window is approximately three times an assumed transit duration, expressed in cadences. This is a configurable physical scale rather than a fixed sample count. A preliminary BLS signal supplies a transit mask for a second detrending pass, protecting likely transit cadences during trend estimation. Before/after figures mark detected transit windows for inspection. A window derived from an inaccurate duration can still distort events; the software cannot replace visual inspection or injection tests over duration and variability regimes.

### 2.2 Periodic search

Astropy BLS searches a duration grid and a period interval beginning at 0.5 days by default. The maximum period is capped at half the observed baseline. Results include period, epoch, duration, depth, depth SNR, power, and the number of distinct observed transit events. Detection requires depth SNR at least seven and at least three observed events. The event count is evaluated from actual in-transit cadences, so a long baseline with gaps does not by itself guarantee adequate coverage.

Half-period and double-period scores are exported as alias diagnostics. They are not treated as interchangeable successful recoveries in completeness calculations. Additional signals are sought by masking windows around detected events and rerunning BLS. This approach can miss overlapping transits and should be followed by a joint analysis for credible multi-planet candidates.

### 2.3 Physical modeling and uncertainty

The physical model uses batman with quadratic limb darkening and numerical exposure integration. Fitted parameters are orbital period, mid-transit epoch, planet-to-star radius ratio, scaled semimajor axis, impact parameter, and a flux offset. A bounded least-squares solution initializes emcee walkers. A shifted time origin improves numerical conditioning, while exported epochs are returned to the original light-curve convention.

Reported marginal intervals use posterior percentiles 16, 50, and 84. When stellar radius and uncertainty are supplied, posterior draws propagate that uncertainty into Earth-radius estimates. Stellar uncertainty is not inferred from the photometry. The likelihood assumes independent Gaussian errors; the orbit is circular; limb darkening is fixed to coefficients 0.3 and 0.2. Those choices make the example reproducible but can underestimate uncertainty or bias geometry for real stars. They must be revisited for publication-quality characterization.

Acceptance fraction and autocorrelation diagnostics accompany every fit. A chain is marked converged only when retained length exceeds fifty estimated autocorrelation times. Failure to estimate reliable autocorrelation time leaves convergence false. A corner plot or a finite numerical interval is therefore not automatically evidence of reliable inference.

### 2.4 Candidate screening

The screening module compares odd and even depths, searches for a secondary dip at phase 0.5, compares depths across individual events, and evaluates a heuristic edge-to-core depth ratio. With a supplied stellar density it also compares duration against a circular central-transit estimate. Missing information remains explicitly unassessed rather than being converted into a passing test.

These diagnostics are deliberately interpretable. They are not a validated false-positive probability model. The depth errors use formal photometric uncertainties and do not fully model uncertainty in the shared baseline or correlated stellar noise. Eccentric binaries can evade a phase-0.5 secondary search. Low signal-to-noise transits can produce unreliable shape ratios. These effects motivate evaluation on labeled samples with both retention and rejection metrics.

### 2.5 Evaluation design

Injection recovery multiplies synthetic batman transits into supplied light curves before preprocessing. Trials vary period, planet radius, and random epoch; the default geometry fixes impact parameter to 0.3. An exact recovery requires a period within 1%, an epoch within one recovered duration, adequate SNR, and at least three observed events. Harmonics are recorded separately. Recovery fractions include Wilson 95% binomial intervals so small or all-success samples do not imply certainty.

Known-target comparisons query the NASA archive and retain published values alongside recovered periods, depths, and durations. Labeled evaluation counts false positives rejected and confirmed planets retained. It reports matched and requested coverage because failed downloads can bias performance. Blind-search results are saved before fetching KOI or TOI labels; crossmatching uses host identifiers and period agreement. The supplied identifier manifest is a demonstration sample, not a probability-selected survey.

## 3. Verified results

Four automated tests passed in the prepared Python 3.12 environment. They check period recovery, posterior production and unconverged-chain labeling, preservation of downward dips during cleaning, rejection of unequal odd/even depths, and rejection of an insufficient baseline.

The end-to-end synthetic demonstration used a 3-day signal, radius ratio 0.08, a 27-day baseline, and Gaussian noise with amplitude 300 parts per million plus slow variability. BLS returned 2.99943475 days, depth approximately 0.006675, and depth SNR approximately 132. The 300-step posterior run produced the expected artifacts, including a model overlay and corner plot. Its convergence flag is false; its narrow intervals must not be used as validated scientific uncertainties.

The first real-data check downloaded Kepler-10 photometry and processed its first available segment. It found a 0.83771248-day signal, demonstrating that archive data can pass through cleaning, search, vetting, and plotting. This smoke test alone does not verify literature agreement across all parameters or a broad range of planet sizes.

A NASA archive comparison retrieved during verification gives Kepler-10 b a period of 0.8374907 days. The smoke-test relative error is 0.02648%, meeting the proposed 1% period criterion for this one signal. The recovered box depth is 187.13 ppm versus a catalog depth of 191.9 ppm; the trial-grid duration is 1.92 hours versus 1.811 hours. Box measurements are not substitutes for uncertainties from a converged physical fit.

The synthetic completeness demonstration ran three random phases in each of nine period-radius bins: periods 1, 3, and 7 days and radii 2, 4, and 8 Earth radii. The 2-Earth-radius signal was recovered in one of three 1-day trials and none of the 3-day or 7-day trials. All 4- and 8-Earth-radius injections were recovered. With three trials, a bin containing three successes still has a Wilson lower 95% bound of about 0.438. These are demonstration measurements on synthetic noise, not an empirical sensitivity claim about bright Kepler or TESS stars.

The authoritative numeric outputs are the JSON/CSV artifacts under `results/`. Archive experiments may add further evidence; incomplete batches and unconverged fits must not be presented as meeting the original research success metrics.

An additional 18 injections used the first available real segment of KIC 8006161, with provisional solar mass and radius assumptions. At 1 day all three phases recovered in every radius bin. At 3 days the recovery counts were 1/3, 2/3, and 3/3 for 2, 4, and 8 Earth radii respectively. Event coverage, formal errors, stellar properties, and variability must be inspected before treating this as quiet-star completeness; the short baseline does not support a general population claim. These products are saved separately in `results/real_injection/`.

## 4. Discussion and remaining scientific validation

The pipeline now supports the proposed analyses in reusable modules, but implementation completeness differs from completion of a scientific study. Acceptance requires successful runs for several known planets, documented literature comparisons with uncertainties, converged posteriors, and a labeled sample large enough to evaluate rejection without excessive loss of real planets. Completeness needs independently inspected quiet stars, measured stellar properties, more random phases, multiple geometries, and sufficient trials per bin. A blind search needs completed target coverage, explicit selection criteria, ranked detections, and post-search catalog reconciliation.

The main methodological risks are detrending loss, aliases, inaccurate formal errors, fixed limb darkening, correlated noise, dilution, and selection bias. The current outputs expose several of these through diagnostic plots, alias scores, incomplete-test indicators, posterior diagnostics, and experiment counts. Further work could add centroid analysis, joint stellar-density constraints, eccentric models, transit timing variation, or TLS. None should be substituted for careful inspection of the present validation evidence.

## References

Kovács, Zucker & Mazeh (2002), BLS; Kreidberg (2015), batman; Foreman-Mackey et al. (2013), emcee. Software/API references: [Lightkurve flatten](https://lightkurve.github.io/lightkurve/reference/api/lightkurve.LightCurve.flatten.html), [batman](https://lkreidberg.github.io/batman/), [NASA archive TAP](https://exoplanetarchive.ipac.caltech.edu/docs/API_resources.html).
