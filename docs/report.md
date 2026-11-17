# A reproducible pipeline for detecting and screening exoplanet transits

**Implementation and verification report · 17 November 2026**

## Abstract

This project implements a modular transit detection workflow for Kepler and TESS photometry. It combines systematics-corrected archive light curves, segment-wise cleaning, Box Least Squares (BLS), exposure-integrated physical transit models, ensemble Markov chain Monte Carlo, and candidate screening. Three Kepler fits and one TESS fit met the retained-chain convergence criterion. The updated vetting retained all four known planets using measured stellar densities and an independent occultation bound for WASP-18. All 360 injections on two low-variability Kepler baselines met the recovery criteria on the tested grid. The blind search crossmatched 22 targets, and the labeled evaluation scored 11 period matches from 20 targets. The results establish performance on the specified data and models. Survey-level completeness and reliability remain unmeasured.

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

The screening module compares odd and even depths, searches for a secondary dip at phase 0.5, compares depths across individual events, and evaluates a heuristic edge-to-core depth ratio. A supplied stellar density enables the circular central-transit duration check. An independent upper bound on a planetary occultation can be supplied for a characterized system. The result retains the significance relative to zero and applies the secondary flag to any excess above the bound. Missing diagnostics remain explicitly unassessed.

Per-event depth errors combine the formal weighted error, the robust out-of-transit scatter (scaled median absolute deviation), and a 5% box-depth term in quadrature. The last term allows for uneven sampling of ingress and egress. Odd/even comparisons use the resulting event weights when both parities have usable events. These are screening heuristics. Shared baselines, correlated noise and stellar activity can still undermine the nominal significance, so the labeled experiment reports retention and rejection with period-match coverage.

### 2.5 Evaluation design

Injection recovery multiplies synthetic batman transits into supplied light curves before preprocessing. Trials vary period, planet radius, and random epoch; the default geometry fixes impact parameter to 0.3. An exact recovery requires a period within 1%, an epoch within one recovered duration, adequate SNR, and at least three observed events. Harmonics are recorded separately. Recovery fractions include Wilson 95% binomial intervals so small or all-success samples do not imply certainty.

Known-target comparisons query the NASA archive and retain published values alongside recovered periods, depths, and durations. Labeled evaluation counts false positives rejected and confirmed planets retained. It reports matched and requested coverage because failed downloads can bias performance. Blind-search results are saved before fetching KOI or TOI labels; crossmatching uses host identifiers and period agreement. The supplied identifier manifest is a demonstration sample, not a probability-selected survey.

## 3. Verified results

Eight automated tests passed in Python 3.12 on Windows. They cover period recovery, posterior production and convergence labeling, the exposure-integrated likelihood, cleaning, eclipse screening, short baselines, and FITS cache integrity. The run is recorded in [test_output.txt](test_output.txt).

The end-to-end synthetic demonstration used a 3-day signal, radius ratio 0.08, a 27-day baseline, and Gaussian noise with amplitude 300 parts per million plus slow variability. BLS returned 2.99943475 days, depth approximately 0.006675, and depth SNR approximately 132. The 300-step posterior run produced the expected artifacts, including a model overlay and corner plot. Its convergence flag is false; its narrow intervals must not be used as validated scientific uncertainties.

The first real-data check downloaded Kepler-10 photometry and processed its first available segment. It found a 0.83771248-day signal, demonstrating that archive data can pass through cleaning, search, vetting, and plotting. This smoke test alone does not verify literature agreement across all parameters or a broad range of planet sizes.

A NASA archive comparison retrieved during verification gives Kepler-10 b a period of 0.8374907 days. The smoke-test relative error is 0.02648%, meeting the proposed 1% period criterion for this one signal. The recovered box depth is 187.13 ppm versus a catalog depth of 191.9 ppm; the trial-grid duration is 1.92 hours versus 1.811 hours. Box measurements are not substitutes for uncertainties from a converged physical fit.

The synthetic completeness demonstration ran three random phases in each of nine period-radius bins: periods 1, 3, and 7 days and radii 2, 4, and 8 Earth radii. The 2-Earth-radius signal was recovered in one of three 1-day trials and none of the 3-day or 7-day trials. All 4- and 8-Earth-radius injections were recovered. With three trials, a bin containing three successes still has a Wilson lower 95% bound of about 0.438. These are demonstration measurements on synthetic noise, not an empirical sensitivity claim about bright Kepler or TESS stars.

The authoritative numeric outputs are the JSON/CSV artifacts under `results/`. Archive experiments may add further evidence; incomplete batches and unconverged fits must not be presented as meeting the original research success metrics.

An additional 18 injections used the first available real segment of KIC 8006161, with provisional solar mass and radius assumptions. At 1 day all three phases recovered in every radius bin. At 3 days the recovery counts were 1/3, 2/3, and 3/3 for 2, 4, and 8 Earth radii respectively. Event coverage, formal errors, stellar properties, and variability must be inspected before treating this as quiet-star completeness; the short baseline does not support a general population claim. These products are saved separately in `results/real_injection/`.

### 3.1 Known targets and depth measurements

Each Kepler fit used two segments and 32 walkers. The initial 10,000-step runs discarded 2,000 steps. Kepler-17 passed the convergence check at that length. Kepler-10 and Kepler-8 had maximum autocorrelation estimates of approximately 417 and 258 steps, both above the 160-step maximum for 8,000 retained samples per walker. The final comparison uses 20,000 steps and 4,000 burn-in for all three and an 80:20 mixture of emcee's differential-evolution and snooker proposals. The saved WASP-18 samples used 10,000 steps, 2,000 burn-in and the stretch move. Its historical convergence correction is recorded beside its result. [emcee proposal moves](https://emcee.readthedocs.io/en/stable/user/moves/)

| Target | BLS period (d) | Published period (d) | Period error | Box / published depth (ppm) | Grid / published duration (h) | Max autocorrelation (steps) | Fit converged |
|---|---:|---:|---:|---:|---:|---:|---|
| Kepler-10 b | 0.83771828 | 0.83749070 | 0.02717% | 153.2 / 191.9 | 1.92 / 1.811 | 222.6 | Yes |
| Kepler-17 b | 1.48574424 | 1.48571080 | 0.00225% | 17,637.9 / 21,801.9 | 1.92 / 2.276 | 56.5 | Yes |
| Kepler-8 b | 3.52346702 | 3.52249910 | 0.02748% | 7,718.0 / 9,145.7 | 2.88 / 3.191 | 76.4 | Yes |
| WASP-18 b | 0.94141493 | 0.94145223 | 0.00396% | 9,746.3 / 10,410.0 | 1.92 / 2.210 | 71.0 | Yes |

| Target | Fitted radius ratio (16th–84th percentile errors) | Archive radius ratio (quoted errors) | Relative difference |
|---|---:|---:|---:|
| Kepler-10 b | 0.012735 (−0.000165 / +0.000264) | 0.012680 (−0.000040 / +0.000040) | +0.4% |
| Kepler-17 b | 0.133512 (−0.000055 / +0.000058) | 0.130310 (−0.000180 / +0.000220) | +2.5% |
| Kepler-8 b | 0.091711 (−0.000668 / +0.000717) | 0.095751 (−0.000230 / +0.000190) | −4.2% |
| WASP-18 b | 0.097490 (−0.000122 / +0.000124) | 0.101800 (−0.001100 / +0.001100) | −4.2% |

The published values are from the NASA Archive composite table. The comparison file also includes Earth-radius estimates using the stellar radii and uncertainties in the known-target manifest. The archive supplies those inputs for Kepler-10, Kepler-8 and WASP-18. For Kepler-17 the fit uses the discovery paper's stellar radius of 1.02 ± 0.03 solar radii. [Known-target comparison](../results/research/known_comparison.csv) · [NASA Archive composite parameters](https://exoplanetarchive.ipac.caltech.edu/docs/API_PSCompPars_columns.html) · [Kepler-17 discovery paper](https://doi.org/10.1088/0067-0049/197/1/14)

The Kepler BLS depths are 15.6–20.2% below the archive depths. To check for attenuation, the fitted transit was divided out of the original PDCSAP flux. The original and divided light curves were passed through the same detrending mask, and their difference was measured with the original box and cadence weights before and after detrending. The change in that differential depth was below 0.1% for each host. Applying the same box to the exposure-integrated transit model made the depth 14.0–16.7% shallower than the model's central depth and reproduced the measured Kepler BLS depths to within 0.5%. The box measurement therefore accounts for most of the gap on these light curves. The remaining radius-ratio differences for Kepler-17, Kepler-8 and WASP-18 exceed the quoted intervals. The experiment isolates detrending response for the fitted profiles, while the radius-ratio comparisons remain sensitive to model assumptions and the source of the archive parameters. [Depth diagnostic](../results/research/depth_diagnostic.csv)

All four pass the revised vetting, with duration assessed using densities of 1.068, 1.421655, 0.65374 and 0.795 g/cm³ in table order. Their depth-consistency p-values are approximately 0.985, 0.247, 0.983 and 0.836. WASP-18's event-weighted odd/even difference is 1.30σ. Its nominal secondary depth in the flattened light curve is 224 ± 10 ppm (about 21.4σ). On normalized PDCSAP flux, a local quadratic baseline on either side of phase 0.5 gave 333 ± 23 ppm across 45 eclipses. This agrees within the combined uncertainties with the published TESS eclipse of 357 ± 14 ppm. The known-target vetter uses 399 ppm, the published depth plus three quoted standard deviations, as the independent upper bound. [WASP-18 TESS eclipse measurement](https://www.nature.com/articles/s41586-023-06230-1)

### 3.2 Injection recovery on low-variability Kepler photometry

The follow-up used quarter 1 PDCSAP flux from KIC 6106415 and KIC 6603624. Visual inspection found low-amplitude variability and occasional gaps. The normalized robust scatter was 29.3 and 34.9 ppm respectively, falling to 25.1 and 30.4 ppm after detrending. Asteroseismic mass and radius inputs were 1.078 solar masses and 1.213 solar radii for KIC 6106415 and 1.09 solar masses and 1.181 solar radii for KIC 6603624. The strongest uninjected BLS boxes had depths of 39.8 and 17.8 ppm. Rescaling their formal SNRs by the ratio of quoted cadence noise to measured robust scatter gave 5.90 and 5.13. That rescaling was used to inspect the baselines. Recovery used the pipeline's ordinary detection criteria. [KIC 6106415 stellar modeling](https://academic.oup.com/mnras/article/479/4/4416/5043232) · [KIC 6603624 stellar modeling](https://doi.org/10.1088/0067-0049/214/2/27)

For each star, 20 independently drawn epochs were tested in each of nine period–radius bins. Periods were 1, 3 and 7 days, radii were 2, 4 and 8 Earth radii, and the impact parameter was fixed at 0.3. Signals were injected into the original PDCSAP flux before cleaning and both BLS passes. All 360 injections met the period, epoch, SNR and observed-event criteria. The Wilson 95% interval for each star/bin is 0.839–1.000. These estimates apply to the two observed quarters and the specified grid. The raw time-series inspection, inputs, trials and recovery tables are saved in [quiet_injection](../results/research/quiet_injection).

### 3.3 Blind search and crossmatch

The frozen search processed 22 of the 24 listed Kepler targets before the KOI crossmatch. KIC 2694497 and KIC 3112389 had no public Kepler long-cadence light curves in either the MAST search or the public Kepler archive listing. Nineteen of the 22 primary peaks met the nominal SNR and observed-event thresholds. Nine of those agreed within 1% of a published KOI period, comprising five confirmed objects and four catalog false positives. The other ten are unvetted signals, several with best periods near the upper edge of the search. The output includes the [crossmatch](../results/research/blind_crossmatch.csv), [unvetted candidates](../results/research/blind_candidates.csv), and [archive misses](../results/research/blind/failures.json). [MAST Kepler data access](https://archive.stsci.edu/kepler/search_retrieve.html) · [NASA cumulative KOI table](https://exoplanetarchive.ipac.caltech.edu/docs/PurposeOfKOITable.html)

### 3.4 Labeled vetting evaluation

The sample used a fixed-seed draw of ten confirmed planets and ten catalog false positives on single-KOI hosts with published periods from 0.55 to 12 days. The first two available Kepler long-cadence segments were searched without using the catalog period or disposition during detection and vetting. The final evaluation requires the recovered primary period to agree with the labeled period within 1%. Results with aliases or other periods contribute to analysis coverage but carry no label for scoring the vetting of the catalog signal. Stellar density and target-specific occultation estimates were omitted for this sample.

| Catalog class | Requested | Analyzed | Period matched | Correctly classified among matches |
|---|---:|---:|---:|---:|
| Confirmed | 10 | 10 | 5 | 5 retained (100%) |
| False positive | 10 | 10 | 6 | 4 rejected (66.7%) |

Overall period-match coverage was 11/20 (55%). Two matched catalog false positives passed the screening flags. [Label sample](../results/research/labels.csv) · [Evaluation metrics](../results/research/vetting_metrics.json)

## 4. Discussion and remaining scientific validation

The sample sizes, search coverage and adopted transit model limit the scientific interpretation. In particular, performance on the labeled set is conditional on primary-period recovery. The injection study tests two quiet Kepler baselines at one impact parameter. Establishing survey-level sensitivity and false-positive rates would require broader sampling of light curves, stellar properties, transit geometries and false-positive populations.

The main methodological risks are detrending loss, aliases, inaccurate formal errors, fixed limb darkening, correlated noise, dilution, and selection bias. The current outputs expose several of these through diagnostic plots, alias scores, incomplete-test indicators, posterior diagnostics, and experiment counts. Further work could add centroid analysis, joint stellar-density constraints, eccentric models, transit timing variation, or TLS. None should be substituted for careful inspection of the present validation evidence.

## References

Kovács, Zucker & Mazeh (2002), BLS; Kreidberg (2015), batman; Foreman-Mackey et al. (2013), emcee. Software/API references: [Lightkurve flatten](https://lightkurve.github.io/lightkurve/reference/api/lightkurve.LightCurve.flatten.html), [batman](https://lkreidberg.github.io/batman/), [NASA archive TAP](https://exoplanetarchive.ipac.caltech.edu/docs/API_resources.html).
