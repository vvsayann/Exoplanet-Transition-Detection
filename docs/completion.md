# Delivery and validation ledger

Updated 17 November 2026. The [report](report.md) describes the measured results and their scope.

| Plan item | Implementation | Validation |
|---|---|---|
| PDCSAP archive access | `src/download.py` with mission/author filtering and a download cache | Processed three known Kepler hosts and WASP-18 using Kepler and TESS photometry |
| Quality cleaning and detrending | `src/preprocess.py`; upward-only clipping, segment handling, protected transit mask | Downward-dip preservation test passed |
| BLS, folding, aliases | `src/search.py` and pipeline plots | All four known periods agree with the archive within 0.03% |
| Physical fit and uncertainties | Exposure-integrated batman, least squares, emcee and corner plots | Kepler-10, Kepler-17 and Kepler-8 pass with 20,000 steps and 4,000 burn-in. The saved WASP-18 fit passes with 10,000 steps and 2,000 burn-in. |
| Radius and depth comparison | Catalog ratio and Earth-radius comparisons, fitted-profile depth diagnostic | Four hosts compared. The box explains most of the Kepler depth gap. The detrending response changes by less than 0.1% in the paired test. Ratios for Kepler-17, Kepler-8 and WASP-18 differ beyond their quoted intervals. |
| Vetting | Event-depth error floor, parity and secondary tests, duration with measured density | All four known planets retained. The WASP-18 test uses a published TESS occultation bound. |
| Labeled vetting | Period-aware matching and coverage by class | 20/20 analyzed, 11/20 period-matched. Retained 5/5 confirmed planets and rejected 4/6 catalog false positives among matches. |
| Multiple signals | Iterative masking and BLS | Available in pipeline; additional detections require inspection |
| Injection recovery | Injections before both search passes, exact recovery, alias tracking and Wilson intervals | 360/360 recovered across two inspected Kepler baselines and nine bins per star, with 20 trials per bin |
| Blind search | Frozen detections and post-search KOI crossmatch | 22/24 analyzed and crossmatched. The other two have no public Kepler long-cadence products. Ten unmatched signals meet the nominal detection criteria and are listed as unvetted. |
| Reproduction | CLI, notebooks, figures and structured outputs | Eight automated tests passed. See [test output](test_output.txt). |

The evidence covers the observed quarters, target sample and tested injection grid. The vetting evaluation measures classification conditional on recovering the catalog period, and the injection results use a fixed impact parameter on two stars. Section 4 of the report discusses the limits on broader scientific claims.

Review files include the [known-target comparison](../results/research/known_comparison.csv), [depth diagnostic](../results/research/depth_diagnostic.csv), [blind candidate list](../results/research/blind_candidates.csv), and [vetting metrics](../results/research/vetting_metrics.json).
