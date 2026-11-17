# Delivery and validation ledger

Updated 22 October 2026. This ledger separates implemented workflows from measured scientific outcomes.

| Plan item | Implementation | Validation |
|---|---|---|
| PDCSAP archive access | `src/download.py`; mission/author filtering and download cache | Kepler-10 archive download succeeded |
| Quality cleaning and detrending | `src/preprocess.py`; upward-only clipping, segment handling, protected transit mask | Downward-dip preservation test passed |
| BLS, folding, aliases | `src/search.py`, pipeline plots | Synthetic 3-day signal recovered at 2.99943475 days; real Kepler-10 smoke test at 0.83771248 days |
| Physical fit and uncertainties | Exposure-integrated batman, least squares, emcee, corner | Posterior generation tested; 300-step demo explicitly fails convergence check |
| Radius conversion | Stellar radius/error arguments; posterior propagation | Formula implemented; measured host-specific inputs must be supplied |
| Multiple known planets | Four-target manifest and batch driver | Kepler-10 completed; broader archive batch in progress; inspect before acceptance |
| Vetting | Odd/even, secondary, duration, shape, event consistency | Synthetic unequal-depth binary rejected; labeled evaluation workflow supplied |
| Multiple signals | Iterative masking and BLS | Available in pipeline; additional detections require inspection |
| Injection recovery | Pre-detrending injections, exact/alias separation, binomial intervals | 27 synthetic trials and 18 real-baseline trials completed; see report for limitations |
| Blind search | 24 identifiers; frozen detections before catalog matching | Manifest supplied; completed coverage and crossmatch must be assessed |
| Reproduction | CLI, six notebooks, dependency lock, figures and structured outputs | Four automated tests passed |

Optional TLS, machine-learning classification, TTV, centroid analysis, and Streamlit are not included. They are extensions in the supplied plan, not prerequisites for the core pipeline.

Acceptance of the full research project requires: completed multi-target archive comparisons; converged physical posteriors and literature agreement; a labeled sample containing both classes with quantified coverage; inspected quiet-star completeness across sufficient trials; and a completed blind search followed by catalog crossmatching. Software implementation and smoke tests alone do not demonstrate these success criteria. See generated output folders and the report for measured evidence.
