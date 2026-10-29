"""Execute archive experiments. Failures remain visible; never synthesize archive results."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pandas as pd
import numpy as np
from src.catalog import catalog, compare_known, crossmatch_blind
from src.experiments import batch, evaluate_labeled, injection_recovery
from src.download import download_light_curve


def main():
    out = Path('results/research'); out.mkdir(parents=True, exist_ok=True)
    known = batch('config/known_targets.csv', output=out/'known', fit=True, max_segments=2, max_period=15.)
    # Select targets independently of dispositions, freeze detections, then query catalogs.
    blind = batch('config/blind_targets.csv', output=out/'blind', fit=False, max_segments=2, max_period=15.)
    planets = catalog('pscomppars', 'data/planets.csv')
    koi = catalog('cumulative', 'data/koi.csv')
    toi = catalog('toi', 'data/toi.csv')
    compare_known(known, planets).to_csv(out/'known_comparison.csv', index=False)
    # Convert saved radius-ratio draws with measured catalog stellar radii.
    physical = []
    rng = np.random.default_rng(42)
    for r in known:
        host = planets[planets.hostname.str.casefold() == r['target'].casefold()]
        if not len(host) or 'fit' not in r:
            continue
        m = host.loc[abs(host.pl_orbper/r['period']-1).idxmin()]
        samples = np.load(out/'known'/r['target']/'posterior.npy')
        error = np.nanmax([abs(m.st_raderr1), abs(m.st_raderr2)])
        if not np.isfinite(m.st_rad) or not np.isfinite(error):
            continue
        radius = samples[:,2]*rng.normal(m.st_rad,error,len(samples))*109.076
        lo, mid, hi = np.percentile(radius,[16,50,84])
        physical.append(dict(target=r['target'], planet=m.pl_name, radius_earth=mid,
            radius_minus=mid-lo, radius_plus=hi-mid, published_radius=m.pl_rade,
            published_radius_plus=m.pl_radeerr1, published_radius_minus=m.pl_radeerr2,
            converged=r['fit']['diagnostics']['converged']))
    pd.DataFrame(physical).to_csv(out/'physical_comparison.csv',index=False)
    crossmatch_blind(blind, koi, toi).to_csv(out/'blind_crossmatch.csv', index=False)
    # One object per host: avoid assigning conflicting candidate labels to a star.
    labels = koi[koi.koi_disposition.isin(['CONFIRMED','FALSE POSITIVE'])].copy()
    labels = labels.groupby('kepid').filter(lambda g: len(g)==1)
    labels = labels.groupby('koi_disposition').head(10)
    labels['target'] = labels.kepid.map(lambda k: f'KIC {int(k)}')
    labels['mission'] = 'Kepler'
    labels[['target','mission']].to_csv(out/'labeled_targets.csv', index=False)
    labels.rename(columns={'koi_disposition':'disposition'})[['target','disposition']].to_csv(out/'labels.csv', index=False)
    batch(out/'labeled_targets.csv', output=out/'labeled', fit=False, max_segments=2, max_period=15.)
    evaluate_labeled(out/'labeled/frozen_results.json', out/'labels.csv', out/'vetting_metrics.json')
    # A named baseline must be independently inspected for variability/known signals.
    # Do not automatically claim a quiet-star completeness measurement.
    print('Archive batches complete. Select and inspect quiet real baselines before injection runs.')


if __name__ == '__main__':
    main()
