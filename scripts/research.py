"""Execute archive experiments. Failures remain visible; never synthesize archive results."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.catalog import catalog, compare_known, crossmatch_blind, rank_unmatched_candidates
from src.experiments import batch, evaluate_labeled


def main():
    out = Path('results/research'); out.mkdir(parents=True, exist_ok=True)
    known = batch('config/known_targets.csv', output=out/'known', fit=True, max_segments=2, max_period=15.,
                  steps=20000, burn=4000, differential_evolution=True)
    # Select targets independently of dispositions, freeze detections, then query catalogs.
    blind = batch('config/blind_targets.csv', output=out/'blind', fit=False, max_segments=2, max_period=15.)
    planets = catalog('pscomppars', 'data/planets.csv')
    koi = catalog('cumulative', 'data/koi.csv')
    toi = catalog('toi', 'data/toi.csv')
    comparison = compare_known(known, planets)
    comparison.to_csv(out/'known_comparison.csv', index=False)
    physical_columns = {
        'target': 'target', 'planet': 'planet', 'fitted_radius_earth': 'radius_earth',
        'fitted_radius_earth_minus': 'radius_minus', 'fitted_radius_earth_plus': 'radius_plus',
        'published_radius_earth': 'published_radius', 'published_radius_earth_plus': 'published_radius_plus',
        'published_radius_earth_minus': 'published_radius_minus', 'fit_converged': 'converged',
    }
    physical = comparison.reindex(columns=physical_columns).dropna(subset=['fitted_radius_earth'])
    physical.rename(columns=physical_columns).to_csv(out/'physical_comparison.csv', index=False)
    crossmatches = crossmatch_blind(blind, koi, toi)
    crossmatches.to_csv(out/'blind_crossmatch.csv', index=False)
    rank_unmatched_candidates(blind, crossmatches).to_csv(out/'blind_candidates.csv', index=False)
    # Select one cataloged signal per host within the search window.
    labels = koi.groupby('kepid').filter(lambda group: len(group) == 1)
    labels = labels[labels.koi_disposition.isin(['CONFIRMED', 'FALSE POSITIVE'])
                    & labels.koi_period.between(.55, 12)]
    labels = labels.sort_values('kepid').groupby('koi_disposition').sample(n=10, random_state=42)
    labels['target'] = labels.kepid.map(lambda k: f'KIC {int(k)}')
    labels['mission'] = 'Kepler'
    labels[['target', 'mission']].to_csv(out/'labeled_targets.csv', index=False)
    labels = labels.rename(columns={'koi_disposition': 'disposition', 'koi_period': 'catalog_period'})
    labels[['target', 'disposition', 'catalog_period', 'kepoi_name']].to_csv(out/'labels.csv', index=False)
    batch(out/'labeled_targets.csv', output=out/'labeled', fit=False, max_segments=2, max_period=15., multiple=False)
    evaluate_labeled(out/'labeled/frozen_results.json', out/'labels.csv', out/'vetting_metrics.json')


if __name__ == '__main__':
    main()
