"""NASA Exoplanet Archive TAP queries, with saved query provenance."""
import json
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive


def catalog(table='pscomppars', output='data/catalog.csv'):
    columns = {
        'pscomppars': 'pl_name,hostname,pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_rade,pl_radeerr1,pl_radeerr2,pl_ratror,pl_ratrorerr1,pl_ratrorerr2,pl_trandep,pl_trandur,st_rad,st_raderr1,st_raderr2,st_dens',
        'cumulative': 'kepid,kepoi_name,koi_disposition,koi_period,koi_depth,koi_duration,koi_prad',
        'toi': 'tid,toi,toipfx,tfopwg_disp,pl_orbper,pl_trandep,pl_trandur'}
    if table not in columns:
        raise ValueError('Unsupported catalog')
    result = NasaExoplanetArchive.query_criteria(table=table, select=columns[table]).to_pandas()
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(path, index=False)
    path.with_suffix('.json').write_text(json.dumps(dict(table=table, columns=columns[table],
        retrieved_utc=datetime.now(timezone.utc).isoformat(), service='https://exoplanetarchive.ipac.caltech.edu/TAP'), indent=2))
    return result


def compare_known(results, planets):
    rows = []
    for r in results:
        matches = planets[planets.hostname.str.casefold() == r['target'].casefold()]
        if len(matches):
            # Match the primary recovered signal to the nearest catalog period;
            # do not misreport every other planet around this host as this detection.
            matches = matches.loc[[abs(matches.pl_orbper/r['period']-1).idxmin()]]
        for _, known in matches.iterrows():
            fit = r.get('fit', {})
            ratio = fit.get('radius_ratio', {})
            radius = fit.get('radius_earth', {})
            rows.append(dict(target=r['target'], planet=known.pl_name, recovered_period=r['period'],
                published_period=known.pl_orbper,
                period_error_percent=100*abs(r['period']/known.pl_orbper-1),
                recovered_depth_ppm=r['depth']*1e6,
                published_depth_ppm=known.pl_trandep*1e4,
                recovered_duration_hours=r['duration']*24,
                published_duration_hours=known.pl_trandur,
                fitted_radius_ratio=ratio.get('median'),
                fitted_radius_ratio_minus=ratio.get('minus'), fitted_radius_ratio_plus=ratio.get('plus'),
                published_radius_ratio=known.get('pl_ratror'),
                published_radius_ratio_minus=abs(known.get('pl_ratrorerr2', float('nan'))),
                published_radius_ratio_plus=known.get('pl_ratrorerr1'),
                fitted_radius_earth=radius.get('median'),
                fitted_radius_earth_minus=radius.get('minus'), fitted_radius_earth_plus=radius.get('plus'),
                published_radius_earth=known.get('pl_rade'),
                published_radius_earth_minus=abs(known.get('pl_radeerr2', float('nan'))),
                published_radius_earth_plus=known.get('pl_radeerr1'),
                fit_converged=fit.get('diagnostics', {}).get('converged')))
    return pd.DataFrame(rows)


def crossmatch_blind(results, koi, toi):
    """Match host identifiers after freezing detections; unmatched means uncataloged signal only."""
    rows = []
    for r in results:
        target = r['target'].upper().strip()
        matches = pd.DataFrame()
        if target.startswith('KIC '):
            matches = koi[pd.to_numeric(koi.kepid, errors='coerce') == int(target.split()[-1])]
            name, period, disposition = 'kepoi_name', 'koi_period', 'koi_disposition'
        elif target.startswith('TIC '):
            matches = toi[pd.to_numeric(toi.tid, errors='coerce') == int(target.split()[-1])]
            name, period, disposition = 'toi', 'pl_orbper', 'tfopwg_disp'
        if len(matches):
            for _, m in matches.iterrows():
                rows.append(dict(target=target, recovered_period=r['period'], catalog_id=str(m[name]),
                    catalog_period=m[period], catalog_disposition=m.get(disposition),
                    period_match=bool(abs(r['period']/m[period]-1)<.01)))
        else:
            rows.append(dict(target=target, recovered_period=r['period'], catalog_id=None,
                             catalog_period=None, catalog_disposition=None, period_match=False))
    return pd.DataFrame(rows)


def rank_unmatched_candidates(results, crossmatches):
    matched = set(crossmatches.loc[crossmatches.period_match, 'target']) if len(crossmatches) else set()
    rows = []
    for r in results:
        if (r['target'].strip().upper() not in matched and r['depth_snr'] >= 7
                and r['n_transits'] >= 3):
            rows.append(dict(target=r['target'], period_days=r['period'], box_depth_ppm=r['depth'] * 1e6,
                             duration_hours=r['duration'] * 24, depth_snr=r['depth_snr'],
                             observed_events=r['n_transits'], status='UNVETTED'))
    columns = ['target', 'period_days', 'box_depth_ppm', 'duration_hours', 'depth_snr', 'observed_events', 'status']
    return pd.DataFrame(rows, columns=columns).sort_values('depth_snr', ascending=False)
