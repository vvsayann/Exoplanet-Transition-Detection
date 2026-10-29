"""NASA Exoplanet Archive TAP queries, with saved query provenance."""
import json
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive


def catalog(table='pscomppars', output='data/catalog.csv'):
    columns = {
        'pscomppars': 'pl_name,hostname,pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_rade,pl_radeerr1,pl_radeerr2,pl_trandep,pl_trandur,st_rad,st_raderr1,st_raderr2,st_dens',
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
            rows.append(dict(target=r['target'], planet=known.pl_name, recovered_period=r['period'],
                published_period=known.pl_orbper,
                period_error_percent=100*abs(r['period']/known.pl_orbper-1),
                recovered_depth_ppm=r['depth']*1e6,
                published_depth_ppm=known.pl_trandep*1e4,
                recovered_duration_hours=r['duration']*24,
                published_duration_hours=known.pl_trandur))
    return pd.DataFrame(rows)


def crossmatch_blind(results, koi, toi):
    """Match host identifiers after freezing detections; unmatched means uncataloged signal only."""
    rows = []
    for r in results:
        target = r['target'].upper().strip()
        matches = pd.DataFrame()
        if target.startswith('KIC '):
            matches = koi[pd.to_numeric(koi.kepid, errors='coerce') == int(target.split()[-1])]
            name, period = 'kepoi_name', 'koi_period'
        elif target.startswith('TIC '):
            matches = toi[pd.to_numeric(toi.tid, errors='coerce') == int(target.split()[-1])]
            name, period = 'toi', 'pl_orbper'
        if len(matches):
            for _, m in matches.iterrows():
                rows.append(dict(target=target, recovered_period=r['period'], catalog_id=str(m[name]),
                    catalog_period=m[period], period_match=bool(abs(r['period']/m[period]-1)<.01)))
        else:
            rows.append(dict(target=target, recovered_period=r['period'], catalog_id=None,
                             catalog_period=None, period_match=False))
    return pd.DataFrame(rows)
