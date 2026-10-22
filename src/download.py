"""PDCSAP archive access; downloads are cached by Lightkurve."""
from pathlib import Path
import time
import warnings
from astropy.io import fits
import lightkurve as lk


def quarantine_incomplete(path):
    """Preserve an unreadable/truncated cached FITS so MAST can download afresh."""
    path = Path(path)
    if not path.is_file():
        return None
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', fits.verify.VerifyWarning)
            with fits.open(path, memmap=False, lazy_load_hdus=False) as hdus:
                for hdu in hdus:
                    info = hdu.fileinfo()
                    if info and path.stat().st_size < info['datLoc'] + info['datSpan']:
                        raise OSError('Truncated FITS data')
                hdus.verify('exception')
    except (OSError, ValueError, fits.verify.VerifyError):
        backup = path.with_name(path.name + f'.incomplete-{time.time_ns()}')
        path.rename(backup)
        warnings.warn(f'Incomplete cached FITS preserved at {backup}; downloading again.', RuntimeWarning)
        return backup
    return None


def download_light_curve(target, mission='Kepler', cache_dir='data', max_segments=4):
    if max_segments < 1:
        raise ValueError('max_segments must be positive')
    author = 'Kepler' if mission == 'Kepler' else 'SPOC'
    exposure = 1800 if mission == 'Kepler' else 120
    result = lk.search_lightcurve(target, mission=mission, author=author, exptime=exposure)
    if len(result) == 0:
        raise ValueError(f'No {mission} {author} light curves at {exposure}s for {target}')
    Path(cache_dir).mkdir(parents=True, exist_ok=True)
    selected = result[:max_segments]
    cache_root = Path(cache_dir).resolve()
    for row in selected.table:
        path = cache_root / 'mastDownload' / str(row['obs_collection']) / str(row['obs_id']) / str(row['productFilename'])
        if not path.resolve().is_relative_to(cache_root):
            raise ValueError('Archive product path escapes the cache directory')
        quarantine_incomplete(path)
    curves = selected.download_all(
        flux_column='pdcsap_flux', quality_bitmask='default', download_dir=cache_dir)
    if curves is None or len(curves) == 0:
        raise RuntimeError(f'Archive download failed for {target}')
    return list(curves)
