import tempfile
import unittest
from pathlib import Path
import warnings
import numpy as np
from astropy.io import fits
from src.download import quarantine_incomplete


class CacheTests(unittest.TestCase):
    def test_valid_fits_retained(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'valid.fits'
            fits.PrimaryHDU(np.ones(2000)).writeto(path)
            self.assertIsNone(quarantine_incomplete(path))
            self.assertTrue(path.exists())

    def test_truncated_fits_preserved_for_redownload(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'partial.fits'
            fits.PrimaryHDU(np.ones(2000)).writeto(path)
            partial = path.read_bytes()[:4000]
            path.write_bytes(partial)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                backup = quarantine_incomplete(path)
            self.assertFalse(path.exists())
            self.assertEqual(backup.read_bytes(), partial)
