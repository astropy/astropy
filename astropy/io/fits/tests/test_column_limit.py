import numpy as np
import pytest
from astropy.io import fits


def test_fits_table_column_limit(tmp_path):
    # 1000 columns -> should hit the guard in _prewriteto
    cols = [
        fits.Column(name=f"c{i}", format="E", array=np.array([1.0]))
        for i in range(1000)
    ]

    hdu = fits.BinTableHDU.from_columns(cols)

    with pytest.raises(ValueError, match="999 columns"):
        hdu.writeto(tmp_path / "too_many_cols.fits")
