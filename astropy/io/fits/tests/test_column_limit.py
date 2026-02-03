import numpy as np
import pytest

from astropy.table import Table


def test_fits_table_column_limit(tmp_path):
    ncols = 1000
    data = {f"c{i}": np.arange(5) for i in range(ncols)}
    t = Table(data)

    with pytest.raises(ValueError, match="999 columns"):
        t.write(tmp_path / "x.fits", format="fits")
