from pathlib import Path

import pytest
from numpy.testing import assert_allclose

from astropy.io import fits


@pytest.fixture
def canonical_data_base_path():
    return Path(__file__).parent / "data"


@pytest.fixture(
    params=(Path(__file__).parent / "data").glob("m13_*.fits"), ids=lambda x: x.name
)
def canonical_int_hdus(request):
    """
    This fixture provides 4 files downloaded from https://fits.gsfc.nasa.gov/registry/tilecompression.html

    Which are used as canonical tests of data not compressed by Astropy.
    """
    with fits.open(request.param) as hdul:
        yield hdul[1]


@pytest.fixture
def original_int_hdu(canonical_data_base_path):
    with fits.open(canonical_data_base_path / "m13.fits") as hdul:
        yield hdul[0]


# pytest-openfiles does not correctly check for open files when the files are
# opened in a fixture, so we skip the check here.
# https://github.com/astropy/pytest-openfiles/issues/32
@pytest.mark.openfiles_ignore
def test_canonical_data(original_int_hdu, canonical_int_hdus):
    assert_allclose(original_int_hdu.data, canonical_int_hdus.data)
