import pytest
import numpy as np
from io import BytesIO
from astropy.table import Table
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning, VerifyError

def test_fits_999_limit_smart_write():
    """Test that writing >999 columns is BLOCKED with a clear error."""
    data = {f'col{i}': [0] for i in range(1000)}
    t = Table(data) 
    
    # Ignore character length warnings during initialization
    with pytest.warns(VerifyWarning, match="is greater than 8 characters"):
        hdu = fits.BinTableHDU(t)
    
    outfile = BytesIO()
    
    # Writing MUST fail with our explicit message from _verify
    with pytest.raises((ValueError, VerifyError), match="more than 999 columns"):
        hdu.writeto(outfile)

def test_fits_large_table_allowed_in_memory():
    """
    Verify that tables > 999 columns can be created and manipulated in memory.
    Ensures no regressions for non-FITS in-memory workflows.
    """
    num_cols = 1200
    data = {f'col{i}': [0] for i in range(num_cols)}
    t = Table(data)
    
    # Table object should be unrestricted
    assert len(t.columns) == num_cols
    
    # HDU creation allowed (with warnings caught)
    with pytest.warns(VerifyWarning):
        hdu = fits.BinTableHDU(t)
    
    # In-memory operations (like adding a column) remain valid
    hdu.columns.add_col(fits.Column(name='EXTRA', format='J', array=[1]))
    
    # SMART FIX: Catch warnings during header update
    with pytest.warns(VerifyWarning, match="is greater than 8 characters"):
        hdu.update_header() 
    
    assert hdu.header['TFIELDS'] == num_cols + 1
    assert len(hdu.columns) == num_cols + 1

def test_fits_999_limit_boundary():
    """Verify exactly 999 columns still pass validation."""
    data = {f'col{i}': [0] for i in range(999)}
    t = Table(data)
    
    outfile = BytesIO()
    # This should pass without error
    t.write(outfile, format='fits')
    assert len(t.columns) == 999