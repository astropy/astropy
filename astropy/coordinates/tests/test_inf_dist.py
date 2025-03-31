import numpy as np
import pytest
import astropy.units as u
from astropy.coordinates import SkyCoord, set_inf_dist, HCRS, ICRS

def test_icrs_to_hcrs_with_missing_distance():
    """Test transformation from ICRS to HCRS, handling missing distance."""

    icrs_coord = SkyCoord(10*u.deg, 20*u.deg, frame=ICRS)  # No distance provided

    with pytest.raises(u.UnitsError, match="ICRS to HCRS transformation requires a distance"):
        icrs_coord.hcrs  # This should raise an error due to missing distance

def test_icrs_to_hcrs_with_set_inf_dist():
    """Test ICRS to HCRS transformation with set_inf_dist enabled."""

    icrs_coord = SkyCoord(10*u.deg, 20*u.deg, frame=ICRS)  # No explicit distance

    with set_inf_dist(True):
        hcrs_coord = icrs_coord.transform_to(HCRS)

    assert hcrs_coord is not None, "Transformation failed with set_inf_dist enabled."
    assert hcrs_coord.distance.unit == u.kpc, "Distance unit is incorrect after transformation."

def test_icrs_to_hcrs_nan_warning():
    """Test for NaN values in ICRS to HCRS transformation."""

    icrs_coord = SkyCoord(0*u.deg, 0*u.deg, 1e12*u.kpc, frame=ICRS)  # Large assumed distance

    hcrs_coord = icrs_coord.transform_to(HCRS)

    assert not (np.isnan(hcrs_coord.ra) or np.isnan(hcrs_coord.dec)), "NaN detected in transformation!"
    print("Test Passed: No NaN detected in ICRS to HCRS transformation.")


def test_icrs_to_gcrs_nan_warning():
    """Test ICRS to GCRS transformation, checking for NaN warnings."""
    
    icrs_coord = SkyCoord(10*u.deg, 20*u.deg, frame='icrs')

    # Transform to GCRS using set_inf_dist()
    with set_inf_dist(True):
        gcrs_coord = icrs_coord.gcrs

    # **Trigger NaN Check**
    if np.isnan(gcrs_coord.ra) or np.isnan(gcrs_coord.dec):
        print("Test Failed: NaN detected in transformed coordinates!")
    else:
        print("Test Passed: No NaN detected in ICRS to GCRS transformation.")

test_icrs_to_gcrs_nan_warning()


# Enable infinite distance assumption
with set_inf_dist(True):
    coord = SkyCoord(0*u.deg, 0*u.deg, frame='icrs')  # No distance provided
    print("ICRS -> HCRS:", coord.hcrs)

