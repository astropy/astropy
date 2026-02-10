import pytest

from astropy.time import Time
from astropy.validation.errors import PhysicalInconsistencyError
from astropy.validation.time import validate_time_scale_consistency


def test_consistent_time_scale():
    t = Time(["2020-01-01", "2020-01-02"], scale="utc")
    assert validate_time_scale_consistency(t, strict=True)


def test_inconsistent_time_scale():
    t1 = Time("2020-01-01", scale="utc")
    t2 = Time("2020-01-02", scale="tt")

    with pytest.raises(PhysicalInconsistencyError):
        validate_time_scale_consistency([t1, t2], strict=True)
