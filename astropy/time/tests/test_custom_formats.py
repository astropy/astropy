import pytest

from astropy.time import Time, TimeFormat


def test_custom_time_format_exception():
    class SpecificException(ValueError):
        pass

    # Note there is no way to remove this after the test!
    class Custom(TimeFormat):
        name = "custom_time_format_test"

        def set_jds(self, val, val2):
            raise SpecificException

    # with pytest.raises(SpecificException):
    #    Time(7., format="custom_time_format_test")
    try:
        Time(7.0, format="custom_time_format_test")
    except SpecificException:
        # This is not what astropy does but it is fine too
        pass
    except ValueError as e:
        assert (hasattr(e, "__cause__")
                and isinstance(e.__cause__, SpecificException))
