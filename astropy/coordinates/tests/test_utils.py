from astropy.time import Time
from astropy.coordinates.builtin_frames.utils import get_polar_motion
from astropy.utils.exceptions import AstropyWarning
import pytest


def test_polar_motion_unsupported_dates():
    msg = r'Tried to get polar motions for times {} IERS.*'

    with pytest.warns(AstropyWarning, match=msg.format('before')):
        get_polar_motion(Time('1900-01-01'))

    with pytest.warns(AstropyWarning, match=msg.format('after')):
        get_polar_motion(Time('2100-01-01'))
