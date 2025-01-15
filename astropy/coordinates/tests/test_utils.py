import pytest

from astropy.coordinates.builtin_frames.utils import (
    get_offset_sun_from_barycenter,
    get_polar_motion,
)
from astropy.coordinates.solar_system import get_body_barycentric_posvel
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time
from astropy.utils.exceptions import AstropyWarning


def test_polar_motion_unsupported_dates():
    msg = r"Tried to get polar motions for times {} IERS.*"

    with pytest.warns(AstropyWarning, match=msg.format("before")):
        get_polar_motion(Time("1900-01-01"))

    with pytest.warns(AstropyWarning, match=msg.format("after")):
        get_polar_motion(Time("2100-01-01"))


def test_sun_from_barycenter_offset():
    time = Time("2020-01-01")
    pos, vel = get_body_barycentric_posvel("sun", time)

    offset = get_offset_sun_from_barycenter(time)
    assert_quantity_allclose(offset.xyz, pos.xyz)
    assert not bool(offset.differentials)

    offset_with_vel = get_offset_sun_from_barycenter(time, include_velocity=True)
    assert_quantity_allclose(offset_with_vel.xyz, pos.xyz)
    assert_quantity_allclose(offset_with_vel.differentials["s"].d_xyz, vel.xyz)

    reverse = get_offset_sun_from_barycenter(time, reverse=True)
    assert_quantity_allclose(reverse.xyz, -pos.xyz)
    assert not bool(reverse.differentials)

    reverse_with_vel = get_offset_sun_from_barycenter(
        time, reverse=True, include_velocity=True
    )
    assert_quantity_allclose(reverse_with_vel.xyz, -pos.xyz)
    assert_quantity_allclose(reverse_with_vel.differentials["s"].d_xyz, -vel.xyz)
