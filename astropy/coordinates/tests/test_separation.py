# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests for `position_angle()`, `separation()` and `separation_3d()` methods. They
are implemented in `BaseCoordinateFrame`, but are also exposed by `SkyCoord`
instances, so they should be tested on both.
"""

from contextlib import nullcontext
from typing import NamedTuple

import pytest

from astropy import units as u
from astropy.coordinates import (
    FK5,
    GCRS,
    ICRS,
    Angle,
    BaseCoordinateFrame,
    Distance,
    Galactic,
    NonRotationTransformationError,
    NonRotationTransformationWarning,
    SkyCoord,
)
from astropy.tests.helper import assert_quantity_allclose


class SeparationExpectation(NamedTuple):
    """
    The coordinates the position angle and separations are relative to
    are different for different tests.
    """

    coord: BaseCoordinateFrame | SkyCoord
    pytest_id: str
    position_angle: u.Quantity
    separation: u.Quantity
    separation_3d: u.Quantity
    reversed_position_angle: u.Quantity

    @property
    def reversed_separation(self) -> u.Quantity:
        return self.separation

    @property
    def reversed_separation_3d(self) -> u.Quantity:
        return self.separation_3d


@pytest.mark.parametrize("coord_class", [SkyCoord, ICRS])
@pytest.mark.parametrize(
    "other_coord",
    [
        SeparationExpectation(
            ICRS(0 * u.deg, 0 * u.deg, 3 * u.pc),
            "no_separation",
            0 * u.deg,
            0 * u.deg,
            0 * u.pc,
            0 * u.deg,
        ),
        SeparationExpectation(
            ICRS(0 * u.deg, 1 * u.deg, 3 * u.pc),
            "dec_separation",
            0 * u.deg,
            1 * u.deg,
            0.05235921 * u.pc,
            180 * u.deg,
        ),
        SeparationExpectation(
            ICRS(1 * u.deg, 0 * u.deg, 3 * u.pc),
            "ra_separation",
            90 * u.deg,
            1 * u.deg,
            0.05235921 * u.pc,
            270 * u.deg,
        ),
        SeparationExpectation(
            ICRS(0 * u.deg, 0 * u.deg, 10 * u.pc),
            "distance_separation",
            0 * u.deg,
            0 * u.deg,
            7 * u.pc,
            0 * u.deg,
        ),
        SeparationExpectation(
            SkyCoord(1 * u.deg, 1 * u.deg, 9 * u.pc),
            "SkyCoord_input",
            44.995636 * u.deg,
            1.4141777 * u.deg,
            6.00137 * u.pc,
            225.00436354 * u.deg,
        ),
    ],
    ids=lambda x: x.pytest_id,
)
@pytest.mark.parametrize("method", ["position_angle", "separation", "separation_3d"])
def test_scalar_coords(coord_class, other_coord, method):
    vernal_equinox = coord_class(0 * u.deg, 0 * u.deg, 3 * u.pc)
    assert_quantity_allclose(
        getattr(vernal_equinox, method)(other_coord.coord), getattr(other_coord, method)
    )
    assert_quantity_allclose(
        getattr(other_coord.coord, method)(vernal_equinox),
        getattr(other_coord, f"reversed_{method}"),
    )


@pytest.mark.parametrize("coord_class", [SkyCoord, ICRS])
@pytest.mark.parametrize(
    "other_coord",
    [
        SeparationExpectation(
            FK5(0 * u.deg, 90 * u.deg, 2 * u.pc),
            "FK5_input",
            245.42603114 * u.deg,
            6.07832112e-06 * u.deg,
            2.12173433e-07 * u.pc,
            65.42602474 * u.deg,
        ),
        SeparationExpectation(
            Galactic(0 * u.deg, 90 * u.deg, 9 * u.pc),
            "Galactic_input",
            347.14052211 * u.deg,
            62.871748 * u.deg,
            8.28158093 * u.pc,
            57.06807474 * u.deg,
        ),
    ],
    ids=lambda x: x.pytest_id,
)
@pytest.mark.parametrize("method", ["position_angle", "separation", "separation_3d"])
def test_scalar_coords_frame_transformation(coord_class, other_coord, method):
    north_pole = coord_class(0 * u.deg, 90 * u.deg, 2 * u.pc)
    assert_quantity_allclose(
        getattr(north_pole, method)(other_coord.coord), getattr(other_coord, method)
    )
    assert_quantity_allclose(
        getattr(other_coord.coord, method)(north_pole),
        getattr(other_coord, f"reversed_{method}"),
    )


@pytest.mark.parametrize("coord_class", [SkyCoord, ICRS])
@pytest.mark.parametrize(
    "other_coord",
    [
        SeparationExpectation(
            ICRS([-1, -2, -3] * u.deg, [0.1, 1.1, 2.1] * u.deg, [11, 13, 17] * u.pc),
            "ICRS_input",
            [275.710887, 272.880921, 271.963607] * u.deg,
            [1.0049871, 2.0021627, 2.9997464] * u.deg,
            [8.000635, 8.004959, 10.016293] * u.pc,
            [95.71001423, 92.84426777, 91.8562681] * u.deg,
        ),
        SeparationExpectation(
            SkyCoord([1, -2, 3] * u.deg, [-0.1, 1.1, 2.1] * u.deg, [11, 13, 17] * u.pc),
            "SkyCoord_input",
            [95.71088692, 272.880921, 88.03639251] * u.deg,
            [1.0049871, 2.0021627, 2.9997464] * u.deg,
            [8.000635, 8.004959, 10.016293] * u.pc,
            [275.71001423, 92.84426777, 268.1437319] * u.deg,
        ),
    ],
    ids=lambda x: x.pytest_id,
)
@pytest.mark.parametrize("method", ["position_angle", "separation", "separation_3d"])
def test_array_coords(coord_class, other_coord, method):
    coord = coord_class(0 * u.deg, [0, 1, 2] * u.deg, [3, 5, 7] * u.pc)
    assert_quantity_allclose(
        getattr(coord, method)(other_coord.coord), getattr(other_coord, method)
    )
    assert_quantity_allclose(
        getattr(other_coord.coord, method)(coord),
        getattr(other_coord, f"reversed_{method}"),
    )


@pytest.mark.parametrize(
    "coord",
    [
        pytest.param(FK5(1 * u.deg, 0 * u.deg, 1 * u.pc), id="FK5"),
        pytest.param(
            SkyCoord(1 * u.deg, 0 * u.deg, 1 * u.pc, frame="fk5"), id="SkyCoord"
        ),
    ],
)
@pytest.mark.parametrize(
    "other_coord",
    [
        SeparationExpectation(
            FK5(1 * u.deg, 0 * u.deg, 1 * u.pc, equinox="B1950"),
            "FK5_B1950",
            66.51310007 * u.deg,
            0.69835342 * u.deg,
            0.01218849 * u.pc,
            246.50823798 * u.deg,
        ),
        SeparationExpectation(
            SkyCoord(1 * u.deg, 0 * u.deg, 10 * u.pc, frame="fk5", equinox="B1950"),
            "SkyCoord_B1950",
            66.51310007 * u.deg,
            0.69835342 * u.deg,
            9.000083 * u.pc,
            246.50823798 * u.deg,
        ),
    ],
    ids=lambda x: x.pytest_id,
)
@pytest.mark.parametrize("method", ["position_angle", "separation", "separation_3d"])
def test_equinox_conversion(coord, other_coord, method):
    """
    Regression test for
        - #868, #891: methods raised errors in case of frames with equinoxes
        - #3106, #15659: unspecified equinox should be interpreted as frame default,
          both in `SkyCoord` methods (#3106) and in `BaseCoordinateFrame` methods with
          `SkyCoord` input (#15659).
        - #5702: `position_angle()` didn't check if equinoxes differed
    """
    assert_quantity_allclose(
        getattr(coord, method)(other_coord.coord), getattr(other_coord, method)
    )
    assert_quantity_allclose(
        getattr(other_coord.coord, method)(coord),
        getattr(other_coord, f"reversed_{method}"),
    )


@pytest.mark.parametrize("other_class", [SkyCoord, ICRS])
@pytest.mark.parametrize("coord_class", [SkyCoord, ICRS])
def test_separation_3d_dimensionless_distance(coord_class, other_class):
    assert_quantity_allclose(
        coord_class(35 * u.deg, 0 * u.deg, 3 * u.one).separation_3d(
            other_class(125 * u.deg, 0 * u.deg, 4 * u.one)
        ),
        5 * u.one,
    )


@pytest.mark.parametrize("dimensionless_class", [SkyCoord, ICRS])
@pytest.mark.parametrize("length_class", [SkyCoord, ICRS])
def test_separation_3d_distance_dimension_mismatch(length_class, dimensionless_class):
    dimensionless_coord = dimensionless_class(1 * u.deg, -2 * u.deg, 14)
    length_coord = length_class(-1 * u.deg, 2 * u.deg, 21 * u.pc)
    error_message = (
        "^Can only apply 'subtract' function to quantities with compatible dimensions$"
    )
    with pytest.raises(u.UnitConversionError, match=error_message):
        dimensionless_coord.separation_3d(length_coord)
    with pytest.raises(u.UnitConversionError, match=error_message):
        length_coord.separation_3d(dimensionless_coord)


@pytest.mark.parametrize("coord_class", [SkyCoord, ICRS])
def test_separation_3d_no_distance(coord_class):
    coord_no_distance = coord_class(0 * u.deg, 0 * u.deg)
    coord_with_distance = ICRS(0 * u.deg, 0 * u.deg, distance=3 * u.pc)
    with pytest.raises(
        ValueError,
        match="^This object does not have a distance; cannot compute 3d separation.$",
    ):
        coord_no_distance.separation_3d(coord_with_distance)
    with pytest.raises(
        u.UnitConversionError,
        match=(
            "^Can only apply 'subtract' function to quantities with compatible"
            " dimensions$"
        ),
    ):
        coord_with_distance.separation_3d(coord_no_distance)


@pytest.mark.parametrize("coord_class", [SkyCoord, ICRS])
@pytest.mark.parametrize(
    "velocity_kwargs",
    [
        pytest.param({}, id="no_velocity"),
        pytest.param({"radial_velocity": -108 * u.km / u.s}, id="radial_velocity"),
        pytest.param(
            {"pm_ra_cosdec": 7 * u.mas / u.s, "pm_dec": -5 * u.mas / u.s},
            id="proper_motion",
        ),
        pytest.param(
            {
                "radial_velocity": -108 * u.km / u.s,
                "pm_ra_cosdec": 7 * u.mas / u.s,
                "pm_dec": -5 * u.mas / u.s,
            },
            id="3d_velocity",
        ),
    ],
)
@pytest.mark.parametrize(
    "other_coord",
    [
        SeparationExpectation(
            ICRS(
                ra=20 * u.deg,
                dec=10 * u.deg,
                distance=8 * u.pc,
                pm_ra_cosdec=14 * u.mas / u.yr,
                pm_dec=-11 * u.mas / u.yr,
                radial_velocity=5 * u.km / u.s,
            ),
            "ICRS_with_velocity",
            134.58168775 * u.deg,
            13.89233851 * u.deg,
            3.36750833 * u.pc,
            317.18593154 * u.deg,
        ),
        SeparationExpectation(
            SkyCoord(ra=20 * u.deg, dec=10 * u.deg, distance=8 * u.pc),
            "SkyCoord_no_velocity",
            134.58168775 * u.deg,
            13.89233851 * u.deg,
            3.36750833 * u.pc,
            317.18593154 * u.deg,
        ),
    ],
    ids=lambda x: x.pytest_id,
)
@pytest.mark.parametrize("method", ["position_angle", "separation", "separation_3d"])
def test_with_velocities(coord_class, velocity_kwargs, other_coord, method):
    coord = coord_class(
        ra=10 * u.deg, dec=20 * u.deg, distance=5 * u.pc, **velocity_kwargs
    )
    assert_quantity_allclose(
        getattr(coord, method)(other_coord.coord), getattr(other_coord, method)
    )
    assert_quantity_allclose(
        getattr(other_coord.coord, method)(coord),
        getattr(other_coord, f"reversed_{method}"),
    )


@pytest.mark.parametrize("coord_class", [SkyCoord, ICRS])
@pytest.mark.parametrize(
    "method,output_type",
    [
        pytest.param(method, type_, id=method)
        for method, type_ in (
            ("position_angle", Angle),
            ("separation", Angle),
            ("separation_3d", Distance),
        )
    ],
)
def test_return_types(coord_class, method, output_type):
    """
    This is especially important for SkyCoord because SkyCoord instances
    expose the methods of their underlying frame at runtime, so they cannot be
    checked statically.
    """
    coord = coord_class(0 * u.deg, 0 * u.deg, 1 * u.pc)
    assert type(getattr(coord, method)(coord)) is output_type


@pytest.mark.parametrize("coord_class", [SkyCoord, ICRS])
@pytest.mark.parametrize(
    "origin_mismatch_kwarg,expectation",
    [
        pytest.param({"origin_mismatch": "ignore"}, nullcontext(), id="ignore"),
        pytest.param(
            {"origin_mismatch": "warn"},
            pytest.warns(
                NonRotationTransformationWarning,
                match="^transforming other coordinates from <GCRS Frame ",
            ),
            id="warn",
        ),
        pytest.param(
            {"origin_mismatch": "error"},
            pytest.raises(
                NonRotationTransformationError,
                match="^refusing to transform other coordinates from <GCRS Frame ",
            ),
            id="error",
        ),
        pytest.param(
            {},
            pytest.warns(
                NonRotationTransformationWarning,
                match="^transforming other coordinates from <GCRS Frame ",
            ),
            id="default",
        ),
        pytest.param(
            {"origin_mismatch": "bad"},
            pytest.raises(
                ValueError,
                match=(
                    r"^origin_mismatch='bad' is invalid\. Allowed values are 'ignore', "
                    r"'warn' or 'error'\.$"
                ),
            ),
            id="invalid",
        ),
    ],
)
def test_separation_origin_mismatch_action(
    coord_class, origin_mismatch_kwarg, expectation
):
    with expectation:
        coord_class(0 * u.deg, 0 * u.deg).separation(
            SkyCoord(0 * u.deg, 0 * u.deg, frame=GCRS), **origin_mismatch_kwarg
        )
