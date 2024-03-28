# Licensed under a 3-clause BSD style license - see LICENSE.rst

import inspect
import sys
from contextlib import nullcontext
from io import StringIO

import numpy as np
import pytest

from astropy import units as u
from astropy.cosmology import core, flrw
from astropy.cosmology.funcs import z_at_value
from astropy.cosmology.funcs.optimize import _z_at_scalar_value
from astropy.cosmology.realizations import (
    WMAP1,
    WMAP3,
    WMAP5,
    WMAP7,
    WMAP9,
    Planck13,
    Planck15,
    Planck18,
)
from astropy.tests.helper import PYTEST_LT_8_0
from astropy.units import allclose
from astropy.utils.compat.optional_deps import HAS_SCIPY
from astropy.utils.exceptions import AstropyUserWarning


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_z_at_value_scalar():
    # These are tests of expected values, and hence have less precision
    # than the roundtrip tests below (test_z_at_value_roundtrip);
    # here we have to worry about the cosmological calculations
    # giving slightly different values on different architectures,
    # there we are checking internal consistency on the same architecture
    # and so can be more demanding
    cosmo = Planck13
    assert allclose(z_at_value(cosmo.age, 2 * u.Gyr), 3.19812268, rtol=1e-6)
    assert allclose(z_at_value(cosmo.lookback_time, 7 * u.Gyr), 0.795198375, rtol=1e-6)
    assert allclose(z_at_value(cosmo.distmod, 46 * u.mag), 1.991389168, rtol=1e-6)
    assert allclose(
        z_at_value(cosmo.luminosity_distance, 1e4 * u.Mpc), 1.36857907, rtol=1e-6
    )
    assert allclose(
        z_at_value(cosmo.luminosity_distance, 26.037193804 * u.Gpc, ztol=1e-10),
        3,
        rtol=1e-9,
    )
    assert allclose(
        z_at_value(cosmo.angular_diameter_distance, 1500 * u.Mpc, zmax=2),
        0.681277696,
        rtol=1e-6,
    )
    assert allclose(
        z_at_value(cosmo.angular_diameter_distance, 1500 * u.Mpc, zmin=2.5),
        3.7914908,
        rtol=1e-6,
    )

    # test behavior when the solution is outside z limits (should
    # raise a CosmologyError)
    with (
        pytest.raises(core.CosmologyError),
        pytest.warns(AstropyUserWarning, match="fval is not bracketed"),
    ):
        z_at_value(cosmo.angular_diameter_distance, 1500 * u.Mpc, zmax=0.5)

    with (
        pytest.raises(core.CosmologyError),
        pytest.warns(AstropyUserWarning, match="fval is not bracketed"),
        np.errstate(over="ignore"),
    ):
        z_at_value(cosmo.angular_diameter_distance, 1500 * u.Mpc, zmin=4.0)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
class Test_ZatValue:
    def setup_class(self):
        self.cosmo = Planck13

    def test_broadcast_arguments(self):
        """Test broadcast of arguments."""
        # broadcasting main argument
        assert allclose(
            z_at_value(self.cosmo.age, [2, 7] * u.Gyr),
            [3.1981206134773115, 0.7562044333305182],
            rtol=1e-6,
        )

        # basic broadcast of secondary arguments
        assert allclose(
            z_at_value(
                self.cosmo.angular_diameter_distance,
                1500 * u.Mpc,
                zmin=[0, 2.5],
                zmax=[2, 4],
            ),
            [0.681277696, 3.7914908],
            rtol=1e-6,
        )

        # more interesting broadcast
        assert allclose(
            z_at_value(
                self.cosmo.angular_diameter_distance,
                1500 * u.Mpc,
                zmin=[[0, 2.5]],
                zmax=[2, 4],
            ),
            [[0.681277696, 3.7914908]],
            rtol=1e-6,
        )

    def test_broadcast_bracket(self):
        """`bracket` has special requirements."""
        # start with an easy one
        assert allclose(
            z_at_value(self.cosmo.age, 2 * u.Gyr, bracket=None),
            3.1981206134773115,
            rtol=1e-6,
        )

        # now actually have a bracket
        assert allclose(
            z_at_value(self.cosmo.age, 2 * u.Gyr, bracket=[0, 4]),
            3.1981206134773115,
            rtol=1e-6,
        )

        # now a bad length
        with pytest.raises(ValueError, match="sequence"):
            z_at_value(self.cosmo.age, 2 * u.Gyr, bracket=[0, 4, 4, 5])

        # now the wrong dtype : an ndarray, but not an object array
        with pytest.raises(TypeError, match="dtype"):
            z_at_value(self.cosmo.age, 2 * u.Gyr, bracket=np.array([0, 4]))

        # now an object array of brackets
        bracket = np.array([[0, 4], [0, 3, 4]], dtype=object)
        assert allclose(
            z_at_value(self.cosmo.age, 2 * u.Gyr, bracket=bracket),
            [3.1981206134773115, 3.1981206134773115],
            rtol=1e-6,
        )

    def test_bad_broadcast(self):
        """Shapes mismatch as expected"""
        with pytest.raises(ValueError, match="broadcast"):
            z_at_value(
                self.cosmo.angular_diameter_distance,
                1500 * u.Mpc,
                zmin=[0, 2.5, 0.1],
                zmax=[2, 4],
            )

    def test_scalar_input_to_output(self):
        """Test scalar input returns a scalar."""
        z = z_at_value(
            self.cosmo.angular_diameter_distance, 1500 * u.Mpc, zmin=0, zmax=2
        )
        assert isinstance(z, u.Quantity)
        assert z.dtype == np.float64
        assert z.shape == ()


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_z_at_value_numpyvectorize():
    """Test that numpy vectorize fails on Quantities.

    If this test starts failing then numpy vectorize can be used instead of
    the home-brewed vectorization. Please submit a PR making the change.
    """
    z_at_value = np.vectorize(
        _z_at_scalar_value, excluded=["func", "method", "verbose"]
    )
    with pytest.raises(u.UnitConversionError, match="dimensionless quantities"):
        z_at_value(Planck15.age, 10 * u.Gyr)


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
def test_z_at_value_verbose(monkeypatch):
    cosmo = Planck13

    # Test the "verbose" flag. Since this uses "print", need to mod stdout
    mock_stdout = StringIO()
    monkeypatch.setattr(sys, "stdout", mock_stdout)

    resx = z_at_value(cosmo.age, 2 * u.Gyr, verbose=True)
    assert str(resx.value) in mock_stdout.getvalue()  # test "verbose" prints res


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
@pytest.mark.parametrize("method", ["Brent", "Golden", "Bounded"])
def test_z_at_value_bracketed(method):
    """
    Test 2 solutions for angular diameter distance by not constraining zmin, zmax,
    but setting `bracket` on the appropriate side of the turning point z.
    Setting zmin / zmax should override `bracket`.
    """
    cosmo = Planck13

    if PYTEST_LT_8_0:
        ctx_fval = nullcontext()
    else:
        ctx_fval = pytest.warns(AstropyUserWarning, match="fval is not bracketed")

    if method == "Bounded":
        with pytest.warns(AstropyUserWarning, match="fval is not bracketed"):
            z = z_at_value(cosmo.angular_diameter_distance, 1500 * u.Mpc, method=method)
        if z > 1.6:
            z = 3.7914908
            bracket = (0.9, 1.5)
        else:
            z = 0.6812777
            bracket = (1.6, 2.0)
        with pytest.warns(UserWarning, match="Option 'bracket' is ignored"), ctx_fval:
            assert allclose(
                z_at_value(
                    cosmo.angular_diameter_distance,
                    1500 * u.Mpc,
                    method=method,
                    bracket=bracket,
                ),
                z,
                rtol=1e-6,
            )
    else:
        assert allclose(
            z_at_value(
                cosmo.angular_diameter_distance,
                1500 * u.Mpc,
                method=method,
                bracket=(0.3, 1.0),
            ),
            0.6812777,
            rtol=1e-6,
        )
        assert allclose(
            z_at_value(
                cosmo.angular_diameter_distance,
                1500 * u.Mpc,
                method=method,
                bracket=(2.0, 4.0),
            ),
            3.7914908,
            rtol=1e-6,
        )
        assert allclose(
            z_at_value(
                cosmo.angular_diameter_distance,
                1500 * u.Mpc,
                method=method,
                bracket=(0.1, 1.5),
            ),
            0.6812777,
            rtol=1e-6,
        )
        assert allclose(
            z_at_value(
                cosmo.angular_diameter_distance,
                1500 * u.Mpc,
                method=method,
                bracket=(0.1, 1.0, 2.0),
            ),
            0.6812777,
            rtol=1e-6,
        )
        with pytest.warns(AstropyUserWarning, match=r"fval is not bracketed"):
            assert allclose(
                z_at_value(
                    cosmo.angular_diameter_distance,
                    1500 * u.Mpc,
                    method=method,
                    bracket=(0.9, 1.5),
                ),
                0.6812777,
                rtol=1e-6,
            )
            assert allclose(
                z_at_value(
                    cosmo.angular_diameter_distance,
                    1500 * u.Mpc,
                    method=method,
                    bracket=(1.6, 2.0),
                ),
                3.7914908,
                rtol=1e-6,
            )
        assert allclose(
            z_at_value(
                cosmo.angular_diameter_distance,
                1500 * u.Mpc,
                method=method,
                bracket=(1.6, 2.0),
                zmax=1.6,
            ),
            0.6812777,
            rtol=1e-6,
        )
        assert allclose(
            z_at_value(
                cosmo.angular_diameter_distance,
                1500 * u.Mpc,
                method=method,
                bracket=(0.9, 1.5),
                zmin=1.5,
            ),
            3.7914908,
            rtol=1e-6,
        )

    if not PYTEST_LT_8_0 and method == "Bounded":
        ctx_bracket = pytest.warns(
            UserWarning, match="Option 'bracket' is ignored by method Bounded"
        )
    else:
        ctx_bracket = nullcontext()

    with (
        pytest.raises(core.CosmologyError),
        pytest.warns(AstropyUserWarning, match="fval is not bracketed"),
        ctx_bracket,
    ):
        z_at_value(
            cosmo.angular_diameter_distance,
            1500 * u.Mpc,
            method=method,
            bracket=(3.9, 5.0),
            zmin=4.0,
        )


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
@pytest.mark.parametrize("method", ["Brent", "Golden", "Bounded"])
def test_z_at_value_unconverged(method):
    """
    Test warnings on non-converged solution when setting `maxfun` to too small iteration number -
    only 'Bounded' returns status value and specific message.
    """
    cosmo = Planck18
    ztol = {"Brent": [1e-4, 1e-4], "Golden": [1e-3, 1e-2], "Bounded": [1e-3, 1e-1]}

    if method == "Bounded":
        ctx = pytest.warns(
            AstropyUserWarning,
            match="Solver returned 1: Maximum number of function calls reached",
        )
    else:
        ctx = pytest.warns(AstropyUserWarning, match="Solver returned None")

    with ctx:
        z0 = z_at_value(
            cosmo.angular_diameter_distance, 1 * u.Gpc, zmax=2, maxfun=13, method=method
        )
    with ctx:
        z1 = z_at_value(
            cosmo.angular_diameter_distance, 1 * u.Gpc, zmin=2, maxfun=13, method=method
        )

    assert allclose(z0, 0.32442, rtol=ztol[method][0])
    assert allclose(z1, 8.18551, rtol=ztol[method][1])


@pytest.mark.skipif(not HAS_SCIPY, reason="test requires scipy")
@pytest.mark.parametrize(
    "cosmo",
    [
        Planck13,
        Planck15,
        Planck18,
        WMAP1,
        WMAP3,
        WMAP5,
        WMAP7,
        WMAP9,
        flrw.LambdaCDM,
        flrw.FlatLambdaCDM,
        flrw.wpwaCDM,
        flrw.w0wzCDM,
        flrw.wCDM,
        flrw.FlatwCDM,
        flrw.w0waCDM,
        flrw.Flatw0waCDM,
    ],
)
def test_z_at_value_roundtrip(cosmo):
    """
    Calculate values from a known redshift, and then check that
    z_at_value returns the right answer.
    """
    z = 0.5

    # Skip Ok, w, de_density_scale because in the Planck cosmologies
    # they are redshift independent and hence uninvertable,
    # *_distance_z1z2 methods take multiple arguments, so require
    # special handling
    # clone is not a redshift-dependent method
    # nu_relative_density is not redshift-dependent in the WMAP cosmologies
    skip = (
        "Ok",
        "Otot",
        "angular_diameter_distance_z1z2",
        "clone",
        "is_equivalent",
        "de_density_scale",
        "w",
    )
    if str(cosmo.name).startswith("WMAP"):
        skip += ("nu_relative_density",)

    methods = inspect.getmembers(cosmo, predicate=inspect.ismethod)

    for name, func in methods:
        if name.startswith("_") or name in skip:
            continue
        fval = func(z)
        # we need zmax here to pick the right solution for
        # angular_diameter_distance and related methods.
        # Be slightly more generous with rtol than the default 1e-8
        # used in z_at_value
        got = z_at_value(func, fval, bracket=[0.3, 1.0], ztol=1e-12)
        assert allclose(got, z, rtol=2e-11), f"Round-trip testing {name} failed"

    # Test distance functions between two redshifts; only for realizations
    if isinstance(getattr(cosmo, "name", None), str):
        z2 = 2.0
        func_z1z2 = [
            lambda z1: cosmo._comoving_distance_z1z2(z1, z2),
            lambda z1: cosmo._comoving_transverse_distance_z1z2(z1, z2),
            lambda z1: cosmo.angular_diameter_distance_z1z2(z1, z2),
        ]
        for func in func_z1z2:
            fval = func(z)
            assert allclose(z, z_at_value(func, fval, zmax=1.5, ztol=1e-12), rtol=2e-11)
