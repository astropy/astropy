# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys
import io

import pytest
import numpy as np

from astropy.cosmology import core, funcs, realizations
from astropy.cosmology.realizations import Planck13, Planck15, Planck18, WMAP5, WMAP7, WMAP9
from astropy.units import allclose
from astropy import units as u
from astropy.utils.exceptions import AstropyUserWarning
from astropy.utils.compat.optional_deps import HAS_SCIPY  # noqa


@pytest.mark.skipif('not HAS_SCIPY')
def test_z_at_value():
    # These are tests of expected values, and hence have less precision
    # than the roundtrip tests below (test_z_at_value_roundtrip);
    # here we have to worry about the cosmological calculations
    # giving slightly different values on different architectures,
    # there we are checking internal consistency on the same architecture
    # and so can be more demanding
    z_at_value = funcs.z_at_value
    cosmo = Planck13
    assert allclose(z_at_value(cosmo.age, 2 * u.Gyr), 3.19812268, rtol=1e-6)
    assert allclose(z_at_value(cosmo.lookback_time, 7 * u.Gyr), 0.795198375, rtol=1e-6)
    assert allclose(z_at_value(cosmo.distmod, 46 * u.mag), 1.991389168, rtol=1e-6)
    assert allclose(z_at_value(cosmo.luminosity_distance, 1e4 * u.Mpc), 1.36857907, rtol=1e-6)
    assert allclose(z_at_value(cosmo.luminosity_distance, 26.037193804 * u.Gpc, ztol=1e-10),
                    3, rtol=1e-9)
    assert allclose(z_at_value(cosmo.angular_diameter_distance, 1500*u.Mpc, zmax=2),
                    0.681277696, rtol=1e-6)
    assert allclose(z_at_value(cosmo.angular_diameter_distance, 1500*u.Mpc, zmin=2.5),
                    3.7914908, rtol=1e-6)

    # test behavior when the solution is outside z limits (should
    # raise a CosmologyError)
    with pytest.raises(core.CosmologyError):
        with pytest.warns(AstropyUserWarning, match=r'fval is not bracketed'):
            z_at_value(cosmo.angular_diameter_distance, 1500*u.Mpc, zmax=0.5)

    with pytest.raises(core.CosmologyError):
        with pytest.warns(AstropyUserWarning, match=r'fval is not bracketed'):
            z_at_value(cosmo.angular_diameter_distance, 1500*u.Mpc, zmin=4.)


@pytest.mark.skipif('not HAS_SCIPY')
def test_z_at_value_verbose(monkeypatch):
    z_at_value = funcs.z_at_value
    cosmo = Planck13

    # Test the "verbose" flag. Since this uses "print", need to mod stdout
    mock_stdout = io.StringIO()
    monkeypatch.setattr(sys, 'stdout', mock_stdout)

    resx = z_at_value(cosmo.age, 2 * u.Gyr, verbose=True)
    assert str(resx) in mock_stdout.getvalue()  # test "verbose" prints res


@pytest.mark.skipif('not HAS_SCIPY')
@pytest.mark.parametrize('method', ['Brent', 'Golden', 'Bounded'])
def test_z_at_value_bracketed(method):
    """
    Test 2 solutions for angular diameter distance by not constraining zmin, zmax,
    but setting `bracket` on the appropriate side of the turning point z.
    Setting zmin / zmax should override `bracket`.
    """
    z_at_value = funcs.z_at_value
    cosmo = Planck13

    if method == 'Bounded':
        with pytest.warns(AstropyUserWarning, match=r'fval is not bracketed'):
            z = z_at_value(cosmo.angular_diameter_distance, 1500*u.Mpc, method=method)
        if z > 1.6:
            z = 3.7914908
            bracket = (0.9, 1.5)
        else:
            z = 0.6812777
            bracket = (1.6, 2.0)
        with pytest.warns(UserWarning, match=r"Option 'bracket' is ignored"):
            assert allclose(z_at_value(cosmo.angular_diameter_distance, 1500*u.Mpc, method=method,
                                       bracket=bracket), z, rtol=1e-6)
    else:
        assert allclose(z_at_value(cosmo.angular_diameter_distance, 1500*u.Mpc, method=method,
                                   bracket=(0.3, 1.0)), 0.6812777, rtol=1e-6)
        assert allclose(z_at_value(cosmo.angular_diameter_distance, 1500*u.Mpc, method=method,
                                   bracket=(2.0, 4.0)), 3.7914908, rtol=1e-6)
        assert allclose(z_at_value(cosmo.angular_diameter_distance, 1500*u.Mpc, method=method,
                                   bracket=(0.1, 1.5)), 0.6812777, rtol=1e-6)
        assert allclose(z_at_value(cosmo.angular_diameter_distance, 1500*u.Mpc, method=method,
                                   bracket=(0.1, 1.0, 2.0)), 0.6812777, rtol=1e-6)
        with pytest.warns(AstropyUserWarning, match=r'fval is not bracketed'):
            assert allclose(z_at_value(cosmo.angular_diameter_distance, 1500*u.Mpc, method=method,
                                       bracket=(0.9, 1.5)), 0.6812777, rtol=1e-6)
            assert allclose(z_at_value(cosmo.angular_diameter_distance, 1500*u.Mpc, method=method,
                                       bracket=(1.6, 2.0)), 3.7914908, rtol=1e-6)
        assert allclose(z_at_value(cosmo.angular_diameter_distance, 1500*u.Mpc, method=method,
                                   bracket=(1.6, 2.0), zmax=1.6), 0.6812777, rtol=1e-6)
        assert allclose(z_at_value(cosmo.angular_diameter_distance, 1500*u.Mpc, method=method,
                                   bracket=(0.9, 1.5), zmin=1.5), 3.7914908, rtol=1e-6)

    with pytest.raises(core.CosmologyError):
        with pytest.warns(AstropyUserWarning, match=r'fval is not bracketed'):
            z_at_value(cosmo.angular_diameter_distance, 1500*u.Mpc, method=method,
                       bracket=(3.9, 5.0), zmin=4.)


@pytest.mark.skipif('not HAS_SCIPY')
@pytest.mark.parametrize('method', ['Brent', 'Golden', 'Bounded'])
def test_z_at_value_unconverged(method):
    """
    Test warnings on non-converged solution when setting `maxfun` to too small iteration number -
    only 'Bounded' returns status value and specific message.
    """
    z_at_value = funcs.z_at_value
    cosmo = Planck18
    ztol = {'Brent': [1e-4, 1e-4], 'Golden': [1e-3, 1e-2], 'Bounded': [1e-3, 1e-1]}

    if method == 'Bounded':
        status, message = 1, 'Maximum number of function calls reached.'
    else:
        status, message = None, 'Unsuccessful'
    diag = rf'Solver returned {status}: {message}'

    with pytest.warns(AstropyUserWarning, match=diag):
        z0 = z_at_value(cosmo.angular_diameter_distance, 1*u.Gpc, zmax=2, maxfun=13, method=method)
    with pytest.warns(AstropyUserWarning, match=diag):
        z1 = z_at_value(cosmo.angular_diameter_distance, 1*u.Gpc, zmin=2, maxfun=13, method=method)

    assert allclose(z0, 0.32442, rtol=ztol[method][0])
    assert allclose(z1, 8.18551, rtol=ztol[method][1])


@pytest.mark.skipif('not HAS_SCIPY')
@pytest.mark.parametrize('cosmo', [Planck13, Planck15, Planck18, WMAP5, WMAP7, WMAP9,
                                   core.LambdaCDM, core.FlatLambdaCDM, core.wpwaCDM, core.w0wzCDM,
                                   core.wCDM, core.FlatwCDM, core.w0waCDM, core.Flatw0waCDM])
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
    skip = ('Ok',
            'angular_diameter_distance_z1z2',
            'clone',
            'de_density_scale', 'w')
    if str(cosmo.name).startswith('WMAP'):
        skip += ('nu_relative_density', )

    import inspect
    methods = inspect.getmembers(cosmo, predicate=inspect.ismethod)

    for name, func in methods:
        if name.startswith('_') or name in skip:
            continue
        print(f'Round-trip testing {name}')
        fval = func(z)
        # we need zmax here to pick the right solution for
        # angular_diameter_distance and related methods.
        # Be slightly more generous with rtol than the default 1e-8
        # used in z_at_value
        assert allclose(z, funcs.z_at_value(func, fval, bracket=[0.3, 1.0], ztol=1e-12), rtol=2e-11)

    # Test distance functions between two redshifts; only for realizations
    if isinstance(cosmo.name, str):
        z2 = 2.0
        func_z1z2 = [
            lambda z1: cosmo._comoving_distance_z1z2(z1, z2),
            lambda z1: cosmo._comoving_transverse_distance_z1z2(z1, z2),
            lambda z1: cosmo.angular_diameter_distance_z1z2(z1, z2)
        ]
        for func in func_z1z2:
            fval = func(z)
            assert allclose(z, funcs.z_at_value(func, fval, zmax=1.5, ztol=1e-12), rtol=2e-11)
