# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io

import pytest

try:
    import matplotlib.pyplot as plt
except ImportError:
    HAS_PLT = False
else:
    HAS_PLT = True

from astropy import units as u
from astropy.coordinates import Angle
from astropy.visualization.units import quantity_support


def teardown_function(function):
    plt.close('all')


@pytest.mark.skipif('not HAS_PLT')
def test_units():
    plt.figure()

    with quantity_support():
        buff = io.BytesIO()

        plt.plot([1, 2, 3] * u.m, [3, 4, 5] * u.kg, label='label')
        plt.plot([105, 210, 315] * u.cm, [3050, 3025, 3010] * u.g)
        plt.legend()
        # Also test fill_between, which requires actual conversion to ndarray
        # with numpy >=1.10 (#4654).
        plt.fill_between([1, 3] * u.m, [3, 5] * u.kg, [3050, 3010] * u.g)
        plt.savefig(buff, format='svg')

        assert plt.gca().xaxis.get_units() == u.m
        assert plt.gca().yaxis.get_units() == u.kg


@pytest.mark.skipif('not HAS_PLT')
def test_units_errbarr():
    pytest.importorskip("matplotlib")
    plt.figure()

    with quantity_support():
        x = [1, 2, 3] * u.s
        y = [1, 2, 3] * u.m
        yerr = [3, 2, 1] * u.cm

        fig, ax = plt.subplots()
        ax.errorbar(x, y, yerr=yerr)

        assert ax.xaxis.get_units() == u.s
        assert ax.yaxis.get_units() == u.m


@pytest.mark.skipif('not HAS_PLT')
def test_incompatible_units():
    # NOTE: minversion check does not work properly for matplotlib dev.
    try:
        # https://github.com/matplotlib/matplotlib/pull/13005
        from matplotlib.units import ConversionError
    except ImportError:
        err_type = u.UnitConversionError
    else:
        err_type = ConversionError

    plt.figure()

    with quantity_support():
        plt.plot([1, 2, 3] * u.m)
        with pytest.raises(err_type):
            plt.plot([105, 210, 315] * u.kg)


@pytest.mark.skipif('not HAS_PLT')
def test_quantity_subclass():
    """Check that subclasses are recognized.

    This sadly is not done by matplotlib.units itself, though
    there is a PR to change it:
    https://github.com/matplotlib/matplotlib/pull/13536
    """
    plt.figure()

    with quantity_support():
        plt.scatter(Angle([1, 2, 3], u.deg), [3, 4, 5] * u.kg)
        plt.scatter([105, 210, 315] * u.arcsec, [3050, 3025, 3010] * u.g)
        plt.plot(Angle([105, 210, 315], u.arcsec), [3050, 3025, 3010] * u.g)

        assert plt.gca().xaxis.get_units() == u.deg
        assert plt.gca().yaxis.get_units() == u.kg


@pytest.mark.skipif('not HAS_PLT')
def test_nested():

    with quantity_support():

        with quantity_support():

            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.scatter(Angle([1, 2, 3], u.deg), [3, 4, 5] * u.kg)

            assert ax.xaxis.get_units() == u.deg
            assert ax.yaxis.get_units() == u.kg

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.scatter(Angle([1, 2, 3], u.arcsec), [3, 4, 5] * u.pc)

        assert ax.xaxis.get_units() == u.arcsec
        assert ax.yaxis.get_units() == u.pc


@pytest.mark.skipif('not HAS_PLT')
def test_empty_hist():

    with quantity_support():
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.hist([1, 2, 3, 4] * u.mmag, bins=100)
        # The second call results in an empty list being passed to the
        # unit converter in matplotlib >= 3.1
        ax.hist([] * u.mmag, bins=100)
