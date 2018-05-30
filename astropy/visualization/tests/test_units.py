# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import io

import pytest

try:
    import matplotlib.pyplot as plt
except ImportError:
    HAS_PLT = False
else:
    HAS_PLT = True

from ... import units as u
from ..units import quantity_support


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

    plt.clf()


@pytest.mark.skipif('not HAS_PLT')
def test_units_errbarr():
    pytest.importorskip("matplotlib", minversion="2.2")
    plt.figure()

    with quantity_support():
        x = [1, 2, 3] * u.s
        y = [1, 2, 3] * u.m
        yerr = [3, 2, 1] * u.cm

        fig, ax = plt.subplots()
        ax.errorbar(x, y, yerr=yerr)

        assert ax.xaxis.get_units() == u.s
        assert ax.yaxis.get_units() == u.m

    plt.clf()


@pytest.mark.skipif('not HAS_PLT')
def test_incompatible_units():
    plt.figure()

    with quantity_support():
        plt.plot([1, 2, 3] * u.m)
        with pytest.raises(u.UnitConversionError):
            plt.plot([105, 210, 315] * u.kg)

    plt.clf()
