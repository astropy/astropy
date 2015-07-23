# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import io

try:
    import matplotlib.pyplot as plt
except ImportError:
    HAS_PLT = False
else:
    HAS_PLT = True


from ...tests.helper import pytest


from ... import units as u
from ..units import quantity_support


@pytest.mark.skipif('not HAS_PLT')
def test_units():
    plt.figure()

    with quantity_support():
        buff = io.BytesIO()

        plt.plot([1, 2, 3] * u.m, [3, 4, 5] * u.kg)
        plt.plot([105, 210, 315] * u.cm, [3050, 3025, 3010] * u.g)
        plt.savefig(buff, format='svg')

        assert plt.gca().xaxis.get_units() == u.m
        assert plt.gca().yaxis.get_units() == u.kg

    plt.clf()


@pytest.mark.skipif('not HAS_PLT')
def test_incompatible_units():
    plt.figure()

    with quantity_support():
        plt.plot([1, 2, 3] * u.m)
        with pytest.raises(u.UnitConversionError):
            plt.plot([105, 210, 315] * u.kg)

    plt.clf()
