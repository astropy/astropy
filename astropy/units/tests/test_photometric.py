# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from ...tests.helper import assert_quantity_allclose

from .. import Magnitude, mgy, nmgy


def test_maggies():
    assert_quantity_allclose(1e-9*mgy, 1*nmgy)
    assert_quantity_allclose(Magnitude((1*nmgy).to(mgy)).value, 22.5)
