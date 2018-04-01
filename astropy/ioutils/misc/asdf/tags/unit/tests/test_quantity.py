# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import io
import pytest

from astropy import units

asdf = pytest.importorskip('asdf', minversion='2.0.0')
from asdf.tests import helpers


def roundtrip_quantity(yaml, quantity):
    buff = helpers.yaml_to_asdf(yaml)
    with asdf.AsdfFile.open(buff) as ff:
        assert (ff.tree['quantity'] == quantity).all()
        buff2 = io.BytesIO()
        ff.write_to(buff2)

    buff2.seek(0)
    with asdf.AsdfFile.open(buff2) as ff:
        assert (ff.tree['quantity'] == quantity).all()

def test_value_scalar(tmpdir):
    testval = 2.71828
    testunit = units.kpc
    yaml = """
quantity: !unit/quantity-1.1.0
    value: {}
    unit: {}
""".format(testval, testunit)

    quantity = units.Quantity(testval, unit=testunit)
    roundtrip_quantity(yaml, quantity)

def test_value_array(tmpdir):
    testval = [3.14159]
    testunit = units.kg
    yaml = """
quantity: !unit/quantity-1.1.0
    value: !core/ndarray-1.0.0 {}
    unit: {}
""".format(testval, testunit)

    quantity = units.Quantity(testval, unit=testunit)
    roundtrip_quantity(yaml, quantity)

def test_value_multiarray(tmpdir):
    testval = [x*2.3081 for x in range(10)]
    testunit = units.ampere
    yaml = """
quantity: !unit/quantity-1.1.0
    value: !core/ndarray-1.0.0 {}
    unit: {}
""".format(testval, testunit)

    quantity = units.Quantity(testval, unit=testunit)
    roundtrip_quantity(yaml, quantity)

def test_value_ndarray(tmpdir):
    from numpy import array, float64
    testval = [[1,2,3],[4,5,6]]
    testunit = units.km
    yaml = """
quantity: !unit/quantity-1.1.0
    value: !core/ndarray-1.0.0
        datatype: float64
        data:
            {}
    unit: {}
""".format(testval, testunit)

    data = array(testval, float64)
    quantity = units.Quantity(data, unit=testunit)
    roundtrip_quantity(yaml, quantity)
