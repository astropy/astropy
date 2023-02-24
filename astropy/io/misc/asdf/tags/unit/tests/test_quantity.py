# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io
import warnings

import pytest

asdf = pytest.importorskip("asdf")


from asdf.exceptions import AsdfDeprecationWarning

with warnings.catch_warnings():
    warnings.filterwarnings(
        "ignore",
        category=AsdfDeprecationWarning,
        message=r"asdf.tests.helpers is deprecated.*",
    )
    from asdf.tests.helpers import yaml_to_asdf

from astropy import units


def roundtrip_quantity(yaml, quantity):
    buff = yaml_to_asdf(yaml)
    with asdf.open(buff) as ff:
        assert (ff.tree["quantity"] == quantity).all()
        buff2 = io.BytesIO()
        ff.write_to(buff2)

    buff2.seek(0)
    with asdf.open(buff2) as ff:
        assert (ff.tree["quantity"] == quantity).all()


def test_value_scalar(tmpdir):
    testval = 2.71828
    testunit = units.kpc
    yaml = f"""
quantity: !unit/quantity-1.1.0
    value: {testval}
    unit: {testunit}
"""

    quantity = units.Quantity(testval, unit=testunit)
    roundtrip_quantity(yaml, quantity)


def test_value_array(tmpdir):
    testval = [3.14159]
    testunit = units.kg
    yaml = f"""
quantity: !unit/quantity-1.1.0
    value: !core/ndarray-1.0.0 {testval}
    unit: {testunit}
"""

    quantity = units.Quantity(testval, unit=testunit)
    roundtrip_quantity(yaml, quantity)


def test_value_multiarray(tmpdir):
    testval = [x * 2.3081 for x in range(10)]
    testunit = units.ampere
    yaml = f"""
quantity: !unit/quantity-1.1.0
    value: !core/ndarray-1.0.0 {testval}
    unit: {testunit}
"""

    quantity = units.Quantity(testval, unit=testunit)
    roundtrip_quantity(yaml, quantity)


def test_value_ndarray(tmpdir):
    from numpy import array, float64

    testval = [[1, 2, 3], [4, 5, 6]]
    testunit = units.km
    yaml = f"""
quantity: !unit/quantity-1.1.0
    value: !core/ndarray-1.0.0
        datatype: float64
        data:
            {testval}
    unit: {testunit}
"""

    data = array(testval, float64)
    quantity = units.Quantity(data, unit=testunit)
    roundtrip_quantity(yaml, quantity)
