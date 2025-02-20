# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test utilities for `astropy.units`.
"""

import textwrap
from importlib.util import module_from_spec, spec_from_file_location

import numpy as np
import pytest

from astropy.units.quantity import Quantity
from astropy.units.utils import quantity_asanyarray


def test_quantity_asanyarray():
    array_of_quantities = [Quantity(1), Quantity(2), Quantity(3)]
    quantity_array = quantity_asanyarray(array_of_quantities)
    assert isinstance(quantity_array, Quantity)
    assert np.issubdtype(quantity_array.dtype, np.inexact)

    array_of_integers = [1, 2, 3]
    np_array = quantity_asanyarray(array_of_integers)
    assert isinstance(np_array, np.ndarray)
    assert np.issubdtype(np_array.dtype, np.integer)

    np_array = quantity_asanyarray(array_of_integers, dtype=np.inexact)
    assert np.issubdtype(np_array.dtype, np.inexact)


def test_composite_unit_definition(tmp_path):
    # Regression test for #17781 - error message was not informative
    module_path = tmp_path / "bad_module.py"
    module_path.write_text(
        textwrap.dedent("""
        from astropy import units as u
        from astropy.units.utils import generate_unit_summary

        km_per_h = u.km / u.h

        __doc__ = generate_unit_summary(globals())
    """)
    )
    spec = spec_from_file_location("bad_module", module_path)
    module = module_from_spec(spec)

    with pytest.raises(
        TypeError, match=r"^'km_per_h' must be defined with 'def_unit\(\)'$"
    ):
        # emulate an import statement
        spec.loader.exec_module(module)
