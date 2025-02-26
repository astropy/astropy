# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test utilities for `astropy.units`.
"""

import textwrap
from importlib.util import module_from_spec, spec_from_file_location

import pytest


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
