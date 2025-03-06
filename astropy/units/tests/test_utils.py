# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test utilities for `astropy.units`.
"""

# ruff: noqa: FLY002

import textwrap
from importlib.util import module_from_spec, spec_from_file_location

import pytest

from astropy.tests.helper import _skip_docstring_tests_with_optimized_python


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


@_skip_docstring_tests_with_optimized_python
def test_unit_definition_module_summary_header():
    from astropy.units import si

    entries = si.__doc__.split("\n\n")
    assert len(entries) == 55
    assert entries[1] == "\n".join(
        (
            ".. list-table:: Available Units",
            "   :header-rows: 1",
            "   :widths: 10 20 20 20 1",
        )
    )
    assert entries[2] == "\n".join(
        (
            "   * - Unit",
            "     - Description",
            "     - Represents",
            "     - Aliases",
            "     - SI Prefixes",
        )
    )


@_skip_docstring_tests_with_optimized_python
def test_unit_definition_module_prefix_only_summary():
    from astropy.units import required_by_vounit

    entries = required_by_vounit.__doc__.split("\n\n")
    assert len(entries) == 6
    assert entries[3] == "\n".join(
        (
            "   * - Prefixes for ``solLum``",
            "     - Solar luminance prefixes",
            "     - :math:`\\mathrm{3.828 \\times 10^{26}\\,W}`",
            "     - ``L_sun``, ``Lsun``",
            "     - Only",
        )
    )
    assert entries[5] == "\n".join(
        (
            "   * - Prefixes for ``solRad``",
            "     - Solar radius prefixes",
            "     - :math:`\\mathrm{6.957 \\times 10^{8}\\,m}`",
            "     - ``R_sun``, ``Rsun``",
            "     - Only",
            "",
        )
    )


@_skip_docstring_tests_with_optimized_python
def test_unit_definition_module_normal_summary():
    from astropy.units import photometric

    entries = photometric.__doc__.split("\n\n")
    assert len(entries) == 10
    assert entries[5] == "\n".join(
        (
            "   * - ``AB``",
            "     - AB magnitude zero flux density (magnitude ``ABmag``).",
            "     - :math:`\\mathrm{3.6307805 \\times 10^{-20}\\,\\frac{erg}{Hz\\,s\\,cm^{2}}}`",
            "     - ``ABflux``",
            "     - No",
        )
    )
    assert entries[8] == "\n".join(
        (
            "   * - ``mgy``",
            "     - Maggies - a linear flux unit that is the flux for a mag=0 object.To tie this onto a specific calibrated unit system, the zero_point_flux equivalency should be used.",
            "     - ",
            "     - ``maggy``",
            "     - Yes",
        )
    )


@_skip_docstring_tests_with_optimized_python
def test_unit_definition_module_function_units_summary():
    from astropy.units.function import units

    entries = units.__doc__.split("\n\n")
    assert len(entries) == 14
    assert entries[8] == "\n".join(
        (
            ".. list-table:: Available Magnitude Units",
            "   :header-rows: 1",
            "   :widths: 10 50 10",
        )
    )
    assert entries[9] == "\n".join(
        (
            "   * - Unit",
            "     - Description",
            "     - Represents",
        )
    )
    assert entries[12] == "\n".join(
        (
            "   * - ``M_bol``",
            "     - Absolute bolometric magnitude: M_bol=0 corresponds to :math:`\\mathrm{3.0128 \\times 10^{28}\\,W}`",
            "     - mag(Bol)",
        )
    )
