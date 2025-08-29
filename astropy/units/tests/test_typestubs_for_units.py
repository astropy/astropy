import sys
import textwrap
from pathlib import Path

from astropy.units.scripts import typestubs_for_units


def test_typestubs_for_units_script(tmp_path: Path, monkeypatch):
    """
    Tests the script that generates module-level .pyi stub files for astropy.units.

    This test runs the script's main function, directing its output to a temporary
    directory. It then performs spot checks on a few representative files
    to ensure key units and their types are correctly generated.
    """
    # monkeypatch command-line arguments to use the temporary directory
    monkeypatch.setattr(
        sys, "argv", ["typestub_for_units", "--output-dir", str(tmp_path)]
    )

    # --- generate the stub files ---
    typestubs_for_units.main()

    # --- spot-check the fundamental si.pyi file ---
    si_stub = tmp_path / "si.pyi"
    check_stub_file(
        si_stub,
        """
        from astropy.units.core import IrreducibleUnit
        from astropy.units.core import PrefixUnit
        from astropy.units.core import Unit

        deg: Unit
        km: PrefixUnit
        m: IrreducibleUnit
        s: IrreducibleUnit
        """,
    )
    # check docstring
    assert '"""kilometer (km)"""' in si_stub.read_text()

    # --- spot-check the special function/units.pyi file ---
    func_stub = tmp_path / "function" / "units.pyi"
    check_stub_file(
        func_stub,
        """
        from astropy.units.function.mixin import RegularFunctionUnit

        mag: RegularFunctionUnit
        dB: RegularFunctionUnit
        dex: IrreducibleFunctionUnit
        """,
    )

    # --- spot-check a representative astrophys.pyi file ---
    astrophys_stub = tmp_path / "astrophys.pyi"
    check_stub_file(
        astrophys_stub,
        """
        from astropy.units.core import Unit

        solMass: Unit
        solRad: Unit
        AU: PrefixUnit
        """,
    )


def check_stub_file(path: Path, expected_content: str):
    """Asserts that a stub file exists and contains all expected lines."""
    assert path.exists(), f"Stub file not found at {path}"
    content = path.read_text()
    assert content.startswith("# This stub file was automatically generated.")
    # textwrap.dedent allows us to write clean multiline strings in the test
    for line in textwrap.dedent(expected_content).strip().splitlines():
        assert line.strip() in content
