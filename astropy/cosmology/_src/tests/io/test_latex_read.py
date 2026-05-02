import pytest

from astropy.cosmology import Cosmology, Planck18


def test_latex_read_roundtrip(tmp_path):
    fp = tmp_path / "test_roundtrip.tex"
    from astropy.cosmology import FlatLambdaCDM

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3, name="TestCosmo")
    cosmo.write(fp, format="ascii.latex")

    # Read back
    cosmo_read = Cosmology.read(fp, format="ascii.latex")

    assert cosmo_read.name == "TestCosmo"
    assert cosmo_read.H0 == cosmo.H0
    assert cosmo_read.Om0 == cosmo.Om0
    assert cosmo_read.__class__ == cosmo.__class__


def test_latex_read_with_units(tmp_path):
    fp = tmp_path / "test_units.tex"
    Planck18.write(fp, format="ascii.latex")

    # This should now work and at least load the name and scalar parameters
    # m_nu might fail as discussed, let's see.
    try:
        cosmo_read = Cosmology.read(fp, format="ascii.latex")
        assert cosmo_read.name == "Planck18"
        assert cosmo_read.H0 == Planck18.H0
    except Exception as e:
        pytest.fail(f"Cosmology.read failed: {e}")


def test_latex_read_no_units(tmp_path):
    # Manually create a LaTeX table without a unit row
    fp = tmp_path / "test_no_units.tex"
    with open(fp, "w") as f:
        f.write(r"""\begin{table}
\begin{tabular}{cccc}
cosmology & name & H0 & Om0 \\
FlatLambdaCDM & Custom & 72.0 & 0.3 \\
\end{tabular}
\end{table}
""")

    cosmo = Cosmology.read(fp, format="ascii.latex")
    assert cosmo.name == "Custom"
    assert cosmo.H0.value == 72.0
    assert cosmo.Om0 == 0.3
