# Licensed under a 3-clause BSD style license - see LICENSE.rst

from unittest import mock

import pytest

from astropy.io.fits import BinTableHDU, HDUList, Header, PrimaryHDU
from astropy.timeseries.io.kepler import kepler_fits_reader
from astropy.units import UnitsWarning
from astropy.utils.data import get_pkg_data_filename


def fake_header(extver, version, timesys, telescop):
    return Header(
        {
            "SIMPLE": "T",
            "BITPIX": 8,
            "NAXIS": 0,
            "EXTVER": extver,
            "VERSION": version,
            "TIMESYS": f"{timesys}",
            "TELESCOP": f"{telescop}",
        }
    )


def fake_hdulist(extver=1, version=2, timesys="TDB", telescop="KEPLER"):
    new_header = fake_header(extver, version, timesys, telescop)
    return [
        HDUList(
            hdus=[
                PrimaryHDU(header=new_header),
                BinTableHDU(header=new_header, name="LIGHTCURVE"),
            ]
        )
    ]


@mock.patch("astropy.io.fits.open", side_effect=fake_hdulist(telescop="MadeUp"))
def test_raise_telescop_wrong(mock_file):
    with pytest.raises(
        NotImplementedError,
        match=(
            r"MadeUp is not implemented, only KEPLER or TESS are supported through this"
            r" reader"
        ),
    ):
        kepler_fits_reader(None)


@mock.patch("astropy.io.fits.open", side_effect=fake_hdulist(extver=2))
def test_raise_extversion_kepler(mock_file):
    with pytest.raises(
        NotImplementedError, match=r"Support for KEPLER v2 files not yet implemented"
    ):
        kepler_fits_reader(None)


@mock.patch("astropy.io.fits.open", side_effect=fake_hdulist(extver=2, telescop="TESS"))
def test_raise_extversion_tess(mock_file):
    with pytest.raises(
        NotImplementedError, match=r"Support for TESS v2 files not yet implemented"
    ):
        kepler_fits_reader(None)


@mock.patch("astropy.io.fits.open", side_effect=fake_hdulist(timesys="TCB"))
def test_raise_timesys_kepler(mock_file):
    with pytest.raises(
        NotImplementedError,
        match=r"Support for TCB time scale not yet implemented in KEPLER reader",
    ):
        kepler_fits_reader(None)


@mock.patch(
    "astropy.io.fits.open", side_effect=fake_hdulist(timesys="TCB", telescop="TESS")
)
def test_raise_timesys_tess(mock_file):
    with pytest.raises(
        NotImplementedError,
        match=r"Support for TCB time scale not yet implemented in TESS reader",
    ):
        kepler_fits_reader(None)


@pytest.mark.remote_data(source="astropy")
def test_kepler_astropy():
    from astropy.units import UnitsWarning

    filename = get_pkg_data_filename("timeseries/kplr010666592-2009131110544_slc.fits")

    with pytest.warns(UnitsWarning):
        timeseries = kepler_fits_reader(filename)

    assert timeseries["time"].format == "isot"
    assert timeseries["time"].scale == "tdb"
    assert timeseries["sap_flux"].unit.to_string() == "electron / s"
    assert len(timeseries) == 14280
    assert len(timeseries.columns) == 20


@pytest.mark.remote_data(source="astropy")
def test_tess_astropy():
    filename = get_pkg_data_filename(
        "timeseries/hlsp_tess-data-alerts_tess_phot_00025155310-s01_tess_v1_lc.fits"
    )
    with pytest.warns((UserWarning, UnitsWarning)) as record:
        timeseries = kepler_fits_reader(filename)

    # we might hit some warnings more than once, but the exact sequence probably
    # does not matter too much, so we'll just try to match the *set* of unique warnings
    unique_warnings = {(wm.category, wm.message.args[0]) for wm in record}
    expected = {
        (UserWarning, "Ignoring 815 rows with NaN times"),
        (
            UnitsWarning,
            (
                "'BJD - 2457000, days' did not parse as fits unit: "
                "At col 0, Unit 'BJD' not supported by the FITS standard.  "
                "If this is meant to be a custom unit, define it with 'u.def_unit'. "
                "To have it recognized inside a file reader or other code, "
                "enable it with 'u.add_enabled_units'. "
                "For details, see https://docs.astropy.org/en/latest/units/combining_and_defining.html"
            ),
        ),
        (
            UnitsWarning,
            (
                "'pixels' did not parse as fits unit: "
                "At col 0, Unit 'pixels' not supported by the FITS standard. "
                "Did you mean pixel? If this is meant to be a custom unit, "
                "define it with 'u.def_unit'. To have it recognized inside a file "
                "reader or other code, enable it with 'u.add_enabled_units'. "
                "For details, see https://docs.astropy.org/en/latest/units/combining_and_defining.html"
            ),
        ),
        (
            UnitsWarning,
            (
                "'e-/s' did not parse as fits unit: "
                "At col 0, Unit 'e' not supported by the FITS standard.  "
                "If this is meant to be a custom unit, define it with 'u.def_unit'. "
                "To have it recognized inside a file reader or other code, "
                "enable it with 'u.add_enabled_units'. "
                "For details, see https://docs.astropy.org/en/latest/units/combining_and_defining.html"
            ),
        ),
    }
    assert unique_warnings == expected, (
        f"Got some unexpected warnings\n{unique_warnings - expected}"
    )
    assert timeseries["time"].format == "isot"
    assert timeseries["time"].scale == "tdb"
    assert timeseries["sap_flux"].unit.to_string() == "electron / s"
    assert len(timeseries) == 19261
    assert len(timeseries.columns) == 20
