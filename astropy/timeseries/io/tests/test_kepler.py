# Licensed under a 3-clause BSD style license - see LICENSE.rst

from unittest import mock

import pytest

from astropy.io.fits import HDUList, Header, PrimaryHDU, BinTableHDU
from astropy.utils.data import get_pkg_data_filename

from astropy.timeseries.io.kepler import kepler_fits_reader


def fake_header(extver, version, timesys, telescop):
    return Header({"SIMPLE": "T",
                   "BITPIX": 8,
                   "NAXIS": 0,
                   "EXTVER": extver,
                   "VERSION": version,
                   'TIMESYS': f"{timesys}",
                   "TELESCOP": f"{telescop}"})


def fake_hdulist(extver=1, version=2, timesys="TDB", telescop="KEPLER"):
    new_header = fake_header(extver, version, timesys, telescop)
    return [HDUList(hdus=[PrimaryHDU(header=new_header),
                    BinTableHDU(header=new_header, name="LIGHTCURVE")])]


@mock.patch("astropy.io.fits.open", side_effect=fake_hdulist(telescop="MadeUp"))
def test_raise_telescop_wrong(mock_file):
    with pytest.raises(NotImplementedError) as exc:
        kepler_fits_reader(None)
    assert exc.value.args[0] == ("MadeUp is not implemented, only KEPLER or TESS are "
                                 "supported through this reader")


@mock.patch("astropy.io.fits.open", side_effect=fake_hdulist(extver=2))
def test_raise_extversion_kepler(mock_file):
    with pytest.raises(NotImplementedError) as exc:
        kepler_fits_reader(None)
    assert exc.value.args[0] == ("Support for KEPLER v2 files not yet "
                                 "implemented")


@mock.patch("astropy.io.fits.open", side_effect=fake_hdulist(extver=2, telescop="TESS"))
def test_raise_extversion_tess(mock_file):
    with pytest.raises(NotImplementedError) as exc:
        kepler_fits_reader(None)
    assert exc.value.args[0] == ("Support for TESS v2 files not yet "
                                 "implemented")


@mock.patch("astropy.io.fits.open", side_effect=fake_hdulist(timesys="TCB"))
def test_raise_timesys_kepler(mock_file):
    with pytest.raises(NotImplementedError) as exc:
        kepler_fits_reader(None)
    assert exc.value.args[0] == ("Support for TCB time scale not yet "
                                 "implemented in KEPLER reader")


@mock.patch("astropy.io.fits.open", side_effect=fake_hdulist(timesys="TCB", telescop="TESS"))
def test_raise_timesys_tess(mock_file):
    with pytest.raises(NotImplementedError) as exc:
        kepler_fits_reader(None)
    assert exc.value.args[0] == ("Support for TCB time scale not yet "
                                 "implemented in TESS reader")


@pytest.mark.remote_data(source='astropy')
def test_kepler_astropy():
    from astropy.units import UnitsWarning

    filename = get_pkg_data_filename('timeseries/kplr010666592-2009131110544_slc.fits')

    with pytest.warns(UnitsWarning):
        timeseries = kepler_fits_reader(filename)

    assert timeseries["time"].format == 'isot'
    assert timeseries["time"].scale == 'tdb'
    assert timeseries["sap_flux"].unit.to_string() == 'electron / s'
    assert len(timeseries) == 14280
    assert len(timeseries.columns) == 20


@pytest.mark.remote_data(source='astropy')
def test_tess_astropy():
    filename = get_pkg_data_filename('timeseries/hlsp_tess-data-alerts_tess_phot_00025155310-s01_tess_v1_lc.fits')
    with pytest.warns(UserWarning, match='Ignoring 815 rows with NaN times'):
        timeseries = kepler_fits_reader(filename)
    assert timeseries["time"].format == 'isot'
    assert timeseries["time"].scale == 'tdb'
    assert timeseries["sap_flux"].unit.to_string() == 'electron / s'
    assert len(timeseries) == 19261
    assert len(timeseries.columns) == 20
