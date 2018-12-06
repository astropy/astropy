from __future__ import print_function, absolute_import, division

from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy import constants
from astropy.tests.helper import pytest
import numpy as np
from numpy.testing import assert_allclose

from ...utils.data import get_pkg_data_filename as data_path
from ..wcs_frame_transformation import (convert_spectral_axis,
                                        determine_ctype_from_vconv,
                                        _cdelt_derivative,
                                        determine_vconv_from_ctype,
                                        get_rest_value_from_wcs, _greisen2006_air_to_vac,
                                        _greisen2006_air_to_vac_deriv, _greisen2006_vac_to_air,
                                        )

def test_cube_wcs_freqtovel():
    header = fits.Header.fromtextfile(data_path('data/cubewcs1.hdr'))
    w1 = wcs.WCS(header)
    # CTYPE3 = 'FREQ'

    newwcs = convert_spectral_axis(w1, 'km/s', 'VRAD',
                                   rest_value=w1.wcs.restfrq*u.Hz)
    assert newwcs.wcs.ctype[2] == 'VRAD'
    assert newwcs.wcs.crval[2] == 305.2461585938794
    assert newwcs.wcs.cunit[2] == u.Unit('km/s')

    newwcs = convert_spectral_axis(w1, 'km/s', 'VRAD')

    assert newwcs.wcs.ctype[2] == 'VRAD'
    assert newwcs.wcs.crval[2] == 305.2461585938794
    assert newwcs.wcs.cunit[2] == u.Unit('km/s')

def test_cube_wcs_freqtovopt():
    header = fits.Header.fromtextfile(data_path('data/cubewcs1.hdr'))
    w1 = wcs.WCS(header)

    w2 = convert_spectral_axis(w1, 'km/s', 'VOPT')

    # TODO: what should w2's values be?  test them

    # these need to be set to zero to test the failure
    w1.wcs.restfrq = 0.0
    w1.wcs.restwav = 0.0

    with pytest.raises(ValueError) as exc:
        convert_spectral_axis(w1, 'km/s', 'VOPT')

    assert exc.value.args[0] == 'If converting from wavelength/frequency to speed, a reference wavelength/frequency is required.'

@pytest.mark.parametrize('wcstype',('Z','W','R','V'))
def test_greisen2006(wcstype):
    # This is the header extracted from Greisen 2006, including many examples
    # of valid transforms.  It should be the gold standard (in principle)
    hdr = fits.Header.fromtextfile(data_path('data/greisen2006.hdr'))

    # We have not implemented frame conversions, so we can only convert bary
    # <-> bary in this case
    wcs0 = wcs.WCS(hdr, key='F')
    wcs1 = wcs.WCS(hdr, key=wcstype)

    if wcstype in ('R','V','Z'):
        if wcs1.wcs.restfrq:
            rest = wcs1.wcs.restfrq*u.Hz
        elif wcs1.wcs.restwav:
            rest = wcs1.wcs.restwav*u.m
    else:
        rest = None

    outunit = u.Unit(wcs1.wcs.cunit[wcs1.wcs.spec])
    out_ctype = wcs1.wcs.ctype[wcs1.wcs.spec]

    wcs2 = convert_spectral_axis(wcs0,
                                 outunit,
                                 out_ctype,
                                 rest_value=rest)
    assert_allclose(wcs2.wcs.cdelt[wcs2.wcs.spec],
                    wcs1.wcs.cdelt[wcs1.wcs.spec], rtol=1.e-3)
    assert_allclose(wcs2.wcs.crval[wcs2.wcs.spec],
                    wcs1.wcs.crval[wcs1.wcs.spec], rtol=1.e-3)
    assert wcs2.wcs.ctype[wcs2.wcs.spec] == wcs1.wcs.ctype[wcs1.wcs.spec]
    assert wcs2.wcs.cunit[wcs2.wcs.spec] == wcs1.wcs.cunit[wcs1.wcs.spec]

    # round trip test:
    inunit = u.Unit(wcs0.wcs.cunit[wcs0.wcs.spec])
    in_ctype = wcs0.wcs.ctype[wcs0.wcs.spec]
    wcs3 = convert_spectral_axis(wcs2,
                                 inunit,
                                 in_ctype,
                                 rest_value=rest)

    assert_allclose(wcs3.wcs.crval[wcs3.wcs.spec],
                    wcs0.wcs.crval[wcs0.wcs.spec], rtol=1.e-3)
    assert_allclose(wcs3.wcs.cdelt[wcs3.wcs.spec],
                    wcs0.wcs.cdelt[wcs0.wcs.spec], rtol=1.e-3)
    assert wcs3.wcs.ctype[wcs3.wcs.spec] == wcs0.wcs.ctype[wcs0.wcs.spec]
    assert wcs3.wcs.cunit[wcs3.wcs.spec] == wcs0.wcs.cunit[wcs0.wcs.spec]

def test_byhand_f2v():
    # VELO-F2V
    CRVAL3F = 1.37847121643E+09
    CDELT3F = 9.764775E+04
    RESTFRQV= 1.420405752E+09
    CRVAL3V = 8.98134229811E+06
    CDELT3V = -2.1217551E+04
    CUNIT3V = 'm/s'
    CUNIT3F = 'Hz'

    crvalf = CRVAL3F * u.Unit(CUNIT3F)
    crvalv = CRVAL3V * u.Unit(CUNIT3V)
    restfreq = RESTFRQV * u.Unit(CUNIT3F)
    cdeltf = CDELT3F * u.Unit(CUNIT3F)
    cdeltv = CDELT3V * u.Unit(CUNIT3V)

    # (Pdb) crval_in,crval_lin1,crval_lin2,crval_out
    # (<Quantity 1378471216.43 Hz>, <Quantity 1378471216.43 Hz>, <Quantity
    # 8981342.29795544 m / s>, <Quantity 8981342.29795544 m / s>) (Pdb)
    # cdelt_in, cdelt_lin1, cdelt_lin2, cdelt_out
    # (<Quantity 97647.75 Hz>, <Quantity 97647.75 Hz>, <Quantity
    # -21217.552294728768 m / s>, <Quantity -21217.552294728768 m / s>)
    crvalv_computed = crvalf.to(CUNIT3V, u.doppler_relativistic(restfreq))
    cdeltv_computed = -4*constants.c*cdeltf*crvalf*restfreq**2 / (crvalf**2+restfreq**2)**2
    cdeltv_computed_byfunction = _cdelt_derivative(crvalf, cdeltf,
                                                   intype='frequency',
                                                   outtype='speed',
                                                   rest=restfreq)
    # this should be EXACT
    assert cdeltv_computed == cdeltv_computed_byfunction

    assert_allclose(crvalv_computed, crvalv, rtol=1.e-3)
    assert_allclose(cdeltv_computed, cdeltv, rtol=1.e-3)

    # round trip
    # (Pdb) crval_in,crval_lin1,crval_lin2,crval_out
    # (<Quantity 8981342.29795544 m / s>, <Quantity 8981342.29795544 m / s>,
    # <Quantity 1377852479.159838 Hz>, <Quantity 1377852479.159838 Hz>)
    # (Pdb) cdelt_in, cdelt_lin1, cdelt_lin2, cdelt_out
    # (<Quantity -21217.552294728768 m / s>, <Quantity -21217.552294728768 m /
    # s>, <Quantity 97647.74999999997 Hz>, <Quantity 97647.74999999997 Hz>)

    crvalf_computed = crvalv_computed.to(CUNIT3F, u.doppler_relativistic(restfreq))
    cdeltf_computed = -(cdeltv_computed * constants.c * restfreq /
                        ((constants.c+crvalv_computed)*(constants.c**2 -
                                               crvalv_computed**2)**0.5))

    assert_allclose(crvalf_computed, crvalf, rtol=1.e-2)
    assert_allclose(cdeltf_computed, cdeltf, rtol=1.e-2)

    cdeltf_computed_byfunction = _cdelt_derivative(crvalv_computed,
                                                   cdeltv_computed,
                                                   intype='speed',
                                                   outtype='frequency',
                                                   rest=restfreq)
    # this should be EXACT
    assert cdeltf_computed == cdeltf_computed_byfunction

def test_byhand_vrad():
    # VRAD
    CRVAL3F = 1.37847121643E+09
    CDELT3F = 9.764775E+04
    RESTFRQR= 1.420405752E+09
    CRVAL3R = 8.85075090419E+06
    CDELT3R = -2.0609645E+04
    CUNIT3R = 'm/s'
    CUNIT3F = 'Hz'

    crvalf = CRVAL3F * u.Unit(CUNIT3F)
    crvalv = CRVAL3R * u.Unit(CUNIT3R)
    restfreq = RESTFRQR * u.Unit(CUNIT3F)
    cdeltf = CDELT3F * u.Unit(CUNIT3F)
    cdeltv = CDELT3R * u.Unit(CUNIT3R)

    # (Pdb) crval_in,crval_lin1,crval_lin2,crval_out
    # (<Quantity 1378471216.43 Hz>, <Quantity 1378471216.43 Hz>, <Quantity 8850750.904040769 m / s>, <Quantity 8850750.904040769 m / s>)
    # (Pdb) cdelt_in, cdelt_lin1, cdelt_lin2, cdelt_out
    # (<Quantity 97647.75 Hz>, <Quantity 97647.75 Hz>, <Quantity -20609.645482954576 m / s>, <Quantity -20609.645482954576 m / s>)
    crvalv_computed = crvalf.to(CUNIT3R, u.doppler_radio(restfreq))
    cdeltv_computed = -(cdeltf / restfreq)*constants.c

    assert_allclose(crvalv_computed, crvalv, rtol=1.e-3)
    assert_allclose(cdeltv_computed, cdeltv, rtol=1.e-3)

    crvalf_computed = crvalv_computed.to(CUNIT3F, u.doppler_radio(restfreq))
    cdeltf_computed = -(cdeltv_computed/constants.c) * restfreq

    assert_allclose(crvalf_computed, crvalf, rtol=1.e-3)
    assert_allclose(cdeltf_computed, cdeltf, rtol=1.e-3)

    # round trip:
    # (Pdb) crval_in,crval_lin1,crval_lin2,crval_out
    # (<Quantity 8850750.904040769 m / s>, <Quantity 8850750.904040769 m / s>, <Quantity 1378471216.43 Hz>, <Quantity 1378471216.43 Hz>)
    # (Pdb) cdelt_in, cdelt_lin1, cdelt_lin2, cdelt_out
    # (<Quantity -20609.645482954576 m / s>, <Quantity -20609.645482954576 m / s>, <Quantity 94888.9338036023 Hz>, <Quantity 94888.9338036023 Hz>)
    # (Pdb) myunit,lin_cunit,out_lin_cunit,outunit
    # WRONG (Unit("m / s"), Unit("m / s"), Unit("Hz"), Unit("Hz"))

def test_byhand_vopt():
    # VOPT: case "Z"
    CRVAL3F = 1.37847121643E+09
    CDELT3F = 9.764775E+04
    CUNIT3F = 'Hz'
    RESTWAVZ= 0.211061139
    #CTYPE3Z = 'VOPT-F2W'
    # This comes from Greisen 2006, but appears to be wrong: CRVAL3Z = 9.120000E+06
    CRVAL3Z = 9.120002206E+06
    CDELT3Z = -2.1882651E+04
    CUNIT3Z = 'm/s'

    crvalf = CRVAL3F * u.Unit(CUNIT3F)
    crvalv = CRVAL3Z * u.Unit(CUNIT3Z)
    restwav = RESTWAVZ * u.m
    cdeltf = CDELT3F * u.Unit(CUNIT3F)
    cdeltv = CDELT3Z * u.Unit(CUNIT3Z)

    # Forward: freq -> vopt
    # crval: (<Quantity 1378471216.43 Hz>, <Quantity 1378471216.43 Hz>, <Quantity 0.2174818410618759 m>, <Quantity 9120002.205689976 m / s>)
    # cdelt: (<Quantity 97647.75 Hz>, <Quantity 97647.75 Hz>, <Quantity -1.540591649098696e-05 m>, <Quantity -21882.652554887027 m / s>)
    #crvalv_computed = crvalf.to(CUNIT3R, u.doppler_radio(restwav))
    crvalw_computed = crvalf.to(u.m, u.spectral())
    crvalw_computed32 = crvalf.astype('float32').to(u.m, u.spectral())
    cdeltw_computed = -(cdeltf / crvalf**2)*constants.c
    cdeltw_computed_byfunction = _cdelt_derivative(crvalf, cdeltf,
                                                   intype='frequency',
                                                   outtype='length',
                                                   rest=None)
    # this should be EXACT
    assert cdeltw_computed == cdeltw_computed_byfunction

    crvalv_computed = crvalw_computed.to(CUNIT3Z, u.doppler_optical(restwav))
    crvalv_computed32 = crvalw_computed32.astype('float32').to(CUNIT3Z, u.doppler_optical(restwav))
    #cdeltv_computed = (cdeltw_computed *
    #                   4*constants.c*crvalw_computed*restwav**2 /
    #                   (restwav**2+crvalw_computed**2)**2)
    cdeltv_computed = (cdeltw_computed / restwav)*constants.c
    cdeltv_computed_byfunction = _cdelt_derivative(crvalw_computed,
                                                   cdeltw_computed,
                                                   intype='length',
                                                   outtype='speed',
                                                   rest=restwav,
                                                   linear=True)

    # Disagreement is 2.5e-7: good, but not really great...
    #assert np.abs((crvalv_computed-crvalv)/crvalv) < 1e-6
    assert_allclose(crvalv_computed, crvalv, rtol=1.e-2)
    assert_allclose(cdeltv_computed, cdeltv, rtol=1.e-2)

    # Round=trip test:
    # from velo_opt -> freq
    # (<Quantity 9120002.205689976 m / s>, <Quantity 0.2174818410618759 m>, <Quantity 1378471216.43 Hz>, <Quantity 1378471216.43 Hz>)
    # (<Quantity -21882.652554887027 m / s>, <Quantity -1.540591649098696e-05 m>, <Quantity 97647.75 Hz>, <Quantity 97647.75 Hz>)

    crvalw_computed = crvalv_computed.to(u.m, u.doppler_optical(restwav))
    cdeltw_computed = (cdeltv_computed/constants.c) * restwav
    cdeltw_computed_byfunction = _cdelt_derivative(crvalv_computed,
                                                   cdeltv_computed,
                                                   intype='speed',
                                                   outtype='length',
                                                   rest=restwav,
                                                   linear=True)
    assert cdeltw_computed == cdeltw_computed_byfunction

    crvalf_computed = crvalw_computed.to(CUNIT3F, u.spectral())
    cdeltf_computed = -cdeltw_computed * constants.c / crvalw_computed**2

    assert_allclose(crvalf_computed, crvalf, rtol=1.e-3)
    assert_allclose(cdeltf_computed, cdeltf, rtol=1.e-3)

    cdeltf_computed_byfunction = _cdelt_derivative(crvalw_computed,
                                                   cdeltw_computed,
                                                   intype='length',
                                                   outtype='frequency',
                                                   rest=None)
    assert cdeltf_computed == cdeltf_computed_byfunction

    # Fails intentionally (but not really worth testing)
    #crvalf_computed = crvalv_computed.to(CUNIT3F, u.spectral()+u.doppler_optical(restwav))
    #cdeltf_computed = -(cdeltv_computed / constants.c) * restwav.to(u.Hz, u.spectral())

    #assert_allclose(crvalf_computed, crvalf, rtol=1.e-3)
    #assert_allclose(cdeltf_computed, cdeltf, rtol=1.e-3)


def test_byhand_f2w():
    CRVAL3F = 1.37847121643E+09
    CDELT3F = 9.764775E+04
    CUNIT3F = 'Hz'
    #CTYPE3W = 'WAVE-F2W'
    CRVAL3W = 0.217481841062
    CDELT3W = -1.5405916E-05
    CUNIT3W = 'm'

    crvalf = CRVAL3F * u.Unit(CUNIT3F)
    crvalw = CRVAL3W * u.Unit(CUNIT3W)
    cdeltf = CDELT3F * u.Unit(CUNIT3F)
    cdeltw = CDELT3W * u.Unit(CUNIT3W)

    crvalf_computed = crvalw.to(CUNIT3F, u.spectral())
    cdeltf_computed = -constants.c * cdeltw / crvalw**2

    assert_allclose(crvalf_computed, crvalf, rtol=0.1)
    assert_allclose(cdeltf_computed, cdeltf, rtol=0.1)

@pytest.mark.parametrize(('ctype','unit','velocity_convention','result'),
                         (('VELO-F2V', "Hz", None, 'FREQ'),
                          ('VELO-F2V', "m", None, 'WAVE-F2W'),
                          ('VOPT', "m", None, 'WAVE'),
                          ('VOPT', "Hz", None, 'FREQ-W2F'),
                          ('VELO', "Hz", None, 'FREQ-V2F'),
                          ('WAVE', "Hz", None, 'FREQ-W2F'),
                          ('FREQ', 'm/s', None, ValueError('A velocity convention must be specified')),
                          ('FREQ', 'm/s', u.doppler_radio, 'VRAD'),
                          ('FREQ', 'm/s', u.doppler_optical, 'VOPT-F2W'),
                          ('FREQ', 'm/s', u.doppler_relativistic, 'VELO-F2V'),
                          ('WAVE', 'm/s', u.doppler_radio, 'VRAD-W2F')))
def test_ctype_determinator(ctype,unit,velocity_convention,result):

    if isinstance(result, Exception):
        with pytest.raises(Exception) as exc:
            determine_ctype_from_vconv(ctype, unit,
                                       velocity_convention=velocity_convention)
        assert exc.value.args[0] == result.args[0]
        assert type(exc.value) == type(result)
    else:
        outctype = determine_ctype_from_vconv(ctype, unit,
                                              velocity_convention=velocity_convention)
        assert outctype == result

@pytest.mark.parametrize(('ctype','vconv'),
                         (('VELO-F2W', u.doppler_optical),
                          ('VELO-F2V', u.doppler_relativistic),
                          ('VRAD', u.doppler_radio),
                          ('VOPT', u.doppler_optical),
                          ('VELO', u.doppler_relativistic),
                          ('WAVE', u.doppler_optical),
                          ('WAVE-F2W', u.doppler_optical),
                          ('WAVE-V2W', u.doppler_optical),
                          ('FREQ', u.doppler_radio),
                          ('FREQ-V2F', u.doppler_radio),
                          ('FREQ-W2F', u.doppler_radio),))
def test_vconv_determinator(ctype, vconv):
    assert determine_vconv_from_ctype(ctype) == vconv

# @pytest.mark.parametrize(('name'),
#                          (('advs'),
#                           ('dvsa'),
#                           ('sdav'),
#                           ('sadv'),
#                           ('vsad'),
#                           ('vad'),
#                           ('adv'),
#                           ))
# def test_vopt_to_freq(name):
#     h = fits.getheader(data_path(name+".fits"))
#     wcs0 = wcs.WCS(h)
#
#     # check to make sure astropy.wcs's "fix" changes VELO-HEL to VOPT
#     assert wcs0.wcs.ctype[wcs0.wcs.spec] == 'VOPT'
#
#     out_ctype = determine_ctype_from_vconv('VOPT', u.Hz)
#
#     wcs1 = convert_spectral_axis(wcs0, u.Hz, out_ctype)
#
#     assert wcs1.wcs.ctype[wcs1.wcs.spec] == 'FREQ-W2F'


@pytest.mark.parametrize('wcstype',('Z','W','R','V','F'))
def test_change_rest_frequency(wcstype):
    # This is the header extracted from Greisen 2006, including many examples
    # of valid transforms.  It should be the gold standard (in principle)
    hdr = fits.Header.fromtextfile(data_path('data/greisen2006.hdr'))

    wcs0 = wcs.WCS(hdr, key=wcstype)

    old_rest = get_rest_value_from_wcs(wcs0)
    if old_rest is None:
        # This test doesn't matter if there was no rest frequency in the first
        # place but I prefer to keep the option open in case we want to try
        # forcing a rest frequency on some of the non-velocity frames at some
        # point
        return
    vconv1 = determine_vconv_from_ctype(hdr['CTYPE3'+wcstype])
    new_rest = (100*u.km/u.s).to(u.Hz, vconv1(old_rest))

    wcs1 = wcs.WCS(hdr, key='V')
    vconv2 = determine_vconv_from_ctype(hdr['CTYPE3V'])

    inunit = u.Unit(wcs0.wcs.cunit[wcs0.wcs.spec])
    outunit = u.Unit(wcs1.wcs.cunit[wcs1.wcs.spec])
    # VELO-F2V
    out_ctype = wcs1.wcs.ctype[wcs1.wcs.spec]

    wcs2 = convert_spectral_axis(wcs0,
                                 outunit,
                                 out_ctype,
                                 rest_value=new_rest)

    sp1 = wcs1.sub([wcs.WCSSUB_SPECTRAL])
    sp2 = wcs2.sub([wcs.WCSSUB_SPECTRAL])

    p_old = sp1.wcs_world2pix([old_rest.to(inunit, vconv1(old_rest)).value,
                               new_rest.to(inunit, vconv1(old_rest)).value],0)
    p_new = sp2.wcs_world2pix([old_rest.to(outunit, vconv2(new_rest)).value,
                               new_rest.to(outunit, vconv2(new_rest)).value],0)

    assert_allclose(p_old, p_new, rtol=1e-3)
    assert_allclose(p_old, p_new, rtol=1e-3)


# from http://classic.sdss.org/dr5/products/spectra/vacwavelength.html
# these aren't accurate enough for my liking, but I can't find a better one readily
air_vac = {
    'H-beta':(4861.363, 4862.721)*u.AA,
    '[O III]':(4958.911, 4960.295)*u.AA,
    '[O III]':(5006.843, 5008.239)*u.AA,
    '[N II]':(6548.05, 6549.86)*u.AA,
    'H-alpha':(6562.801, 6564.614)*u.AA,
    '[N II]':(6583.45, 6585.27)*u.AA,
    '[S II]':(6716.44, 6718.29)*u.AA,
    '[S II]':(6730.82, 6732.68)*u.AA,
}


@pytest.mark.parametrize(('air','vac'), air_vac.values())
def test__greisen2006_air_to_vac(air, vac):
    # This is the accuracy provided by the line list we have.
    # I'm not sure if the formula are incorrect or if the reference wavelengths
    # are, but this is an accuracy of only 6 km/s, which is *very bad* for
    # astrophysical applications.
    assert np.abs((_greisen2006_air_to_vac(air)- vac)) < 0.15*u.AA
    assert np.abs((_greisen2006_vac_to_air(vac)- air)) < 0.15*u.AA

    assert np.abs((_greisen2006_air_to_vac(air)- vac)/vac) < 2e-5
    assert np.abs((_greisen2006_vac_to_air(vac)- air)/air) < 2e-5

    # round tripping
    assert np.abs((_greisen2006_vac_to_air(_greisen2006_air_to_vac(air))-air))/air < 1e-8
    assert np.abs((_greisen2006_air_to_vac(_greisen2006_vac_to_air(vac))-vac))/vac < 1e-8

def test_byhand_awav2vel():
    # AWAV
    CRVAL3A = (6560*u.AA).to(u.m).value
    CDELT3A = (1.0*u.AA).to(u.m).value
    CUNIT3A = 'm'
    CRPIX3A = 1.0
    # restwav MUST be vacuum
    restwl = _greisen2006_air_to_vac(6562.81*u.AA)
    RESTWAV = restwl.to(u.m).value
    CRVAL3V = (CRVAL3A*u.m).to(u.m/u.s,
                               u.doppler_optical(restwl)).value
    CDELT3V = (CDELT3A*u.m*_greisen2006_air_to_vac_deriv(CRVAL3A*u.m)/restwl) * constants.c
    CUNIT3V = 'm/s'

    mywcs = wcs.WCS(naxis=1)
    mywcs.wcs.ctype[0] = 'AWAV'
    mywcs.wcs.crval[0] = CRVAL3A
    mywcs.wcs.crpix[0] = CRPIX3A
    mywcs.wcs.cunit[0] = CUNIT3A
    mywcs.wcs.cdelt[0] = CDELT3A
    mywcs.wcs.restwav = RESTWAV
    mywcs.wcs.set()


    newwcs = convert_spectral_axis(mywcs, u.km/u.s,
                                   determine_ctype_from_vconv(mywcs.wcs.ctype[0],
                                                              u.km/u.s,
                                                              'optical'))

    newwcs.wcs.set()
    assert newwcs.wcs.cunit[0] == 'm / s'
    np.testing.assert_almost_equal(newwcs.wcs.crval,
                                   _greisen2006_air_to_vac(CRVAL3A*u.m)
                                   .to(u.m/u.s,
                                       u.doppler_optical(restwl)).value)
    # Check that the cdelts match the expected cdelt, 1 angstrom / rest
    # wavelength (vac)
    np.testing.assert_almost_equal(newwcs.wcs.cdelt, CDELT3V.to(u.m/u.s).value)
    # Check that the reference wavelength is 2.81 angstroms up
    np.testing.assert_almost_equal(newwcs.wcs_pix2world((2.81,), 0), 0.0, decimal=3)


    # Go through a full-on sanity check:
    vline = 100*u.km/u.s
    wave_line_vac = vline.to(u.AA, u.doppler_optical(restwl))
    wave_line_air = _greisen2006_vac_to_air(wave_line_vac)

    pix_line_input = mywcs.wcs_world2pix((wave_line_air.to(u.m).value,), 0)
    pix_line_output = newwcs.wcs_world2pix((vline.to(u.m/u.s).value,), 0)

    np.testing.assert_almost_equal(pix_line_output, pix_line_input, decimal=4)

def test_byhand_awav2wav():
    # AWAV
    CRVAL3A = (6560*u.AA).to(u.m).value
    CDELT3A = (1.0*u.AA).to(u.m).value
    CUNIT3A = 'm'
    CRPIX3A = 1.0

    mywcs = wcs.WCS(naxis=1)
    mywcs.wcs.ctype[0] = 'AWAV'
    mywcs.wcs.crval[0] = CRVAL3A
    mywcs.wcs.crpix[0] = CRPIX3A
    mywcs.wcs.cunit[0] = CUNIT3A
    mywcs.wcs.cdelt[0] = CDELT3A
    mywcs.wcs.set()


    newwcs = convert_spectral_axis(mywcs, u.AA, 'WAVE')
    newwcs.wcs.set()

    np.testing.assert_almost_equal(newwcs.wcs_pix2world((0,),0),
                                   _greisen2006_air_to_vac(mywcs.wcs_pix2world((0,),0)*u.m).value)

    np.testing.assert_almost_equal(newwcs.wcs_pix2world((10,),0),
                                   _greisen2006_air_to_vac(mywcs.wcs_pix2world((10,),0)*u.m).value)

    # At least one of the components MUST change
    assert not (mywcs.wcs.crval[0] == newwcs.wcs.crval[0]
                and mywcs.wcs.crpix[0] == newwcs.wcs.crpix[0])

class test_nir_sinfoni_base(object):
    def setup_method(self, method):
        CD3_3   = 0.000245000002905726 # CD rotation matrix
        CTYPE3  = 'WAVE    '           # wavelength axis in microns
        CRPIX3  =                1109. # Reference pixel in z
        CRVAL3  =     2.20000004768372 # central wavelength
        CDELT3  = 0.000245000002905726 # microns per pixel
        CUNIT3  = 'um      '           # spectral unit
        SPECSYS = 'TOPOCENT'           # Coordinate reference frame

        self.rest_wavelength = 2.1218*u.um

        self.mywcs = wcs.WCS(naxis=1)
        self.mywcs.wcs.ctype[0] = CTYPE3
        self.mywcs.wcs.crval[0] = CRVAL3
        self.mywcs.wcs.crpix[0] = CRPIX3
        self.mywcs.wcs.cunit[0] = CUNIT3
        self.mywcs.wcs.cdelt[0] = CDELT3
        self.mywcs.wcs.cd = [[CD3_3]]
        self.mywcs.wcs.specsys = SPECSYS
        self.mywcs.wcs.set()

        self.wavelengths = np.array([[2.12160005e-06,   2.12184505e-06,   2.12209005e-06]])

        np.testing.assert_almost_equal(self.mywcs.wcs_pix2world([788,789,790], 0),
                                       self.wavelengths)

    def test_nir_sinfoni_example_optical(self):

        mywcs = self.mywcs.copy()

        velocities_opt = ((self.wavelengths*u.m-self.rest_wavelength)/(self.wavelengths*u.m) * constants.c).to(u.km/u.s)

        newwcs_opt = convert_spectral_axis(mywcs, u.km/u.s, 'VOPT',
                                           rest_value=self.rest_wavelength)
        assert newwcs_opt.wcs.cunit[0] == u.km/u.s
        newwcs_opt.wcs.set()
        worldpix_opt = newwcs_opt.wcs_pix2world([788,789,790], 0)
        assert newwcs_opt.wcs.cunit[0] == u.m/u.s

        np.testing.assert_almost_equal(worldpix_opt,
                                       velocities_opt.to(newwcs_opt.wcs.cunit[0]).value)

    def test_nir_sinfoni_example_radio(self):

        mywcs = self.mywcs.copy()

        velocities_rad = ((self.wavelengths*u.m-self.rest_wavelength)/(self.rest_wavelength) * constants.c).to(u.km/u.s)

        newwcs_rad = convert_spectral_axis(mywcs, u.km/u.s, 'VRAD',
                                           rest_value=self.rest_wavelength)
        assert newwcs_rad.wcs.cunit[0] == u.km/u.s
        newwcs_rad.wcs.set()
        worldpix_rad = newwcs_rad.wcs_pix2world([788,789,790], 0)
        assert newwcs_rad.wcs.cunit[0] == u.m/u.s

        np.testing.assert_almost_equal(worldpix_rad,
                                       velocities_rad.to(newwcs_rad.wcs.cunit[0]).value)
