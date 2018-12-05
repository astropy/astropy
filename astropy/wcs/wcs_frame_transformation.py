from __future__ import print_function, absolute_import, division

from astropy import units as u
from astropy import constants
import warnings

from .utils import ExperimentalImplementationWarning

def _parse_velocity_convention(vc):
    if vc in (u.doppler_radio, 'radio', 'RADIO', 'VRAD', 'F', 'FREQ'):
        return u.doppler_radio
    elif vc in (u.doppler_optical, 'optical', 'OPTICAL', 'VOPT', 'W', 'WAVE'):
        return u.doppler_optical
    elif vc in (u.doppler_relativistic, 'relativistic', 'RELATIVE', 'VREL',
                'speed', 'V', 'VELO'):
        return u.doppler_relativistic

# These are the only linear transformations allowed
LINEAR_CTYPES = {u.doppler_optical: 'VOPT', u.doppler_radio: 'VRAD',
                 u.doppler_relativistic: 'VELO'}
LINEAR_CTYPE_CHARS = {u.doppler_optical: 'W', u.doppler_radio: 'F',
                      u.doppler_relativistic: 'V'}

ALL_CTYPES = {'speed': LINEAR_CTYPES,
              'frequency': 'FREQ',
              'length': 'WAVE'}

CTYPE_TO_PHYSICALTYPE = {'WAVE': 'length',
                         'AIR': 'air wavelength',
                         'AWAV': 'air wavelength',
                         'FREQ': 'frequency',
                         'VELO': 'speed',
                         'VRAD': 'speed',
                         'VOPT': 'speed',
                         }

CTYPE_CHAR_TO_PHYSICALTYPE = {'W': 'length',
                              'A': 'air wavelength',
                              'F': 'frequency',
                              'V': 'speed'}

CTYPE_TO_PHYSICALTYPE.update(CTYPE_CHAR_TO_PHYSICALTYPE)

PHYSICAL_TYPE_TO_CTYPE = dict([(v,k) for k,v in
                               CTYPE_CHAR_TO_PHYSICALTYPE.items()])
PHYSICAL_TYPE_TO_CHAR = {'speed': 'V',
                         'frequency': 'F',
                         'length': 'W'}

# Used to indicate the intial / final sampling system
WCS_UNIT_DICT = {'F': u.Hz, 'W': u.m, 'V': u.m/u.s}
PHYS_UNIT_DICT = {'length': u.m, 'frequency': u.Hz, 'speed': u.m/u.s}

LINEAR_CUNIT_DICT = {'VRAD': u.Hz, 'VOPT': u.m, 'FREQ': u.Hz, 'WAVE': u.m,
                     'VELO': u.m/u.s, 'AWAV': u.m}
LINEAR_CUNIT_DICT.update(WCS_UNIT_DICT)


def get_rest_value_from_wcs(mywcs):
    if mywcs.wcs.restfrq:
        ref_value = mywcs.wcs.restfrq*u.Hz
        return ref_value
    elif mywcs.wcs.restwav:
        ref_value = mywcs.wcs.restwav*u.m
        return ref_value

def convert_spectral_axis(mywcs, outunit, out_ctype, rest_value=None):
    """
    Convert a spectral axis from its unit to a specified out unit with a given
    output ctype

    Only VACUUM units are supported (not air)

    Process:
        1. Convert the input unit to its equivalent linear unit
        2. Convert the input linear unit to the output linear unit
        3. Convert the output linear unit to the output unit

    Parameters
    ----------
    mywcs : `astropy.wcs.WCS`
        The FITS WCS object to transform
    outunit : `astropy.units.Unit`
        The output unit, which must be a spectral unit (wavelength/frequency
        or velocity equivalent)
    rest_value : `astropy.units.Quantity`
        The rest wavelength/frequency to use when converting to/from velocity
    """

    # If the WCS includes a rest frequency/wavelength, convert it to frequency
    # or wavelength first.  This allows the possibility of changing the rest
    # frequency
    wcs_rv = get_rest_value_from_wcs(mywcs)
    inunit = u.Unit(mywcs.wcs.cunit[mywcs.wcs.spec])
    outunit = u.Unit(outunit)

    # If wcs_rv is set and speed -> speed, then we're changing the reference
    # location and we need to convert to meters or Hz first
    if ((inunit.physical_type == 'speed' and
         outunit.physical_type == 'speed' and
         wcs_rv is not None)):
        mywcs = convert_spectral_axis(mywcs, wcs_rv.unit,
                                      ALL_CTYPES[wcs_rv.unit.physical_type],
                                      rest_value=wcs_rv)
        inunit = u.Unit(mywcs.wcs.cunit[mywcs.wcs.spec])
    elif (inunit.physical_type == 'speed' and outunit.physical_type == 'speed'
          and wcs_rv is None):
        # If there is no reference change, we want an identical WCS, since
        # WCS doesn't know about units *at all*
        newwcs = mywcs.deepcopy()
        return newwcs
        #crval_out = (mywcs.wcs.crval[mywcs.wcs.spec] * inunit).to(outunit)
        #cdelt_out = (mywcs.wcs.cdelt[mywcs.wcs.spec] * inunit).to(outunit)
        #newwcs.wcs.cdelt[newwcs.wcs.spec] = cdelt_out.value
        #newwcs.wcs.cunit[newwcs.wcs.spec] = cdelt_out.unit.to_string(format='fits')
        #newwcs.wcs.crval[newwcs.wcs.spec] = crval_out.value
        #newwcs.wcs.ctype[newwcs.wcs.spec] = out_ctype
        #return newwcs

    in_spec_ctype = mywcs.wcs.ctype[mywcs.wcs.spec]

    # Check whether we need to convert the rest value first
    ref_value = None
    if outunit.physical_type == 'speed':
        if rest_value is None:
            rest_value = wcs_rv
            if rest_value is None:
                raise ValueError("If converting from wavelength/frequency to speed, "
                                 "a reference wavelength/frequency is required.")
        ref_value = rest_value.to(u.Hz, u.spectral())
    elif inunit.physical_type == 'speed':
        # The rest frequency and wavelength should be equivalent
        if rest_value is not None:
            ref_value = rest_value
        elif wcs_rv is not None:
            ref_value = wcs_rv
        else:
            raise ValueError("If converting from speed to wavelength/frequency, "
                             "a reference wavelength/frequency is required.")


    # If the input unit is not linearly sampled, its linear equivalent will be
    # the 8th character in the ctype, and the linearly-sampled ctype will be
    # the 6th character
    # e.g.: VOPT-F2V
    lin_ctype = (in_spec_ctype[7] if len(in_spec_ctype) > 4 else in_spec_ctype[:4])
    lin_cunit = (LINEAR_CUNIT_DICT[lin_ctype] if lin_ctype in LINEAR_CUNIT_DICT
                 else mywcs.wcs.cunit[mywcs.wcs.spec])
    in_vcequiv = _parse_velocity_convention(in_spec_ctype[:4])

    out_ctype_conv = out_ctype[7] if len(out_ctype) > 4 else out_ctype[:4]
    if CTYPE_TO_PHYSICALTYPE[out_ctype_conv] == 'air wavelength':
        raise NotImplementedError("Conversion to air wavelength is not supported.")
    out_lin_cunit = (LINEAR_CUNIT_DICT[out_ctype_conv] if out_ctype_conv in
                     LINEAR_CUNIT_DICT else outunit)
    out_vcequiv = _parse_velocity_convention(out_ctype_conv)

    # Load the input values
    crval_in = (mywcs.wcs.crval[mywcs.wcs.spec] * inunit)
    # the cdelt matrix may not be correctly populated: need to account for cd,
    # cdelt, and pc
    cdelt_in = (mywcs.pixel_scale_matrix[mywcs.wcs.spec, mywcs.wcs.spec] *
                inunit)

    if in_spec_ctype == 'AWAV':
        warnings.warn("Support for air wavelengths is experimental and only "
                      "works in the forward direction (air->vac, not vac->air).",
                      ExperimentalImplementationWarning
                     )
        cdelt_in = _greisen2006_air_to_vac_deriv(crval_in) * cdelt_in
        crval_in = _greisen2006_air_to_vac(crval_in)
        in_spec_ctype = 'WAVE'

    # 1. Convert input to input, linear
    if in_vcequiv is not None and ref_value is not None:
        crval_lin1 = crval_in.to(lin_cunit, u.spectral() + in_vcequiv(ref_value))
    else:
        crval_lin1 = crval_in.to(lin_cunit, u.spectral())
    cdelt_lin1 = _cdelt_derivative(crval_in,
                                   cdelt_in,
                                   # equivalent: inunit.physical_type
                                   intype=CTYPE_TO_PHYSICALTYPE[in_spec_ctype[:4]],
                                   outtype=lin_cunit.physical_type,
                                   rest=ref_value,
                                   linear=True
                                   )

    # 2. Convert input, linear to output, linear
    if ref_value is None:
        if in_vcequiv is not None:
            pass # consider raising a ValueError here; not clear if this is valid
        crval_lin2 = crval_lin1.to(out_lin_cunit, u.spectral())
    else:
        # at this stage, the transition can ONLY be relativistic, because the V
        # frame (as a linear frame) is only defined as "apparent velocity"
        crval_lin2 = crval_lin1.to(out_lin_cunit, u.spectral() +
                                   u.doppler_relativistic(ref_value))

    # For cases like VRAD <-> FREQ and VOPT <-> WAVE, this will be linear too:
    linear_middle = in_vcequiv == out_vcequiv

    cdelt_lin2 = _cdelt_derivative(crval_lin1, cdelt_lin1,
                                   intype=lin_cunit.physical_type,
                                   outtype=CTYPE_TO_PHYSICALTYPE[out_ctype_conv],
                                   rest=ref_value,
                                   linear=linear_middle)

    # 3. Convert output, linear to output
    if out_vcequiv is not None and ref_value is not None:
        crval_out = crval_lin2.to(outunit, out_vcequiv(ref_value) + u.spectral())
        #cdelt_out = cdelt_lin2.to(outunit, out_vcequiv(ref_value) + u.spectral())
        cdelt_out = _cdelt_derivative(crval_lin2,
                                      cdelt_lin2,
                                      intype=CTYPE_TO_PHYSICALTYPE[out_ctype_conv],
                                      outtype=outunit.physical_type,
                                      rest=ref_value,
                                      linear=True
                                      ).to(outunit)
    else:
        crval_out = crval_lin2.to(outunit, u.spectral())
        cdelt_out = cdelt_lin2.to(outunit, u.spectral())


    if crval_out.unit != cdelt_out.unit:
        # this should not be possible, but it's a sanity check
        raise ValueError("Conversion failed: the units of cdelt and crval don't match.")

    # A cdelt of 0 would be meaningless
    if cdelt_out.value == 0:
        raise ValueError("Conversion failed: the output CDELT would be 0.")

    newwcs = mywcs.deepcopy()
    if hasattr(newwcs.wcs,'cd'):
        newwcs.wcs.cd[newwcs.wcs.spec, newwcs.wcs.spec] = cdelt_out.value
        # todo: would be nice to have an assertion here that no off-diagonal
        # values for the spectral WCS are nonzero, but this is a nontrivial
        # check
    else:
        newwcs.wcs.cdelt[newwcs.wcs.spec] = cdelt_out.value
    newwcs.wcs.cunit[newwcs.wcs.spec] = cdelt_out.unit.to_string(format='fits')
    newwcs.wcs.crval[newwcs.wcs.spec] = crval_out.value
    newwcs.wcs.ctype[newwcs.wcs.spec] = out_ctype
    if rest_value is not None:
        if rest_value.unit.physical_type == 'frequency':
            newwcs.wcs.restfrq = rest_value.to(u.Hz).value
        elif rest_value.unit.physical_type == 'length':
            newwcs.wcs.restwav = rest_value.to(u.m).value
        else:
            raise ValueError("Rest Value was specified, but not in frequency or length units")

    return newwcs

def _cdelt_derivative(crval, cdelt, intype, outtype, linear=False, rest=None):
    """
    Determine the local derivative (``cdelt``) at a specified location
    (``crval``) going from one frame (``intype``) to another frame
    (``outtype``).

    This is not intended to be used outside of `convert_spectral_axis`, since
    the ``linear`` keyword needs to be predetermined.

    Parameters
    ----------
    crval : float
        The reference location in the input frame
    cdelt : float
        The cdelt (spacing between pixels) in the input frame defined at the
        input crval
    intype : string
    outtype : string
        The input and output frame types (i.e., length, frequency, speed)
    linear : bool
        Use a linear transformation?
    rest : None or quantity
        The rest value at which the velocity is defined
    """
    if intype == outtype:
        return cdelt
    elif set((outtype,intype)) == set(('length','frequency')):
        # Symmetric equations!
        return (-constants.c / crval**2 * cdelt).to(PHYS_UNIT_DICT[outtype])
    elif outtype in ('frequency','length') and intype == 'speed':
        if linear:
            numer = cdelt * rest.to(PHYS_UNIT_DICT[outtype], u.spectral())
            denom = constants.c
        else:
            numer = cdelt * constants.c * rest.to(PHYS_UNIT_DICT[outtype], u.spectral())
            denom = (constants.c + crval)*(constants.c**2 - crval**2)**0.5
        if outtype == 'frequency':
            return (-numer/denom).to(PHYS_UNIT_DICT[outtype], u.spectral())
        else:
            return (numer/denom).to(PHYS_UNIT_DICT[outtype], u.spectral())
    elif outtype == 'speed' and intype in ('frequency','length'):

        if linear:
            numer = cdelt * constants.c
            denom = rest.to(PHYS_UNIT_DICT[intype], u.spectral())
        else:
            numer = 4 * constants.c * crval * rest.to(crval.unit, u.spectral())**2 * cdelt
            denom = (crval**2 + rest.to(crval.unit, u.spectral())**2)**2
        if intype == 'frequency':
            return (-numer/denom).to(PHYS_UNIT_DICT[outtype], u.spectral())
        else:
            return (numer/denom).to(PHYS_UNIT_DICT[outtype], u.spectral())
    elif intype == 'air wavelength':
        raise TypeError("Air wavelength should be converted to vacuum earlier.")
    elif outtype == 'air wavelength':
        raise TypeError("Conversion to air wavelength not supported.")
    else:
        raise ValueError("Invalid in/out frames")


def _greisen2006_air_to_vac(wavelength):
    """
    Implements the air to vacuum wavelength conversion described in eqn 65 of
    Griesen 2006
    """
    wlum = wavelength.to(u.um).value
    return (1+1e-6*(287.6155+1.62887/wlum**2+0.01360/wlum**4)) * wavelength

def _greisen2006_vac_to_air(wavelength):
    """
    Griesen 2006 reports that the error in naively inverting Eqn 65 is less
    than 10^-9 and therefore acceptable.  This is therefore eqn 67
    """
    wlum = wavelength.to(u.um).value
    nl = (1+1e-6*(287.6155+1.62887/wlum**2+0.01360/wlum**4))
    return wavelength/nl

def _greisen2006_air_to_vac_deriv(wavelength):
    """
    Eqn 66
    """
    wlum = wavelength.to(u.um).value
    return (1+1e-6*(287.6155 - 1.62887/wlum**2 - 0.04080/wlum**4))
