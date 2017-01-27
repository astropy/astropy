from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from collections import OrderedDict

import numpy as np
from .core import Fittable1DModel
from .parameters import Parameter
from .. import units as u

__all__ = ['BlackBody1D']


class BlackBody1D(Fittable1DModel):
    """
    One dimensional blackbody model.

    Parameters
    ----------
    temperature : :class:`~astropy.units.Quantity`
        Blackbody temperature
    bolometric_flux : :class:`~astropy.units.Quantity`
        The bolometric flux of the blackbody (i.e. the integral over the
        spectral axis)
    """

    # We parametrize this model with a temperature and a bolometric flux. The
    # bolometric flux is the integral of the model over the spectral axis. This
    # is more useful than simply having an amplitude parameter.
    temperature = Parameter(default=5000, min=0, unit=u.K)
    bolometric_flux = Parameter(default=1, unit=u.erg / u.cm ** 2 / u.s)

    # We allow values without units to be passed when evaluating the model, and
    # in this case the input x values are assumed to be frequencies in Hz.
    input_units_allow_dimensionless = True

    def evaluate(self, x, temperature, bolometric_flux):
        """Evaluate the model.

        Parameters
        ----------
        x : float or `~numpy.ndarray` or `~astropy.units.Quantity`
            Frquency at which to compute the blackbody. If no units are given,
            this defaults to Hz.

        temperature : float or `~numpy.ndarray` or `~astropy.units.Quantity`
            Temperature of the blackbody. If no units are given, this defaults
            to Kelvin.

        bolometric_flux : float or `~numpy.ndarray` or `~astropy.units.Quantity`
            Desired integral for the blackbody.

        Returns
        -------
        y : number or ndarray
            Blackbody spectrum. The units are determined from the units of
            ``bolometric_flux``.
        """

        from ..constants import sigma_sb
        from ..analytic_functions.blackbody import blackbody_nu

        # We need to make sure that we attach units to the temperature if it
        # doens't have any units. We do this because even though blackbody_nu
        # can take temperature values without units, the / temperature ** 4
        # factor needs units to be defined.
        temperature = u.Quantity(temperature, u.K)

        # We normalize the returned blackbody so that the integral would be
        # unity, and we then multiply by the bolometric flux. A normalized
        # blackbody has f_nu = pi * B_nu / (sigma * T^4), which is what we
        # calculate here. We convert to 1/Hz to make sure the units are
        # simplified as much as possible, then we multiply by the bolometric
        # flux to get the normalization right.
        fnu = (np.pi * u.sr * blackbody_nu(x, temperature) /
               sigma_sb / temperature ** 4).to(1 / u.Hz) * bolometric_flux

        # If the bolometric_flux parameter has no unit, we should drop the /Hz
        # and return a unitless value. This occurs for instance during fitting,
        # since we drop the units temporarily.
        if hasattr(bolometric_flux, 'unit'):
            return fnu
        else:
            return fnu.value

    @property
    def input_units(self):
        # The input units are those of the 'x' value, which should always be
        # Hz. Because we do this, and because input_units_allow_dimensionless
        # is set to True, dimensionless values are assumed to be in Hz.
        return {'x': u.Hz}

    def _parameter_units_for_data_units(self, xunit, yunit, zunit=None):
        return OrderedDict([('temperature', u.K),
                            ('bolometric_flux', yunit * u.Hz)])
