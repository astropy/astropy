# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from numpy.lib.recfunctions import structured_to_unstructured

from astropy.cosmology.utils import is_structured
from astropy.modeling import Parameter as ModelParameter
from astropy.table import Column

FULLQUALNAME_SUBSTITUTIONS = {
    "astropy.cosmology.flrw.base.FLRW": "astropy.cosmology.flrw.FLRW",
    "astropy.cosmology.flrw.lambdacdm.LambdaCDM": "astropy.cosmology.flrw.LambdaCDM",
    "astropy.cosmology.flrw.lambdacdm.FlatLambdaCDM": "astropy.cosmology.flrw.FlatLambdaCDM",
    "astropy.cosmology.flrw.w0wacdm.w0waCDM": "astropy.cosmology.flrw.w0waCDM",
    "astropy.cosmology.flrw.w0wacdm.Flatw0waCDM": "astropy.cosmology.flrw.Flatw0waCDM",
    "astropy.cosmology.flrw.w0wzcdm.w0wzCDM": "astropy.cosmology.flrw.w0wzCDM",
    "astropy.cosmology.flrw.w0cdm.wCDM": "astropy.cosmology.flrw.wCDM",
    "astropy.cosmology.flrw.w0cdm.FlatwCDM": "astropy.cosmology.flrw.FlatwCDM",
    "astropy.cosmology.flrw.wpwazpcdm.wpwaCDM": "astropy.cosmology.flrw.wpwaCDM",
}
"""Substitutions mapping the actual qualitative name to its preferred value."""


def convert_parameter_to_column(parameter, value, meta=None, unstructure=False):
    """Convert a |Cosmology| Parameter to a Table |Column|.

    Parameters
    ----------
    parameter : `astropy.cosmology.parameter.Parameter`
    value : Any
    meta : dict or None, optional
        Information from the Cosmology's metadata.
    unstructure : bool (optional, keyword-only)
        Whether to convert parameters with structured dtypes to plain,
        unstructured arrays. If `True` (default) the column shape will be
        (number of rows, length of dtype).
        See :func:`~numpy.lib.recfunctions.structured_to_unstructured` for
        details.

    Returns
    -------
    `astropy.table.Column`
    """
    format = None if value is None else parameter.format_spec
    shape = (1,) + np.shape(value)  # minimum of 1d

    value = np.reshape(value, shape)
    if unstructure and is_structured(value):
        value = structured_to_unstructured(value)

    col = Column(data=value,
                 name=parameter.name,
                 dtype=None,  # inferred from the data
                 description=parameter.__doc__,
                 format=format,
                 meta=meta)

    return col


def convert_parameter_to_model_parameter(parameter, value, meta=None):
    """Convert a Cosmology Parameter to a Model Parameter.

    Parameters
    ----------
    parameter : `astropy.cosmology.parameter.Parameter`
    value : Any
    meta : dict or None, optional
        Information from the Cosmology's metadata.
        This function will use any of: 'getter', 'setter', 'fixed', 'tied',
        'min', 'max', 'bounds', 'prior', 'posterior'.

    Returns
    -------
    `astropy.modeling.Parameter`
    """
    if is_structured(value):
        value = structured_to_unstructured(value)

    # Get from meta information relavant to Model
    extra = {k: v for k, v in (meta or {}).items()
             if k in ('getter', 'setter', 'fixed', 'tied', 'min', 'max',
                      'bounds', 'prior', 'posterior')}

    return ModelParameter(description=parameter.__doc__,
                          default=value,
                          unit=getattr(value, "unit", None),
                          **extra)
