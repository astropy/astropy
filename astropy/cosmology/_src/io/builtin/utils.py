# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

from astropy.modeling import Parameter as ModelParameter
from astropy.table import Column

FULLQUALNAME_SUBSTITUTIONS = {
    "astropy.cosmology._src.flrw.base.FLRW": "astropy.cosmology.FLRW",
    "astropy.cosmology._src.flrw.lambdacdm.LambdaCDM": "astropy.cosmology.LambdaCDM",
    "astropy.cosmology._src.flrw.lambdacdm.FlatLambdaCDM": (
        "astropy.cosmology.FlatLambdaCDM"
    ),
    "astropy.cosmology._src.flrw.w0wacdm.w0waCDM": "astropy.cosmology.w0waCDM",
    "astropy.cosmology._src.flrw.w0wacdm.Flatw0waCDM": "astropy.cosmology.Flatw0waCDM",
    "astropy.cosmology._src.flrw.w0wzcdm.w0wzCDM": "astropy.cosmology.w0wzCDM",
    "astropy.cosmology._src.flrw.w0cdm.wCDM": "astropy.cosmology.wCDM",
    "astropy.cosmology._src.flrw.w0cdm.FlatwCDM": "astropy.cosmology.FlatwCDM",
    "astropy.cosmology._src.flrw.wpwazpcdm.wpwaCDM": "astropy.cosmology.wpwaCDM",
    # ===== deprecated paths =====
    "astropy.cosmology.flrw.base.FLRW": "astropy.cosmology.FLRW",
    "astropy.cosmology.flrw.lambdacdm.LambdaCDM": "astropy.cosmology.LambdaCDM",
    "astropy.cosmology.flrw.lambdacdm.FlatLambdaCDM": (
        "astropy.cosmology.flrw.FlatLambdaCDM"
    ),
    "astropy.cosmology.flrw.w0wacdm.w0waCDM": "astropy.cosmology.w0waCDM",
    "astropy.cosmology.flrw.w0wacdm.Flatw0waCDM": "astropy.cosmology.Flatw0waCDM",
    "astropy.cosmology.flrw.w0wzcdm.w0wzCDM": "astropy.cosmology.w0wzCDM",
    "astropy.cosmology.flrw.w0cdm.wCDM": "astropy.cosmology.wCDM",
    "astropy.cosmology.flrw.w0cdm.FlatwCDM": "astropy.cosmology.FlatwCDM",
    "astropy.cosmology.flrw.wpwazpcdm.wpwaCDM": "astropy.cosmology.wpwaCDM",
}
"""Substitutions mapping the actual qualified name to its preferred value."""


def convert_parameter_to_column(parameter, value, meta=None):
    """Convert a |Cosmology| Parameter to a Table |Column|.

    Parameters
    ----------
    parameter : `astropy.cosmology._src.parameter.Parameter`
    value : Any
    meta : dict or None, optional
        Information from the Cosmology's metadata.

    Returns
    -------
    `astropy.table.Column`
    """
    shape = (1,) + np.shape(value)  # minimum of 1d

    col = Column(
        data=np.reshape(value, shape),
        name=parameter.name,
        dtype=None,  # inferred from the data
        description=parameter.__doc__,
        format=None,
        meta=meta,
    )

    return col


def convert_parameter_to_model_parameter(parameter, value, meta=None):
    """Convert a Cosmology Parameter to a Model Parameter.

    Parameters
    ----------
    parameter : `astropy.cosmology._src.parameter.Parameter`
    value : Any
    meta : dict or None, optional
        Information from the Cosmology's metadata.
        This function will use any of: 'getter', 'setter', 'fixed', 'tied',
        'min', 'max', 'bounds', 'prior', 'posterior'.

    Returns
    -------
    `astropy.modeling.Parameter`
    """
    # Get from meta information relevant to Model
    attrs = (
        "getter",
        "setter",
        "fixed",
        "tied",
        "min",
        "max",
        "bounds",
        "prior",
        "posterior",
    )
    extra = {k: v for k, v in (meta or {}).items() if k in attrs}

    return ModelParameter(
        description=parameter.__doc__,
        default=value,
        unit=getattr(value, "unit", None),
        **extra,
    )
