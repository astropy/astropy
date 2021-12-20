# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

from astropy.cosmology.parameter import Parameter
from astropy.table import Column
from astropy.modeling import Parameter as ModelParameter


def _convert_Parameter_to_Column(parameter, value, param_meta=None):
    """Convert a Cosmology Parameter to a Table Column.

    Parameters
    ----------
    parameter : `astropy.cosmology.parameter.Parameter`
    value : Any
    meta : dict or None, optional
        Information from the Cosmology's metadata.

    Returns
    -------
    `astropy.table.Column`
    """
    if value is None:  # Create None column
        data = [None]
        dtype = np.object_
        format = None
    else:  # Create empty column. Assign later.
        data = None  # assigned later
        dtype = None
        format = parameter.format_spec

    col = Column(data=data,
                 name=parameter.name,
                 dtype=dtype,
                 shape=np.shape(value),  # minimum of 1d
                 length=1,  # Cosmology is scalar
                 description=parameter.__doc__,
                 unit=parameter.unit,
                 format=format,
                 meta=param_meta,
                 copy=False,
                 copy_indices=True)
    if value is not None:
        col[:] = value  # Assign value in-place

    return col


def _convert_Parameter_to_Model_Parameter(parameter, value, meta=None):
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
    # Get from meta information relavant to Model
    extra = {k: v for k, v in (meta or {}).items()
             if k in ('getter', 'setter', 'fixed', 'tied', 'min', 'max',
                      'bounds', 'prior', 'posterior')}

    return ModelParameter(description=parameter.__doc__,
                          default=value,
                          unit=getattr(value, "unit", None),
                          **extra)
