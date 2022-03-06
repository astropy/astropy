# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Module to test utility mixins
"""

import pytest

from ... import units as u
from .. import models
from ..functional_models import _InverseTrigonometric1D

models_with_amplitudes = [cls for cls in models.__dict__.values()
                          if hasattr(cls, 'amplitude')
                          and not issubclass(cls, _InverseTrigonometric1D)]


@pytest.mark.parametrize('cls', models_with_amplitudes)
def test_return_units(cls):
    model = cls(amplitude=4 * u.micron)
    return_units = model.return_units
    assert len(return_units) == model.n_outputs
    assert all(unit == u.micron for unit in return_units.values())
