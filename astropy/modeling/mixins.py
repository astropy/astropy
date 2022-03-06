# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module provides utility mixins for the models package.
"""


class AmplitudeReturnUnitMixin:

    @property
    def return_units(self):
        if self.amplitude.unit is None:
            return None
        return dict.fromkeys(self.outputs, self.amplitude.unit)
