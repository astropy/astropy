import numpy as np
import astropy.units as u

__all__ = [
    "QuantityOperatorsMixin",
]


class QuantityOperatorsMixin(np.lib.mixins.NDArrayOperatorsMixin):

    def __mul__(self, other):
        if isinstance(other, (u.UnitBase, str)):
            try:
                other = u.Quantity(1, other)
            except Exception:
                return NotImplemented
        return super().__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, (u.UnitBase, str)):
            try:
                other = u.Quantity(1, other)
            except Exception:
                return NotImplemented
        return super().__truediv__(other)

    def __lshift__(self, other):
        if isinstance(other, (u.UnitBase, str)):
            try:
                other = u.Quantity(1, other)
            except Exception:
                return NotImplemented
        return super().__lshift__(other)
