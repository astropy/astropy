import pytest
import dataclasses
import numpy as np
import astropy.units as u
import astropy.units.mixins


@dataclasses.dataclass
class DuckQuantity(astropy.units.mixins.QuantityOperatorsMixin):

    data: float | np.ndarray | u.Quantity

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        inputs = [
            inp.data if isinstance(inp, DuckQuantity) else inp
            for inp in inputs
        ]

        return DuckQuantity(
            getattr(ufunc, method)(*inputs, **kwargs)
        )


@pytest.mark.parametrize(
    argnames='operand_left',
    argvalues=[
        DuckQuantity(3),
        DuckQuantity(np.array([3, 4])),
        DuckQuantity(3 * u.mm),
        DuckQuantity([3, 4] * u.mm),
    ]
)
@pytest.mark.parametrize('operand_right', [u.mm, u.s])
class TestQuantityOperatorsMixin:

    def _unwrap(self, operand):
        if isinstance(operand, DuckQuantity):
            return operand.data
        elif isinstance(operand, str):
            return u.Unit(operand)
        else:
            return operand

    def test__mul__(
            self,
            operand_left: str | u.UnitBase | u.FunctionUnitBase | DuckQuantity,
            operand_right: str | u.UnitBase | u.FunctionUnitBase | DuckQuantity,
    ):

        result = operand_left * operand_right

        operand_left_unwrapped = self._unwrap(operand_left)
        operand_right_unwrapped = self._unwrap(operand_right)
        result_expected = operand_left_unwrapped * operand_right_unwrapped

        assert isinstance(result, DuckQuantity)
        assert np.all(result.data == result_expected)

    def test__mul__reversed(
            self,
            operand_left: str | u.UnitBase | u.FunctionUnitBase | DuckQuantity,
            operand_right: str | u.UnitBase | u.FunctionUnitBase | DuckQuantity,
    ):
        self.test__mul__(
            operand_left=operand_right,
            operand_right=operand_left
        )

    def test__truediv__(
            self,
            operand_left: str | u.UnitBase | u.FunctionUnitBase | DuckQuantity,
            operand_right: str | u.UnitBase | u.FunctionUnitBase | DuckQuantity,
    ):

        result = operand_left / operand_right

        operand_left_unwrapped = self._unwrap(operand_left)
        operand_right_unwrapped = self._unwrap(operand_right)
        result_expected = operand_left_unwrapped / operand_right_unwrapped

        assert isinstance(result, DuckQuantity)
        assert np.all(result.data == result_expected)

    def test__truediv__reversed(
            self,
            operand_left: str | u.UnitBase | u.FunctionUnitBase | DuckQuantity,
            operand_right: str | u.UnitBase | u.FunctionUnitBase | DuckQuantity,
    ):
        self.test__truediv__(
            operand_left=operand_right,
            operand_right=operand_left,
        )

    def test__lshift__(
            self,
            operand_left: str | u.UnitBase | u.FunctionUnitBase | DuckQuantity,
            operand_right: str | u.UnitBase | u.FunctionUnitBase | DuckQuantity,
    ):

        result = operand_left << operand_right

        operand_left_unwrapped = self._unwrap(operand_left)
        operand_right_unwrapped = self._unwrap(operand_right)
        result_expected = operand_left_unwrapped << operand_right_unwrapped

        assert isinstance(result, DuckQuantity)
        assert np.all(result.data == result_expected)
