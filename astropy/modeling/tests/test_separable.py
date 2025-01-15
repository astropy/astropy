# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test separability of models.

"""
import numpy as np

# pylint: disable=invalid-name
import pytest
from numpy.testing import assert_allclose

from astropy.modeling import custom_model, models
from astropy.modeling.core import ModelDefinitionError
from astropy.modeling.models import Mapping
from astropy.modeling.separable import (
    _arith_oper,
    _cdot,
    _coord_matrix,
    _cstack,
    is_separable,
    separability_matrix,
)

sh1 = models.Shift(1, name="shift1")
sh2 = models.Shift(2, name="sh2")
scl1 = models.Scale(1, name="scl1")
scl2 = models.Scale(2, name="scl2")
map1 = Mapping((0, 1, 0, 1), name="map1")
map2 = Mapping((0, 0, 1), name="map2")
map3 = Mapping((0, 0), name="map3")
rot = models.Rotation2D(2, name="rotation")
p2 = models.Polynomial2D(1, name="p2")
p22 = models.Polynomial2D(2, name="p22")
p1 = models.Polynomial1D(1, name="p1")


cm_4d_expected = (
    np.array([False, False, True, True]),
    np.array(
        [
            [True, True, False, False],
            [True, True, False, False],
            [False, False, True, False],
            [False, False, False, True],
        ]
    ),
)


compound_models = {
    "cm1": (
        map3 & sh1 | rot & sh1 | sh1 & sh2 & sh1,
        (
            np.array([False, False, True]),
            np.array([[True, False], [True, False], [False, True]]),
        ),
    ),
    "cm2": (
        sh1 & sh2 | rot | map1 | p2 & p22,
        (np.array([False, False]), np.array([[True, True], [True, True]])),
    ),
    "cm3": (
        map2 | rot & scl1,
        (
            np.array([False, False, True]),
            np.array([[True, False], [True, False], [False, True]]),
        ),
    ),
    "cm4": (
        sh1 & sh2 | map2 | rot & scl1,
        (
            np.array([False, False, True]),
            np.array([[True, False], [True, False], [False, True]]),
        ),
    ),
    "cm5": (
        map3 | sh1 & sh2 | scl1 & scl2,
        (np.array([False, False]), np.array([[True], [True]])),
    ),
    "cm7": (
        map2 | p2 & sh1,
        (np.array([False, True]), np.array([[True, False], [False, True]])),
    ),
    "cm8": (rot & (sh1 & sh2), cm_4d_expected),
    "cm9": (rot & sh1 & sh2, cm_4d_expected),
    "cm10": ((rot & sh1) & sh2, cm_4d_expected),
    "cm11": (
        rot & sh1 & (scl1 & scl2),
        (
            np.array([False, False, True, True, True]),
            np.array(
                [
                    [True, True, False, False, False],
                    [True, True, False, False, False],
                    [False, False, True, False, False],
                    [False, False, False, True, False],
                    [False, False, False, False, True],
                ]
            ),
        ),
    ),
}


def test_coord_matrix():
    c = _coord_matrix(p2, "left", 2)
    assert_allclose(np.array([[1, 1], [0, 0]]), c)
    c = _coord_matrix(p2, "right", 2)
    assert_allclose(np.array([[0, 0], [1, 1]]), c)
    c = _coord_matrix(p1, "left", 2)
    assert_allclose(np.array([[1], [0]]), c)
    c = _coord_matrix(p1, "left", 1)
    assert_allclose(np.array([[1]]), c)
    c = _coord_matrix(sh1, "left", 2)
    assert_allclose(np.array([[1], [0]]), c)
    c = _coord_matrix(sh1, "right", 2)
    assert_allclose(np.array([[0], [1]]), c)
    c = _coord_matrix(sh1, "right", 3)
    assert_allclose(np.array([[0], [0], [1]]), c)
    c = _coord_matrix(map3, "left", 2)
    assert_allclose(np.array([[1], [1]]), c)
    c = _coord_matrix(map3, "left", 3)
    assert_allclose(np.array([[1], [1], [0]]), c)


def test_cdot():
    result = _cdot(sh1, scl1)
    assert_allclose(result, np.array([[1]]))

    result = _cdot(rot, p2)
    assert_allclose(result, np.array([[2, 2]]))

    result = _cdot(rot, rot)
    assert_allclose(result, np.array([[2, 2], [2, 2]]))

    result = _cdot(Mapping((0, 0)), rot)
    assert_allclose(result, np.array([[2], [2]]))

    with pytest.raises(
        ModelDefinitionError,
        match=r"Models cannot be combined with the \"|\" operator; .*",
    ):
        _cdot(sh1, map1)


def test_cstack():
    result = _cstack(sh1, scl1)
    assert_allclose(result, np.array([[1, 0], [0, 1]]))

    result = _cstack(sh1, rot)
    assert_allclose(result, np.array([[1, 0, 0], [0, 1, 1], [0, 1, 1]]))
    result = _cstack(rot, sh1)
    assert_allclose(result, np.array([[1, 1, 0], [1, 1, 0], [0, 0, 1]]))


def test_arith_oper():
    # Models as inputs
    result = _arith_oper(sh1, scl1)
    assert_allclose(result, np.array([[1]]))
    result = _arith_oper(rot, rot)
    assert_allclose(result, np.array([[1, 1], [1, 1]]))

    # ndarray
    result = _arith_oper(np.array([[1, 2], [3, 4]]), np.array([[1, 2], [3, 4]]))
    assert_allclose(result, np.array([[1, 1], [1, 1]]))

    # Error
    with pytest.raises(
        ModelDefinitionError, match=r"Unsupported operands for arithmetic operator: .*"
    ):
        _arith_oper(sh1, map1)


@pytest.mark.parametrize(("compound_model", "result"), compound_models.values())
def test_separable(compound_model, result):
    assert_allclose(is_separable(compound_model), result[0])
    assert_allclose(separability_matrix(compound_model), result[1])


def test_custom_model_separable():
    @custom_model
    def model_a(x):
        return x

    assert model_a().separable

    @custom_model
    def model_c(x, y):
        return x + y

    assert not model_c().separable
    assert np.all(separability_matrix(model_c()) == [True, True])
