# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test separability of models.

"""
# pylint: disable=invalid-name
import pytest
import numpy as np
from numpy.testing import assert_allclose

from astropy.modeling import custom_model, models
from astropy.modeling.models import Mapping
from astropy.modeling.separable import (_coord_matrix, is_separable, _cdot,
                                        _cstack, _arith_oper, separability_matrix)
from astropy.modeling.core import ModelDefinitionError


sh1 = models.Shift(1, name='shift1')
sh2 = models.Shift(2, name='sh2')
scl1 = models.Scale(1, name='scl1')
scl2 = models.Scale(2, name='scl2')
map1 = Mapping((0, 1, 0, 1), name='map1')
map2 = Mapping((0, 0, 1), name='map2')
map3 = Mapping((0, 0), name='map3')
rot = models.Rotation2D(2, name='rotation')
p2 = models.Polynomial2D(1, name='p2')
p22 = models.Polynomial2D(2, name='p22')
p1 = models.Polynomial1D(1, name='p1')


compound_models = {
    'cm1': (map3 & sh1 | rot & sh1 | sh1 & sh2 & sh1,
            (np.array([False, False, True]),
             np.array([[True, False], [True, False], [False, True]]))
            ),
    'cm2': (sh1 & sh2 | rot | map1 | p2 & p22,
            (np.array([False, False]),
             np.array([[True, True], [True, True]]))
            ),
    'cm3': (map2 | rot & scl1,
            (np.array([False, False, True]),
             np.array([[True, False], [True, False], [False, True]]))
            ),
    'cm4': (sh1 & sh2 | map2 | rot & scl1,
            (np.array([False, False, True]),
             np.array([[True, False], [True, False], [False, True]]))
            ),
    'cm5': (map3 | sh1 & sh2 | scl1 & scl2,
            (np.array([False, False]),
             np.array([[True], [True]]))
            ),
    'cm7': (map2 | p2 & sh1,
            (np.array([False, True]),
             np.array([[True, False], [False, True]]))
            )
}


def test_coord_matrix():
    c = _coord_matrix(p2, 'left', 2)
    assert_allclose(np.array([[1, 1], [0, 0]]), c)
    c = _coord_matrix(p2, 'right', 2)
    assert_allclose(np.array([[0, 0], [1, 1]]), c)
    c = _coord_matrix(p1, 'left', 2)
    assert_allclose(np.array([[1], [0]]), c)
    c = _coord_matrix(p1, 'left', 1)
    assert_allclose(np.array([[1]]), c)
    c = _coord_matrix(sh1, 'left', 2)
    assert_allclose(np.array([[1], [0]]), c)
    c = _coord_matrix(sh1, 'right', 2)
    assert_allclose(np.array([[0], [1]]), c)
    c = _coord_matrix(sh1, 'right', 3)
    assert_allclose(np.array([[0], [0], [1]]), c)
    c = _coord_matrix(map3, 'left', 2)
    assert_allclose(np.array([[1], [1]]), c)
    c = _coord_matrix(map3, 'left', 3)
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

    with pytest.raises(ModelDefinitionError,
                       match=r"Models cannot be combined with the \"|\" operator; .*"):
        _cdot(sh1, map1)


def test_cstack():
    result = _cstack(sh1, scl1)
    assert_allclose(result, np.array([[1, 0], [0, 1]]))

    result = _cstack(sh1, rot)
    assert_allclose(result,
                    np.array([[1, 0, 0],
                              [0, 1, 1],
                              [0, 1, 1]])
                    )
    result = _cstack(rot, sh1)
    assert_allclose(result,
                    np.array([[1, 1, 0],
                              [1, 1, 0],
                              [0, 0, 1]])
                    )


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
    with pytest.raises(ModelDefinitionError, match=r"Unsupported operands for arithmetic operator: .*"):
        _arith_oper(sh1, map1)


@pytest.mark.parametrize(('compound_model', 'result'), compound_models.values())
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


def test_nested_compound_model_separability():
    """
    Test that nested CompoundModels maintain correct separability.

    Regression test for issue where nested CompoundModels would incorrectly
    show inputs/outputs as non-separable when they should be separable.
    """
    # Create a simple separable compound model
    cm = models.Linear1D(10) & models.Linear1D(5)

    # Verify the simple compound model is separable
    expected_simple = np.array([
        [True, False],
        [False, True]
    ])
    assert_allclose(separability_matrix(cm), expected_simple)

    # Create a non-nested version for comparison
    direct = models.Pix2Sky_TAN() & models.Linear1D(10) & models.Linear1D(5)
    expected_direct = np.array([
        [True, True, False, False],
        [True, True, False, False],
        [False, False, True, False],
        [False, False, False, True]
    ])
    assert_allclose(separability_matrix(direct), expected_direct)

    # Now test the nested version - this should give the same result as direct
    nested = models.Pix2Sky_TAN() & cm
    expected_nested = np.array([
        [True, True, False, False],
        [True, True, False, False],
        [False, False, True, False],
        [False, False, False, True]
    ])
    assert_allclose(separability_matrix(nested), expected_nested)

    # Additional test: nested on the left side
    cm_left = models.Linear1D(10) & models.Linear1D(5)
    nested_left = cm_left & models.Pix2Sky_TAN()
    expected_nested_left = np.array([
        [True, False, False, False],
        [False, True, False, False],
        [False, False, True, True],
        [False, False, True, True]
    ])
    assert_allclose(separability_matrix(nested_left), expected_nested_left)

    # Test with multiple levels of nesting
    cm1 = models.Linear1D(1) & models.Linear1D(2)
    cm2 = models.Linear1D(3) & models.Linear1D(4)
    deeply_nested = cm1 & cm2
    expected_deeply_nested = np.array([
        [True, False, False, False],
        [False, True, False, False],
        [False, False, True, False],
        [False, False, False, True]
    ])
    assert_allclose(separability_matrix(deeply_nested), expected_deeply_nested)
