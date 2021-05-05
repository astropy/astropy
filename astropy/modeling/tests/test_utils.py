# Licensed under a 3-clause BSD style license - see LICENSE.rst
# pylint: disable=invalid-name
import operator

import numpy as np
import pytest
import unittest.mock as mk

from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.modeling.utils import ExpressionTree as ET, ellipse_extent
from astropy.modeling.models import Ellipse2D, Gaussian1D, Gaussian2D

from astropy.modeling.utils import (_SpecialOperatorsDict,
                                    ComplexBoundingBox, _BoundingBox)


def test_traverse_postorder_duplicate_subtrees():
    """
    Regression test for a bug in `ExpressionTree.traverse_postorder`
    where given an expression like ``(1 + 2) + (1 + 2)`` where the two proper
    subtrees are actually the same object.
    """
    with pytest.warns(AstropyDeprecationWarning):
        subtree = ET('+', ET(1), ET(2))
        tree = ET('+', subtree, subtree)
    traversal = [n.value for n in tree.traverse_postorder()]
    assert traversal == [1, 2, '+', 1, 2, '+', '+']


def test_tree_evaluate_subexpression():
    """Test evaluating a subexpression from an expression tree."""

    operators = {'+': operator.add, '-': operator.sub, '*': operator.mul,
                 '/': operator.truediv, '**': operator.pow}
    # The full expression represented by this tree is:
    # 1.0 + 2 - 3 * 4 / 5 ** 6 (= 2.999232 if you must know)
    with pytest.warns(AstropyDeprecationWarning):
        tree = ET('+', ET(1.0), ET('-', ET(2.0),
                                   ET('*', ET(3.0), ET('/', ET(4.0),
                                                       ET('**', ET(5.0), ET(6.0))))))

    def test_slice(start, stop, expected):
        assert np.allclose(tree.evaluate(operators, start=start, stop=stop),
                           expected)

    assert tree.evaluate(operators) == (1.0 + 2.0 - 3.0 * 4.0 / 5.0 ** 6.0)
    test_slice(0, 5, (1.0 + 2.0 - 3.0 * 4.0 / 5.0))
    test_slice(0, 4, (1.0 + 2.0 - 3.0 * 4.0))
    test_slice(0, 3, (1.0 + 2.0 - 3.0))
    test_slice(0, 2, (1.0 + 2.0))
    test_slice(0, 1, 1.0)

    test_slice(1, 6, (2.0 - 3.0 * 4.0 / 5.0 ** 6.0))
    test_slice(1, 5, (2.0 - 3.0 * 4.0 / 5.0))
    test_slice(1, 4, (2.0 - 3.0 * 4.0))
    test_slice(1, 3, (2.0 - 3.0))
    test_slice(1, 2, 2.0)

    test_slice(2, 6, (3.0 * 4.0 / 5.0 ** 6.0))
    test_slice(2, 5, (3.0 * 4.0 / 5.0))
    test_slice(2, 4, (3.0 * 4.0))
    test_slice(2, 3, 3.0)

    test_slice(3, 6, (4.0 / 5.0 ** 6.0))
    test_slice(3, 5, (4.0 / 5.0))
    test_slice(3, 4, 4.0)

    test_slice(4, 6, (5.0 ** 6.0))
    test_slice(4, 5, 5.0)

    test_slice(5, 6, 6.0)


def test_ellipse_extent():
    # Test this properly bounds the ellipse

    imshape = (100, 100)
    coords = y, x = np.indices(imshape)

    amplitude = 1
    x0 = 50
    y0 = 50
    a = 30
    b = 10
    theta = np.pi / 4

    model = Ellipse2D(amplitude, x0, y0, a, b, theta)

    dx, dy = ellipse_extent(a, b, theta)

    limits = ((y0 - dy, y0 + dy), (x0 - dx, x0 + dx))

    model.bounding_box = limits

    actual = model.render(coords=coords)

    expected = model(x, y)

    # Check that the full ellipse is captured
    np.testing.assert_allclose(expected, actual, atol=0, rtol=1)

    # Check the bounding_box isn't too large
    limits = np.array(limits).flatten()
    for i in [0, 1]:
        s = actual.sum(axis=i)
        diff = np.abs(limits[2 * i] - np.where(s > 0)[0][0])
        assert diff < 1


def test_ComplexBoundingBox__init__():
    bbox = {1: (-1, 0), 2: (0, 1)}
    bounding_box = ComplexBoundingBox(bbox)

    assert bounding_box == bbox
    assert bounding_box._model is None
    assert bounding_box._slice_arg is None

    bounding_box = ComplexBoundingBox(bbox, Gaussian1D(), 'x')
    assert bounding_box == bbox
    assert (bounding_box._model.parameters == Gaussian1D().parameters).all()
    assert bounding_box._slice_arg == 0


def test_ComplexBoundingBox__get_arg_index():
    bounding_box = ComplexBoundingBox({}, Gaussian2D())

    assert bounding_box._get_arg_index(0) == 0
    assert bounding_box._get_arg_index(1) == 1
    with pytest.raises(ValueError):
        bounding_box._get_arg_index(2)

    assert bounding_box._get_arg_index('x') == 0
    assert bounding_box._get_arg_index('y') == 1
    with pytest.raises(ValueError):
        bounding_box._get_arg_index('z')


def test_ComplexBoundingBox_validate():
    bbox = {1: (-1, 0), 2: (0, 1)}
    model = Gaussian1D()
    bounding_box = ComplexBoundingBox.validate(model, bbox, 'x')

    assert bounding_box == bbox
    assert bounding_box._model == model
    assert bounding_box._slice_arg == 0
    for slice_box in bounding_box.values():
        assert isinstance(slice_box, _BoundingBox)

    model = Gaussian2D()
    bounding_box = ComplexBoundingBox.validate(model, bbox, 'x',
                                               remove_slice_arg=True)
    assert bounding_box == bbox
    assert bounding_box._model == model
    assert bounding_box._slice_arg == 0
    for slice_box in bounding_box.values():
        assert isinstance(slice_box, _BoundingBox)


def test_ComplexBoundingBox_set_slice_arg():
    bounding_box = ComplexBoundingBox((), slice_arg='arg')
    assert bounding_box._slice_arg == 'arg'

    bounding_box.set_slice_arg(None)
    assert bounding_box._slice_arg is None

    bounding_box._model = Gaussian1D()
    with pytest.raises(ValueError):
        bounding_box.set_slice_arg('arg')

    with pytest.raises(ValueError):
        bounding_box.set_slice_arg(('x', 'y'))

    bounding_box.set_slice_arg('x')
    assert bounding_box._slice_arg == 0
    bounding_box.set_slice_arg(0)
    assert bounding_box._slice_arg == 0

    bounding_box._model = Gaussian2D()
    bounding_box.set_slice_arg(('x', 'y'))
    assert bounding_box._slice_arg == (0, 1)
    bounding_box.set_slice_arg((0, 1))
    assert bounding_box._slice_arg == (0, 1)


def test_ComplexBoundingBox__get_slice_index():
    bounding_box = ComplexBoundingBox({}, Gaussian2D())

    inputs = [mk.MagicMock(), mk.MagicMock(), mk.MagicMock()]
    assert bounding_box._get_slice_index(inputs, 'x') == inputs[0]
    assert bounding_box._get_slice_index(inputs, 'y') == inputs[1]
    with pytest.raises(RuntimeError):
        bounding_box._get_slice_index(None, 'x')

    inputs = [np.array(1), np.array(2), np.array(3)]
    assert bounding_box._get_slice_index(inputs, 'x') == 1
    assert bounding_box._get_slice_index(inputs, 'y') == 2
    with pytest.raises(RuntimeError):
        bounding_box._get_slice_index(None, 'x')


def test_ComplexBoundingBox_get_bounding_box():
    inputs = [np.array(1), np.array(2), np.array(3)]

    bounding_box = ComplexBoundingBox({}, Gaussian2D())
    assert bounding_box._slice_arg is None
    assert bounding_box.get_bounding_box(inputs) is None
    with pytest.raises(RuntimeError):
        bounding_box.get_bounding_box(inputs, slice_index=mk.MagicMock())

    bbox = {(1, 2): mk.MagicMock(), (3, 4): mk.MagicMock()}
    bounding_box = ComplexBoundingBox(bbox, Gaussian2D(), ('x', 'y'))
    assert bounding_box.get_bounding_box(inputs) == bbox[(1, 2)]
    assert bounding_box.get_bounding_box(inputs, slice_index=(3, 4)) == bbox[(3, 4)]
    with pytest.raises(RuntimeError):
        bounding_box.get_bounding_box([np.array(4), np.array(5)])
    bounding_box = ComplexBoundingBox(bbox, Gaussian2D(), (0, 1))

    assert bounding_box.get_bounding_box(inputs) == bbox[(1, 2)]
    assert bounding_box.get_bounding_box(inputs, slice_index=(3, 4)) == bbox[(3, 4)]
    with pytest.raises(RuntimeError):
        bounding_box.get_bounding_box([np.array(4), np.array(5)])

    bbox = {1: mk.MagicMock(), 2: mk.MagicMock()}
    bounding_box = ComplexBoundingBox(bbox, Gaussian2D(), 'x')
    assert bounding_box.get_bounding_box(inputs) == bbox[1]
    assert bounding_box.get_bounding_box(inputs, slice_index=2) == bbox[2]
    with pytest.raises(RuntimeError):
        bounding_box.get_bounding_box([np.array(3)])

    bounding_box = ComplexBoundingBox(bbox, Gaussian2D(), 0)
    assert bounding_box.get_bounding_box(inputs) == bbox[1]
    assert bounding_box.get_bounding_box(inputs, slice_index=2) == bbox[2]
    with pytest.raises(RuntimeError):
        bounding_box.get_bounding_box([np.array(3)])


def test__SpecialOperatorsDict__set_value():
    key = 'test'
    val = 'value'

    special_operators = _SpecialOperatorsDict()
    assert key not in special_operators

    special_operators._set_value(key, val)
    assert key in special_operators
    assert special_operators[key] == val

    with pytest.raises(ValueError, match='Special operator "test" already exists'):
        special_operators._set_value(key, val)


def test__SpecialOperatorsDict___setitem__():
    key = 'test'
    val = 'value'

    special_operators = _SpecialOperatorsDict()
    assert key not in special_operators

    with pytest.deprecated_call():
        special_operators[key] = val
    assert key in special_operators
    assert special_operators[key] == val


def test__SpecialOperatorsDict__get_unique_id():
    special_operators = _SpecialOperatorsDict()
    assert special_operators._unique_id == 0

    assert special_operators._get_unique_id() == 1
    assert special_operators._unique_id == 1

    assert special_operators._get_unique_id() == 2
    assert special_operators._unique_id == 2

    assert special_operators._get_unique_id() == 3
    assert special_operators._unique_id == 3


def test__SpecialOperatorsDict_add():
    special_operators = _SpecialOperatorsDict()

    operator_name = 'test'
    operator = 'operator'

    key0 = special_operators.add(operator_name, operator)
    assert key0 == (operator_name, special_operators._unique_id)
    assert key0 in special_operators
    assert special_operators[key0] == operator

    key1 = special_operators.add(operator_name, operator)
    assert key1 == (operator_name, special_operators._unique_id)
    assert key1 in special_operators
    assert special_operators[key1] == operator

    assert key0 != key1
