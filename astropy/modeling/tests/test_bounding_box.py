# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.modeling.models import Gaussian1D, Gaussian2D
from astropy.modeling.bounding_box import (BoundingBox, CompoundBoundingBox)

import numpy as np
import pytest
import unittest.mock as mk


def test_CompoundBoundingBox__init__():
    bbox = {1: (-1, 0), 2: (0, 1)}
    bounding_box = CompoundBoundingBox(bbox)

    assert bounding_box == bbox
    assert bounding_box._model is None
    assert bounding_box._slice_arg is None

    bounding_box = CompoundBoundingBox(bbox, Gaussian1D(), 'x')
    assert bounding_box == bbox
    assert (bounding_box._model.parameters == Gaussian1D().parameters).all()
    assert bounding_box._slice_arg == 0


def test_CompoundBoundingBox__get_arg_index():
    bounding_box = CompoundBoundingBox({}, Gaussian2D())

    assert bounding_box._get_arg_index(0) == 0
    assert bounding_box._get_arg_index(1) == 1
    with pytest.raises(ValueError):
        bounding_box._get_arg_index(2)

    assert bounding_box._get_arg_index('x') == 0
    assert bounding_box._get_arg_index('y') == 1
    with pytest.raises(ValueError):
        bounding_box._get_arg_index('z')


def test_CompoundBoundingBox_validate():
    bbox = {1: (-1, 0), 2: (0, 1)}
    model = Gaussian1D()
    bounding_box = CompoundBoundingBox.validate(model, bbox, 'x')

    assert bounding_box == bbox
    assert bounding_box._model == model
    assert bounding_box._slice_arg == 0
    for slice_box in bounding_box.values():
        assert isinstance(slice_box, BoundingBox)

    model = Gaussian2D()
    bounding_box = CompoundBoundingBox.validate(model, bbox, 'x',
                                               remove_slice_arg=True)
    assert bounding_box == bbox
    assert bounding_box._model == model
    assert bounding_box._slice_arg == 0
    for slice_box in bounding_box.values():
        assert isinstance(slice_box, BoundingBox)


def test_CompoundBoundingBox_set_slice_arg():
    bounding_box = CompoundBoundingBox((), slice_arg='arg')
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


def test_CompoundBoundingBox__get_slice_index():
    bounding_box = CompoundBoundingBox({}, Gaussian2D())

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


def test_CompoundBoundingBox_get_bounding_box():
    inputs = [np.array(1), np.array(2), np.array(3)]

    bounding_box = CompoundBoundingBox({}, Gaussian2D())
    assert bounding_box._slice_arg is None
    assert bounding_box.get_bounding_box(inputs) is None
    with pytest.raises(RuntimeError):
        bounding_box.get_bounding_box(inputs, slice_index=mk.MagicMock())

    bbox = {(1, 2): mk.MagicMock(), (3, 4): mk.MagicMock()}
    bounding_box = CompoundBoundingBox(bbox, Gaussian2D(), ('x', 'y'))
    assert bounding_box.get_bounding_box(inputs) == bbox[(1, 2)]
    assert bounding_box.get_bounding_box(inputs, slice_index=(3, 4)) == bbox[(3, 4)]
    with pytest.raises(RuntimeError):
        bounding_box.get_bounding_box([np.array(4), np.array(5)])
    bounding_box = CompoundBoundingBox(bbox, Gaussian2D(), (0, 1))

    assert bounding_box.get_bounding_box(inputs) == bbox[(1, 2)]
    assert bounding_box.get_bounding_box(inputs, slice_index=(3, 4)) == bbox[(3, 4)]
    with pytest.raises(RuntimeError):
        bounding_box.get_bounding_box([np.array(4), np.array(5)])

    bbox = {1: mk.MagicMock(), 2: mk.MagicMock()}
    bounding_box = CompoundBoundingBox(bbox, Gaussian2D(), 'x')
    assert bounding_box.get_bounding_box(inputs) == bbox[1]
    assert bounding_box.get_bounding_box(inputs, slice_index=2) == bbox[2]
    with pytest.raises(RuntimeError):
        bounding_box.get_bounding_box([np.array(3)])

    bounding_box = CompoundBoundingBox(bbox, Gaussian2D(), 0)
    assert bounding_box.get_bounding_box(inputs) == bbox[1]
    assert bounding_box.get_bounding_box(inputs, slice_index=2) == bbox[2]
    with pytest.raises(RuntimeError):
        bounding_box.get_bounding_box([np.array(3)])
