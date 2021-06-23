# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.modeling.models import Gaussian1D, Gaussian2D
from astropy.modeling.bounding_box import (BoundingBox,
                                           _BaseModelArgument, ModelArgument,
                                           ModelArguments, CompoundBoundingBox)

import numpy as np
import pytest
import unittest.mock as mk


class TestModelArgument:
    def test_create(self):
        argument = ModelArgument('test', True, 7)
        assert isinstance(argument, tuple)
        assert isinstance(argument, _BaseModelArgument)
        assert argument.name == 'test' == argument[0]
        assert argument.remove == True == argument[1]
        assert argument.index == 7 == argument[2]
        assert argument == ('test', True, 7)

    def test__get_index(self):
        argument = ModelArgument('test', True, 7)

        # name is None
        assert argument._get_index(None, mk.MagicMock()) is None

        # model is None
        assert argument._get_index(mk.MagicMock(), None) is None
        assert argument._get_index(mk.MagicMock()) is None

        # neither is None
        model = mk.MagicMock()
        model.inputs = ['test', 'name']
        assert argument._get_index('test', model) == 0
        assert argument._get_index('name', model) == 1

        # Error
        with pytest.raises(ValueError) as err:
            argument._get_index('other', model)
        assert str(err.value) == \
            "other is not an input of your model inputs: ['test', 'name']."

    def test_validate(self):
        model = mk.MagicMock()
        model.inputs = [mk.MagicMock(), 'name']

        # Success with full data
        assert ModelArgument.validate(model, 'name', True, 1) == ('name', 1, True)

        # Fail with full data
        with pytest.raises(IndexError) as err:
            ModelArgument.validate(model, 'name', index=2)
        assert str(err.value) == \
            "Index should be 1, but was given 2."

        # Success with missing inputs
        assert ModelArgument.validate(model, 'name')  == ('name', False, 1)
        assert ModelArgument.validate(model, 'name', index=1) == ('name', False, 1)
        assert ModelArgument.validate(name='name', index=1)  == ('name', False, 1)

        # Fail with missing inputs
        with pytest.raises(ValueError) as err:
            ModelArgument.validate(model, index=1)
        assert str(err.value) == \
            "Enough information must be given so that both name and index can be determined."

        with pytest.raises(ValueError) as err:
            ModelArgument.validate(model)
        assert str(err.value) == \
            "Enough information must be given so that both name and index can be determined."

        with pytest.raises(ValueError) as err:
            ModelArgument.validate(name='name')
        assert str(err.value) == \
            "Enough information must be given so that both name and index can be determined."

        with pytest.raises(ValueError) as err:
            ModelArgument.validate(index=1)
        assert str(err.value) == \
            "Enough information must be given so that both name and index can be determined."

        with pytest.raises(ValueError) as err:
            ModelArgument.validate()
        assert str(err.value) == \
            "Enough information must be given so that both name and index can be determined."

    def test_get_slice(self):
        argument = ModelArgument('test', True, 7)
        value_kwarg = mk.MagicMock()
        kwargs = {'test': value_kwarg, 'other': mk.MagicMock()}

        # Test success
        kwargs = {'test': value_kwarg, 'other': mk.MagicMock()}
        assert argument.get_slice(**kwargs) == value_kwarg

        # Test fail
        kwargs = {'thing': value_kwarg, 'other': mk.MagicMock()}
        with pytest.raises(ValueError) as err:
            argument.get_slice(**kwargs)
        assert str(err.value) == \
            f"Cannot find a valid input corresponding to test in: {kwargs}."

    def test__removed_bounding_box(self):
        argument = ModelArgument('test', True, 7)
        bbox = argument._removed_bounding_box()

        assert isinstance(bbox, BoundingBox)
        assert bbox == (-np.inf, np.inf)

        for num in [0, -1, -100, 3.14, 10^100, -10^100]:
            assert not num < bbox[0]
            assert num > bbox[0]

            assert not num > bbox[1]
            assert num < bbox[1]

    def test__add_bounding_box(self):
        argument = ModelArgument('test', True, 2)

        # bounding_box is dim 1
        bounding_box = BoundingBox(((0, 1), (2, 3)))
        bbox = argument._add_bounding_box(bounding_box)
        assert isinstance(bbox, BoundingBox)
        assert bbox == ((0, 1), (2, 3), (-np.inf, np.inf))

        # bounding_box is only a tuple
        bounding_box = ((0, 1), (2, 3))
        bbox = argument._add_bounding_box(bounding_box)
        assert isinstance(bbox, BoundingBox)
        assert bbox == ((0, 1), (2, 3), (-np.inf, np.inf))

        # bounding_box is dim > 1
        base_box = [(0, 1), (2, 3)]
        true_box = [(0, 1), (2, 3), (-np.inf, np.inf)]
        for _ in range(3):
            base_box.append((0, 1))
            true_box.append((0, 1))
            bounding_box = BoundingBox(tuple(base_box))

            bbox = argument._add_bounding_box(bounding_box)
            assert isinstance(bbox, BoundingBox)
            assert bbox == tuple(true_box)

    def test_add_bounding_box(self):
        bounding_box = mk.MagicMock()
        with mk.patch.object(ModelArgument, '_add_bounding_box',
                             autospec=True) as mkAdd:
            # No Remove
            argument = ModelArgument('test', False, 7)
            assert argument.add_bounding_box(bounding_box) == bounding_box
            assert mkAdd.call_args_list == []

            # Remove
            argument = ModelArgument('test', True, 7)
            assert argument.add_bounding_box(bounding_box) == mkAdd.return_value
            assert mkAdd.call_args_list == [mk.call(argument, bounding_box)]

    def test_add_removed_axis(self):
        # Test no change, due to not removed
        argument = ModelArgument('test', False, 2)
        assert (argument.add_removed_axis(np.array([0, 1])) == np.array([0, 1])).all()

        # Test no change, due to removed and present
        argument = ModelArgument('test', True, 2)
        assert (argument.add_removed_axis(np.array([0, 1, 2])) == np.array([0, 1, 2])).all()

        # Add change
        argument = ModelArgument('test', True, 2)
        assert (argument.add_removed_axis(np.array([0, 1])) == np.array([0, 1, 2])).all()


class TestModelArguments:
    def test_properties(self):
        args = [ModelArgument(mk.MagicMock(), mk.MagicMock(), mk.MagicMock())
                for _ in range(3)]
        arguments = ModelArguments(args)

        assert arguments._arguments == args
        assert arguments.arguments == args

        assert arguments.names == [arg[0] for arg in args]
        assert arguments.indices == {arg[0]: arg[2] for arg in args}

    def test_sorted(self):
        arg0 = ModelArgument(mk.MagicMock(), mk.MagicMock(), 0)
        arg1 = ModelArgument(mk.MagicMock(), mk.MagicMock(), 1)
        arg2 = ModelArgument(mk.MagicMock(), mk.MagicMock(), 2)
        arg3 = ModelArgument(mk.MagicMock(), mk.MagicMock(), 3)

        args = [arg2, arg0, arg3, arg1]
        arguments = ModelArguments(args)

        new = arguments.sorted
        assert isinstance(new, ModelArguments)
        assert new.arguments == [arg0, arg1, arg2, arg3]

    def test_removed(self):
        arg0 = ModelArgument(mk.MagicMock(), True, 0)
        arg1 = ModelArgument(mk.MagicMock(), False, 1)
        arg2 = ModelArgument(mk.MagicMock(), True, 2)
        arg3 = ModelArgument(mk.MagicMock(), False, 3)

        args = [arg2, arg0, arg3, arg1]
        arguments = ModelArguments(args)

        assert arguments.removed == [arg2, arg0]

    def test_validate(self):
        names = [f"name{idx}" for idx in range(3)]
        model = mk.MagicMock()
        model.inputs = names

        # just names
        arguments = ModelArguments.validate(model, names)
        assert isinstance(arguments, ModelArguments)
        for idx, arg in enumerate(arguments.arguments):
            assert isinstance(arg, ModelArgument)
            assert arg == (names[idx], False, idx)

        # Names and remove
        args = [(name, True) for name in names]
        arguments = ModelArguments.validate(model, args)
        assert isinstance(arguments, ModelArguments)
        for idx, arg in enumerate(arguments.arguments):
            assert isinstance(arg, ModelArgument)
            assert arg == (names[idx], True, idx)

        # Names, remove, and index
        args = [(name, True, idx) for idx, name in enumerate(names)]
        arguments = ModelArguments.validate(model, args)
        assert isinstance(arguments, ModelArguments)
        for idx, arg in enumerate(arguments.arguments):
            assert isinstance(arg, ModelArgument)
            assert arg == (names[idx], True, idx)

        # Other ModelArguments
        args = ModelArguments([ModelArgument(name, True, idx) for idx, name in enumerate(names)])
        arguments = ModelArguments.validate(model, args)
        assert isinstance(arguments, ModelArguments)
        for idx, arg in enumerate(arguments.arguments):
            assert isinstance(arg, ModelArgument)
            assert arg == (names[idx], True, idx)
        assert id(args) != id(arguments)

        # No arguments
        arguments = ModelArguments.validate(model)
        assert isinstance(arguments, ModelArguments)
        assert arguments.arguments == []

        # Invalid
        args = [(name, True, idx + 1) for idx, name in enumerate(names)]
        with pytest.raises(IndexError):
            arguments = ModelArguments.validate(model, args)

    def test__eq__(self):
        args = [ModelArgument(f"name{idx}", idx, mk.MagicMock())
                for idx in range(3)]
        arguments0 = ModelArguments(args.copy())
        arguments1 = ModelArguments(args.copy())

        # Equal
        assert arguments0 == arguments1
        assert id(arguments0) != id(arguments1)

        # Not equal do to argument differences
        arguments1._arguments.append(mk.MagicMock())
        assert not arguments0 == arguments1

        # Not equal do to different type
        assert not arguments0 == mk.MagicMock()

    def test_get_slice(self):
        args = [ModelArgument(f"name{idx}", mk.MagicMock(), idx)
                for idx in range(3)]
        kwargs = {arg.name: mk.MagicMock() for arg in args}

        # Full tuple
        arguments = ModelArguments(args)
        assert arguments.get_slice(**kwargs) == tuple(kwargs.values())

        # Single argument
        arguments = ModelArguments([args[1]])
        assert arguments.get_slice(**kwargs) == kwargs[args[1].name]

    def test_add_bounding_box(self):
        args = [ModelArgument(f"name{idx}", idx, mk.MagicMock())
                for idx in range(3)]
        arguments = ModelArguments(args)
        bounding_box = mk.MagicMock()

        new_bbox = [mk.MagicMock() for _ in range(3)]
        with mk.patch.object(ModelArgument, 'add_bounding_box',
                             autospec=True, side_effect=new_bbox) as mkAdd:
            assert arguments.add_bounding_box(bounding_box) == new_bbox[2]
            assert mkAdd.call_args_list == [
                mk.call(args[0], bounding_box),
                mk.call(args[1], new_bbox[0]),
                mk.call(args[2], new_bbox[1]),
            ]

    def test_add_removed_axes(self):
        args = [ModelArgument(f"name{idx}", idx, mk.MagicMock())
                for idx in range(3)]
        arguments = ModelArguments(args)
        axes_ind = mk.MagicMock()

        new_axes_ind = [mk.MagicMock() for _ in range(3)]
        with mk.patch.object(ModelArgument, 'add_removed_axis',
                             autospec=True, side_effect=new_axes_ind) as mkAdd:
            assert arguments.add_removed_axes(axes_ind) == new_axes_ind[2]
            assert mkAdd.call_args_list == [
                mk.call(args[0], axes_ind),
                mk.call(args[1], new_axes_ind[0]),
                mk.call(args[2], new_axes_ind[1]),
            ]


class TestCompoundBoundingBox:
    def test___init__(self):
        create = mk.MagicMock()
        bbox = {1: (-1, 0), 2: (0, 1)}

        bounding_box = CompoundBoundingBox(bbox)
        assert bounding_box == bbox
        assert bounding_box._model is None
        assert bounding_box._slice_args == ModelArguments([])
        assert bounding_box._create_slice is None

        bounding_box = CompoundBoundingBox(bbox, Gaussian1D(), 'x', create)
        assert bounding_box == bbox
        assert (bounding_box._model.parameters == Gaussian1D().parameters).all()
        assert bounding_box._slice_args == ModelArguments([ModelArgument('x', False, 0)])
        assert bounding_box._create_slice == create

        bounding_box = CompoundBoundingBox(bbox, Gaussian1D(), [('x', True)], create)
        assert bounding_box == bbox
        assert (bounding_box._model.parameters == Gaussian1D().parameters).all()
        assert bounding_box._slice_args == ModelArguments([ModelArgument('x', True, 0)])
        assert bounding_box._create_slice == create

    def test_slice_args(self):
        bbox = {1: (-1, 0), 2: (0, 1)}
        bounding_box = CompoundBoundingBox(bbox, slice_args=[('x', True, 0)])

        assert bounding_box._slice_args == ModelArguments([ModelArgument('x', True, 0)])
        assert bounding_box.slice_args == ModelArguments([ModelArgument('x', True, 0)])

    def test_slice_names(self):
        bbox = {1: (-1, 0), 2: (0, 1)}
        bounding_box = CompoundBoundingBox(bbox, Gaussian2D(), slice_args=['x', 'y'])

        assert bounding_box.slice_names == ['x', 'y']

    def test_slice_indicies(self):
        bbox = {1: (-1, 0), 2: (0, 1)}
        bounding_box = CompoundBoundingBox(bbox, Gaussian2D(), slice_args=['x', 'y'])

        assert bounding_box.slice_indicies == {'x': 0, 'y': 1}

    def test_validate(self):
        # Main validate
        bbox = {1: (-1, 0), 2: (0, 1)}
        bounding_box = CompoundBoundingBox.validate(bbox, slice_args=['x', ('y', True)], model=Gaussian2D())
        assert bounding_box == bbox
        for bbox_slice in bounding_box.values():
            assert isinstance(bbox_slice, BoundingBox)
        assert (bounding_box._model.parameters == Gaussian2D().parameters).all()
        assert bounding_box._slice_args == ModelArguments([ModelArgument('x', False, 0),
                                                           ModelArgument('y', True, 1)])

        # Re validate
        new_bounding_box = CompoundBoundingBox.validate(bounding_box, model=Gaussian2D())
        assert new_bounding_box == bounding_box
        assert new_bounding_box._slice_args == bounding_box._slice_args

        # Simple validate
        bbox = {1: (-1, 0), 2: (0, 1)}
        bounding_box = CompoundBoundingBox.validate(bbox, slice_args=[('y', True)], model=Gaussian2D())
        assert bounding_box == bbox
        for bbox_slice in bounding_box.values():
            assert isinstance(bbox_slice, BoundingBox)
        assert (bounding_box._model.parameters == Gaussian2D().parameters).all()
        assert bounding_box._slice_args == ModelArguments([ModelArgument('y', True, 1)])

    def test__get_slice(self):
        bbox = {1: (-1, 0), 2: (0, 1)}
        bounding_box = CompoundBoundingBox(bbox, Gaussian2D(), slice_args=['x', 'y'])

        # Normal Get
        assert bounding_box._get_slice(1) == (-1, 0)
        assert bounding_box._get_slice(2) == (0, 1)

        # Error
        with pytest.raises(RuntimeError) as err:
            bounding_box._get_slice(3)
        assert str(err.value) == \
            "No bounding_box is defined for slice: 3!"

        # New box
        model = mk.MagicMock()
        create = mk.MagicMock()
        bounding_box._model = model
        bounding_box._create_slice = create
        assert 3 not in bounding_box
        assert bounding_box._get_slice(3) == create.return_value
        assert 3 in bounding_box
        assert bounding_box[3] == create.return_value

        assert create.call_args_list == [mk.call(model, 3)]

    def test__get_bounding_box(self):
        bbox = {1: (-1, 0), 2: (0, 1)}
        bounding_box = CompoundBoundingBox(bbox, Gaussian2D(), slice_args=['y'])

        assert bounding_box._get_bounding_box(y=1) == (-1, 0)
        assert bounding_box._get_bounding_box(y=2) == (0, 1)

    def test__add_bounding_box(self):
        origin_bbox = BoundingBox((1, 2))
        bbox = {1: (-1, 0), 2: (0, 1)}

        bounding_box = CompoundBoundingBox(bbox, Gaussian2D(), slice_args=['x', 'y'])
        assert bounding_box._add_bounding_box(origin_bbox) == (1, 2)

        bounding_box = CompoundBoundingBox(bbox, Gaussian2D(), slice_args=['x', ('y', True)])
        assert bounding_box._add_bounding_box(origin_bbox) == ((1, 2), (-np.inf, np.inf))

        bounding_box = CompoundBoundingBox(bbox, Gaussian2D(), slice_args=[('x', True), ('y', True)])
        assert bounding_box._add_bounding_box(origin_bbox) == ((-np.inf, np.inf), (-np.inf, np.inf), (1, 2))

        bounding_box = CompoundBoundingBox(bbox, Gaussian2D(), slice_args=['y', 'x'])
        assert bounding_box._add_bounding_box(origin_bbox) == (1, 2)

        bounding_box = CompoundBoundingBox(bbox, Gaussian2D(), slice_args=['y', ('x', True)])
        assert bounding_box._add_bounding_box(origin_bbox) == ((-np.inf, np.inf), (1, 2))

        bounding_box = CompoundBoundingBox(bbox, Gaussian2D(), slice_args=[('y', True), ('x', True)])
        assert bounding_box._add_bounding_box(origin_bbox) == ((-np.inf, np.inf), (-np.inf, np.inf), (1, 2))

    def test_get_bounding_box(self):
        bbox = {1: (-1, 0), 2: (0, 1)}
        bounding_box = CompoundBoundingBox(bbox, Gaussian2D(), slice_args=['y'])
        kwargs = {'test': mk.MagicMock()}

        with mk.patch.object(CompoundBoundingBox, '_add_bounding_box',
                             autospec=True) as mkAdd:
            with mk.patch.object(CompoundBoundingBox, '_get_bounding_box',
                                 autospec=True) as mkGet:
                assert bounding_box.get_bounding_box(**kwargs) == mkAdd.return_value
                assert mkAdd.call_args_list == [mk.call(bounding_box, mkGet.return_value)]
                assert mkGet.call_args_list == [mk.call(bounding_box, **kwargs)]

    def test_add_removed_axes(self):
        bbox = {1: (-1, 0), 2: (0, 1)}
        bounding_box = CompoundBoundingBox(bbox, Gaussian2D(), slice_args=['y'])
        axes_ind = mk.MagicMock()

        with mk.patch.object(ModelArguments, 'add_removed_axes',
                             autospec=True) as mkAdd:
            with mk.patch.object(np, 'argsort', autospec=True) as mkSort:
                assert bounding_box.add_removed_axes(axes_ind) == mkSort.return_value
                assert mkSort.call_args_list == [mk.call(mkAdd.return_value)]
                assert mkAdd.call_args_list == [mk.call(bounding_box.slice_args, axes_ind)]


# def test_CompoundBoundingBox__init__():
#     bbox = {1: (-1, 0), 2: (0, 1)}
#     bounding_box = CompoundBoundingBox(bbox)

#     assert bounding_box == bbox
#     assert bounding_box._model is None
#     assert bounding_box._slice_arg is None

#     bounding_box = CompoundBoundingBox(bbox, Gaussian1D(), 'x')
#     assert bounding_box == bbox
#     assert (bounding_box._model.parameters == Gaussian1D().parameters).all()
#     assert bounding_box._slice_arg == 0


# def test_CompoundBoundingBox__get_arg_index():
#     bounding_box = CompoundBoundingBox({}, Gaussian2D())

#     assert bounding_box._get_arg_index(0) == 0
#     assert bounding_box._get_arg_index(1) == 1
#     with pytest.raises(ValueError):
#         bounding_box._get_arg_index(2)

#     assert bounding_box._get_arg_index('x') == 0
#     assert bounding_box._get_arg_index('y') == 1
#     with pytest.raises(ValueError):
#         bounding_box._get_arg_index('z')


# def test_CompoundBoundingBox_validate():
#     bbox = {1: (-1, 0), 2: (0, 1)}
#     model = Gaussian1D()
#     bounding_box = CompoundBoundingBox.validate(model, bbox, 'x')

#     assert bounding_box == bbox
#     assert bounding_box._model == model
#     assert bounding_box._slice_arg == 0
#     for slice_box in bounding_box.values():
#         assert isinstance(slice_box, BoundingBox)

#     model = Gaussian2D()
#     bounding_box = CompoundBoundingBox.validate(model, bbox, 'x',
#                                                remove_slice_arg=True)
#     assert bounding_box == bbox
#     assert bounding_box._model == model
#     assert bounding_box._slice_arg == 0
#     for slice_box in bounding_box.values():
#         assert isinstance(slice_box, BoundingBox)


# def test_CompoundBoundingBox_set_slice_arg():
#     bounding_box = CompoundBoundingBox((), slice_arg='arg')
#     assert bounding_box._slice_arg == 'arg'

#     bounding_box.set_slice_arg(None)
#     assert bounding_box._slice_arg is None

#     bounding_box._model = Gaussian1D()
#     with pytest.raises(ValueError):
#         bounding_box.set_slice_arg('arg')

#     with pytest.raises(ValueError):
#         bounding_box.set_slice_arg(('x', 'y'))

#     bounding_box.set_slice_arg('x')
#     assert bounding_box._slice_arg == 0
#     bounding_box.set_slice_arg(0)
#     assert bounding_box._slice_arg == 0

#     bounding_box._model = Gaussian2D()
#     bounding_box.set_slice_arg(('x', 'y'))
#     assert bounding_box._slice_arg == (0, 1)
#     bounding_box.set_slice_arg((0, 1))
#     assert bounding_box._slice_arg == (0, 1)


# def test_CompoundBoundingBox__get_slice_index():
#     bounding_box = CompoundBoundingBox({}, Gaussian2D())

#     inputs = [mk.MagicMock(), mk.MagicMock(), mk.MagicMock()]
#     assert bounding_box._get_slice_index(inputs, 'x') == inputs[0]
#     assert bounding_box._get_slice_index(inputs, 'y') == inputs[1]
#     with pytest.raises(RuntimeError):
#         bounding_box._get_slice_index(None, 'x')

#     inputs = [np.array(1), np.array(2), np.array(3)]
#     assert bounding_box._get_slice_index(inputs, 'x') == 1
#     assert bounding_box._get_slice_index(inputs, 'y') == 2
#     with pytest.raises(RuntimeError):
#         bounding_box._get_slice_index(None, 'x')


# def test_CompoundBoundingBox_get_bounding_box():
#     inputs = [np.array(1), np.array(2), np.array(3)]

#     bounding_box = CompoundBoundingBox({}, Gaussian2D())
#     assert bounding_box._slice_arg is None
#     assert bounding_box.get_bounding_box(inputs) is None
#     with pytest.raises(RuntimeError):
#         bounding_box.get_bounding_box(inputs, slice_index=mk.MagicMock())

#     bbox = {(1, 2): mk.MagicMock(), (3, 4): mk.MagicMock()}
#     bounding_box = CompoundBoundingBox(bbox, Gaussian2D(), ('x', 'y'))
#     assert bounding_box.get_bounding_box(inputs) == bbox[(1, 2)]
#     assert bounding_box.get_bounding_box(inputs, slice_index=(3, 4)) == bbox[(3, 4)]
#     with pytest.raises(RuntimeError):
#         bounding_box.get_bounding_box([np.array(4), np.array(5)])
#     bounding_box = CompoundBoundingBox(bbox, Gaussian2D(), (0, 1))

#     assert bounding_box.get_bounding_box(inputs) == bbox[(1, 2)]
#     assert bounding_box.get_bounding_box(inputs, slice_index=(3, 4)) == bbox[(3, 4)]
#     with pytest.raises(RuntimeError):
#         bounding_box.get_bounding_box([np.array(4), np.array(5)])

#     bbox = {1: mk.MagicMock(), 2: mk.MagicMock()}
#     bounding_box = CompoundBoundingBox(bbox, Gaussian2D(), 'x')
#     assert bounding_box.get_bounding_box(inputs) == bbox[1]
#     assert bounding_box.get_bounding_box(inputs, slice_index=2) == bbox[2]
#     with pytest.raises(RuntimeError):
#         bounding_box.get_bounding_box([np.array(3)])

#     bounding_box = CompoundBoundingBox(bbox, Gaussian2D(), 0)
#     assert bounding_box.get_bounding_box(inputs) == bbox[1]
#     assert bounding_box.get_bounding_box(inputs, slice_index=2) == bbox[2]
#     with pytest.raises(RuntimeError):
#         bounding_box.get_bounding_box([np.array(3)])
