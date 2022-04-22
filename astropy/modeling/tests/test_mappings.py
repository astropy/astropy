# Licensed under a 3-clause BSD style license - see LICENSE.rst
# pylint: disable=invalid-name
import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from astropy import units as u
from astropy.modeling.fitting import DogBoxLSQFitter, LevMarLSQFitter, LMLSQFitter, TRFLSQFitter
from astropy.modeling.models import Gaussian1D, Identity, Mapping, Rotation2D, Shift, UnitsMapping
from astropy.utils import NumpyRNGContext
from astropy.utils.compat.optional_deps import HAS_SCIPY  # noqa: F401

fitters = [LevMarLSQFitter, TRFLSQFitter, LMLSQFitter, DogBoxLSQFitter]


def test_swap_axes():
    x = np.zeros((2, 3))
    y = np.ones((2, 3))

    mapping = Mapping((1, 0))
    assert mapping(1, 2) == (2.0, 1.0)
    assert mapping.inverse(2, 1) == (1, 2)
    assert_array_equal(mapping(x, y), (y, x))
    assert_array_equal(mapping.inverse(y, x), (x, y))


def test_duplicate_axes():
    mapping = Mapping((0, 1, 0, 1))
    assert mapping(1, 2) == (1.0, 2., 1., 2)
    assert mapping.inverse(1, 2, 1, 2) == (1, 2)
    assert mapping.inverse.n_inputs == 4
    assert mapping.inverse.n_outputs == 2


def test_drop_axes_1():
    mapping = Mapping((0,), n_inputs=2)
    assert mapping(1, 2) == (1.)


def test_drop_axes_2():
    mapping = Mapping((1, ))
    assert mapping(1, 2) == (2.)
    with pytest.raises(NotImplementedError):
        mapping.inverse


def test_drop_axes_3():
    mapping = Mapping((1,), n_inputs=2)
    assert mapping.n_inputs == 2
    rotation = Rotation2D(60)
    model = rotation | mapping
    assert_allclose(model(1, 2), 1.86602540378)


@pytest.mark.parametrize('name', [None, 'test_name'])
def test_bad_inputs(name):
    mapping = Mapping((1, 0), name=name)

    if name is None:
        name = "Mapping"

    x = [np.ones((2, 3))*idx for idx in range(5)]
    for idx in range(1, 6):
        if idx == 2:
            continue

        with pytest.raises(TypeError) as err:
            mapping.evaluate(*x[:idx])
        assert str(err.value) == f"{name} expects 2 inputs; got {idx}"


def test_identity():
    x = np.zeros((2, 3))
    y = np.ones((2, 3))

    ident1 = Identity(1)
    shift = Shift(1)
    rotation = Rotation2D(angle=60)
    model = ident1 & shift | rotation
    assert_allclose(model(1, 2), (-2.098076211353316, 2.3660254037844393))
    res_x, res_y = model(x, y)
    assert_allclose((res_x, res_y),
                    (np.array([[-1.73205081, -1.73205081, -1.73205081],
                               [-1.73205081, -1.73205081, -1.73205081]]),
                     np.array([[1., 1., 1.],
                               [1., 1., 1.]])))
    assert_allclose(model.inverse(res_x, res_y), (x, y), atol=1.e-10)


# https://github.com/astropy/astropy/pull/6018
@pytest.mark.skipif('not HAS_SCIPY')
@pytest.mark.parametrize('fitter', fitters)
def test_fittable_compound(fitter):
    fitter = fitter()

    m = Identity(1) | Mapping((0, )) | Gaussian1D(1, 5, 4)
    x = np.arange(10)
    y_real = m(x)
    dy = 0.005
    with NumpyRNGContext(1234567):
        n = np.random.normal(0., dy, x.shape)
    y_noisy = y_real + n
    new_model = fitter(m, x, y_noisy)
    y_fit = new_model(x)
    assert_allclose(y_fit, y_real, atol=dy)


def test_identity_repr():
    m = Identity(1, name='foo')
    assert repr(m) == "<Identity(1, name='foo')>"

    m = Identity(1)
    assert repr(m) == "<Identity(1)>"


def test_mapping_repr():
    m = Mapping([0, 1], name='foo')
    assert repr(m) == "<Mapping([0, 1], name='foo')>"

    m = Mapping([0, 1])
    assert repr(m) == "<Mapping([0, 1])>"


class TestUnitsMapping:
    def test___init__(self):
        # Set values
        model = UnitsMapping(((u.m, None),),
                             input_units_equivalencies='test_eqiv',
                             input_units_allow_dimensionless=True,
                             name='test')
        assert model._mapping == ((u.m, None),)
        assert model._input_units_strict == {'x': True}
        assert model.input_units_equivalencies == 'test_eqiv'
        assert model.input_units_allow_dimensionless == {'x': True}
        assert model.name == 'test'
        assert model._input_units == {'x': u.m}

        # Default values
        model = UnitsMapping(((u.K, None),))
        assert model._mapping == ((u.K, None),)
        assert model._input_units_strict == {'x': True}
        assert model.input_units_equivalencies is None
        assert model.input_units_allow_dimensionless == {'x': False}
        assert model.name is None
        assert model._input_units == {'x': u.K}

        # Error
        with pytest.raises(ValueError) as err:
            UnitsMapping(((u.m, None), (u.m, u.K)))
        assert str(err.value) == "If one return unit is None, then all must be None"

    def test_evaluate(self):
        model = UnitsMapping(((u.m, None),))
        assert model(10*u.m) == 10

        model = UnitsMapping(((u.m, u.K),))
        assert model(10*u.m) == 10 * u.K

        model = UnitsMapping(((u.m, None), (u.K, None)),)
        assert model(10*u.m, 20*u.K) == (10, 20)

        model = UnitsMapping(((u.m, u.K), (u.K, u.m)),)
        assert model(10*u.m, 20*u.K) == (10*u.K, 20*u.m)

    def test_repr(self):
        model = UnitsMapping(((u.m, None),), name='foo')
        assert repr(model) == f"<UnitsMapping((({repr(u.m)}, None),), name='foo')>"

        model = UnitsMapping(((u.m, None),))
        assert repr(model) == f"<UnitsMapping((({repr(u.m)}, None),))>"
