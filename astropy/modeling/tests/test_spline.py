# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Tests for spline models and fitters"""
import unittest.mock as mk

import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy.modeling.core import FittableModel, ModelDefinitionError
from astropy.modeling.fitting import (
    SplineExactKnotsFitter, SplineInterpolateFitter, SplineSmoothingFitter, SplineSplrepFitter)
from astropy.modeling.parameters import Parameter
from astropy.modeling.spline import Spline1D, _Spline, _SplineFitter
from astropy.utils.compat.optional_deps import HAS_SCIPY  # noqa: F401
# pylint: disable=invalid-name
from astropy.utils.exceptions import AstropyUserWarning

npts = 50
nknots = 10
np.random.seed(42)
test_w = np.random.rand(npts)
test_t = [-1, 0, 1]
noise = np.random.randn(npts)

degree_tests = [1, 2, 3, 4, 5]
wieght_tests = [None, test_w]
smoothing_tests = [None, 0.01]


class TestSpline:
    def setup_class(self):
        self.num_opt = 3
        self.optional_inputs = {f'test{i}': mk.MagicMock() for i in range(self.num_opt)}
        self.extra_kwargs = {f'new{i}': mk.MagicMock() for i in range(self.num_opt)}

        class Spline(_Spline):
            optional_inputs = {'test': 'test'}

            def _init_parameters(self):
                super()._init_parameters()

            def _init_data(self, knots, coeffs, bounds=None):
                super()._init_data(knots, coeffs, bounds=bounds)

        self.Spline = Spline

    def test___init__(self):
        # empty spline
        spl = self.Spline()
        assert spl._t is None
        assert spl._c is None
        assert spl._user_knots is False
        assert spl._degree is None
        assert spl._test is None

        assert not hasattr(spl, 'degree')

        # Call _init_spline
        with mk.patch.object(_Spline, '_init_spline',
                             autospec=True) as mkInit:
            # No call (knots=None)
            spl = self.Spline()
            assert mkInit.call_args_list == []

            knots = mk.MagicMock()
            coeffs = mk.MagicMock()
            bounds = mk.MagicMock()
            spl = self.Spline(knots=knots, coeffs=coeffs, bounds=bounds)
            assert mkInit.call_args_list == [mk.call(spl, knots, coeffs, bounds)]

            assert spl._t is None
            assert spl._c is None
            assert spl._user_knots is False
            assert spl._degree is None
            assert spl._test is None

        # Coeffs but no knots
        with pytest.raises(ValueError) as err:
            self.Spline(coeffs=mk.MagicMock())
        assert str(err.value) == "If one passes a coeffs vector one needs to also pass knots!"

    def test_param_names(self):
        # no parameters
        spl = self.Spline()
        assert spl.param_names == ()

        knot_names = tuple([mk.MagicMock() for _ in range(3)])
        spl._knot_names = knot_names
        assert spl.param_names == knot_names

        coeff_names = tuple([mk.MagicMock() for _ in range(3)])
        spl._coeff_names = coeff_names
        assert spl.param_names == knot_names + coeff_names

    def test__optional_arg(self):

        spl = self.Spline()
        assert spl._optional_arg('test') == '_test'

    def test__create_optional_inputs(self):
        class Spline(self.Spline):
            optional_inputs = self.optional_inputs

            def __init__(self):
                self._create_optional_inputs()

        spl = Spline()
        for arg in self.optional_inputs:
            attribute = spl._optional_arg(arg)
            assert hasattr(spl, attribute)
            assert getattr(spl, attribute) is None

        with pytest.raises(ValueError,
                           match=r"Optional argument .* already exists in this class!"):
            spl._create_optional_inputs()

    def test__intercept_optional_inputs(self):
        class Spline(self.Spline):
            optional_inputs = self.optional_inputs

            def __init__(self):
                self._create_optional_inputs()

        spl = Spline()
        new_kwargs = spl._intercept_optional_inputs(**self.extra_kwargs)
        for arg, value in self.optional_inputs.items():
            attribute = spl._optional_arg(arg)
            assert getattr(spl, attribute) is None
        assert new_kwargs == self.extra_kwargs

        kwargs = self.extra_kwargs.copy()
        for arg in self.optional_inputs:
            kwargs[arg] = mk.MagicMock()
        new_kwargs = spl._intercept_optional_inputs(**kwargs)
        for arg, value in self.optional_inputs.items():
            attribute = spl._optional_arg(arg)
            assert getattr(spl, attribute) is not None
            assert getattr(spl, attribute) == kwargs[arg]
            assert getattr(spl, attribute) != value
            assert arg not in new_kwargs
        assert new_kwargs == self.extra_kwargs
        assert kwargs != self.extra_kwargs

        with pytest.raises(RuntimeError,
                           match=r".* has already been set, something has gone wrong!"):
            spl._intercept_optional_inputs(**kwargs)

    def test_evaluate(self):
        class Spline(self.Spline):
            optional_inputs = self.optional_inputs

        spl = Spline()

        # No options passed in and No options set
        new_kwargs = spl.evaluate(**self.extra_kwargs)
        for arg, value in self.optional_inputs.items():
            assert new_kwargs[arg] == value
        for arg, value in self.extra_kwargs.items():
            assert new_kwargs[arg] == value
        assert len(new_kwargs) == (len(self.optional_inputs) + len(self.extra_kwargs))

        # No options passed in and Options set
        kwargs = self.extra_kwargs.copy()
        for arg in self.optional_inputs:
            kwargs[arg] = mk.MagicMock()
        spl._intercept_optional_inputs(**kwargs)
        new_kwargs = spl.evaluate(**self.extra_kwargs)
        assert new_kwargs == kwargs
        for arg in self.optional_inputs:
            attribute = spl._optional_arg(arg)
            assert getattr(spl, attribute) is None

        # Options passed in
        set_kwargs = self.extra_kwargs.copy()
        for arg in self.optional_inputs:
            kwargs[arg] = mk.MagicMock()
        spl._intercept_optional_inputs(**set_kwargs)
        kwargs = self.extra_kwargs.copy()
        for arg in self.optional_inputs:
            kwargs[arg] = mk.MagicMock()
        assert set_kwargs != kwargs
        new_kwargs = spl.evaluate(**kwargs)
        assert new_kwargs == kwargs

    def test___call__(self):
        spl = self.Spline()

        args = tuple([mk.MagicMock() for _ in range(3)])
        kwargs = {f"test{idx}": mk.MagicMock() for idx in range(3)}
        new_kwargs = {f"new_test{idx}": mk.MagicMock() for idx in range(3)}
        with mk.patch.object(_Spline, "_intercept_optional_inputs",
                             autospec=True, return_value=new_kwargs) as mkIntercept:
            with mk.patch.object(FittableModel, "__call__",
                                 autospec=True) as mkCall:
                assert mkCall.return_value == spl(*args, **kwargs)
                assert mkCall.call_args_list == [mk.call(spl, *args, **new_kwargs)]
                assert mkIntercept.call_args_list == [mk.call(spl, **kwargs)]

    def test__create_parameter(self):
        np.random.seed(37)
        base_vec = np.random.random(20)
        test = base_vec.copy()
        fixed_test = base_vec.copy()

        class Spline(self.Spline):
            @property
            def test(self):
                return test

            @property
            def fixed_test(self):
                return fixed_test

        spl = Spline()
        assert (spl.test == test).all()
        assert (spl.fixed_test == fixed_test).all()

        for index in range(20):
            name = f"test_name{index}"
            spl._create_parameter(name, index, 'test')
            assert hasattr(spl, name)
            param = getattr(spl, name)
            assert isinstance(param, Parameter)
            assert param.model == spl
            assert param.fixed is False
            assert param.value == test[index] == spl.test[index] == base_vec[index]
            new_set = np.random.random()
            param.value = new_set
            assert spl.test[index] == new_set
            assert spl.test[index] != base_vec[index]
            new_get = np.random.random()
            spl.test[index] = new_get
            assert param.value == new_get
            assert param.value != new_set

        for index in range(20):
            name = f"fixed_test_name{index}"
            spl._create_parameter(name, index, 'fixed_test', True)
            assert hasattr(spl, name)
            param = getattr(spl, name)
            assert isinstance(param, Parameter)
            assert param.model == spl
            assert param.fixed is True
            assert param.value == fixed_test[index] == spl.fixed_test[index] == base_vec[index]
            new_set = np.random.random()
            param.value = new_set
            assert spl.fixed_test[index] == new_set
            assert spl.fixed_test[index] != base_vec[index]
            new_get = np.random.random()
            spl.fixed_test[index] = new_get
            assert param.value == new_get
            assert param.value != new_set

    def test__create_parameters(self):
        np.random.seed(37)
        test = np.random.random(20)

        class Spline(self.Spline):
            @property
            def test(self):
                return test

        spl = Spline()

        fixed = mk.MagicMock()
        with mk.patch.object(_Spline, '_create_parameter',
                             autospec=True) as mkCreate:
            params = spl._create_parameters("test_param", "test", fixed)
            assert params == tuple([f"test_param{idx}" for idx in range(20)])
            assert mkCreate.call_args_list == [
                mk.call(spl, f"test_param{idx}", idx, 'test', fixed) for idx in range(20)
            ]

    def test__init_parameters(self):
        spl = self.Spline()

        with pytest.raises(NotImplementedError) as err:
            spl._init_parameters()
        assert str(err.value) == "This needs to be implemented"

    def test__init_data(self):
        spl = self.Spline()

        with pytest.raises(NotImplementedError) as err:
            spl._init_data(mk.MagicMock(), mk.MagicMock(), mk.MagicMock())
        assert str(err.value) == "This needs to be implemented"

        with pytest.raises(NotImplementedError) as err:
            spl._init_data(mk.MagicMock(), mk.MagicMock())
        assert str(err.value) == "This needs to be implemented"

    def test__init_spline(self):
        spl = self.Spline()

        knots = mk.MagicMock()
        coeffs = mk.MagicMock()
        bounds = mk.MagicMock()
        with mk.patch.object(_Spline, "_init_parameters",
                             autospec=True) as mkParameters:
            with mk.patch.object(_Spline, "_init_data",
                                 autospec=True) as mkData:
                main = mk.MagicMock()
                main.attach_mock(mkParameters, 'parameters')
                main.attach_mock(mkData, 'data')

                spl._init_spline(knots, coeffs, bounds)
                assert main.mock_calls == [
                    mk.call.data(spl, knots, coeffs, bounds=bounds),
                    mk.call.parameters(spl)
                ]

    def test__init_tck(self):
        spl = self.Spline()
        assert spl._c is None
        assert spl._t is None
        assert spl._degree is None

        spl = self.Spline(degree=4)
        assert spl._c is None
        assert spl._t is None
        assert spl._degree == 4


@pytest.mark.skipif('not HAS_SCIPY')
class TestSpline1D:
    def setup_class(self):
        def func(x, noise=0):
            return np.exp(-x**2) + 0.1*noise

        self.x = np.linspace(-3, 3, npts)
        self.y = func(self.x, noise)
        self.truth = func(self.x)

        arg_sort = np.argsort(self.x)
        np.random.shuffle(arg_sort)

        self.x_s = self.x[arg_sort]
        self.y_s = func(self.x_s, noise[arg_sort])

        self.npts_out = 1000
        self.xs = np.linspace(-3, 3, self.npts_out)

        self.t = np.linspace(-3, 3, nknots)[1:-1]

    def check_parameter(self, spl, base_name, name, index, value, fixed):
        assert base_name in name
        assert index == int(name.split(base_name)[-1])
        knot_name = f"{base_name}{index}"
        assert knot_name == name
        assert hasattr(spl, name)
        param = getattr(spl, name)
        assert isinstance(param, Parameter)
        assert param.name == name
        assert param.value == value(index)
        assert param.model == spl
        assert param.fixed is fixed

    def check_parameters(self, spl, params, base_name, value, fixed):
        for idx, name in enumerate(params):
            self.check_parameter(spl, base_name, name, idx, value, fixed)

    def update_parameters(self, spl, knots, value):
        for name in knots:
            param = getattr(spl, name)
            param.value = value
            assert param.value == value

    def test___init__with_no_knot_information(self):
        spl = Spline1D()
        assert spl._degree == 3
        assert spl._user_knots is False
        assert spl._t is None
        assert spl._c is None
        assert spl._nu is None

        # Check no parameters created
        assert len(spl._knot_names) == 0
        assert len(spl._coeff_names) == 0

    def test___init__with_number_of_knots(self):
        spl = Spline1D(knots=10)

        # Check baseline data
        assert spl._degree == 3
        assert spl._user_knots is False
        assert spl._nu is None

        # Check vector data
        assert len(spl._t) == 18
        t = np.zeros(18)
        t[-4:] = 1
        assert (spl._t == t).all()
        assert len(spl._c) == 18
        assert (spl._c == np.zeros(18)).all()

        # Check all parameter names created:
        assert len(spl._knot_names) == 18
        assert len(spl._coeff_names) == 18

        # Check knot values:
        def value0(idx):
            if idx < 18 - 4:
                return 0
            else:
                return 1
        self.check_parameters(spl, spl._knot_names, "knot", value0, True)

        # Check coeff values:
        def value1(idx):
            return 0
        self.check_parameters(spl, spl._coeff_names, "coeff", value1, False)

    def test___init__with_full_custom_knots(self):
        t = 17*np.arange(20) - 32
        spl = Spline1D(knots=t)

        # Check baseline data
        assert spl._degree == 3
        assert spl._user_knots is True
        assert spl._nu is None

        # Check vector data
        assert (spl._t == t).all()
        assert len(spl._c) == 20
        assert (spl._c == np.zeros(20)).all()

        # Check all parameter names created
        assert len(spl._knot_names) == 20
        assert len(spl._coeff_names) == 20

        # Check knot values:
        def value0(idx):
            return t[idx]
        self.check_parameters(spl, spl._knot_names, "knot", value0, True)

        # Check coeff values
        def value1(idx):
            return 0
        self.check_parameters(spl, spl._coeff_names, "coeff", value1, False)

    def test___init__with_interior_custom_knots(self):
        t = np.arange(1, 20)
        spl = Spline1D(knots=t, bounds=[0, 20])
        # Check baseline data
        assert spl._degree == 3
        assert spl._user_knots is True
        assert spl._nu is None

        # Check vector data
        assert len(spl._t) == 27
        assert (spl._t[4:-4] == t).all()
        assert (spl._t[:4] == 0).all()
        assert (spl._t[-4:] == 20).all()

        assert len(spl._c) == 27
        assert (spl._c == np.zeros(27)).all()

        # Check knot values:
        def value0(idx):
            if idx < 4:
                return 0
            elif idx >= 19 + 4:
                return 20
            else:
                return t[idx-4]
        self.check_parameters(spl, spl._knot_names, "knot", value0, True)

        # Check coeff values
        def value1(idx):
            return 0
        self.check_parameters(spl, spl._coeff_names, "coeff", value1, False)

    def test___init__with_user_knots_and_coefficients(self):
        t = 17*np.arange(20) - 32
        c = np.linspace(-1, 1, 20)
        spl = Spline1D(knots=t, coeffs=c)

        # Check baseline data
        assert spl._degree == 3
        assert spl._user_knots is True
        assert spl._nu is None

        # Check vector data
        assert (spl._t == t).all()
        assert len(spl._c) == 20
        assert (spl._c == c).all()

        # Check all parameter names created
        assert len(spl._knot_names) == 20
        assert len(spl._coeff_names) == 20

        # Check knot values:
        def value0(idx):
            return t[idx]
        self.check_parameters(spl, spl._knot_names, "knot", value0, True)

        # Check coeff values
        def value1(idx):
            return c[idx]
        self.check_parameters(spl, spl._coeff_names, "coeff", value1, False)

    def test___init__errors(self):
        # Bad knot type
        knots = 3.5
        with pytest.raises(ValueError) as err:
            Spline1D(knots=knots)
        assert str(err.value) == f"Knots: {knots} must be iterable or value"

        # Not enough knots
        for idx in range(8):
            with pytest.raises(ValueError) as err:
                Spline1D(knots=np.arange(idx))
            assert str(err.value) == "Must have at least 8 knots."

        # Bad scipy spline
        t = np.arange(20)[::-1]
        with pytest.raises(ValueError):
            Spline1D(knots=t)

    def test_parameter_array_link(self):
        spl = Spline1D(10)

        # Check knot base values
        def value0(idx):
            if idx < 18 - 4:
                return 0
            else:
                return 1
        self.check_parameters(spl, spl._knot_names, "knot", value0, True)

        # Check knot vector -> knot parameter link
        t = np.arange(18)
        spl._t = t.copy()

        def value1(idx):
            return t[idx]
        self.check_parameters(spl, spl._knot_names, "knot", value1, True)

        # Check knot parameter -> knot vector link
        self.update_parameters(spl, spl._knot_names, 3)
        assert (spl._t[:] == 3).all()

        # Check coeff base values
        def value2(idx):
            return 0
        self.check_parameters(spl, spl._coeff_names, "coeff", value2, False)

        # Check coeff vector -> coeff parameter link
        c = 5 * np.arange(18) + 18
        spl._c = c.copy()

        def value3(idx):
            return c[idx]
        self.check_parameters(spl, spl._coeff_names, "coeff", value3, False)

        # Check coeff parameter -> coeff vector link
        self.update_parameters(spl, spl._coeff_names, 4)
        assert (spl._c[:] == 4).all()

    def test_two_splines(self):
        spl0 = Spline1D(knots=10)
        spl1 = Spline1D(knots=15, degree=2)

        assert spl0._degree == 3
        assert len(spl0._t) == 18
        t = np.zeros(18)
        t[-4:] = 1
        assert (spl0._t == t).all()
        assert len(spl0._c) == 18
        assert (spl0._c == np.zeros(18)).all()
        assert spl1._degree == 2
        assert len(spl1._t) == 21
        t = np.zeros(21)
        t[-3:] = 1
        assert (spl1._t == t).all()
        assert len(spl1._c) == 21
        assert (spl1._c == np.zeros(21)).all()

        # Check all knot names created
        assert len(spl0._knot_names) == 18
        assert len(spl1._knot_names) == 21

        # Check knot base values
        def value0(idx):
            if idx < 18 - 4:
                return 0
            else:
                return 1
        self.check_parameters(spl0, spl0._knot_names, "knot", value0, True)

        def value1(idx):
            if idx < 21 - 3:
                return 0
            else:
                return 1
        self.check_parameters(spl1, spl1._knot_names, "knot", value1, True)

        # Check knot vector -> knot parameter link
        t0 = 7 * np.arange(18) + 27
        t1 = 11 * np.arange(21) + 19
        spl0._t[:] = t0.copy()
        spl1._t[:] = t1.copy()

        def value2(idx):
            return t0[idx]
        self.check_parameters(spl0, spl0._knot_names, "knot", value2, True)

        def value3(idx):
            return t1[idx]
        self.check_parameters(spl1, spl1._knot_names, "knot", value3, True)

        # Check knot parameter -> knot vector link
        self.update_parameters(spl0, spl0._knot_names, 3)
        self.update_parameters(spl1, spl1._knot_names, 4)
        assert (spl0._t[:] == 3).all()
        assert (spl1._t[:] == 4).all()

        # Check all coeff names created
        assert len(spl0._coeff_names) == 18
        assert len(spl1._coeff_names) == 21

        # Check coeff base values
        def value4(idx):
            return 0
        self.check_parameters(spl0, spl0._coeff_names, "coeff", value4, False)
        self.check_parameters(spl1, spl1._coeff_names, "coeff", value4, False)

        # Check coeff vector -> coeff parameter link
        c0 = 17 * np.arange(18) + 14
        c1 = 37 * np.arange(21) + 47
        spl0._c[:] = c0.copy()
        spl1._c[:] = c1.copy()

        def value5(idx):
            return c0[idx]
        self.check_parameters(spl0, spl0._coeff_names, "coeff", value5, False)

        def value6(idx):
            return c1[idx]
        self.check_parameters(spl1, spl1._coeff_names, "coeff", value6, False)

        # Check coeff parameter -> coeff vector link
        self.update_parameters(spl0, spl0._coeff_names, 5)
        self.update_parameters(spl1, spl1._coeff_names, 6)
        assert (spl0._t[:] == 3).all()
        assert (spl1._t[:] == 4).all()
        assert (spl0._c[:] == 5).all()
        assert (spl1._c[:] == 6).all()

    def test__knot_names(self):
        # no parameters
        spl = Spline1D()
        assert spl._knot_names == ()

        # some parameters
        knot_names = [f"knot{idx}" for idx in range(18)]

        spl = Spline1D(10)
        assert spl._knot_names == tuple(knot_names)

    def test__coeff_names(self):
        # no parameters
        spl = Spline1D()
        assert spl._coeff_names == ()

        # some parameters
        coeff_names = [f"coeff{idx}" for idx in range(18)]

        spl = Spline1D(10)
        assert spl._coeff_names == tuple(coeff_names)

    def test_param_names(self):
        # no parameters
        spl = Spline1D()
        assert spl.param_names == ()

        # some parameters
        knot_names = [f"knot{idx}" for idx in range(18)]
        coeff_names = [f"coeff{idx}" for idx in range(18)]
        param_names = knot_names + coeff_names

        spl = Spline1D(10)
        assert spl.param_names == tuple(param_names)

    def test_t(self):
        # no parameters
        spl = Spline1D()
        # test get
        assert spl._t is None
        assert (spl.t == [0, 0, 0, 0, 1, 1, 1, 1]).all()
        # test set
        with pytest.raises(ValueError) as err:
            spl.t = mk.MagicMock()
        assert str(err.value) == "The model parameters must be initialized before setting knots."

        # with parameters
        spl = Spline1D(10)
        # test get
        t = np.zeros(18)
        t[-4:] = 1
        assert (spl._t == t).all()
        assert (spl.t == t).all()
        # test set
        spl.t = (np.arange(18) + 15)
        assert (spl._t == (np.arange(18) + 15)).all()
        assert (spl.t == (np.arange(18) + 15)).all()
        assert (spl.t != t).all()
        # set error
        for idx in range(30):
            if idx == 18:
                continue
            with pytest.raises(ValueError) as err:
                spl.t = np.arange(idx)
            assert str(err.value) == "There must be exactly as many knots as previously defined."

    def test_c(self):
        # no parameters
        spl = Spline1D()
        # test get
        assert spl._c is None
        assert (spl.c == [0, 0, 0, 0, 0, 0, 0, 0]).all()
        # test set
        with pytest.raises(ValueError) as err:
            spl.c = mk.MagicMock()
        assert str(err.value) == "The model parameters must be initialized before setting coeffs."

        # with parameters
        spl = Spline1D(10)
        # test get
        assert (spl._c == np.zeros(18)).all()
        assert (spl.c == np.zeros(18)).all()
        # test set
        spl.c = (np.arange(18) + 15)
        assert (spl._c == (np.arange(18) + 15)).all()
        assert (spl.c == (np.arange(18) + 15)).all()
        assert (spl.c != np.zeros(18)).all()
        # set error
        for idx in range(30):
            if idx == 18:
                continue
            with pytest.raises(ValueError) as err:
                spl.c = np.arange(idx)
            assert str(err.value) == "There must be exactly as many coeffs as previously defined."

    def test_degree(self):
        # default degree
        spl = Spline1D()
        # test get
        assert spl._degree == 3
        assert spl.degree == 3
        # test set

        # non-default degree
        spl = Spline1D(degree=2)
        # test get
        assert spl._degree == 2
        assert spl.degree == 2

    def test__initialized(self):
        # no parameters
        spl = Spline1D()
        assert spl._initialized is False

        # with parameters
        spl = Spline1D(knots=10, degree=2)
        assert spl._initialized is True

    def test_tck(self):
        # no parameters
        spl = Spline1D()
        # test get
        assert (spl.t == [0, 0, 0, 0, 1, 1, 1, 1]).all()
        assert (spl.c == [0, 0, 0, 0, 0, 0, 0, 0]).all()
        assert spl.degree == 3
        tck = spl.tck
        assert (tck[0] == spl.t).all()
        assert (tck[1] == spl.c).all()
        assert tck[2] == spl.degree
        # test set
        assert spl._t is None
        assert spl._c is None
        assert spl._knot_names == ()
        assert spl._coeff_names == ()
        t = np.array([0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5])
        np.random.seed(619)
        c = np.random.random(12)
        k = 3
        spl.tck = (t, c, k)
        assert (spl._t == t).all()
        assert (spl._c == c).all()
        assert spl.degree == k

        def value0(idx):
            return t[idx]
        self.check_parameters(spl, spl._knot_names, "knot", value0, True)

        def value1(idx):
            return c[idx]
        self.check_parameters(spl, spl._coeff_names, "coeff", value1, False)

        # with parameters
        spl = Spline1D(knots=10, degree=2)
        # test get
        t = np.zeros(16)
        t[-3:] = 1
        assert (spl.t == t).all()
        assert (spl.c == np.zeros(16)).all()
        assert spl.degree == 2
        tck = spl.tck
        assert (tck[0] == spl.t).all()
        assert (tck[1] == spl.c).all()
        assert tck[2] == spl.degree
        # test set
        t = 5*np.arange(16) + 11
        c = 7*np.arange(16) + 13
        k = 2
        spl.tck = (t, c, k)
        assert (spl.t == t).all()
        assert (spl.c == c).all()
        assert spl.degree == k
        tck = spl.tck
        assert (tck[0] == spl.t).all()
        assert (tck[1] == spl.c).all()
        assert tck[2] == spl.degree

        # Error
        with pytest.raises(ValueError) as err:
            spl.tck = (t, c, 4)
        assert str(err.value) == "tck has incompatible degree!"

    def test_bspline(self):
        from scipy.interpolate import BSpline

        # no parameters
        spl = Spline1D()
        bspline = spl.bspline

        assert isinstance(bspline, BSpline)
        assert (bspline.tck[0] == spl.tck[0]).all()
        assert (bspline.tck[1] == spl.tck[1]).all()
        assert bspline.tck[2] == spl.tck[2]

        t = np.array([0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5])
        np.random.seed(619)
        c = np.random.random(12)
        k = 3

        def value0(idx):
            return t[idx]

        def value1(idx):
            return c[idx]

        # set (bspline)
        spl = Spline1D()
        assert spl._t is None
        assert spl._c is None
        assert spl._knot_names == ()
        assert spl._coeff_names == ()
        bspline = BSpline(t, c, k)
        spl.bspline = bspline
        assert (spl._t == t).all()
        assert (spl._c == c).all()
        assert spl.degree == k
        self.check_parameters(spl, spl._knot_names, "knot", value0, True)
        self.check_parameters(spl, spl._coeff_names, "coeff", value1, False)

        # set (tuple spline)
        spl = Spline1D()
        assert spl._t is None
        assert spl._c is None
        assert spl._knot_names == ()
        assert spl._coeff_names == ()
        spl.bspline = (t, c, k)
        assert (spl._t == t).all()
        assert (spl._c == c).all()
        assert spl.degree == k
        self.check_parameters(spl, spl._knot_names, "knot", value0, True)
        self.check_parameters(spl, spl._coeff_names, "coeff", value1, False)

        # with parameters
        spl = Spline1D(knots=10, degree=2)
        bspline = spl.bspline

        assert isinstance(bspline, BSpline)
        assert (bspline.tck[0] == spl.tck[0]).all()
        assert (bspline.tck[1] == spl.tck[1]).all()
        assert bspline.tck[2] == spl.tck[2]

    def test_knots(self):
        # no parameters
        spl = Spline1D()
        assert spl.knots == []

        # with parameters
        spl = Spline1D(10)
        knots = spl.knots
        assert len(knots) == 18

        for knot in knots:
            assert isinstance(knot, Parameter)
            assert hasattr(spl, knot.name)
            assert getattr(spl, knot.name) == knot

    def test_coeffs(self):
        # no parameters
        spl = Spline1D()
        assert spl.coeffs == []

        # with parameters
        spl = Spline1D(10)
        coeffs = spl.coeffs
        assert len(coeffs) == 18

        for coeff in coeffs:
            assert isinstance(coeff, Parameter)
            assert hasattr(spl, coeff.name)
            assert getattr(spl, coeff.name) == coeff

    def test__init_parameters(self):
        spl = Spline1D()

        with mk.patch.object(Spline1D, '_create_parameters',
                             autospec=True) as mkCreate:
            spl._init_parameters()
            assert mkCreate.call_args_list == [
                mk.call(spl, "knot", "t", fixed=True),
                mk.call(spl, "coeff", "c")
            ]

    def test__init_bounds(self):
        spl = Spline1D()

        has_bounds, lower, upper = spl._init_bounds()
        assert has_bounds is False
        assert (lower == [0, 0, 0, 0]).all()
        assert (upper == [1, 1, 1, 1]).all()
        assert spl._user_bounding_box is None

        has_bounds, lower, upper = spl._init_bounds((-5, 5))
        assert has_bounds is True
        assert (lower == [-5, -5, -5, -5]).all()
        assert (upper == [5, 5, 5, 5]).all()
        assert spl._user_bounding_box == (-5, 5)

    def test__init_knots(self):
        np.random.seed(19)
        lower = np.random.random(4)
        upper = np.random.random(4)

        # Integer
        with mk.patch.object(Spline1D, "bspline",
                             new_callable=mk.PropertyMock) as mkBspline:
            spl = Spline1D()
            assert spl._t is None
            spl._init_knots(10, mk.MagicMock(), lower, upper)
            t = np.concatenate((lower, np.zeros(10), upper))
            assert (spl._t == t).all()
            assert mkBspline.call_args_list == [mk.call()]

        # vector with bounds
        with mk.patch.object(Spline1D, "bspline",
                             new_callable=mk.PropertyMock) as mkBspline:
            knots = np.random.random(10)
            spl = Spline1D()
            assert spl._t is None
            spl._init_knots(knots, True, lower, upper)
            t = np.concatenate((lower, knots, upper))
            assert (spl._t == t).all()
            assert mkBspline.call_args_list == [mk.call()]

        # vector with no bounds
        with mk.patch.object(Spline1D, "bspline",
                             new_callable=mk.PropertyMock) as mkBspline:
            knots = np.random.random(10)
            spl = Spline1D()
            assert spl._t is None
            spl._init_knots(knots, False, lower, upper)
            assert (spl._t == knots).all()
            assert mkBspline.call_args_list == [mk.call()]

            # error
            for num in range(8):
                knots = np.random.random(num)
                spl = Spline1D()
                assert spl._t is None
                with pytest.raises(ValueError) as err:
                    spl._init_knots(knots, False, lower, upper)
                assert str(err.value) == "Must have at least 8 knots."

        # Error
        spl = Spline1D()
        assert spl._t is None
        with pytest.raises(ValueError) as err:
            spl._init_knots(0.5, False, lower, upper)
        assert str(err.value) == "Knots: 0.5 must be iterable or value"

    def test__init_coeffs(self):
        np.random.seed(492)
        # No coeffs
        with mk.patch.object(Spline1D, "bspline",
                             new_callable=mk.PropertyMock) as mkBspline:
            spl = Spline1D()
            assert spl._c is None
            spl._t = [1, 2, 3, 4]
            spl._init_coeffs()
            assert (spl._c == [0, 0, 0, 0]).all()
            assert mkBspline.call_args_list == [mk.call()]

        # Some coeffs
        with mk.patch.object(Spline1D, "bspline",
                             new_callable=mk.PropertyMock) as mkBspline:
            coeffs = np.random.random(10)
            spl = Spline1D()
            assert spl._c is None
            spl._init_coeffs(coeffs)
            assert (spl._c == coeffs).all()
            assert mkBspline.call_args_list == [mk.call()]

    def test__init_data(self):
        spl = Spline1D()

        knots = mk.MagicMock()
        coeffs = mk.MagicMock()
        bounds = mk.MagicMock()
        has_bounds = mk.MagicMock()
        lower = mk.MagicMock()
        upper = mk.MagicMock()
        with mk.patch.object(Spline1D, '_init_bounds', autospec=True,
                             return_value=(has_bounds, lower, upper)) as mkBounds:
            with mk.patch.object(Spline1D, '_init_knots',
                                 autospec=True) as mkKnots:
                with mk.patch.object(Spline1D, '_init_coeffs',
                                     autospec=True) as mkCoeffs:
                    main = mk.MagicMock()
                    main.attach_mock(mkBounds, 'bounds')
                    main.attach_mock(mkKnots, 'knots')
                    main.attach_mock(mkCoeffs, 'coeffs')

                    spl._init_data(knots, coeffs, bounds)
                    assert main.mock_calls == [
                        mk.call.bounds(spl, bounds),
                        mk.call.knots(spl, knots, has_bounds, lower, upper),
                        mk.call.coeffs(spl, coeffs)
                    ]

    def test_evaluate(self):
        spl = Spline1D()

        args = tuple([mk.MagicMock() for _ in range(3)])
        kwargs = {f"test{idx}": mk.MagicMock() for idx in range(3)}
        new_kwargs = {f"new_test{idx}": mk.MagicMock() for idx in range(3)}

        with mk.patch.object(_Spline, 'evaluate', autospec=True,
                             return_value=new_kwargs) as mkEval:
            with mk.patch.object(Spline1D, "bspline",
                                 new_callable=mk.PropertyMock) as mkBspline:
                assert mkBspline.return_value.return_value == spl.evaluate(*args, **kwargs)
                assert mkBspline.return_value.call_args_list == [mk.call(args[0], **new_kwargs)]
                assert mkBspline.call_args_list == [mk.call()]
                assert mkEval.call_args_list == [mk.call(spl, *args, **kwargs)]

        # Error
        for idx in range(5, 8):
            with mk.patch.object(_Spline, 'evaluate', autospec=True,
                                 return_value={'nu': idx}):
                with pytest.raises(RuntimeError) as err:
                    spl.evaluate(*args, **kwargs)
                assert str(err.value) == "Cannot evaluate a derivative of order higher than 4"

    def check_knots_created(self, spl, k):
        def value0(idx):
            return self.x[0]

        def value1(idx):
            return self.x[-1]

        for idx in range(k + 1):
            name = f"knot{idx}"
            self.check_parameter(spl, "knot", name, idx, value0, True)

            index = len(spl.t) - (k + 1) + idx
            name = f"knot{index}"
            self.check_parameter(spl, "knot", name, index, value1, True)

        def value3(idx):
            return spl.t[idx]

        assert len(spl._knot_names) == len(spl.t)
        for idx, name in enumerate(spl._knot_names):
            assert name == f"knot{idx}"
            self.check_parameter(spl, "knot", name, idx, value3, True)

    def check_coeffs_created(self, spl):
        def value(idx):
            return spl.c[idx]

        assert len(spl._coeff_names) == len(spl.c)
        for idx, name in enumerate(spl._coeff_names):
            assert name == f"coeff{idx}"
            self.check_parameter(spl, "coeff", name, idx, value, False)

    @staticmethod
    def check_base_spline(spl, t, c, k):
        """Check the base spline form"""
        if t is None:
            assert spl._t is None
        else:
            assert_allclose(spl._t, t)

        if c is None:
            assert spl._c is None
        else:
            assert_allclose(spl._c, c)

        assert spl.degree == k
        assert spl._bounding_box is None

    def check_spline_fit(self, fit_spl, spline, fitter, atol_fit, atol_truth):
        """Check the spline fit"""
        assert_allclose(fit_spl.t, spline._eval_args[0])
        assert_allclose(fit_spl.c, spline._eval_args[1])
        assert_allclose(fitter.fit_info['spline']._eval_args[0], spline._eval_args[0])
        assert_allclose(fitter.fit_info['spline']._eval_args[1],  spline._eval_args[1])

        # check that _parameters are correct
        assert len(fit_spl._parameters) == len(fit_spl.t) + len(fit_spl.c)
        assert_allclose(fit_spl._parameters[:len(fit_spl.t)], fit_spl.t)
        assert_allclose(fit_spl._parameters[len(fit_spl.t):], fit_spl.c)

        # check that parameters are correct
        assert len(fit_spl.parameters) == len(fit_spl.t) + len(fit_spl.c)
        assert_allclose(fit_spl.parameters[:len(fit_spl.t)], fit_spl.t)
        assert_allclose(fit_spl.parameters[len(fit_spl.t):], fit_spl.c)

        assert_allclose(spline.get_residual(), fitter.fit_info['resid'])

        assert_allclose(fit_spl(self.x), spline(self.x))
        assert_allclose(fit_spl(self.x), fitter.fit_info['spline'](self.x))

        assert_allclose(fit_spl(self.x), self.y, atol=atol_fit)
        assert_allclose(fit_spl(self.x), self.truth, atol=atol_truth)

    def check_bbox(self, spl, fit_spl, fitter, w, **kwargs):
        """Check the spline fit with bbox option"""
        bbox = [self.x[0], self.x[-1]]
        bbox_spl = fitter(spl, self.x, self.y, weights=w, bbox=bbox, **kwargs)
        assert bbox_spl.bounding_box == tuple(bbox)
        assert_allclose(fit_spl.t, bbox_spl.t)
        assert_allclose(fit_spl.c, bbox_spl.c)

    def check_knots_warning(self, fitter, knots, k, w, **kwargs):
        """Check that the knots warning is raised"""
        spl = Spline1D(knots=knots, degree=k)
        with pytest.warns(AstropyUserWarning):
            fitter(spl, self.x, self.y, weights=w, **kwargs)

    @pytest.mark.parametrize('w', wieght_tests)
    @pytest.mark.parametrize('k', degree_tests)
    def test_interpolate_fitter(self, w, k):
        fitter = SplineInterpolateFitter()
        assert fitter.fit_info == {'resid': None, 'spline': None}

        spl = Spline1D(degree=k)
        self.check_base_spline(spl, None, None, k)

        fit_spl = fitter(spl, self.x, self.y, weights=w)
        self.check_base_spline(spl, None, None, k)

        assert len(fit_spl.t) == (len(self.x) + k + 1) == len(fit_spl._knot_names)
        self.check_knots_created(fit_spl, k)
        self.check_coeffs_created(fit_spl)
        assert fit_spl._bounding_box is None

        from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline
        spline = InterpolatedUnivariateSpline(self.x, self.y, w=w, k=k)
        assert isinstance(fitter.fit_info['spline'], UnivariateSpline)

        assert spline.get_residual() == 0
        self.check_spline_fit(fit_spl, spline, fitter, 0, 1)
        self.check_bbox(spl, fit_spl, fitter, w)

        knots = np.linspace(self.x[0], self.x[-1], len(self.x) + k + 1)
        self.check_knots_warning(fitter, knots, k, w)

    @pytest.mark.parametrize('w', wieght_tests)
    @pytest.mark.parametrize('k', degree_tests)
    @pytest.mark.parametrize('s', smoothing_tests)
    def test_smoothing_fitter(self, w, k, s):
        fitter = SplineSmoothingFitter()
        assert fitter.fit_info == {'resid': None, 'spline': None}

        spl = Spline1D(degree=k)
        self.check_base_spline(spl, None, None, k)

        fit_spl = fitter(spl, self.x, self.y, s=s, weights=w)
        self.check_base_spline(spl, None, None, k)

        self.check_knots_created(fit_spl, k)
        self.check_coeffs_created(fit_spl)
        assert fit_spl._bounding_box is None

        from scipy.interpolate import UnivariateSpline
        spline = UnivariateSpline(self.x, self.y, w=w, k=k, s=s)
        assert isinstance(fitter.fit_info['spline'], UnivariateSpline)

        self.check_spline_fit(fit_spl, spline, fitter, 1, 1)
        self.check_bbox(spl, fit_spl, fitter, w, s=s)

        # test warning
        knots = fit_spl.t.copy()
        self.check_knots_warning(fitter, knots, k, w, s=s)

    @pytest.mark.parametrize('w', wieght_tests)
    @pytest.mark.parametrize('k', degree_tests)
    def test_exact_knots_fitter(self, w, k):
        fitter = SplineExactKnotsFitter()
        assert fitter.fit_info == {'resid': None, 'spline': None}

        knots = [-1, 0, 1]
        t = np.concatenate(([self.x[0]]*(k + 1), knots, [self.x[-1]]*(k + 1)))
        c = np.zeros(len(t))

        # With knots preset
        spl = Spline1D(knots=knots, degree=k, bounds=[self.x[0], self.x[-1]])
        self.check_base_spline(spl, t, c, k)
        assert (spl.t_interior == knots).all()

        fit_spl = fitter(spl, self.x, self.y, weights=w)
        self.check_base_spline(spl, t, c, k)
        assert (spl.t_interior == knots).all()

        assert len(fit_spl.t) == len(t) == len(fit_spl._knot_names)
        self.check_knots_created(fit_spl, k)
        self.check_coeffs_created(fit_spl)
        assert fit_spl._bounding_box is None

        from scipy.interpolate import LSQUnivariateSpline, UnivariateSpline
        spline = LSQUnivariateSpline(self.x, self.y, knots, w=w, k=k)
        assert isinstance(fitter.fit_info['spline'], UnivariateSpline)

        assert_allclose(spline.get_residual(), 0.1, atol=1)
        assert_allclose(fitter.fit_info['spline'].get_residual(), 0.1, atol=1)
        self.check_spline_fit(fit_spl, spline, fitter, 1, 1)
        self.check_bbox(spl, fit_spl, fitter, w)

        # Pass knots via fitter function
        with pytest.warns(AstropyUserWarning):
            fitter(spl, self.x, self.y, t=knots, weights=w)

        # pass no knots
        spl = Spline1D(degree=k)
        with pytest.raises(RuntimeError) as err:
            fitter(spl, self.x, self.y, weights=w)
        assert str(err.value) == "No knots have been provided"

    @pytest.mark.parametrize('w', wieght_tests)
    @pytest.mark.parametrize('k', degree_tests)
    @pytest.mark.parametrize('s', smoothing_tests)
    def test_splrep_fitter_no_knots(self, w, k, s):
        fitter = SplineSplrepFitter()
        assert fitter.fit_info == {'fp': None, 'ier': None, 'msg': None}

        spl = Spline1D(degree=k)
        self.check_base_spline(spl, None, None, k)

        fit_spl = fitter(spl, self.x, self.y, s=s, weights=w)
        self.check_base_spline(spl, None, None, k)

        self.check_knots_created(fit_spl, k)
        self.check_coeffs_created(fit_spl)
        assert fit_spl._bounding_box is None

        from scipy.interpolate import BSpline, splrep
        tck, spline_fp, spline_ier, spline_msg = splrep(self.x, self.y,
                                                        w=w, k=k, s=s, full_output=1)
        assert_allclose(fit_spl.t, tck[0])
        assert_allclose(fit_spl.c, tck[1])

        assert fitter.fit_info['fp'] == spline_fp
        assert fitter.fit_info['ier'] == spline_ier
        assert fitter.fit_info['msg'] == spline_msg

        spline = BSpline(*tck)
        assert_allclose(fit_spl(self.x), spline(self.x))

        assert_allclose(fit_spl(self.x), self.y, atol=1)
        assert_allclose(fit_spl(self.x), self.truth, atol=1)

        self.check_bbox(spl, fit_spl, fitter, w, s=s)

    @pytest.mark.parametrize('w', wieght_tests)
    @pytest.mark.parametrize('k', degree_tests)
    def test_splrep_fitter_with_knots(self, w, k):
        fitter = SplineSplrepFitter()
        assert fitter.fit_info == {'fp': None, 'ier': None, 'msg': None}

        knots = [-1, 0, 1]
        t = np.concatenate(([self.x[0]]*(k + 1), knots, [self.x[-1]]*(k + 1)))
        c = np.zeros(len(t))

        # With knots preset
        spl = Spline1D(knots=knots, degree=k, bounds=[self.x[0], self.x[-1]])
        self.check_base_spline(spl, t, c, k)
        assert (spl.t_interior == knots).all()

        fit_spl = fitter(spl, self.x, self.y, weights=w)
        self.check_base_spline(spl, t, c, k)
        assert (spl.t_interior == knots).all()

        self.check_knots_created(fit_spl, k)
        self.check_coeffs_created(fit_spl)
        assert fit_spl._bounding_box is None

        from scipy.interpolate import BSpline, splrep
        tck, spline_fp, spline_ier, spline_msg = splrep(self.x, self.y,
                                                        w=w, k=k, t=knots, full_output=1)
        assert_allclose(fit_spl.t, tck[0])
        assert_allclose(fit_spl.c, tck[1])

        assert fitter.fit_info['fp'] == spline_fp
        assert fitter.fit_info['ier'] == spline_ier
        assert fitter.fit_info['msg'] == spline_msg

        spline = BSpline(*tck)
        assert_allclose(fit_spl(self.x), spline(self.x))

        assert_allclose(fit_spl(self.x), self.y, atol=1)
        assert_allclose(fit_spl(self.x), self.truth, atol=1)

        self.check_bbox(spl, fit_spl, fitter, w)

        # test warning
        with pytest.warns(AstropyUserWarning):
            fitter(spl, self.x, self.y, t=knots, weights=w)

        # With no knots present
        spl = Spline1D(degree=k)
        self.check_base_spline(spl, None, None, k)

        fit_spl = fitter(spl, self.x, self.y, t=knots, weights=w)
        self.check_base_spline(spl, None, None, k)

        self.check_knots_created(fit_spl, k)
        self.check_coeffs_created(fit_spl)
        assert fit_spl._bounding_box is None

        from scipy.interpolate import BSpline, splrep
        tck = splrep(self.x, self.y, w=w, k=k, t=knots)
        assert_allclose(fit_spl.t, tck[0])
        assert_allclose(fit_spl.c, tck[1])

        spline = BSpline(*tck)
        assert_allclose(fit_spl(self.x), spline(self.x))

        assert_allclose(fit_spl(self.x), self.y, atol=1)
        assert_allclose(fit_spl(self.x), self.truth, atol=1)

        self.check_bbox(spl, fit_spl, fitter, w, t=knots)

    def generate_spline(self, w=None, bbox=[None]*2, k=None, s=None, t=None):
        if k is None:
            k = 3

        from scipy.interpolate import BSpline, splrep

        tck = splrep(self.x, self.y, w=w, xb=bbox[0], xe=bbox[1],
                     k=k, s=s, t=t)

        return BSpline(*tck)

    def test_derivative(self):
        bspline = self.generate_spline()

        spl = Spline1D()
        spl.bspline = bspline
        assert_allclose(spl.t, bspline.t)
        assert_allclose(spl.c, bspline.c)
        assert spl.degree == bspline.k

        # 1st derivative
        d_bspline = bspline.derivative(nu=1)
        assert_allclose(d_bspline(self.xs),       bspline(self.xs, nu=1))
        assert_allclose(d_bspline(self.xs, nu=1), bspline(self.xs, nu=2))
        assert_allclose(d_bspline(self.xs, nu=2), bspline(self.xs, nu=3))
        assert_allclose(d_bspline(self.xs, nu=3), bspline(self.xs, nu=4))

        der = spl.derivative()
        assert_allclose(der.t, d_bspline.t)
        assert_allclose(der.c, d_bspline.c)
        assert der.degree == d_bspline.k == 2
        assert_allclose(der.evaluate(self.xs),       spl.evaluate(self.xs, nu=1))
        assert_allclose(der.evaluate(self.xs, nu=1), spl.evaluate(self.xs, nu=2))
        assert_allclose(der.evaluate(self.xs, nu=2), spl.evaluate(self.xs, nu=3))
        assert_allclose(der.evaluate(self.xs, nu=3), spl.evaluate(self.xs, nu=4))

        # 2nd derivative
        d_bspline = bspline.derivative(nu=2)
        assert_allclose(d_bspline(self.xs),       bspline(self.xs, nu=2))
        assert_allclose(d_bspline(self.xs, nu=1), bspline(self.xs, nu=3))
        assert_allclose(d_bspline(self.xs, nu=2), bspline(self.xs, nu=4))

        der = spl.derivative(nu=2)
        assert_allclose(der.t, d_bspline.t)
        assert_allclose(der.c, d_bspline.c)
        assert der.degree == d_bspline.k == 1
        assert_allclose(der.evaluate(self.xs),       spl.evaluate(self.xs, nu=2))
        assert_allclose(der.evaluate(self.xs, nu=1), spl.evaluate(self.xs, nu=3))
        assert_allclose(der.evaluate(self.xs, nu=2), spl.evaluate(self.xs, nu=4))

        # 3rd derivative
        d_bspline = bspline.derivative(nu=3)
        assert_allclose(d_bspline(self.xs),       bspline(self.xs, nu=3))
        assert_allclose(d_bspline(self.xs, nu=1), bspline(self.xs, nu=4))

        der = spl.derivative(nu=3)
        assert_allclose(der.t, d_bspline.t)
        assert_allclose(der.c, d_bspline.c)
        assert der.degree == d_bspline.k == 0
        assert_allclose(der.evaluate(self.xs),       spl.evaluate(self.xs, nu=3))
        assert_allclose(der.evaluate(self.xs, nu=1), spl.evaluate(self.xs, nu=4))

        # Too many derivatives
        for nu in range(4, 9):
            with pytest.raises(ValueError) as err:
                spl.derivative(nu=nu)
            assert str(err.value) == "Must have nu <= 3"

    def test_antiderivative(self):
        bspline = self.generate_spline()

        spl = Spline1D()
        spl.bspline = bspline

        # 1st antiderivative
        a_bspline = bspline.antiderivative(nu=1)
        assert_allclose(bspline(self.xs),       a_bspline(self.xs, nu=1))
        assert_allclose(bspline(self.xs, nu=1), a_bspline(self.xs, nu=2))
        assert_allclose(bspline(self.xs, nu=2), a_bspline(self.xs, nu=3))
        assert_allclose(bspline(self.xs, nu=3), a_bspline(self.xs, nu=4))
        assert_allclose(bspline(self.xs, nu=4), a_bspline(self.xs, nu=5))

        anti = spl.antiderivative()
        assert_allclose(anti.t, a_bspline.t)
        assert_allclose(anti.c, a_bspline.c)
        assert anti.degree == a_bspline.k == 4
        assert_allclose(spl.evaluate(self.xs),       anti.evaluate(self.xs, nu=1))
        assert_allclose(spl.evaluate(self.xs, nu=1), anti.evaluate(self.xs, nu=2))
        assert_allclose(spl.evaluate(self.xs, nu=2), anti.evaluate(self.xs, nu=3))
        assert_allclose(spl.evaluate(self.xs, nu=3), anti.evaluate(self.xs, nu=4))
        assert_allclose(spl.evaluate(self.xs, nu=4), anti.evaluate(self.xs, nu=5))

        # 2nd antiderivative
        a_bspline = bspline.antiderivative(nu=2)
        assert_allclose(bspline(self.xs),       a_bspline(self.xs, nu=2))
        assert_allclose(bspline(self.xs, nu=1), a_bspline(self.xs, nu=3))
        assert_allclose(bspline(self.xs, nu=2), a_bspline(self.xs, nu=4))
        assert_allclose(bspline(self.xs, nu=3), a_bspline(self.xs, nu=5))
        assert_allclose(bspline(self.xs, nu=4), a_bspline(self.xs, nu=6))

        anti = spl.antiderivative(nu=2)
        assert_allclose(anti.t, a_bspline.t)
        assert_allclose(anti.c, a_bspline.c)
        assert anti.degree == a_bspline.k == 5
        assert_allclose(spl.evaluate(self.xs),       anti.evaluate(self.xs, nu=2))
        assert_allclose(spl.evaluate(self.xs, nu=1), anti.evaluate(self.xs, nu=3))
        assert_allclose(spl.evaluate(self.xs, nu=2), anti.evaluate(self.xs, nu=4))
        assert_allclose(spl.evaluate(self.xs, nu=3), anti.evaluate(self.xs, nu=5))
        assert_allclose(spl.evaluate(self.xs, nu=4), anti.evaluate(self.xs, nu=6))

        # Too many anti derivatives
        for nu in range(3, 9):
            with pytest.raises(ValueError) as err:
                spl.antiderivative(nu=nu)
            assert str(err.value) == ("Supported splines can have max degree 5, "
                                      f"antiderivative degree will be {nu + 3}")

    def test__SplineFitter_error(self):
        spl = Spline1D()

        class SplineFitter(_SplineFitter):
            def _fit_method(self, model, x, y, **kwargs):
                super()._fit_method(model, x, y, **kwargs)

        fitter = SplineFitter()

        with pytest.raises(ValueError) as err:
            fitter(spl, mk.MagicMock(), mk.MagicMock(), mk.MagicMock())
        assert str(err.value) == "1D model can only have 2 data points."

        with pytest.raises(ModelDefinitionError) as err:
            fitter(mk.MagicMock(), mk.MagicMock(), mk.MagicMock())
        assert str(err.value) == "Only spline models are compatible with this fitter."

        with pytest.raises(NotImplementedError) as err:
            fitter(spl, mk.MagicMock(), mk.MagicMock())
        assert str(err.value) == "This has not been implemented for _SplineFitter."
