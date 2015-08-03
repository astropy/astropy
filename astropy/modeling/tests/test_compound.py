# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import inspect

import numpy as np

from numpy.testing.utils import (assert_allclose, assert_array_equal,
                                 assert_almost_equal)

from ...extern.six.moves import cPickle as pickle
from ...tests.helper import pytest

from ..core import Model, ModelDefinitionError
from ..parameters import Parameter
from ..models import (Const1D, Shift, Scale, Rotation2D, Gaussian1D,
                      Gaussian2D, Polynomial1D, Polynomial2D,
                      Chebyshev2D, Legendre2D, Chebyshev1D, Legendre1D,
                      AffineTransformation2D, Identity, Mapping)


@pytest.mark.parametrize(('expr', 'result'),
                         [(lambda x, y: x + y, 5.0),
                          (lambda x, y: x - y, -1.0),
                          (lambda x, y: x * y, 6.0),
                          (lambda x, y: x / y, 2.0 / 3.0),
                          (lambda x, y: x ** y, 8.0)])
def test_two_model_class_arithmetic_1d(expr, result):
    # Const1D is perhaps the simplest model to test basic arithmetic with.
    # TODO: Should define more tests later on for more complicated
    # combinations of models

    S = expr(Const1D, Const1D)

    assert issubclass(S, Model)
    assert S.n_inputs == 1
    assert S.n_outputs == 1

    # Initialize an instance of the model, providing values for the two
    # "amplitude" parameters
    s = S(2, 3)

    # It shouldn't matter what input we evaluate on since this is a constant
    # function
    out = s(0)
    assert out == result
    assert isinstance(out, float)


@pytest.mark.parametrize(('expr', 'result'),
                         [(lambda x, y: x + y, 5.0),
                          (lambda x, y: x - y, -1.0),
                          (lambda x, y: x * y, 6.0),
                          (lambda x, y: x / y, 2.0 / 3.0),
                          (lambda x, y: x ** y, 8.0)])
def test_two_model_instance_arithmetic_1d(expr, result):
    """
    Like test_two_model_class_arithmetic_1d, but creates a new model from two
    model *instances* with fixed parameters.
    """

    s = expr(Const1D(2), Const1D(3))

    assert isinstance(s, Model)
    assert s.n_inputs == 1
    assert s.n_outputs == 1

    out = s(0)
    assert out == result
    assert isinstance(out, float)


@pytest.mark.parametrize(('expr', 'result'),
                         [(lambda x, y: x + y, 5.0),
                          (lambda x, y: x - y, -1.0),
                          (lambda x, y: x * y, 6.0),
                          (lambda x, y: x / y, 2.0 / 3.0),
                          (lambda x, y: x ** y, 8.0)])
def test_two_model_mixed_arithmetic_1d(expr, result):
    """
    Like test_two_model_class_arithmetic_1d, but creates a new model from an
    expression of one model class with one model instance (and vice-versa).
    """

    S1 = expr(Const1D, Const1D(3))
    S2 = expr(Const1D(2), Const1D)

    for cls in (S1, S2):
        assert issubclass(cls, Model)
        assert cls.n_inputs == 1
        assert cls.n_outputs == 1

    # Requires values for both amplitudes even though one of them them has a
    # default
    # TODO: We may wish to fix that eventually, so that if a parameter has a
    # default it doesn't *have* to be given in the init
    s1 = S1(2, 3)
    s2 = S2(2, 3)

    for out in (s1(0), s2(0)):
        assert out == result
        assert isinstance(out, float)


def test_simple_two_model_class_compose_1d():
    """
    Shift and Scale are two of the simplest models to test model composition
    with.
    """

    S1 = Shift | Scale  # First shift then scale
    assert issubclass(S1, Model)
    assert S1.n_inputs == 1
    assert S1.n_outputs == 1

    s1 = S1(2, 3)  # Shift by 2 and scale by 3
    assert s1(1) == 9.0

    S2 = Scale | Shift  # First scale then shift
    assert issubclass(S2, Model)
    assert S2.n_inputs == 1
    assert S2.n_outputs == 1

    s2 = S2(2, 3)  # Scale by 2 then shift by 3
    assert s2(1) == 5.0

    # Test with array inputs
    assert_array_equal(s2([1, 2, 3]), [5.0, 7.0, 9.0])


def test_simple_two_model_class_compose_2d():
    """
    A simple example consisting of two rotations.
    """

    R = Rotation2D | Rotation2D
    assert issubclass(R, Model)
    assert R.n_inputs == 2
    assert R.n_outputs == 2

    r1 = R(45, 45)  # Rotate twice by 45 degrees
    assert_allclose(r1(0, 1), (-1, 0), atol=1e-10)

    r2 = R(90, 90)  # Rotate twice by 90 degrees
    assert_allclose(r2(0, 1), (0, -1), atol=1e-10)

    # Compose R with itself to produce 4 rotations
    R2 = R | R

    r3 = R2(45, 45, 45, 45)
    assert_allclose(r3(0, 1), (0, -1), atol=1e-10)


class TestCompositeLegacy(object):
    """
    Tests inspired by the original _CompositeModel tests in test_core.py,
    this implements the equivalent tests implemented in the new framework.

    Note: These aren't *exactly* the same as the original tests, as they used
    overly trivial models (polynomials with all coeffs 0).
    """

    def setup_class(self):
        self.y, self.x = np.mgrid[:5, :5]

    def test_single_array_input(self):
        p1 = Polynomial1D(3, c0=1, c1=2, c2=3, c3=4)
        p2 = Polynomial1D(3, c0=2, c1=3, c2=4, c3=5)
        m = p1 | p2
        assert_almost_equal(p2(p1(self.x)), m(self.x))

    def test_labeledinput_1(self):
        # Note: No actual use of LabeledInput in this test; this just uses the
        # same name for symmetry with the old tests
        p1 = Polynomial1D(3, c0=1, c1=2, c2=3, c3=4)
        p2 = Polynomial2D(3, c0_0=1, c2_0=2, c0_1=3, c2_1=4)
        m = p2 | p1
        assert_almost_equal(p1(p2(self.x, self.y)), m(self.x, self.y))

    def test_labledinput_2(self):
        rot = Rotation2D(angle=23.4)
        offx = Shift(-2)
        offy = Shift(1.2)
        m = rot | (offx & Identity(1)) | (Identity(1) & offy)

        x, y = rot(self.x, self.y)
        x = offx(x)
        y = offy(y)

        assert_almost_equal(x, m(self.x, self.y)[0])
        assert_almost_equal(y, m(self.x, self.y)[1])

        a = np.deg2rad(23.4)
        # For kicks
        matrix = [[np.cos(a), -np.sin(a)],
                  [np.sin(a), np.cos(a)]]
        x, y = AffineTransformation2D(matrix, [-2, 1.2])(self.x, self.y)
        assert_almost_equal(x, m(self.x, self.y)[0])
        assert_almost_equal(y, m(self.x, self.y)[1])

    def test_multiple_input(self):
        """
        Despite the name, this actually tests inverting composite models,
        which is not yet supported in the new framework (but should be).
        """

        rot = Rotation2D(-60)
        m = rot | rot
        xx, yy = m(self.x, self.y)
        x0, y0 = m.inverse(xx, yy)
        assert_almost_equal(x0, self.x)
        assert_almost_equal(y0, self.y)


def test_expression_formatting():
    """
    Test that the expression strings from compound models are formatted
    correctly.
    """

    # For the purposes of this test it doesn't matter a great deal what
    # model(s) are used in the expression, I don't think
    G = Gaussian1D
    G2 = Gaussian2D

    M = G + G
    assert M._format_expression() == '[0] + [1]'

    M = G + G + G
    assert M._format_expression() == '[0] + [1] + [2]'

    M = G + G * G
    assert M._format_expression() == '[0] + [1] * [2]'

    M = G * G + G
    assert M._format_expression() == '[0] * [1] + [2]'

    M = G + G * G + G
    assert M._format_expression() == '[0] + [1] * [2] + [3]'

    M = (G + G) * (G + G)
    assert M._format_expression() == '([0] + [1]) * ([2] + [3])'

    # This example uses parentheses in the expression, but those won't be
    # preserved in the expression formatting since they technically aren't
    # necessary, and there's no way to know that they were originally
    # parenthesized (short of some deep, and probably not worthwhile
    # introspection)
    M = (G * G) + (G * G)
    assert M._format_expression() == '[0] * [1] + [2] * [3]'

    M = G ** G
    assert M._format_expression() == '[0] ** [1]'

    M = G + G ** G
    assert M._format_expression() == '[0] + [1] ** [2]'

    M = (G + G) ** G
    assert M._format_expression() == '([0] + [1]) ** [2]'

    M = G + G | G
    assert M._format_expression() == '[0] + [1] | [2]'

    M = G + (G | G)
    assert M._format_expression() == '[0] + ([1] | [2])'

    M = G & G | G2
    assert M._format_expression() == '[0] & [1] | [2]'

    M = G & (G | G)
    assert M._format_expression() == '[0] & ([1] | [2])'


def test_indexing_on_class():
    """
    Test indexing on compound model class objects, including cases where the
    submodels are classes, as well as instances, or both.
    """

    g = Gaussian1D(1, 2, 3, name='g')
    p = Polynomial1D(2, name='p')

    M = Gaussian1D + Const1D
    assert M[0] is Gaussian1D
    assert M[1] is Const1D
    assert M['Gaussian1D'] is M[0]
    assert M['Const1D'] is M[1]

    M = Gaussian1D + p
    assert M[0] is Gaussian1D
    assert M[1] is p
    assert M['Gaussian1D'] is M[0]
    assert M['p'] is M[1]

    m = g + p
    assert isinstance(m[0], Gaussian1D)
    assert isinstance(m[1], Polynomial1D)
    assert isinstance(m['g'], Gaussian1D)
    assert isinstance(m['p'], Polynomial1D)

    # Test negative indexing
    assert isinstance(m[-1], Polynomial1D)
    assert isinstance(m[-2], Gaussian1D)

    with pytest.raises(IndexError):
        m[42]

    with pytest.raises(IndexError):
        m['foobar']


# TODO: It would be good if there were an easier way to interrogate a compound
# model class for what expression it represents.  Not sure what that would look
# like though.
def test_slicing_on_class():
    """
    Test slicing a simple compound model class using integers.
    """

    A = Const1D.rename('A')
    B = Const1D.rename('B')
    C = Const1D.rename('C')
    D = Const1D.rename('D')
    E = Const1D.rename('E')
    F = Const1D.rename('F')

    M = A + B - C * D / E ** F

    assert M[0:1] is A
    # This test will also check that the correct parameter names are generated
    # for each slice (fairly trivial in this case since all the submodels have
    # the same parameter, but if any corner cases are found that aren't covered
    # by this test we can do something different...)
    assert M[0:1].param_names == ('amplitude',)
    # This looks goofy but if you slice by name to the sub-model of the same
    # name it should just return that model, logically.
    assert M['A':'A'] is A
    assert M['A':'A'].param_names == ('amplitude',)
    assert M[5:6] is F
    assert M[5:6].param_names == ('amplitude',)
    assert M['F':'F'] is F
    assert M['F':'F'].param_names == ('amplitude',)

    # 1 + 2
    assert M[:2](1, 2)(0) == 3
    assert M[:2].param_names == ('amplitude_0', 'amplitude_1')
    assert M[:'B'](1, 2)(0) == 3
    assert M[:'B'].param_names == ('amplitude_0', 'amplitude_1')
    # 2 - 3
    assert M[1:3](2, 3)(0) == -1
    assert M[1:3].param_names == ('amplitude_1', 'amplitude_2')
    assert M['B':'C'](2, 3)(0) == -1
    assert M['B':'C'].param_names == ('amplitude_1', 'amplitude_2')
    # 3 * 4
    assert M[2:4](3, 4)(0) == 12
    assert M[2:4].param_names == ('amplitude_2', 'amplitude_3')
    assert M['C':'D'](3, 4)(0) == 12
    assert M['C':'D'].param_names == ('amplitude_2', 'amplitude_3')
    # 4 / 5
    assert M[3:5](4, 5)(0) == 0.8
    assert M[3:5].param_names == ('amplitude_3', 'amplitude_4')
    assert M['D':'E'](4, 5)(0) == 0.8
    assert M['D':'E'].param_names == ('amplitude_3', 'amplitude_4')
    # 5 ** 6
    assert M[4:6](5, 6)(0) == 15625
    assert M[4:6].param_names == ('amplitude_4', 'amplitude_5')
    assert M['E':'F'](5, 6)(0) == 15625
    assert M['E':'F'].param_names == ('amplitude_4', 'amplitude_5')


def test_slicing_on_instance():
    """
    Test slicing a simple compound model class using integers.
    """

    A = Const1D.rename('A')
    B = Const1D.rename('B')
    C = Const1D.rename('C')
    D = Const1D.rename('D')
    E = Const1D.rename('E')
    F = Const1D.rename('F')

    M = A + B - C * D / E ** F
    m = M(1, 2, 3, 4, 5, 6)

    assert isinstance(m[0:1], A)
    assert isinstance(m['A':'A'], A)
    assert isinstance(m[5:6], F)
    assert isinstance(m['F':'F'], F)

    # 1 + 2
    assert m[:'B'](0) == 3
    assert m[:'B'].param_names == ('amplitude_0', 'amplitude_1')
    assert np.all(m[:'B'].parameters == [1, 2])
    # 2 - 3
    assert m['B':'C'](0) == -1
    assert m['B':'C'].param_names == ('amplitude_1', 'amplitude_2')
    assert np.all(m['B':'C'].parameters == [2, 3])
    # 3 * 4
    assert m['C':'D'](0) == 12
    assert m['C':'D'].param_names == ('amplitude_2', 'amplitude_3')
    assert np.all(m['C':'D'].parameters == [3, 4])
    # 4 / 5
    assert m['D':'E'](0) == 0.8
    assert m['D':'E'].param_names == ('amplitude_3', 'amplitude_4')
    assert np.all(m['D':'E'].parameters == [4, 5])
    # 5 ** 6
    assert m['E':'F'](0) == 15625
    assert m['E':'F'].param_names == ('amplitude_4', 'amplitude_5')
    assert np.all(m['E':'F'].parameters == [5, 6])


def test_indexing_on_instance():
    """Test indexing on compound model instances."""

    M = Gaussian1D + Const1D
    m = M(1, 0, 0.1, 2)
    assert isinstance(m[0], Gaussian1D)
    assert isinstance(m[1], Const1D)
    assert isinstance(m['Gaussian1D'], Gaussian1D)
    assert isinstance(m['Const1D'], Const1D)

    # Test parameter equivalence
    assert m[0].amplitude == 1 == m.amplitude_0
    assert m[0].mean == 0 == m.mean_0
    assert m[0].stddev == 0.1 == m.stddev_0
    assert m[1].amplitude == 2 == m.amplitude_1

    # Test that parameter value updates are symmetric between the compound
    # model and the submodel returned by indexing
    const = m[1]
    m.amplitude_1 = 42
    assert const.amplitude == 42
    const.amplitude = 137
    assert m.amplitude_1 == 137


    # Similar couple of tests, but now where the compound model was created
    # from model instances
    g = Gaussian1D(1, 2, 3, name='g')
    p = Polynomial1D(2, name='p')
    m = g + p
    assert m[0].name == 'g'
    assert m[1].name == 'p'
    assert m['g'].name == 'g'
    assert m['p'].name == 'p'

    poly = m[1]
    m.c0_1 = 12345
    assert poly.c0 == 12345
    poly.c1 = 6789
    assert m.c1_1 == 6789

    # Ensure this did *not* modify the original models we used as templates
    assert p.c0 == 0
    assert p.c1 == 0

    # Test negative indexing
    assert isinstance(m[-1], Polynomial1D)
    assert isinstance(m[-2], Gaussian1D)

    with pytest.raises(IndexError):
        m[42]

    with pytest.raises(IndexError):
        m['foobar']


def test_basic_compound_inverse():
    """
    Test basic inversion of compound models in the limited sense supported for
    models made from compositions and joins only.
    """

    t = (Shift(2) & Shift(3)) | (Scale(2) & Scale(3)) | Rotation2D(90)
    assert_allclose(t.inverse(*t(0, 1)), (0, 1))


@pytest.mark.parametrize('model', [
    Shift(0) + Shift(0) | Shift(0),
    Shift(0) - Shift(0) | Shift(0),
    Shift(0) * Shift(0) | Shift(0),
    Shift(0) / Shift(0) | Shift(0),
    Shift(0) ** Shift(0) | Shift(0),
    Gaussian1D(1, 2, 3) | Gaussian1D(4, 5, 6)])
def test_compound_unsupported_inverse(model):
    """
    Ensure inverses aren't supported in cases where it shouldn't be.
    """

    with pytest.raises(NotImplementedError):
        model.inverse


def test_mapping_basic_permutations():
    """
    Tests a couple basic examples of the Mapping model--specifically examples
    that merely permute the outputs.
    """

    x, y = Rotation2D(90)(1, 2)

    RS = Rotation2D | Mapping((1, 0))
    x_prime, y_prime = RS(90)(1, 2)
    assert_allclose((x, y), (y_prime, x_prime))

    # A more complicated permutation
    M = Rotation2D & Scale
    m = M(90, 2)
    x, y, z = m(1, 2, 3)

    MS = M | Mapping((2, 0, 1))
    ms = MS(90, 2)
    x_prime, y_prime, z_prime = ms(1, 2, 3)
    assert_allclose((x, y, z), (y_prime, z_prime, x_prime))


def test_mapping_inverse():
    """Tests inverting a compound model that includes a `Mapping`."""

    RS = Rotation2D & Scale

    # Rotates 2 of the coordinates and scales the third--then rotates on a
    # different axis and scales on the axis of rotation.  No physical meaning
    # here just a simple test
    M = RS | Mapping([2, 0, 1]) | RS

    m = M(12.1, 13.2, 14.3, 15.4)

    assert_allclose((0, 1, 2), m.inverse(*m(0, 1, 2)), atol=1e-08)


def test_identity_input():
    """
    Test a case where an Identity (or Mapping) model is the first in a chain
    of composite models and thus is responsible for handling input broadcasting
    properly.

    Regression test for https://github.com/astropy/astropy/pull/3362
    """

    ident1 = Identity(1)
    shift = Shift(1)
    rotation = Rotation2D(angle=90)
    model = ident1 & shift | rotation
    assert_allclose(model(1, 2), [-3.0, 1.0])

    # Same test case but using class composition
    TestModel = ident1 & Shift | Rotation2D
    model = TestModel(offset_1=1, angle_2=90)
    assert_allclose(model(1, 2), [-3.0, 1.0])


def test_slicing_on_instances_2():
    """
    More slicing tests.

    Regression test for https://github.com/embray/astropy/pull/10
    """

    model_a = Shift(1, name='a')
    model_b = Shift(2, name='b')
    model_c = Rotation2D(3, name='c')
    model_d = Scale(2, name='d')
    model_e = Scale(3, name='e')

    m = (model_a & model_b) | model_c | (model_d & model_e)

    with pytest.raises(ModelDefinitionError):
        # The slice can't actually be taken since the resulting model cannot be
        # evaluated
        assert m[1:].submodel_names == ('b', 'c', 'd', 'e')

    assert m[:].submodel_names == ('a', 'b', 'c', 'd', 'e')
    assert m['a':].submodel_names == ('a', 'b', 'c', 'd', 'e')

    with pytest.raises(ModelDefinitionError):
        assert m['c':'d'].submodel_names == ('c', 'd')

    assert m[1:2].name == 'b'
    assert m[2:7].submodel_names == ('c', 'd', 'e')
    with pytest.raises(IndexError):
        m['x']
    with pytest.raises(IndexError):
        m['a' : 'r']

    with pytest.raises(ModelDefinitionError):
        assert m[-4:4].submodel_names == ('b', 'c', 'd')

    with pytest.raises(ModelDefinitionError):
        assert m[-4:-2].submodel_names == ('b', 'c')


def test_slicing_on_instances_3():
    """
    Like `test_slicing_on_instances_2` but uses a compound model that does not
    have any invalid slices due to the resulting model being invalid
    (originally test_slicing_on_instances_2 passed without any
    ModelDefinitionErrors being raised, but that was before we prevented
    invalid models from being created).
    """

    model_a = Shift(1, name='a')
    model_b = Shift(2, name='b')
    model_c = Gaussian1D(3, 0, 0.1, name='c')
    model_d = Scale(2, name='d')
    model_e = Scale(3, name='e')

    m = (model_a + model_b) | model_c | (model_d + model_e)

    assert m[1:].submodel_names == ('b', 'c', 'd', 'e')
    assert m[:].submodel_names == ('a', 'b', 'c', 'd', 'e')
    assert m['a':].submodel_names == ('a', 'b', 'c', 'd', 'e')
    assert m['c':'d'].submodel_names == ('c', 'd')
    assert m[1:2].name == 'b'
    assert m[2:7].submodel_names == ('c', 'd', 'e')
    with pytest.raises(IndexError):
        m['x']
    with pytest.raises(IndexError):
        m['a' : 'r']
    assert m[-4:4].submodel_names == ('b', 'c', 'd')
    assert m[-4:-2].submodel_names == ('b', 'c')


def test_slicing_on_instance_with_parameterless_model():
    """
    Regression test to fix an issue where the indices attached to parameter
    names on a compound model were not handled properly when one or more
    submodels have no parameters.  This was especially evident in slicing.
    """

    p2 = Polynomial2D(1, c0_0=1, c1_0=2, c0_1=3)
    p1 = Polynomial2D(1, c0_0=1, c1_0=2, c0_1=3)
    mapping = Mapping((0, 1, 0, 1))
    offx = Shift(-2, name='x_translation')
    offy = Shift(-1, name='y_translation')
    aff = AffineTransformation2D(matrix=[[1, 2], [3, 4]], name='rotation')
    model = mapping | (p1 & p2) | (offx & offy) | aff

    assert model.param_names == ('c0_0_1', 'c1_0_1', 'c0_1_1',
                                 'c0_0_2', 'c1_0_2', 'c0_1_2',
                                 'offset_3', 'offset_4',
                                 'matrix_5', 'translation_5')
    assert model(1, 2) == (23.0, 53.0)

    m = model[3:]
    assert m.param_names == ('offset_3', 'offset_4', 'matrix_5',
                             'translation_5')
    assert m(1, 2) == (1.0, 1.0)


def test_compound_model_with_nonstandard_broadcasting():
    """
    Ensure that the ``standard_broadcasting`` flag is properly propagated when
    creating compound models.

    See the commit message for the commit in which this was added for more
    details.
    """

    offx = Shift(1)
    offy = Shift(2)
    rot = AffineTransformation2D([[0, -1], [1, 0]])
    m = (offx & offy) | rot

    x, y = m(0, 0)
    assert x == -2
    assert y == 1

    # make sure conversion back to scalars is working properly
    assert isinstance(x, float)
    assert isinstance(y, float)

    x, y = m([0, 1, 2], [0, 1, 2])
    assert np.all(x == [-2, -3, -4])
    assert np.all(y == [1, 2, 3])


def test_compound_model_classify_attributes():
    """
    Regression test for an issue raised here:
    https://github.com/astropy/astropy/pull/3231#discussion_r22221123

    The issue is that part of the `help` implementation calls a utility
    function called `inspect.classify_class_attrs`, which was leading to an
    infinite recursion.

    This is a useful test in its own right just in that it tests that compound
    models can be introspected in some useful way without crashing--this works
    as sort of a test of its somewhat complicated internal state management.

    This test does not check any of the results of
    `~inspect.classify_class_attrs`, though it might be useful to at some
    point.
    """

    inspect.classify_class_attrs(Gaussian1D + Gaussian1D)


def test_invalid_operands():
    """
    Test that certain operators do not work with models whose inputs/outputs do
    not match up correctly.
    """

    with pytest.raises(ModelDefinitionError):
        Rotation2D | Gaussian1D

    with pytest.raises(ModelDefinitionError):
        Rotation2D(90) | Gaussian1D(1, 0, 0.1)

    with pytest.raises(ModelDefinitionError):
        Rotation2D + Gaussian1D

    with pytest.raises(ModelDefinitionError):
        Rotation2D(90) + Gaussian1D(1, 0, 0.1)


class _ConstraintsTestA(Model):
    stddev = Parameter(default=0, min=0, max=0.3)
    mean = Parameter(default=0, fixed=True)

    @staticmethod
    def evaluate(stddev, mean):
        return stddev, mean


class _ConstraintsTestB(Model):
    mean = Parameter(default=0, fixed=True)

    @staticmethod
    def evaluate(mean):
        return mean


@pytest.mark.parametrize('model',
        [Gaussian1D(bounds={'stddev': (0, 0.3)}, fixed={'mean': True}) +
            Gaussian1D(fixed={'mean': True}),
         (_ConstraintsTestA + _ConstraintsTestB)()])
def test_inherit_constraints(model):
    """
    Various tests for copying of constraint values between compound models and
    their members.

    There are two versions of this test: One where a compound model is created
    from two model instances, and another where a compound model is created
    from two model classes that have default constraints set on some of their
    parameters.

    Regression test for https://github.com/astropy/astropy/issues/3481
    """

    # Lots of assertions in this test as there are multiple interfaces to
    # parameter constraints

    assert 'stddev_0' in model.bounds
    assert model.bounds['stddev_0'] == (0, 0.3)
    assert model.stddev_0.bounds == (0, 0.3)
    assert 'mean_0' in model.fixed
    assert model.fixed['mean_0'] is True
    assert model.mean_0.fixed is True
    assert 'mean_1' in model.fixed
    assert model.fixed['mean_1'] is True
    assert model.mean_1.fixed is True

    # Great, all the constraints were inherited properly
    # Now what about if we update them through the sub-models?
    model[0].stddev.bounds = (0, 0.4)
    assert model.bounds['stddev_0'] == (0, 0.4)
    assert model.stddev_0.bounds == (0, 0.4)
    assert model[0].stddev.bounds == (0, 0.4)
    assert model[0].bounds['stddev'] == (0, 0.4)

    model[0].bounds['stddev'] = (0.1, 0.5)
    assert model.bounds['stddev_0'] == (0.1, 0.5)
    assert model.stddev_0.bounds == (0.1, 0.5)
    assert model[0].stddev.bounds == (0.1, 0.5)
    assert model[0].bounds['stddev'] == (0.1, 0.5)

    model[1].mean.fixed = False
    assert model.fixed['mean_1'] is False
    assert model.mean_1.fixed is False
    assert model[1].mean.fixed is False
    assert model[1].fixed['mean'] is False

    model[1].fixed['mean'] = True
    assert model.fixed['mean_1'] is True
    assert model.mean_1.fixed is True
    assert model[1].mean.fixed is True
    assert model[1].fixed['mean'] is True


def test_compound_custom_inverse():
    """
    Test that a compound model with a custom inverse has that inverse applied
    when the inverse of another model, of which it is a component, is computed.
    Regression test for https://github.com/astropy/astropy/issues/3542
    """

    poly = Polynomial1D(1, c0=1, c1=2)
    scale = Scale(1)
    shift = Shift(1)

    model1 = poly | scale
    model1.inverse = poly

    # model1 now has a custom inverse (the polynomial itself, ignoring the
    # trivial scale factor)
    model2 = shift | model1

    assert_allclose(model2.inverse(1), (poly | shift.inverse)(1))

    # Make sure an inverse is not allowed if the models were combined with the
    # wrong operator, or if one of the models doesn't have an inverse defined
    with pytest.raises(NotImplementedError):
        (shift + model1).inverse

    with pytest.raises(NotImplementedError):
        (model1 & poly).inverse


@pytest.mark.parametrize('poly',
                         [Chebyshev1D(5), Legendre1D(5), Polynomial1D(5)])
def test_compound_with_polynomials_1d(poly):
    """
    Tests that polynomials are scaled when used in compound models.
    Issue #3699
    """

    poly.parameters = [1, 2, 3, 4, 1, 2]
    shift = Shift(3)
    model = poly | shift
    x = np.arange(20)
    result_compound = model(x)
    result = shift(poly(x))
    assert_allclose(result, result_compound)


@pytest.mark.parametrize('poly',
                         [Chebyshev2D(1, 2), Polynomial2D(2),
                          Legendre2D(1, 2)])
def test_compound_with_polynomials_2d(poly):
    """
    Tests that polynomials are scaled when used in compound models.
    Issue #3699
    """
    poly.parameters = [1, 2, 3, 4, 1, 2]
    shift = Shift(3)
    model = poly | shift
    x, y = np.mgrid[:20, :37]
    result_compound = model(x, y)
    result = shift(poly(x, y))
    assert_allclose(result, result_compound)


# has to be defined at module level since pickling doesn't work right (in
# general) for classes defined in functions
class _TestPickleModel(Gaussian1D + Gaussian1D):
    pass


@pytest.mark.skipif(str("sys.version_info < (2, 7, 3)"))
def test_pickle_compound():
    """
    Regression test for
    https://github.com/astropy/astropy/issues/3867#issuecomment-114547228
    """

    # Test pickling a compound model class
    GG = Gaussian1D + Gaussian1D
    GG2 = pickle.loads(pickle.dumps(GG))
    assert GG.param_names == GG2.param_names
    assert GG.__name__ == GG2.__name__
    # Test that it works, or at least evaluates successfully
    assert GG()(0.12345) == GG2()(0.12345)

    # Test pickling a compound model instance
    g1 = Gaussian1D(1.0, 0.0, 0.1)
    g2 = Gaussian1D([2.0, 3.0], [0.0, 0.0], [0.2, 0.3])
    m = g1 + g2
    m2 = pickle.loads(pickle.dumps(m))
    assert m.param_names == m2.param_names
    assert m.__class__.__name__ == m2.__class__.__name__
    assert np.all(m.parameters == m2.parameters)
    assert np.all(m(0) == m2(0))

    # Test picling a concrete class
    p = pickle.dumps(_TestPickleModel, protocol=0)
    # Note: This is very dependent on the specific protocol, but the point of
    # this test is that the "concrete" model is pickled in a very simple way
    # that only specifies the module and class name, and is unpickled by
    # re-importing the class from the module in which it was defined.  This
    # should still work for concrete subclasses of compound model classes that
    # were dynamically generated through an expression
    exp = b'castropy.modeling.tests.test_compound\n_TestPickleModel\np0\n.'
    # When testing against the expected value we drop the memo length field
    # at the end, which may differ between runs
    assert p[:p.rfind(b'p')] == exp[:exp.rfind(b'p')]
    assert pickle.loads(p) is _TestPickleModel


@pytest.mark.skipif(str("sys.version_info >= (2, 7, 3)"))
def test_pickle_compound_fallback():
    """
    Test fallback for pickling compound model on old versions of Python
    affected by http://bugs.python.org/issue7689
    """

    gg = (Gaussian1D + Gaussian1D)()
    with pytest.raises(RuntimeError):
        pickle.dumps(gg)
