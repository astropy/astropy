# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for model evaluation.
Compare the results of some models with other programs.
"""
from __future__ import division
import numpy as np
from numpy.testing.utils import assert_equal, assert_almost_equal
from ...tests.helper import pytest
from .. import models
from .. import LabeledInput, SCompositeModel, PCompositeModel

def test_LabeledInput():
    x = np.ones((2,3))
    y = np.arange(5)
    li = LabeledInput([x, y], ['x', 'y'])
    assert_equal(li.x, x)
    assert_equal(li['y'], y)

    li.add('a', 10)
    assert li.a == 10

    li.add(b=42)
    assert li.b == 42

    li.add('c', 43, d=44)
    assert li.c == 43
    assert li.d == 44

    del li.x
    with pytest.raises(AttributeError) as exc:
        li.x
    assert exc.value.args[0] == "'LabeledInput' object has no attribute '_x'"


def test_SCompositeModel_example():
    """Test the example given in the SCompositeModel docstring"""
    # from astropy.modeling import models, LabeledInput, SCompositeModel
    # Set up the serial composite model:
    # 2D rotation followed by a shift in x and y
    rotation = models.MatrixRotation2D(angle=90)
    shift_x = models.ShiftModel(2)
    shift_y = models.ShiftModel(5)
    model = SCompositeModel([rotation, shift_x, shift_y],
                            inmap=[['x', 'y'], ['x'], ['y']],
                            outmap=[['x', 'y'], ['x'], ['y']])
    # Evaluate the model
    input_pos = LabeledInput([0, 1], ["x", "y"])
    output_pos = model(input_pos)
    # 90 deg clockwise rotation: [0, 1] -> [1, 0]
    # x shift by 2:              [1, 0] -> [3, 0]
    # y shift by 5:              [3, 0] -> [3, 5]
    assert_almost_equal([output_pos.x, output_pos.y], [3, 5])


class TestSComposite(object):
    """
    Test composite models evaluation in series
    """
    def setup_class(self):
        self.x, self.y = np.mgrid[:5, :5]
        self.p1= models.Poly1DModel(3)
        self.p11= models.Poly1DModel(3)
        self.p2 = models.Poly2DModel(3)
        
    def test_single_array_input(self):
        model = SCompositeModel([self.p1, self.p11])
        sresult = model(self.x)
        xx = self.p11(self.p1(self.x))
        assert_almost_equal(xx, sresult)
        
    def test_labeledinput(self):
        labeled_input = LabeledInput([self.x, self.y], ['x', 'y'])
        model = SCompositeModel([self.p2, self.p1], [['x', 'y'], ['z']], [['z'], ['z']])
        sresult = model(labeled_input)
        z = self.p2(self.x, self.y)
        z1 = self.p1(z)
        assert_almost_equal(z1, sresult.z)
        
    def test_multiple_arrays(self):
        model = SCompositeModel([self.p2, self.p1], [['x', 'y'], ['z']], [['z'], ['z']])
        sresult = model(self.x, self.y)
        z = self.p2(self.x, self.y)
        z1 = self.p1(z)
        assert_almost_equal(z1, sresult)

    def test_add_model(self):
        # Redo the test_multiple_arrays example, but build the
        # composite model in steps using add_model
        model = SCompositeModel([self.p2], [['x', 'y']], ['z'])
        model.add_model(self.p1, 'z', 'z')
        sresult = model(self.x, self.y)
        z = self.p2(self.x, self.y)
        z1 = self.p1(z)
        assert_almost_equal(z1, sresult)

    def test_str(self):
        model = SCompositeModel([self.p2, self.p1], [['x', 'y'], ['z']], [['z'], ['z']])
        # TODO: do we care about the format to stay the same e.g. to be able to parse it?
        # If so we should add checks here. 
        str(model)
        repr(model)
        
class TestPComposite(object):
    """
    Test composite models evaluation in parallel
    """
    def setup_class(self):
        self.x, self.y = np.mgrid[:5, :5]
        self.p1= models.Poly1DModel(3)
        self.p11= models.Poly1DModel(3)
        self.p2 = models.Poly2DModel(3)
        
    def test_single_array_input(self):
        model = PCompositeModel([self.p1, self.p11])
        presult = model(self.x)
        delta11 = self.p11(self.x) - self.x
        delta1 = self.p1(self.x) - self.x
        xx = self.x + delta1 + delta11
        assert_almost_equal(xx, presult)
        
    def test_labeledinput(self):
        labeled_input = LabeledInput([self.x, self.y], ['x', 'y'])
        model = PCompositeModel([self.p1, self.p11], inmap=['x'], outmap=['x'])
        presult = model(labeled_input)
        delta11 = self.p11(self.x) - self.x
        delta1 = self.p1(self.x) - self.x
        xx = self.x + delta1 + delta11
        assert_almost_equal(xx, presult.x)

    def test_str(self):
        model = PCompositeModel([self.p1, self.p11])
        # TODO: do we care about the format to stay the same e.g. to be able to parse it?
        # If so we should add checks here. 
        str(model)
        repr(model)
