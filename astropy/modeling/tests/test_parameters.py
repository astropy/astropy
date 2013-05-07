# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests models.parameters
"""
from .. import models, fitting
from . import irafutil
from ..utils import InputParameterError
import numpy as np
from numpy.testing import utils
from ...utils.data import get_pkg_data_filename
from ...tests.helper import pytest
from .. import ParametricModel, Parameter

class TestParModel(ParametricModel):
    """
    A toy model to test parameters machinery
    """
    param_names = ['coeff', 'e']

    def __init__(self, coeff, e, param_dim=1): 
        self._coeff = Parameter(name='coeff', val=coeff, mclass=self, param_dim=param_dim)
        self._e = Parameter(name='e', val=e, mclass=self, param_dim=param_dim)
        ParametricModel.__init__(self, self.param_names, ndim=1, outdim=1, param_dim=param_dim)

    def __call__(self):
        pass
        
class TestParameters(object):
    def setup_class(self):
        """
        Unit tests for parameters
        
        Read an iraf database file created by onedspec.identify.
        Use the information to create a 1D Chebyshev model and 
        perform the same fit.
        Create also a gausian model.
        """
        test_file = get_pkg_data_filename('data/idcompspec.fits')
        f = open(test_file)
        lines = f.read()
        reclist = lines.split("begin")
        f.close()
        record = irafutil.IdentifyRecord(reclist[1])
        self.icoeff = record.coeff
        order = int(record.fields['order'])
        self.model = models.Chebyshev1DModel(order-1)
        self.gmodel = models.Gaussian1DModel(2, mean=3, stddev=4)
        self.linear_fitter = fitting.LinearLSQFitter(self.model)
        self.x = record.x
        self.y = record.z
        self.yy = np.array([record.z, record.z])
        
    def test_set_slice(self):
        """
        Tests updating the parameters attribute with a slice.
        This is what fitters internally do.
        """
        self.model.parameters[:] = np.array([3, 4, 5, 6, 7])
        assert(self.model.parameters == [3., 4., 5., 6., 7.])
       
    def test_set_parameters_as_list(self):
        """
        Tests updating parameters using a list.
        """
        self.model.parameters = [30, 40, 50, 60, 70]
        assert(self.model.parameters == [30., 40., 50., 60, 70])
        
    def test_set_parameters_as_array(self):
        """
        Tests updating parameters using an array.
        """
        self.model.parameters = np.array([3, 4, 5, 6, 7])
        assert(self.model.parameters == [3., 4., 5., 6., 7.])
        
    def test_set_as_list(self):
        """
        Parameters can be reset only by using a list or an array
        """
        with pytest.raises(TypeError):
            self.model.parameters = (1, 2, 3, 4, 5)
        
    def test_set_model_attr_seq(self):
        """
        Tests updating the parameters attribute when a model's
        parameter (in this case coeff) is updated.
        """
        self.model.parameters=[0, 0., 0., 0, 0]
        self.model.c0 = 7
        assert(self.model.parameters==[7, 0., 0., 0, 0])
    
    def test_set_model_attr_num(self):
        """
        Update the parameter list when a model's parameter is updated.
        """
        self.gmodel.amplitude = 7
        assert(self.gmodel.parameters==[7, 3, 4])
        
    def test_set_item(self):
        """
        Update the parameters using indexing.
        """
        self.model.parameters = [1, 2, 3, 4, 5]
        self.model.parameters[0] = 10.
        assert(self.model.parameters == [10, 2, 3, 4, 5])
        assert(self.model.c0 == 10)
        
    def test_wrong_size1(self):
        """
        Tests raising an error when attempting to reset the parameters
        using a list of a different size.
        """
        with pytest.raises(InputParameterError):
            self.model.parameters = [1, 2, 3]
     
    def test_wrong_size2(self):
        """
        Tests raising an exception when attemppting to update a model's
        parameter (in this case coeff) with a sequence of the wrong size.
        """
        with pytest.raises(InputParameterError):
            self.model.c0 = [1, 2, 3]

    def test_wrong_shape(self):
        """
        Tests raising an exception when attemppting to update a model's
        parameter and the new value has the wrong shape.
        """
        with pytest.raises(InputParameterError):
            self.gmodel.amplitude = [1, 2]
    
    def test_par_against_iraf(self):
        """
        Test the fitter modifies model.parameters.
        
        Uses an iraf example.
        """
        self.linear_fitter(self.x, self.y)
        utils.assert_allclose(self.model.parameters, 
                              np.array([4826.1066602783685,                                                                   952.8943813407858,
                                        12.641236013982386,
                                        -1.7910672553339604,
                                        0.90252884366711317]), 
                              rtol=10**(-2))
                            
    def testPoly1D(self):
        d={'c0':11, 'c1':12,'c2':13, 'c3':14}
        p1=models.Poly1DModel(3, **d)
        utils.assert_equal(p1.parameters, [11, 12, 13, 14])
        
    def test_poly1d_multiple_sets(self):
        p1=models.Poly1DModel(3, param_dim=3) 
        utils.assert_equal(p1.parameters, [0.0, 0.0, 0.0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0])
        utils.assert_equal(p1.c0, [0., 0, 0])
        p1.c0 = [10,10,10]
        utils.assert_equal(p1.parameters, [10.0, 10.0, 10.0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0])
    
    def test_par_slicing(self):
        """
        Test assigning to a parameter slice  
        """
        p1=models.Poly1DModel(3, param_dim=3)
        p1.c0[:2] = [10,10]
        utils.assert_equal(p1.parameters, [10.0, 10.0, 0.0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0])
        
    def test_poly2d(self):
        p2=models.Poly2DModel(degree=3)
        p2.c0_0=5
        utils.assert_equal(p2.parameters, [5, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        
    def test_poly2d_multiple_sets(self):
        kw={'c0_0':[2,3], 'c1_0':[1,2], 'c2_0':[4,5],
            'c0_1':[1,1],'c0_2':[2,2], 'c1_1':[5,5]}
        p2=models.Poly2DModel(2, **kw)
        utils.assert_equal(p2.parameters, [2, 3, 1, 2, 4, 5,
                                           1, 1, 2, 2, 5, 5])
        
    def test_non_fittable_model_parameters1d(self):
        sh1 = models.ShiftModel(2)
        sh1.offsets  = 3
        assert(sh1.offsets[0]== 3)

    def test_non_fittable_model_parametersnd(self):
        sc1 = models.ScaleModel([2,2])
        sc1.factors  = [3,3]
        assert(sc1.factors == [3,3])
    
    def test_non_fittable_model_parameters_wrong_shape(self):
        sh1 = models.ShiftModel(2)
        with pytest.raises(InputParameterError):
            sh1.offsets  = [3, 3]

class TestMultipleParameterSets(object):
    def setup_class(self):
        self.x1 = np.arange(1, 10, .1)
        self.x, self.y = np.mgrid[:10, :7]
        self.x11 = np.array([self.x1, self.x1]).T
        self.gmodel = models.Gaussian1DModel([12, 10], [3.5, 5.2], stddev=[.4, .7])
        
    def test_change_par(self):
        """
        Test that a change to one parameter as a set propagates
        to param_sets.
        """
        self.gmodel.amplitude = [1, 10]
        utils.assert_almost_equal(self.gmodel.param_sets, np.array([[1., 10], [3.5, 5.2], [0.4, 0.7]]))
        utils.assert_almost_equal(self.gmodel.parameters, [1.0, 10.0, 3.5, 5.2, 0.4, 0.7])
        
    def test_change_par2(self):
        """
        Test that a change to one single parameter in a set propagates
        to param_sets.
        """
        self.gmodel.amplitude[0] = 11
        utils.assert_almost_equal(self.gmodel.param_sets, np.array([[11., 10], [3.5, 5.2], [0.4, 0.7]]))
        utils.assert_almost_equal(self.gmodel.parameters, [11.0, 10.0, 3.5, 5.2, 0.4, 0.7])
        
    def test_change_parameters(self):
        self.gmodel.parameters=[13, 10, 9, 5.2, 0.4, 0.7]
        utils.assert_almost_equal(self.gmodel.amplitude, [13., 10.])
        utils.assert_almost_equal(self.gmodel.mean, [9., 5.2])
        
    def test_object_pars(self):
        l2 = TestParModel(coeff=[[1,2], [3,4]], e=(2,3),param_dim=2)
        utils.assert_almost_equal(l2.parameters, [1.0, 2.0, 3.0, 4.0, 2.0, 3.0])
        #utils.assert_almost_equal(l2.param_sets, np.array([[[1,2.],[3., 4.]],
         #                                               [2., 3.]], dtype=np.object))
        
    def test_wrong_number_of_pars(self):
        with pytest.raises(InputParameterError):
            l2 = TestParModel(coeff=[[1,2],[3,4]], e=(2,3,4),param_dim=2)
            
    def test_wrong_number_of_pars2(self):
        with pytest.raises(InputParameterError):
            l2 = TestParModel(coeff=[[1, 2], [3, 4]], e=4, param_dim=2)
