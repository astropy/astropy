# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Module to test fitting routines
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)


import numpy as np

from numpy.testing.utils import assert_allclose
from ...tests.helper import pytest

try:
    from sherpa.data import DataSimulFit
    HAS_SHERPA = True
except ImportError:
    HAS_SHERPA = False

import warnings

from astropy.modeling.astro_sherpa import SherpaFitter, Dataset, ConvertedModel
from sherpa.stats import Chi2
from astropy.modeling.models import Gaussian1D, Gaussian2D

_RANDOM_SEED = 0x1337
np.random.seed(_RANDOM_SEED)


class TestSherpaFitter(object):

    def setup_class(self):
        # make data and models to use later!
        err = 0.1

        self.x1 = np.arange(1, 10, .1)
        self.dx1 = np.ones(self.x1.shape) * (.1 / 2.)
        self.x2 = np.arange(1, 10, .05)
        self.dx2 = np.ones(self.x2.shape) * (.05 / 2.)

        self.model1d = Gaussian1D(mean=5, amplitude=10, stddev=0.8)
        self.model1d_2 = Gaussian1D(mean=4, amplitude=5, stddev=0.2)

        self.tmodel1d = self.model1d.copy()
        self.tmodel1d_2 = self.model1d_2.copy()

        self.y1 = self.model1d(self.x1)
        self.y1 += err * np.random.uniform(-1., 1., size=self.y1.size)
        self.dy1 = err * np.random.uniform(0.5, 1., size=self.y1.size)

        self.y2 = self.model1d_2(self.x2)
        self.y2 += err * np.random.uniform(-1., 1., size=self.y2.size)
        self.dy2 = err * np.random.uniform(0.5, 1., size=self.y2.size)

        self.model1d.mean = 4
        self.model1d.amplitude = 6
        self.model1d.stddev = 0.5

        self.model1d_2.mean = 5
        self.model1d_2.amplitude = 10
        self.model1d_2.stddev = 0.3

        self.xx2, self.xx1 = np.mgrid[1:10:.1, 1:10:.05]
        self.shape = self.xx2.shape
        self.xx1 = self.xx1.flatten()
        self.xx2 = self.xx2.flatten()

        self.model2d = Gaussian2D(amplitude=10, x_mean=5, y_mean=6, x_stddev=0.8, y_stddev=1.1)
        self.model2d.theta.fixed = True

        self.yy = self.model2d(self.xx1, self.xx2)
        self.dxx1 = err * np.random.uniform(0.5, 1., size=self.xx1.size)
        self.dxx2 = err * np.random.uniform(0.5, 1., size=self.xx2.size)
        self.dyy = err * np.random.uniform(0.5, 1., size=self.yy.size)

        self.tmodel2d = self.model2d.copy()
        self.model2d.amplitude = 5
        self.model2d.x_mean = 6
        self.model2d.y_mean = 5
        self.model2d.x_stddev = 0.2
        self.model2d.y_stddev = 0.7

        # to stop stddev going negitive and getting div by zero error
        self.model1d.stddev.min = 1e-99
        self.model1d_2.stddev.min = 1e-99
        self.model2d.x_stddev.min = 1e-99
        self.model2d.y_stddev.min = 1e-99

        self.fitter = SherpaFitter(statistic="Chi2")

    def test_make_datasets_single(self):
        """
        Makes 1d datasets with and without errors to check it's working correctly.
        """
        data_noerr = Dataset(1, x=self.x1, y=self.y1).data
        assert_allclose(data_noerr.get_x(), self.x1)
        assert_allclose(data_noerr.get_y(), self.y1)

        data_xerr = Dataset(1, x=self.x1, y=self.y1, xerr=self.dx1).data
        assert_allclose(data_xerr.get_x(), self.x1)
        assert_allclose(data_xerr.get_y(), self.y1)
        assert_allclose(data_xerr.get_xerr(), 2 * self.dx1)

        data_yerr = Dataset(1, x=self.x1, y=self.y1, yerr=self.dy1).data
        assert_allclose(data_yerr.get_x(), self.x1)
        assert_allclose(data_yerr.get_y(), self.y1)
        assert_allclose(data_yerr.get_yerr(), self.dy1)

        data_botherr = Dataset(1, x=self.x1, y=self.y1,yerr=self.dy1, xerr=self.dx1).data
        assert_allclose(data_botherr.get_x(), self.x1)
        assert_allclose(data_botherr.get_y(), self.y1)
        assert_allclose(data_botherr.get_yerr(), self.dy1)
        assert_allclose(data_botherr.get_xerr(), 2 * self.dx1)

    def test_make_datasets_2d(self):
        """
        Makes 2d datasets with and without errors to check it's working correctly.
        """
        data_noerr = Dataset(2, x=self.xx1, y=self.xx2, z=self.yy).data
        assert_allclose(data_noerr.get_x0(), self.xx1)
        assert_allclose(data_noerr.get_x1(), self.xx2)
        assert_allclose(data_noerr.get_y(), self.yy)

        data_xyerr = Dataset(2, x=self.xx1, y=self.xx2, xerr=self.dxx1, yerr=self.dxx2, z=self.yy).data
        assert_allclose(data_xyerr.get_indep(), np.vstack([self.xx1 - self.dxx1, self.xx2 - self.dxx2, self.xx1 + self.dxx1, self.xx2 + self.dxx2]))
        assert_allclose(data_xyerr.get_y(), self.yy)

        data_xyzerr = Dataset(2, x=self.xx1, y=self.xx2, xerr=self.dxx1, yerr=self.dxx2, z=self.yy, zerr=self.dyy).data
        assert_allclose(data_xyerr.get_indep(), np.vstack([self.xx1 - self.dxx1, self.xx2 - self.dxx2, self.xx1 + self.dxx1, self.xx2 + self.dxx2]))
        assert_allclose(data_xyzerr.get_y(), self.yy)
        assert_allclose(data_xyzerr.get_yerr(), self.dyy)

    def test_make_datasets_multi(self):
        """
        Makes DataSimulFit datasets with and without errors.
        """

        data_noerr = Dataset(1, x=[self.x1, self.x2], y=[self.y1, self.y2]).data
        assert len(data_noerr.datasets) == 2
        assert isinstance(data_noerr, DataSimulFit)
        map(lambda x, y: assert_allclose(x, y), map(lambda x: x.get_x(), data_noerr.datasets), [self.x1, self.x2])
        map(lambda x, y: assert_allclose(x, y), map(lambda x: x.get_y(), data_noerr.datasets), [self.y1, self.y2])

        data_xerr = Dataset(1, x=[self.x1, self.x2], y=[self.y1, self.y2], xerr=[self.dx1, self.dx2]).data
        assert len(data_xerr.datasets) == 2
        assert isinstance(data_xerr, DataSimulFit)
        map(lambda x, y: assert_allclose(x, y), map(lambda x: x.get_x(), data_xerr.datasets), [self.x1, self.x2])
        map(lambda x, y: assert_allclose(x, y), map(lambda x: x.get_y(), data_xerr.datasets), [self.y1, self.y2])
        map(lambda x, y: assert_allclose(x, y), map(lambda x: x.get_xerr(), data_xerr.datasets), [2 * self.dx1, 2 * self.dx2])

        data_yerr = Dataset(1, x=[self.x1, self.x2], y=[self.y1, self.y2], yerr=[self.dy1, self.dy2]).data
        assert len(data_yerr.datasets) == 2
        assert isinstance(data_yerr, DataSimulFit)
        map(lambda x, y: assert_allclose(x, y), map(lambda x: x.get_x(), data_yerr.datasets), [self.x1, self.x2])
        map(lambda x, y: assert_allclose(x, y), map(lambda x: x.get_y(), data_yerr.datasets), [self.y1, self.y2])
        map(lambda x, y: assert_allclose(x, y), map(lambda x: x.get_yerr(), data_yerr.datasets), [self.dy1, self.dy2])

        data_xyerr = Dataset(1, x=[self.x1, self.x2], y=[self.y1, self.y2], yerr=[self.dy1, self.dy2], xerr=[self.dx1, self.dx2]).data
        assert len(data_xyerr.datasets) == 2
        assert isinstance(data_xyerr, DataSimulFit)
        map(lambda x, y: assert_allclose(x, y), map(lambda x: x.get_x(), data_xyerr.datasets), [self.x1, self.x2])
        map(lambda x, y: assert_allclose(x, y), map(lambda x: x.get_y(), data_xyerr.datasets), [self.y1, self.y2])
        map(lambda x, y: assert_allclose(x, y), map(lambda x: x.get_yerr(), data_xyerr.datasets), [self.dy1, self.dy2])
        map(lambda x, y: assert_allclose(x, y), map(lambda x: x.get_xerr(), data_xyerr.datasets), [2 * self.dx1, 2 * self.dx2])
    def test_convert_model_1d(self):
        """
        Check that the parameter constraints (Tied,Fixed and bounds)
        are correctly converted to the sherpa models.
        """

        def tiefunc(mo):
            return mo.stddev

        amodel = self.model1d.copy()
        amodel.amplitude.fixed = True
        amodel.stddev.max = 10.0
        amodel.mean.tied = tiefunc
        smodel = ConvertedModel(amodel).sherpa_model

        assert smodel.amplitude.frozen
        assert smodel.stddev.max == 10.0
        assert smodel.mean.link == smodel.stddev
        amodel.mean.value = amodel.stddev.value
        assert_allclose(smodel(self.x1), amodel(self.x1))

    def test_convert_model_2d(self):
        """
        Check that the parameter constraints (Tied,Fixed and bounds)
        are correctly converted to the sherpa models.
        """
        def tiefunc(mo):
            return mo.x_stddev

        amodel = self.model2d.copy()
        amodel.amplitude.fixed = True
        amodel.x_stddev.max = 10.0
        amodel.y_stddev.tied = tiefunc
        smodel = ConvertedModel(amodel).sherpa_model

        assert smodel.amplitude.frozen
        assert smodel.x_stddev.max == 10.0
        assert smodel.y_stddev.link == smodel.x_stddev

        amodel.y_stddev.value = amodel.x_stddev.value
        assert_allclose(smodel(self.xx1, self.xx2), amodel(self.xx1, self.xx2))

    def test_single_dataset_single_model(self):
        """Test a single model with a single dataset."""
        fmod = self.fitter(self.model1d.copy(), self.x1, self.y1, yerr=self.dy1)
        for pp in fmod.param_names:
            assert_allclose(getattr(fmod, pp), getattr(self.tmodel1d, pp), rtol=0.05)

    def test_single_dataset_two_models(self):
        """Test a two models with a single dataset."""
        fmod = self.fitter([self.model1d.copy(), self.model1d.copy()], self.x1, self.y1, yerr=self.dy1)
        for ff in fmod:
            for nn, pp in enumerate(ff.param_names):
                assert_allclose(getattr(ff, pp), getattr(self.tmodel1d, pp), rtol=0.05)

    def test_two_dataset_single_model(self):
        """Tests two datasets with a single model."""
        fmod = self.fitter(self.model1d.copy(), [self.x1, self.x1], [
                           self.y1, self.y1], yerr=[self.dy1, self.dy1])
        for ff in fmod:
            for nn, pp in enumerate(ff.param_names):
                assert_allclose(getattr(ff, pp), getattr(
                    self.tmodel1d, pp), rtol=0.05)

    def test_two_dataset_two_models(self):
        """Tests two models with two datasets simultaneously."""
        fmod = self.fitter([self.model1d.copy(), self.model1d_2.copy()], [self.x1, self.x2], [self.y1, self.y2], yerr=[self.dy1, self.dy2])
        for ff, mm in zip(fmod, [self.tmodel1d, self.tmodel1d_2]):
            for nn, pp in enumerate(ff.param_names):
                assert_allclose(getattr(ff, pp), getattr(mm, pp), rtol=0.05)

    def test_2d_fit(self):
        """
        Fits a 2d dataset.
        """
        fmod = self.fitter(self.model2d.copy(), self.xx1, self.xx2, self.yy, zerr=self.dyy)
        for pp in fmod.param_names:
            if getattr(self.tmodel2d, pp).fixed is False:
                assert_allclose(getattr(fmod, pp), getattr(self.tmodel2d, pp), rtol=0.05)

    def test_condition_fit(self):
        """
        Some deliberatly bad fits to check that contraints are working!
        """
        def tiefunc(mo):
            return mo.stddev

        amodel = self.model1d.copy()
        amodel.amplitude.fixed = True
        amodel.stddev.max = 0.7
        amodel.stddev.value = 0.6
        amodel.stddev.min = 0.5
        amodel.mean.tied = tiefunc
        with warnings.catch_warnings():
            # Sherpa uses something which throws this warning
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            fmod = self.fitter(amodel, self.x1, self.y1, yerr=self.dy1)

        assert fmod.amplitude.value == self.model1d.amplitude.value
        assert round(fmod.stddev.value, 1) == 0.7 or round(
            fmod.stddev.value, 1) == 0.5
        assert fmod.mean.value == fmod.stddev.value

    def test_error_methods(self):
        pass

    def test_opt_methods(self):
        pass

    def test_bkg_doesnt_explode(self):
        from astropy.modeling.models import Polynomial1D
        m = Polynomial1D(2)

        x = np.arange(0, 10, 0.1)
        y = 2 + 3 * x**2 + 0.5 * x
        bkg = x

        sfit = SherpaFitter(statistic="cash", estmethod='covariance')
        sfit(m, x, y, bkg=bkg)
        #TODO: Make this better!
