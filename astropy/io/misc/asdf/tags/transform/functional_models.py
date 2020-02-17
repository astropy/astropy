# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from numpy.testing import assert_array_equal

from asdf import yamlutil

from astropy.modeling import functional_models
from .basic import TransformType
from . import _parameter_to_value


__all__ = ['Gaussian1DType', 'Gaussian2DType']


class Gaussian1DType(TransformType):
    name = "transform/gaussian1d"
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.Gaussian1D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.Gaussian1D(amplitude=node['amplitude'],
                                            mean=node['mean'],
                                            stddev=node['stddev'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'mean': _parameter_to_value(model.mean),
                'stddev': _parameter_to_value(model.stddev)}
        return yamlutil.custom_tree_to_tagged_tree(node, ctx)

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, functional_models.Gaussian1D) and
                isinstance(b, functional_models.Gaussian1D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.mean, b.mean)
        assert_array_equal(a.stddev, b.stddev)


class Gaussian2DType(TransformType):
    name = "transform/gaussian2d"
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.Gaussian2D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.Gaussian2D(amplitude=node['amplitude'],
                                            x_mean=node['x_mean'],
                                            y_mean=node['y_mean'],
                                            x_stddev=node['x_stddev'],
                                            y_stddev=node['y_stddev'],
                                            theta=node['theta'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'x_mean': _parameter_to_value(model.x_mean),
                'y_mean': _parameter_to_value(model.y_mean),
                'x_stddev': _parameter_to_value(model.x_stddev),
                'y_stddev': _parameter_to_value(model.y_stddev),
                'theta': _parameter_to_value(model.theta)
                }

        return yamlutil.custom_tree_to_tagged_tree(node, ctx)

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, functional_models.Gaussian2D) and
                isinstance(b, functional_models.Gaussian2D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.x_mean, b.x_mean)
        assert_array_equal(a.y_mean, b.y_mean)
        assert_array_equal(a.x_stddev, b.x_stddev)
        assert_array_equal(a.y_stddev, b.y_stddev)
        assert_array_equal(a.theta, b.theta)
