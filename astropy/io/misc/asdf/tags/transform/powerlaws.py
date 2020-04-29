# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from numpy.testing import assert_array_equal

from astropy.modeling import powerlaws
from .basic import TransformType
from . import _parameter_to_value


__all__ = ['PowerLaw1DType', 'BrokenPowerLaw1DType',
           'SmoothlyBrokenPowerLaw1DType', 'ExponentialCutoffPowerLaw1DType',
           'LogParabola1DType']


class PowerLaw1DType(TransformType):
    name = 'transform/power_law1d'
    version = '1.0.0'
    types = ['astropy.modeling.powerlaws.PowerLaw1D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return powerlaws.PowerLaw1D(amplitude=node['amplitude'],
                                    x_0=node['x_0'],
                                    alpha=node['alpha'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'x_0': _parameter_to_value(model.x_0),
                'alpha': _parameter_to_value(model.alpha)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, powerlaws.PowerLaw1D) and
                isinstance(b, powerlaws.PowerLaw1D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.x_0, b.x_0)
        assert_array_equal(a.alpha, b.alpha)


class BrokenPowerLaw1DType(TransformType):
    name = 'transform/broken_power_law1d'
    version = '1.0.0'
    types = ['astropy.modeling.powerlaws.BrokenPowerLaw1D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return powerlaws.BrokenPowerLaw1D(amplitude=node['amplitude'],
                                          x_break=node['x_break'],
                                          alpha_1=node['alpha_1'],
                                          alpha_2=node['alpha_2'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'x_break': _parameter_to_value(model.x_break),
                'alpha_1': _parameter_to_value(model.alpha_1),
                'alpha_2': _parameter_to_value(model.alpha_2)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, powerlaws.BrokenPowerLaw1D) and
                isinstance(b, powerlaws.BrokenPowerLaw1D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.x_break, b.x_break)
        assert_array_equal(a.alpha_1, b.alpha_1)
        assert_array_equal(a.alpha_2, b.alpha_2)


class SmoothlyBrokenPowerLaw1DType(TransformType):
    name = 'transform/smoothly_broken_power_law1d'
    version = '1.0.0'
    types = ['astropy.modeling.powerlaws.SmoothlyBrokenPowerLaw1D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return powerlaws.SmoothlyBrokenPowerLaw1D(amplitude=node['amplitude'],
                                                  x_break=node['x_break'],
                                                  alpha_1=node['alpha_1'],
                                                  alpha_2=node['alpha_2'],
                                                  delta=node['delta'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'x_break': _parameter_to_value(model.x_break),
                'alpha_1': _parameter_to_value(model.alpha_1),
                'alpha_2': _parameter_to_value(model.alpha_2),
                'delta': _parameter_to_value(model.delta)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, powerlaws.SmoothlyBrokenPowerLaw1D) and
                isinstance(b, powerlaws.SmoothlyBrokenPowerLaw1D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.x_break, b.x_break)
        assert_array_equal(a.alpha_1, b.alpha_1)
        assert_array_equal(a.alpha_2, b.alpha_2)
        assert_array_equal(a.delta, b.delta)


class ExponentialCutoffPowerLaw1DType(TransformType):
    name = 'transform/exponential_cutoff_power_law1d'
    version = '1.0.0'
    types = ['astropy.modeling.powerlaws.ExponentialCutoffPowerLaw1D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return powerlaws.ExponentialCutoffPowerLaw1D(amplitude=node['amplitude'],
                                                     x_0=node['x_0'],
                                                     alpha=node['alpha'],
                                                     x_cutoff=node['x_cutoff'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'x_0': _parameter_to_value(model.x_0),
                'alpha': _parameter_to_value(model.alpha),
                'x_cutoff': _parameter_to_value(model.x_cutoff)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, powerlaws.ExponentialCutoffPowerLaw1D) and
                isinstance(b, powerlaws.ExponentialCutoffPowerLaw1D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.x_0, b.x_0)
        assert_array_equal(a.alpha, b.alpha)
        assert_array_equal(a.x_cutoff, b.x_cutoff)


class LogParabola1DType(TransformType):
    name = 'transform/log_parabola1d'
    version = '1.0.0'
    types = ['astropy.modeling.powerlaws.LogParabola1D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return powerlaws.LogParabola1D(amplitude=node['amplitude'],
                                       x_0=node['x_0'],
                                       alpha=node['alpha'],
                                       beta=node['beta'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'x_0': _parameter_to_value(model.x_0),
                'alpha': _parameter_to_value(model.alpha),
                'beta': _parameter_to_value(model.beta)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, powerlaws.LogParabola1D) and
                isinstance(b, powerlaws.LogParabola1D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.x_0, b.x_0)
        assert_array_equal(a.alpha, b.alpha)
        assert_array_equal(a.beta, b.beta)
