# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from numpy.testing import assert_array_equal


from astropy.modeling import functional_models, physical_models
from .basic import TransformType
from . import _parameter_to_value


__all__ = ['BlackBody', 'Drude1DType', 'Plummer1DType']


class BlackBody(TransformType):
    name = 'transform/blackbody'
    version = '1.0.0'
    types = ['astropy.modeling.physical_models.BlackBody']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return physical_models.BlackBody(scale=node['scale'],
                                         temperature=node['temperature'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'scale': _parameter_to_value(model.scale),
                'temperature': _parameter_to_value(model.temperature)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, physical_models.BlackBody) and
                isinstance(b, physical_models.BlackBody))
        assert_array_equal(a.scale, b.scale)
        assert_array_equal(a.temperature, b.temperature)


class Drude1DType(TransformType):
    name = 'transform/drude1d'
    version = '1.0.0'
    types = ['astropy.modeling.physical_models.Drude1D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return physical_models.Drude1D(amplitude=node['amplitude'],
                                       x_0=node['x_0'],
                                       fwhm=node['fwhm'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'x_0': _parameter_to_value(model.x_0),
                'fwhm': _parameter_to_value(model.fwhm)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, physical_models.Drude1D) and
                isinstance(b, physical_models.Drude1D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.x_0, b.x_0)
        assert_array_equal(a.fwhm, b.fwhm)


class Plummer1DType(TransformType):
    name = 'transform/plummer1d'
    version = '1.0.0'
    types = ['astropy.modeling.physical_models.Plummer1D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return physical_models.Plummer1D(mass=node['mass'],
                                         r_plum=node['r_plum'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'mass': _parameter_to_value(model.mass),
                'r_plum': _parameter_to_value(model.r_plum)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, physical_models.Plummer1D) and
                isinstance(b, physical_models.Plummer1D))
        assert_array_equal(a.mass, b.mass)
        assert_array_equal(a.r_plum, b.r_plum)
