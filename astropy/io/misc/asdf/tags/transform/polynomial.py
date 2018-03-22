# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import numpy as np
from numpy.testing import assert_array_equal

from asdf import yamlutil

from astropy import modeling
from .basic import TransformType


__all__ = ['ShiftType', 'ScaleType', 'PolynomialType']


class ShiftType(TransformType):
    name = "transform/shift"
    types = ['astropy.modeling.models.Shift']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        offset = node['offset']
        if not np.isscalar(offset):
            raise NotImplementedError(
                "Asdf currently only supports scalar inputs to Shift transform.")

        return modeling.models.Shift(offset)

    @classmethod
    def to_tree_transform(cls, model, ctx):
        return {'offset': model.offset.value}
        #return yamlutil.custom_tree_to_tagged_tree(node, ctx)

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, modeling.models.Shift) and
                isinstance(b, modeling.models.Shift))
        assert_array_equal(a.offset.value, b.offset.value)


class ScaleType(TransformType):
    name = "transform/scale"
    types = ['astropy.modeling.models.Scale']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        factor = node['factor']
        if not np.isscalar(factor):
            raise NotImplementedError(
                "Asdf currently only supports scalar inputs to Scale transform.")

        return modeling.models.Scale(factor)

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'factor': model.factor.value}
        return yamlutil.custom_tree_to_tagged_tree(node, ctx)

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, modeling.models.Scale) and
                isinstance(b, modeling.models.Scale))
        assert_array_equal(a.factor, b.factor)


class PolynomialType(TransformType):
    name = "transform/polynomial"
    types = ['astropy.modeling.models.Polynomial1D',
             'astropy.modeling.models.Polynomial2D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        coefficients = np.asarray(node['coefficients'])
        n_dim = coefficients.ndim

        if n_dim == 1:
            model = modeling.models.Polynomial1D(coefficients.size - 1)
            model.parameters = coefficients
        elif n_dim == 2:
            shape = coefficients.shape
            degree = shape[0] - 1
            if shape[0] != shape[1]:
                raise TypeError("Coefficients must be an (n+1, n+1) matrix")

            coeffs = {}
            for i in range(shape[0]):
                for j in range(shape[0]):
                    if i + j < degree + 1:
                        name = 'c' + str(i) + '_' +str(j)
                        coeffs[name] = coefficients[i, j]
            model = modeling.models.Polynomial2D(degree, **coeffs)
        else:
            raise NotImplementedError(
                "Asdf currently only supports 1D or 2D polynomial transform.")
        return model

    @classmethod
    def to_tree_transform(cls, model, ctx):
        if isinstance(model, modeling.models.Polynomial1D):
            coefficients = np.array(model.parameters)
        elif isinstance(model, modeling.models.Polynomial2D):
            degree = model.degree
            coefficients = np.zeros((degree + 1, degree + 1))
            for i in range(degree + 1):
                for j in range(degree + 1):
                    if i + j < degree + 1:
                        name = 'c' + str(i) + '_' +str(j)
                        coefficients[i, j] = getattr(model, name).value
        node = {'coefficients': coefficients}
        return yamlutil.custom_tree_to_tagged_tree(node, ctx)

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, (modeling.models.Polynomial1D, modeling.models.Polynomial2D)) and
                isinstance(b, (modeling.models.Polynomial1D, modeling.models.Polynomial2D)))
        assert_array_equal(a.parameters, b.parameters)
