# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from numpy.testing import assert_array_equal

from asdf.versioning import AsdfVersion

import astropy.units as u
from astropy import modeling
from astropy.io.misc.asdf.tags.transform.basic import TransformType
from . import _parameter_to_value

__all__ = ['ShiftType', 'ScaleType', 'Linear1DType']


class ShiftType(TransformType):
    name = "transform/shift"
    version = '1.2.0'
    types = ['astropy.modeling.models.Shift']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        offset = node['offset']
        if not isinstance(offset, u.Quantity) and not np.isscalar(offset):
            raise NotImplementedError(
                "Asdf currently only supports scalar inputs to Shift transform.")

        return modeling.models.Shift(offset)

    @classmethod
    def to_tree_transform(cls, model, ctx):
        offset = model.offset
        return {'offset': _parameter_to_value(offset)}

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, modeling.models.Shift) and
                isinstance(b, modeling.models.Shift))
        assert_array_equal(a.offset.value, b.offset.value)


class ScaleType(TransformType):
    name = "transform/scale"
    version = '1.2.0'
    types = ['astropy.modeling.models.Scale']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        factor = node['factor']
        if not isinstance(factor, u.Quantity) and not np.isscalar(factor):
            raise NotImplementedError(
                "Asdf currently only supports scalar inputs to Scale transform.")

        return modeling.models.Scale(factor)

    @classmethod
    def to_tree_transform(cls, model, ctx):
        factor = model.factor
        return {'factor': _parameter_to_value(factor)}

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, modeling.models.Scale) and
                isinstance(b, modeling.models.Scale))
        assert_array_equal(a.factor, b.factor)


class MultiplyType(TransformType):
    name = "transform/multiplyscale"
    version = '1.0.0'
    types = ['astropy.modeling.models.Multiply']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        factor = node['factor']
        return modeling.models.Multiply(factor)

    @classmethod
    def to_tree_transform(cls, model, ctx):
        factor = model.factor
        return {'factor': _parameter_to_value(factor)}

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, modeling.models.Multiply) and
                isinstance(b, modeling.models.Multiply))
        assert_array_equal(a.factor, b.factor)


class PolynomialTypeBase(TransformType):
    DOMAIN_WINDOW_MIN_VERSION = AsdfVersion("1.2.0")

    name = "transform/polynomial"
    types = ['astropy.modeling.models.Polynomial1D',
             'astropy.modeling.models.Polynomial2D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        coefficients = np.asarray(node['coefficients'])
        n_dim = coefficients.ndim

        if n_dim == 1:
            domain = node.get('domain', None)
            window = node.get('window', None)

            model = modeling.models.Polynomial1D(coefficients.size - 1,
                                                 domain=domain, window=window)
            model.parameters = coefficients
        elif n_dim == 2:
            x_domain, y_domain = tuple(node.get('domain', (None, None)))
            x_window, y_window = tuple(node.get('window', (None, None)))
            shape = coefficients.shape
            degree = shape[0] - 1
            if shape[0] != shape[1]:
                raise TypeError("Coefficients must be an (n+1, n+1) matrix")

            coeffs = {}
            for i in range(shape[0]):
                for j in range(shape[0]):
                    if i + j < degree + 1:
                        name = 'c' + str(i) + '_' + str(j)
                        coeffs[name] = coefficients[i, j]
            model = modeling.models.Polynomial2D(degree,
                                                 x_domain=x_domain,
                                                 y_domain=y_domain,
                                                 x_window=x_window,
                                                 y_window=y_window,
                                                 **coeffs)
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
                        name = 'c' + str(i) + '_' + str(j)
                        coefficients[i, j] = getattr(model, name).value
        node = {'coefficients': coefficients}
        typeindex = cls.types.index(model.__class__)
        ndim = (typeindex % 2) + 1

        if cls.version >= PolynomialTypeBase.DOMAIN_WINDOW_MIN_VERSION:
            # Schema versions prior to 1.2 included an unrelated "domain"
            # property.  We can't serialize the new domain values with those
            # versions because they don't validate.
            if ndim == 1:
                if model.domain is not None:
                    node['domain'] = model.domain
                if model.window is not None:
                    node['window'] = model.window
            else:
                if model.x_domain or model.y_domain is not None:
                    node['domain'] = (model.x_domain, model.y_domain)
                if model.x_window or model.y_window is not None:
                    node['window'] = (model.x_window, model.y_window)

        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, (modeling.models.Polynomial1D, modeling.models.Polynomial2D)) and
                isinstance(b, (modeling.models.Polynomial1D, modeling.models.Polynomial2D)))
        assert_array_equal(a.parameters, b.parameters)

        if cls.version > PolynomialTypeBase.DOMAIN_WINDOW_MIN_VERSION:
            # Schema versions prior to 1.2 are known not to serialize
            # domain or window.
            if isinstance(a, modeling.models.Polynomial1D):
                assert a.domain == b.domain
                assert a.window == b.window
            else:
                assert a.x_domain == b.x_domain
                assert a.x_window == b.x_window
                assert a.y_domain == b.y_domain
                assert a.y_window == b.y_window


class PolynomialType1_0(PolynomialTypeBase):
    version = "1.0.0"


class PolynomialType1_1(PolynomialTypeBase):
    version = "1.1.0"


class PolynomialType1_2(PolynomialTypeBase):
    version = "1.2.0"


class OrthoPolynomialType(TransformType):
    name = "transform/ortho_polynomial"
    types = ['astropy.modeling.models.Legendre1D',
             'astropy.modeling.models.Legendre2D',
             'astropy.modeling.models.Chebyshev1D',
             'astropy.modeling.models.Chebyshev2D',
             'astropy.modeling.models.Hermite1D',
             'astropy.modeling.models.Hermite2D']
    typemap = {
        'legendre': 0,
        'chebyshev': 2,
        'hermite': 4,
    }

    invtypemap = dict([[v, k] for k, v in typemap.items()])

    version = "1.0.0"

    @classmethod
    def from_tree_transform(cls, node, ctx):
        coefficients = np.asarray(node['coefficients'])
        n_dim = coefficients.ndim
        poly_type = node['polynomial_type']
        if n_dim == 1:
            domain = node.get('domain', None)
            window = node.get('window', None)
            model = cls.types[cls.typemap[poly_type]](coefficients.size - 1,
                                                      domain=domain, window=window)
            model.parameters = coefficients
        elif n_dim == 2:
            x_domain, y_domain = tuple(node.get('domain', (None, None)))
            x_window, y_window = tuple(node.get('window', (None, None)))
            coeffs = {}
            shape = coefficients.shape
            x_degree = shape[0] - 1
            y_degree = shape[1] - 1
            for i in range(x_degree + 1):
                for j in range(y_degree + 1):
                    name = f'c{i}_{j}'
                    coeffs[name] = coefficients[i, j]
            model = cls.types[cls.typemap[poly_type]+1](x_degree, y_degree,
                                                        x_domain=x_domain,
                                                        y_domain=y_domain,
                                                        x_window=x_window,
                                                        y_window=y_window,
                                                        **coeffs)
        else:
            raise NotImplementedError(
                "Asdf currently only supports 1D or 2D polynomial transforms.")
        return model

    @classmethod
    def to_tree_transform(cls, model, ctx):
        typeindex = cls.types.index(model.__class__)
        poly_type = cls.invtypemap[int(typeindex/2)*2]
        ndim = (typeindex % 2) + 1
        if ndim == 1:
            coefficients = np.array(model.parameters)
        else:
            coefficients = np.zeros((model.x_degree + 1, model.y_degree + 1))
            for i in range(model.x_degree + 1):
                for j in range(model.y_degree + 1):
                    name = f'c{i}_{j}'
                    coefficients[i, j] = getattr(model, name).value
        node = {'polynomial_type': poly_type, 'coefficients': coefficients}
        if ndim == 1:
            if model.domain is not None:
                node['domain'] = model.domain
            if model.window is not None:
                node['window'] = model.window
        else:
            if model.x_domain or model.y_domain is not None:
                node['domain'] = (model.x_domain, model.y_domain)
            if model.x_window or model.y_window is not None:
                node['window'] = (model.x_window, model.y_window)
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        # There should be a more elegant way of doing this
        TransformType.assert_equal(a, b)
        assert ((isinstance(a, (modeling.models.Legendre1D,   modeling.models.Legendre2D)) and
                 isinstance(b, (modeling.models.Legendre1D,   modeling.models.Legendre2D))) or
                (isinstance(a, (modeling.models.Chebyshev1D,  modeling.models.Chebyshev2D)) and
                 isinstance(b, (modeling.models.Chebyshev1D,  modeling.models.Chebyshev2D))) or
                (isinstance(a, (modeling.models.Hermite1D,    modeling.models.Hermite2D)) and
                 isinstance(b, (modeling.models.Hermite1D,    modeling.models.Hermite2D))))
        assert_array_equal(a.parameters, b.parameters)


class Linear1DType(TransformType):
    name = "transform/linear1d"
    version = '1.0.0'
    types = ['astropy.modeling.models.Linear1D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        slope = node.get('slope', None)
        intercept = node.get('intercept', None)

        return modeling.models.Linear1D(slope=slope, intercept=intercept)

    @classmethod
    def to_tree_transform(cls, model, ctx):
        return {
            'slope': _parameter_to_value(model.slope),
            'intercept': _parameter_to_value(model.intercept),
        }

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, modeling.models.Linear1D) and
                isinstance(b, modeling.models.Linear1D))
        assert_array_equal(a.slope, b.slope)
        assert_array_equal(a.intercept, b.intercept)
