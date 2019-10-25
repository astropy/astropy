# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from numpy.testing import assert_array_equal

from asdf import yamlutil
from astropy.tests.helper import assert_quantity_allclose
from astropy import modeling
from astropy.modeling import models
from .basic import TransformType


__all__ = ['GratingEquationType']


class GratingEquationType(TransformType):
    name = "transform/grating_equation"
    version = '1.0.0'
    types = ['modeling.models.AnglesFromGratingEquation',
             'modeling.models.WavelengthFromGratingEquation']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        groove_density = node['groove_density']
        order = node['order']
        output = node['output']
        if output == "wavelength":
            model = WavelengthFromGratingEquation(groove_density=groove_density,
                                                  spectral_order=order)
        elif output == "angle":
            model = AngleFromGratingEquation(groove_density=groove_density,
                                             spectral_order=order)
        else:
            raise ValueError("Can't create a GratingEquation model with "
                             "output {0}".format(output))
        return model

    @classmethod
    def to_tree_transform(cls, model, ctx):
        if model.groove_density.unit is not None:
            groove_density = u.Quantity(model.groove_density.value,
                                        unit=model.groove_density.unit)
        else:
            groove_density = model.groove_density.value
        node = {'order': model.spectral_order.value,
                'groove_density': groove_density
                }
        if isinstance(model, models.AnglesFromGratingEquation3D):
            node['output'] = 'angle'
        elif isinstance(model, models.WavelengthFromGratingEquation):
            node['output'] = 'wavelength'
        else:
            raise TypeError("Can't serialize an instance of {0}"
                            .format(model.__class__.__name__))
        return yamlutil.custom_tree_to_tagged_tree(node, ctx)

    @classmethod
    def assert_equal(cls, a, b):
        TransformType.assert_equal(a, b)
        if isinstance(a, modeling.models.AnglesFromGratingEquation3D):
            assert isinstance(b, modeling.models.AnglesFromGratingEquation3D)
        elif isinstance(a, modeling.models.WavelengthFromGratingEquation):
            assert isinstance(b, modeling.models.WavelengthFromGratingEquation)
        assert_quantity_allclose(a.groove_density, b.groove_density)
        assert a.spectral_order.value == b.spectral_order.value
