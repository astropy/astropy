# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from numpy.testing import assert_array_equal

from asdf import yamlutil

from astropy import modeling
from .basic import TransformType
from . import _parameter_to_value


__all__ = ['AffineType', 'Rotate2DType', 'Rotate3DType']


class AffineType(TransformType):
    name = "transform/affine"
    version = '1.2.0'
    types = ['astropy.modeling.projections.AffineTransformation2D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        matrix = node['matrix']
        translation = node['translation']
        if matrix.shape != (2, 2):
            raise NotImplementedError(
                "asdf currently only supports 2x2 (2D) rotation transformation "
                "matrices")
        if translation.shape != (2,):
            raise NotImplementedError(
                "asdf currently only supports 2D translation transformations.")

        return modeling.projections.AffineTransformation2D(
            matrix=matrix, translation=translation)

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'matrix': _parameter_to_value(model.matrix),
                'translation': _parameter_to_value(model.translation)}
        return yamlutil.custom_tree_to_tagged_tree(node, ctx)

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (a.__class__ == b.__class__)
        assert_array_equal(a.matrix, b.matrix)
        assert_array_equal(a.translation, b.translation)


class Rotate2DType(TransformType):
    name = "transform/rotate2d"
    version = '1.2.0'
    types = ['astropy.modeling.rotations.Rotation2D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return modeling.rotations.Rotation2D(node['angle'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'angle': _parameter_to_value(model.angle)}
        return yamlutil.custom_tree_to_tagged_tree(node, ctx)

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, modeling.rotations.Rotation2D) and
                isinstance(b, modeling.rotations.Rotation2D))
        assert_array_equal(a.angle, b.angle)


class Rotate3DType(TransformType):
    name = "transform/rotate3d"
    version = '1.2.0'
    types = ['astropy.modeling.rotations.RotateNative2Celestial',
             'astropy.modeling.rotations.RotateCelestial2Native',
             'astropy.modeling.rotations.EulerAngleRotation']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        if node['direction'] == 'native2celestial':
            return modeling.rotations.RotateNative2Celestial(node["phi"],
                                                             node["theta"],
                                                             node["psi"])
        elif node['direction'] == 'celestial2native':
            return modeling.rotations.RotateCelestial2Native(node["phi"],
                                                             node["theta"],
                                                             node["psi"])
        else:
            return modeling.rotations.EulerAngleRotation(node["phi"],
                                                         node["theta"],
                                                         node["psi"],
                                                         axes_order=node["direction"])


    @classmethod
    def to_tree_transform(cls, model, ctx):
        if isinstance(model, modeling.rotations.RotateNative2Celestial):
            try:
                node = {"phi": _parameter_to_value(model.lon),
                        "theta": _parameter_to_value(model.lat),
                        "psi": _parameter_to_value(model.lon_pole),
                        "direction": "native2celestial"
                        }
            except AttributeError:
                node = {"phi": model.lon,
                        "theta": model.lat,
                        "psi": model.lon_pole,
                        "direction": "native2celestial"
                        }
        elif isinstance(model, modeling.rotations.RotateCelestial2Native):
            try:
                node = {"phi": _parameter_to_value(model.lon),
                        "theta": _parameter_to_value(model.lat),
                        "psi": _parameter_to_value(model.lon_pole),
                        "direction": "celestial2native"
                        }
            except AttributeError:
                node = {"phi": model.lon,
                        "theta": model.lat,
                        "psi": model.lon_pole,
                        "direction": "celestial2native"
                        }
        else:
            node = {"phi": _parameter_to_value(model.phi),
                    "theta": _parameter_to_value(model.theta),
                    "psi": _parameter_to_value(model.psi),
                    "direction": model.axes_order
                    }

        return yamlutil.custom_tree_to_tagged_tree(node, ctx)

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert a.__class__ == b.__class__
        if a.__class__.__name__ == "EulerAngleRotation":
            assert_array_equal(a.phi, b.phi)
            assert_array_equal(a.psi, b.psi)
            assert_array_equal(a.theta, b.theta)
        else:
            assert_array_equal(a.lon, b.lon)
            assert_array_equal(a.lat, b.lat)
            assert_array_equal(a.lon_pole, b.lon_pole)


class GenericProjectionType(TransformType):
    @classmethod
    def from_tree_transform(cls, node, ctx):
        args = []
        for param_name, default in cls.params:
            args.append(node.get(param_name, default))

        if node['direction'] == 'pix2sky':
            return cls.types[0](*args)
        else:
            return cls.types[1](*args)

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {}
        if isinstance(model, cls.types[0]):
            node['direction'] = 'pix2sky'
        else:
            node['direction'] = 'sky2pix'
        for param_name, default in cls.params:
            val = getattr(model, param_name).value
            if val != default:
                node[param_name] = val
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert a.__class__ == b.__class__


_generic_projections = {
    'zenithal_perspective': ('ZenithalPerspective', (('mu', 0.0), ('gamma', 0.0)), '1.2.0'),
    'gnomonic': ('Gnomonic', (), None),
    'stereographic': ('Stereographic', (), None),
    'slant_orthographic': ('SlantOrthographic', (('xi', 0.0), ('eta', 0.0)), None),
    'zenithal_equidistant': ('ZenithalEquidistant', (), None),
    'zenithal_equal_area': ('ZenithalEqualArea', (), None),
    'airy': ('Airy', (('theta_b', 90.0),), '1.2.0'),
    'cylindrical_perspective': ('CylindricalPerspective', (('mu', 0.0), ('lam', 0.0)), '1.2.0'),
    'cylindrical_equal_area': ('CylindricalEqualArea', (('lam', 0.0),), '1.2.0'),
    'plate_carree': ('PlateCarree', (), None),
    'mercator': ('Mercator', (), None),
    'sanson_flamsteed': ('SansonFlamsteed', (), None),
    'parabolic': ('Parabolic', (), None),
    'molleweide': ('Molleweide', (), None),
    'hammer_aitoff': ('HammerAitoff', (), None),
    'conic_perspective': ('ConicPerspective', (('sigma', 0.0), ('delta', 0.0)), '1.2.0'),
    'conic_equal_area': ('ConicEqualArea', (('sigma', 0.0), ('delta', 0.0)), '1.2.0'),
    'conic_equidistant': ('ConicEquidistant', (('sigma', 0.0), ('delta', 0.0)), '1.2.0'),
    'conic_orthomorphic': ('ConicOrthomorphic', (('sigma', 0.0), ('delta', 0.0)), '1.2.0'),
    'bonne_equal_area': ('BonneEqualArea', (('theta1', 0.0),), '1.2.0'),
    'polyconic': ('Polyconic', (), None),
    'tangential_spherical_cube': ('TangentialSphericalCube', (), None),
    'cobe_quad_spherical_cube': ('COBEQuadSphericalCube', (), None),
    'quad_spherical_cube': ('QuadSphericalCube', (), None),
    'healpix': ('HEALPix', (('H', 4.0), ('X', 3.0)), None),
    'healpix_polar': ('HEALPixPolar', (), None)
}


def make_projection_types():
    for tag_name, (name, params, version) in _generic_projections.items():
        class_name = '{0}Type'.format(name)
        types = ['astropy.modeling.projections.Pix2Sky_{0}'.format(name),
                 'astropy.modeling.projections.Sky2Pix_{0}'.format(name)]

        members = {'name': 'transform/{0}'.format(tag_name),
                   'types': types,
                   'params': params}
        if version:
            members['version'] = version

        globals()[class_name] = type(
            str(class_name),
            (GenericProjectionType,),
            members)

        __all__.append(class_name)

make_projection_types()
