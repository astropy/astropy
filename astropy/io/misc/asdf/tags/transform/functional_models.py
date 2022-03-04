# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from numpy.testing import assert_array_equal

from astropy.modeling import functional_models
from astropy.io.misc.asdf.tags.transform.basic import TransformType
from . import _parameter_to_value


__all__ = ['AiryDisk2DType', 'Box1DType', 'Box2DType',
           'Disk2DType', 'Ellipse2DType', 'Exponential1DType',
           'Gaussian1DType', 'Gaussian2DType', 'KingProjectedAnalytic1DType',
           'Logarithmic1DType', 'Lorentz1DType', 'Moffat1DType',
           'Moffat2DType', 'Planar2D', 'RedshiftScaleFactorType',
           'RickerWavelet1DType', 'RickerWavelet2DType', 'Ring2DType',
           'Sersic1DType', 'Sersic2DType',
           'Sine1DType', 'Cosine1DType', 'Tangent1DType',
           'ArcSine1DType', 'ArcCosine1DType', 'ArcTangent1DType',
           'Trapezoid1DType', 'TrapezoidDisk2DType', 'Voigt1DType']


class AiryDisk2DType(TransformType):
    name = 'transform/airy_disk2d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.AiryDisk2D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.AiryDisk2D(amplitude=node['amplitude'],
                                            x_0=node['x_0'],
                                            y_0=node['y_0'],
                                            radius=node['radius'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'x_0': _parameter_to_value(model.x_0),
                'y_0': _parameter_to_value(model.y_0),
                'radius': _parameter_to_value(model.radius)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, functional_models.AiryDisk2D) and
                isinstance(b, functional_models.AiryDisk2D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.x_0, b.x_0)
        assert_array_equal(a.y_0, b.y_0)
        assert_array_equal(a.radius, b.radius)


class Box1DType(TransformType):
    name = 'transform/box1d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.Box1D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.Box1D(amplitude=node['amplitude'],
                                       x_0=node['x_0'],
                                       width=node['width'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'x_0': _parameter_to_value(model.x_0),
                'width': _parameter_to_value(model.width)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, functional_models.Box1D) and
                isinstance(b, functional_models.Box1D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.x_0, b.x_0)
        assert_array_equal(a.width, b.width)


class Box2DType(TransformType):
    name = 'transform/box2d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.Box2D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.Box2D(amplitude=node['amplitude'],
                                       x_0=node['x_0'],
                                       x_width=node['x_width'],
                                       y_0=node['y_0'],
                                       y_width=node['y_width'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'x_0': _parameter_to_value(model.x_0),
                'x_width': _parameter_to_value(model.x_width),
                'y_0': _parameter_to_value(model.y_0),
                'y_width': _parameter_to_value(model.y_width)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, functional_models.Box2D) and
                isinstance(b, functional_models.Box2D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.x_0, b.x_0)
        assert_array_equal(a.x_width, b.x_width)
        assert_array_equal(a.y_0, b.y_0)
        assert_array_equal(a.y_width, b.y_width)


class Disk2DType(TransformType):
    name = 'transform/disk2d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.Disk2D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.Disk2D(amplitude=node['amplitude'],
                                        x_0=node['x_0'],
                                        y_0=node['y_0'],
                                        R_0=node['R_0'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'x_0': _parameter_to_value(model.x_0),
                'y_0': _parameter_to_value(model.y_0),
                'R_0': _parameter_to_value(model.R_0)}

        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, functional_models.Disk2D) and
                isinstance(b, functional_models.Disk2D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.x_0, b.x_0)
        assert_array_equal(a.y_0, b.y_0)
        assert_array_equal(a.R_0, b.R_0)


class Ellipse2DType(TransformType):
    name = 'transform/ellipse2d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.Ellipse2D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.Ellipse2D(amplitude=node['amplitude'],
                                           x_0=node['x_0'],
                                           y_0=node['y_0'],
                                           a=node['a'],
                                           b=node['b'],
                                           theta=node['theta'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'x_0': _parameter_to_value(model.x_0),
                'y_0': _parameter_to_value(model.y_0),
                'a': _parameter_to_value(model.a),
                'b': _parameter_to_value(model.b),
                'theta': _parameter_to_value(model.theta)}

        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, functional_models.Ellipse2D) and
                isinstance(b, functional_models.Ellipse2D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.x_0, b.x_0)
        assert_array_equal(a.y_0, b.y_0)
        assert_array_equal(a.a, b.a)
        assert_array_equal(a.b, b.b)
        assert_array_equal(a.theta, b.theta)


class Exponential1DType(TransformType):
    name = 'transform/exponential1d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.Exponential1D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.Exponential1D(amplitude=node['amplitude'],
                                               tau=node['tau'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'tau': _parameter_to_value(model.tau)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, functional_models.Exponential1D) and
                isinstance(b, functional_models.Exponential1D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.tau, b.tau)


class Gaussian1DType(TransformType):
    name = 'transform/gaussian1d'
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
        return node

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
    name = 'transform/gaussian2d'
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
                'theta': _parameter_to_value(model.theta)}

        return node

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


class KingProjectedAnalytic1DType(TransformType):
    name = 'transform/king_projected_analytic1d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.KingProjectedAnalytic1D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.KingProjectedAnalytic1D(
                                            amplitude=node['amplitude'],
                                            r_core=node['r_core'],
                                            r_tide=node['r_tide'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'r_core': _parameter_to_value(model.r_core),
                'r_tide': _parameter_to_value(model.r_tide)}

        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, functional_models.KingProjectedAnalytic1D) and
                isinstance(b, functional_models.KingProjectedAnalytic1D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.r_core, b.r_core)
        assert_array_equal(a.r_tide, b.r_tide)


class Logarithmic1DType(TransformType):
    name = 'transform/logarithmic1d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.Logarithmic1D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.Logarithmic1D(amplitude=node['amplitude'],
                                               tau=node['tau'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'tau': _parameter_to_value(model.tau)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, functional_models.Logarithmic1D) and
                isinstance(b, functional_models.Logarithmic1D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.tau, b.tau)


class Lorentz1DType(TransformType):
    name = 'transform/lorentz1d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.Lorentz1D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.Lorentz1D(amplitude=node['amplitude'],
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
        assert (isinstance(a, functional_models.Lorentz1D) and
                isinstance(b, functional_models.Lorentz1D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.x_0, b.x_0)
        assert_array_equal(a.fwhm, b.fwhm)


class Moffat1DType(TransformType):
    name = 'transform/moffat1d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.Moffat1D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.Moffat1D(amplitude=node['amplitude'],
                                          x_0=node['x_0'],
                                          gamma=node['gamma'],
                                          alpha=node['alpha'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'x_0': _parameter_to_value(model.x_0),
                'gamma': _parameter_to_value(model.gamma),
                'alpha': _parameter_to_value(model.alpha)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, functional_models.Moffat1D) and
                isinstance(b, functional_models.Moffat1D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.x_0, b.x_0)
        assert_array_equal(a.gamma, b.gamma)
        assert_array_equal(a.alpha, b.alpha)


class Moffat2DType(TransformType):
    name = 'transform/moffat2d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.Moffat2D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.Moffat2D(amplitude=node['amplitude'],
                                          x_0=node['x_0'],
                                          y_0=node['y_0'],
                                          gamma=node['gamma'],
                                          alpha=node['alpha'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'x_0': _parameter_to_value(model.x_0),
                'y_0': _parameter_to_value(model.y_0),
                'gamma': _parameter_to_value(model.gamma),
                'alpha': _parameter_to_value(model.alpha)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, functional_models.Moffat2D) and
                isinstance(b, functional_models.Moffat2D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.x_0, b.x_0)
        assert_array_equal(a.y_0, b.y_0)
        assert_array_equal(a.gamma, b.gamma)
        assert_array_equal(a.alpha, b.alpha)


class Planar2D(TransformType):
    name = 'transform/planar2d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.Planar2D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.Planar2D(slope_x=node['slope_x'],
                                          slope_y=node['slope_y'],
                                          intercept=node['intercept'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'slope_x': _parameter_to_value(model.slope_x),
                'slope_y': _parameter_to_value(model.slope_y),
                'intercept': _parameter_to_value(model.intercept)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, functional_models.Planar2D) and
                isinstance(b, functional_models.Planar2D))
        assert_array_equal(a.slope_x, b.slope_x)
        assert_array_equal(a.slope_y, b.slope_y)
        assert_array_equal(a.intercept, b.intercept)


class RedshiftScaleFactorType(TransformType):
    name = 'transform/redshift_scale_factor'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.RedshiftScaleFactor']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.RedshiftScaleFactor(z=node['z'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'z': _parameter_to_value(model.z)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, functional_models.RedshiftScaleFactor) and
                isinstance(b, functional_models.RedshiftScaleFactor))
        assert_array_equal(a.z, b.z)


class RickerWavelet1DType(TransformType):
    name = 'transform/ricker_wavelet1d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.RickerWavelet1D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.RickerWavelet1D(amplitude=node['amplitude'],
                                                 x_0=node['x_0'],
                                                 sigma=node['sigma'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'x_0': _parameter_to_value(model.x_0),
                'sigma': _parameter_to_value(model.sigma)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, functional_models.RickerWavelet1D) and
                isinstance(b, functional_models.RickerWavelet1D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.x_0, b.x_0)
        assert_array_equal(a.sigma, b.sigma)


class RickerWavelet2DType(TransformType):
    name = 'transform/ricker_wavelet2d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.RickerWavelet2D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.RickerWavelet2D(amplitude=node['amplitude'],
                                                 x_0=node['x_0'],
                                                 y_0=node['y_0'],
                                                 sigma=node['sigma'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'x_0': _parameter_to_value(model.x_0),
                'y_0': _parameter_to_value(model.y_0),
                'sigma': _parameter_to_value(model.sigma)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, functional_models.RickerWavelet2D) and
                isinstance(b, functional_models.RickerWavelet2D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.x_0, b.x_0)
        assert_array_equal(a.y_0, b.y_0)
        assert_array_equal(a.sigma, b.sigma)


class Ring2DType(TransformType):
    name = 'transform/ring2d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.Ring2D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.Ring2D(amplitude=node['amplitude'],
                                        x_0=node['x_0'],
                                        y_0=node['y_0'],
                                        r_in=node['r_in'],
                                        width=node['width'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'x_0': _parameter_to_value(model.x_0),
                'y_0': _parameter_to_value(model.y_0),
                'r_in': _parameter_to_value(model.r_in),
                'width': _parameter_to_value(model.width)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, functional_models.Ring2D) and
                isinstance(b, functional_models.Ring2D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.x_0, b.x_0)
        assert_array_equal(a.y_0, b.y_0)
        assert_array_equal(a.r_in, b.r_in)
        assert_array_equal(a.width, b.width)


class Sersic1DType(TransformType):
    name = 'transform/sersic1d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.Sersic1D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.Sersic1D(amplitude=node['amplitude'],
                                          r_eff=node['r_eff'],
                                          n=node['n'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'r_eff': _parameter_to_value(model.r_eff),
                'n': _parameter_to_value(model.n)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, functional_models.Sersic1D) and
                isinstance(b, functional_models.Sersic1D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.r_eff, b.r_eff)
        assert_array_equal(a.n, b.n)


class Sersic2DType(TransformType):
    name = 'transform/sersic2d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.Sersic2D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.Sersic2D(amplitude=node['amplitude'],
                                          r_eff=node['r_eff'],
                                          n=node['n'],
                                          x_0=node['x_0'],
                                          y_0=node['y_0'],
                                          ellip=node['ellip'],
                                          theta=node['theta'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'r_eff': _parameter_to_value(model.r_eff),
                'n': _parameter_to_value(model.n),
                'x_0': _parameter_to_value(model.x_0),
                'y_0': _parameter_to_value(model.y_0),
                'ellip': _parameter_to_value(model.ellip),
                'theta': _parameter_to_value(model.theta)

                }
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, functional_models.Sersic2D) and
                isinstance(b, functional_models.Sersic2D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.r_eff, b.r_eff)
        assert_array_equal(a.n, b.n)
        assert_array_equal(a.x_0, b.x_0)
        assert_array_equal(a.y_0, b.y_0)
        assert_array_equal(a.ellip, b.ellip)
        assert_array_equal(a.theta, b.theta)


class Trigonometric1DType(TransformType):
    _model = None

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return cls._model(amplitude=node['amplitude'],
                          frequency=node['frequency'],
                          phase=node['phase'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'frequency': _parameter_to_value(model.frequency),
                'phase': _parameter_to_value(model.phase)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, cls._model) and
                isinstance(b, cls._model))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.frequency, b.frequency)
        assert_array_equal(a.phase, b.phase)


class Sine1DType(Trigonometric1DType):
    name = 'transform/sine1d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.Sine1D']

    _model = functional_models.Sine1D


class Cosine1DType(Trigonometric1DType):
    name = 'transform/cosine1d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.Cosine1D']

    _model = functional_models.Cosine1D


class Tangent1DType(Trigonometric1DType):
    name = 'transform/tangent1d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.Tangent1D']

    _model = functional_models.Tangent1D


class ArcSine1DType(Trigonometric1DType):
    name = 'transform/arcsine1d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.ArcSine1D']

    _model = functional_models.ArcSine1D


class ArcCosine1DType(Trigonometric1DType):
    name = 'transform/arccosine1d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.ArcCosine1D']

    _model = functional_models.ArcCosine1D


class ArcTangent1DType(Trigonometric1DType):
    name = 'transform/arctangent1d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.ArcTangent1D']

    _model = functional_models.ArcTangent1D


class Trapezoid1DType(TransformType):
    name = 'transform/trapezoid1d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.Trapezoid1D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.Trapezoid1D(amplitude=node['amplitude'],
                                             x_0=node['x_0'],
                                             width=node['width'],
                                             slope=node['slope'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'x_0': _parameter_to_value(model.x_0),
                'width': _parameter_to_value(model.width),
                'slope': _parameter_to_value(model.slope)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, functional_models.Trapezoid1D) and
                isinstance(b, functional_models.Trapezoid1D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.x_0, b.x_0)
        assert_array_equal(a.width, b.width)
        assert_array_equal(a.slope, b.slope)


class TrapezoidDisk2DType(TransformType):
    name = 'transform/trapezoid_disk2d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.TrapezoidDisk2D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.TrapezoidDisk2D(amplitude=node['amplitude'],
                                                 x_0=node['x_0'],
                                                 y_0=node['y_0'],
                                                 R_0=node['R_0'],
                                                 slope=node['slope'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'amplitude': _parameter_to_value(model.amplitude),
                'x_0': _parameter_to_value(model.x_0),
                'y_0': _parameter_to_value(model.y_0),
                'R_0': _parameter_to_value(model.R_0),
                'slope': _parameter_to_value(model.slope)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, functional_models.TrapezoidDisk2D) and
                isinstance(b, functional_models.TrapezoidDisk2D))
        assert_array_equal(a.amplitude, b.amplitude)
        assert_array_equal(a.x_0, b.x_0)
        assert_array_equal(a.y_0, b.y_0)
        assert_array_equal(a.R_0, b.R_0)
        assert_array_equal(a.slope, b.slope)


class Voigt1DType(TransformType):
    name = 'transform/voigt1d'
    version = '1.0.0'
    types = ['astropy.modeling.functional_models.Voigt1D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.Voigt1D(x_0=node['x_0'],
                                         amplitude_L=node['amplitude_L'],
                                         fwhm_L=node['fwhm_L'],
                                         fwhm_G=node['fwhm_G'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'x_0': _parameter_to_value(model.x_0),
                'amplitude_L': _parameter_to_value(model.amplitude_L),
                'fwhm_L': _parameter_to_value(model.fwhm_L),
                'fwhm_G': _parameter_to_value(model.fwhm_G)}
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, functional_models.Voigt1D) and
                isinstance(b, functional_models.Voigt1D))
        assert_array_equal(a.x_0, b.x_0)
        assert_array_equal(a.amplitude_L, b.amplitude_L)
        assert_array_equal(a.fwhm_L, b.fwhm_L)
        assert_array_equal(a.fwhm_G, b.fwhm_G)
