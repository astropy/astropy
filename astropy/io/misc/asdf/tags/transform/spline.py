from astropy.modeling.models import Spline1D
from astropy.io.misc.asdf.tags.transform.basic import TransformType


__all__ = ['SplineType']


class SplineType(TransformType):
    name = 'transform/spline1d'
    version = '1.0.0'
    types = ['astropy.modeling.spline.Spline1D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return Spline1D(knots=node['knots'],
                        coeffs=node['coefficients'],
                        degree=node['degree'])

    @classmethod
    def to_tree_transform(cls, model, ctx):
        return {
            "knots": model.t,
            "coefficients": model.c,
            "degree": model.degree
        }
