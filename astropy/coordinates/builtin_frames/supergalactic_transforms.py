# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.matrix_utilities import matrix_transpose, rotation_matrix
from astropy.coordinates.transformations import StaticMatrixTransform

from .galactic import Galactic
from .supergalactic import Supergalactic


@frame_transform_graph.transform(StaticMatrixTransform, Galactic, Supergalactic)
def gal_to_supergal():
    return (
        rotation_matrix(90, "z")
        @ rotation_matrix(90 - Supergalactic._nsgp_gal.b.degree, "y")
        @ rotation_matrix(Supergalactic._nsgp_gal.l.degree, "z")
    )


@frame_transform_graph.transform(StaticMatrixTransform, Supergalactic, Galactic)
def supergal_to_gal():
    return matrix_transpose(gal_to_supergal())
