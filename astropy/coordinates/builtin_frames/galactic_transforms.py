# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.coordinates.matrix_utilities import (rotation_matrix,
                                                  matrix_product,
                                                  matrix_transpose)
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.transformations import DynamicMatrixTransform

from .fk5 import FK5
from .fk4 import FK4NoETerms
from .utils import EQUINOX_B1950, EQUINOX_J2000
from .galactic import Galactic


# Galactic to/from FK4/FK5 ----------------------->
# can't be static because the equinox is needed
@frame_transform_graph.transform(DynamicMatrixTransform, FK5, Galactic)
def fk5_to_gal(fk5coord, galframe):
    # need precess to J2000 first
    pmat = fk5coord._precession_matrix(fk5coord.equinox, EQUINOX_J2000)
    mat1 = rotation_matrix(180 - Galactic._lon0_J2000.degree, 'z')
    mat2 = rotation_matrix(90 - Galactic._ngp_J2000.dec.degree, 'y')
    mat3 = rotation_matrix(Galactic._ngp_J2000.ra.degree, 'z')

    return matrix_product(mat1, mat2, mat3, pmat)


@frame_transform_graph.transform(DynamicMatrixTransform, Galactic, FK5)
def _gal_to_fk5(galcoord, fk5frame):
    return matrix_transpose(fk5_to_gal(fk5frame, galcoord))


@frame_transform_graph.transform(DynamicMatrixTransform, FK4NoETerms, Galactic)
def fk4_to_gal(fk4coords, galframe):
    mat1 = rotation_matrix(180 - Galactic._lon0_B1950.degree, 'z')
    mat2 = rotation_matrix(90 - Galactic._ngp_B1950.dec.degree, 'y')
    mat3 = rotation_matrix(Galactic._ngp_B1950.ra.degree, 'z')
    matprec = fk4coords._precession_matrix(fk4coords.equinox, EQUINOX_B1950)

    return matrix_product(mat1, mat2, mat3, matprec)


@frame_transform_graph.transform(DynamicMatrixTransform, Galactic, FK4NoETerms)
def gal_to_fk4(galcoords, fk4frame):
    return matrix_transpose(fk4_to_gal(fk4frame, galcoords))
