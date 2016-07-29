# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)



import numpy as np

from ..baseframe import frame_transform_graph
from ..transformations import DynamicMatrixTransform
from ..matrix_utilities import matrix_product, matrix_transpose


from .fk4 import FK4NoETerms
from .fk5 import FK5
from .utils import EQUINOX_B1950, EQUINOX_J2000


# FK5 to/from FK4 ------------------->
# B1950->J2000 matrix from Murray 1989 A&A 218,325 eqn 28
_B1950_TO_J2000_M = np.array(
    [[0.9999256794956877, -0.0111814832204662, -0.0048590038153592],
     [0.0111814832391717,  0.9999374848933135, -0.0000271625947142],
     [0.0048590037723143, -0.0000271702937440,  0.9999881946023742]])

_FK4_CORR = np.array(
    [[-0.0026455262, -1.1539918689, +2.1111346190],
     [+1.1540628161, -0.0129042997, +0.0236021478],
     [-2.1112979048, -0.0056024448, +0.0102587734]]) * 1.e-6

def _fk4_B_matrix(obstime):
    """
    This is a correction term in the FK4 transformations because FK4 is a
    rotating system - see Murray 89 eqn 29
    """
    # Note this is *julian century*, not besselian
    T = (obstime.jyear - 1950.) / 100.
    if getattr(T, 'shape', ()):
        # Ensure we broadcast possibly arrays of times properly.
        T.shape += (1, 1)
    return _B1950_TO_J2000_M + _FK4_CORR * T


# This transformation can't be static because the observation date is needed.
@frame_transform_graph.transform(DynamicMatrixTransform, FK4NoETerms, FK5)
def fk4_no_e_to_fk5(fk4noecoord, fk5frame):
    # Correction terms for FK4 being a rotating system
    B = _fk4_B_matrix(fk4noecoord.obstime)

    # construct both precession matricies - if the equinoxes are B1950 and
    # J2000, these are just identity matricies
    pmat1 = fk4noecoord._precession_matrix(fk4noecoord.equinox, EQUINOX_B1950)
    pmat2 = fk5frame._precession_matrix(EQUINOX_J2000, fk5frame.equinox)

    return matrix_product(pmat2, B, pmat1)


# This transformation can't be static because the observation date is needed.
@frame_transform_graph.transform(DynamicMatrixTransform, FK5, FK4NoETerms)
def fk5_to_fk4_no_e(fk5coord, fk4noeframe):
    # Get transposed version of the rotating correction terms... so with the
    # transpose this takes us from FK5/J200 to FK4/B1950
    B = matrix_transpose(_fk4_B_matrix(fk4noeframe.obstime))

    # construct both precession matricies - if the equinoxes are B1950 and
    # J2000, these are just identity matricies
    pmat1 = fk5coord._precession_matrix(fk5coord.equinox, EQUINOX_J2000)
    pmat2 = fk4noeframe._precession_matrix(EQUINOX_B1950, fk4noeframe.equinox)

    return matrix_product(pmat2, B, pmat1)
