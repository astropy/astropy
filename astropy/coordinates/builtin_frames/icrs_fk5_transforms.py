# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.matrix_utilities import matrix_transpose, rotation_matrix
from astropy.coordinates.transformations import DynamicMatrixTransform

from .fk5 import FK5
from .icrs import ICRS
from .utils import EQUINOX_J2000


def _icrs_to_fk5_matrix():
    """
    B-matrix from USNO circular 179.  Used by the ICRS->FK5 transformation
    functions.
    """
    eta0 = -19.9 / 3600000.0
    xi0 = 9.1 / 3600000.0
    da0 = -22.9 / 3600000.0

    return (
        rotation_matrix(-eta0, "x")
        @ rotation_matrix(xi0, "y")
        @ rotation_matrix(da0, "z")
    )


# define this here because it only needs to be computed once
_ICRS_TO_FK5_J2000_MAT = _icrs_to_fk5_matrix()


@frame_transform_graph.transform(DynamicMatrixTransform, ICRS, FK5)
def icrs_to_fk5(icrscoord, fk5frame):
    # ICRS is by design very close to J2000 equinox
    pmat = fk5frame._precession_matrix(EQUINOX_J2000, fk5frame.equinox)
    return pmat @ _ICRS_TO_FK5_J2000_MAT


# can't be static because the equinox is needed
@frame_transform_graph.transform(DynamicMatrixTransform, FK5, ICRS)
def fk5_to_icrs(fk5coord, icrsframe):
    # ICRS is by design very close to J2000 equinox
    pmat = fk5coord._precession_matrix(fk5coord.equinox, EQUINOX_J2000)
    return matrix_transpose(_ICRS_TO_FK5_J2000_MAT) @ pmat
