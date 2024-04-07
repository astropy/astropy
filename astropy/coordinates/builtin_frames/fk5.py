# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.coordinates import earth_orientation as earth
from astropy.coordinates.attributes import TimeAttribute
from astropy.coordinates.baseframe import base_doc, frame_transform_graph
from astropy.coordinates.transformations import DynamicMatrixTransform
from astropy.utils.decorators import format_doc

from .baseradec import BaseRADecFrame, doc_components
from .utils import EQUINOX_J2000

__all__ = ["FK5"]


doc_footer = """
    Other parameters
    ----------------
    equinox : `~astropy.time.Time`
        The equinox of this frame.
"""


@format_doc(base_doc, components=doc_components, footer=doc_footer)
class FK5(BaseRADecFrame):
    """
    A coordinate or frame in the FK5 system.

    Note that this is a barycentric version of FK5 - that is, the origin for
    this frame is the Solar System Barycenter, *not* the Earth geocenter.

    The frame attributes are listed under **Other Parameters**.
    """

    equinox = TimeAttribute(default=EQUINOX_J2000, doc="The equinox time")

    @staticmethod
    def _precession_matrix(oldequinox, newequinox):
        """
        Compute and return the precession matrix for FK5 based on Capitaine et
        al. 2003/IAU2006.  Used inside some of the transformation functions.

        Parameters
        ----------
        oldequinox : `~astropy.time.Time`
            The equinox to precess from.
        newequinox : `~astropy.time.Time`
            The equinox to precess to.

        Returns
        -------
        newcoord : array
            The precession matrix to transform to the new equinox
        """
        return earth.precession_matrix_Capitaine(oldequinox, newequinox)


# This is the "self-transform".  Defined at module level because the decorator
#  needs a reference to the FK5 class


@frame_transform_graph.transform(DynamicMatrixTransform, FK5, FK5)
def fk5_to_fk5(fk5coord1, fk5frame2):
    return fk5coord1._precession_matrix(fk5coord1.equinox, fk5frame2.equinox)
