# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Helpers functions for different kinds of WCSAxes instances.
"""

import numpy as np
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredEllipse, AnchoredSizeBar

import astropy.units as u
from astropy.wcs.utils import proj_plane_pixel_scales

__all__ = ["add_beam", "add_scalebar"]

CORNERS = {
    "top right": 1,
    "top left": 2,
    "bottom left": 3,
    "bottom right": 4,
    "right": 5,
    "left": 6,
    "bottom": 8,
    "top": 9,
}


def add_beam(
    ax,
    header=None,
    major=None,
    minor=None,
    angle=None,
    corner="bottom left",
    frame=False,
    borderpad=0.4,
    pad=0.5,
    **kwargs,
):
    """
    Display the beam shape and size.

    Parameters
    ----------
    ax : :class:`~astropy.visualization.wcsaxes.WCSAxes`
        WCSAxes instance in which the beam shape and size is displayed. The WCS
        must be celestial.
    header : :class:`~astropy.io.fits.Header`, optional
        Header containing the beam parameters. If specified, the ``BMAJ``,
        ``BMIN``, and ``BPA`` keywords will be searched in the FITS header
        to set the major and minor axes and the position angle on the sky.
    major : float or :class:`~astropy.units.Quantity`, optional
        Major axis of the beam in degrees or an angular quantity.
    minor : float, or :class:`~astropy.units.Quantity`, optional
        Minor axis of the beam in degrees or an angular quantity.
    angle : float or :class:`~astropy.units.Quantity`, optional
        Position angle of the beam on the sky in degrees or an angular
        quantity in the anticlockwise direction.
    corner : str, optional
        The beam location. Acceptable values are ``'left'``, ``'right'``,
        ``'top'``, 'bottom', ``'top left'``, ``'top right'``, ``'bottom left'``
        (default), and ``'bottom right'``.
    frame : bool, optional
        Whether to display a frame behind the beam (default is ``False``).
    borderpad : float, optional
        Border padding, in fraction of the font size. Default is 0.4.
    pad : float, optional
        Padding around the beam, in fraction of the font size. Default is 0.5.
    kwargs
        Additional arguments are passed to :class:`matplotlib.patches.Ellipse`.

    Notes
    -----
    This function may be inaccurate when:

    - The pixel scales at the reference pixel are different from the pixel scales
      within the image extent (e.g., when the reference pixel is well outside of
      the image extent and the projection is non-linear)
    - The pixel scales in the two directions are very different from each other
      (e.g., rectangular pixels)

    """
    if header and major:
        raise ValueError(
            "Either header or major/minor/angle must be specified, not both."
        )

    if header:
        major = header["BMAJ"]
        minor = header["BMIN"]
        angle = header["BPA"]

    if isinstance(major, u.Quantity):
        major = major.to(u.degree).value

    if isinstance(minor, u.Quantity):
        minor = minor.to(u.degree).value

    if isinstance(angle, u.Quantity):
        angle = angle.to(u.degree).value

    if ax.wcs.is_celestial:
        pix_scale = proj_plane_pixel_scales(ax.wcs)
        sx = pix_scale[0]
        sy = pix_scale[1]
        degrees_per_pixel = np.sqrt(sx * sy)
    else:
        raise ValueError("Cannot show beam when WCS is not celestial")

    minor /= degrees_per_pixel
    major /= degrees_per_pixel

    corner = CORNERS[corner]

    beam = AnchoredEllipse(
        ax.transData,
        width=minor,
        height=major,
        angle=angle,
        loc=corner,
        pad=pad,
        borderpad=borderpad,
        frameon=frame,
    )
    beam.ellipse.set(**kwargs)

    ax.add_artist(beam)


def add_scalebar(
    ax,
    length,
    label=None,
    corner="bottom right",
    frame=False,
    borderpad=0.4,
    pad=0.5,
    **kwargs,
):
    """Add a scale bar.

    Parameters
    ----------
    ax : :class:`~astropy.visualization.wcsaxes.WCSAxes`
        WCSAxes instance in which the scale bar is displayed. The WCS must be
        celestial.
    length : float or :class:`~astropy.units.Quantity`
        The length of the scalebar in degrees or an angular quantity
    label : str, optional
        Label to place below the scale bar
    corner : str, optional
        Where to place the scale bar. Acceptable values are:, ``'left'``,
        ``'right'``, ``'top'``, ``'bottom'``, ``'top left'``, ``'top right'``,
        ``'bottom left'`` and ``'bottom right'`` (default)
    frame : bool, optional
        Whether to display a frame behind the scale bar (default is ``False``)
    borderpad : float, optional
        Border padding, in fraction of the font size. Default is 0.4.
    pad : float, optional
        Padding around the scale bar, in fraction of the font size. Default is 0.5.
    kwargs
        Additional arguments are passed to
        :class:`mpl_toolkits.axes_grid1.anchored_artists.AnchoredSizeBar`.

    Notes
    -----
    This function may be inaccurate when:

    - The pixel scales at the reference pixel are different from the pixel scales
      within the image extent (e.g., when the reference pixel is well outside of
      the image extent and the projection is non-linear)
    - The pixel scales in the two directions are very different from each other
      (e.g., rectangular pixels)

    """
    if isinstance(length, u.Quantity):
        length = length.to(u.degree).value

    if ax.wcs.is_celestial:
        pix_scale = proj_plane_pixel_scales(ax.wcs)
        sx = pix_scale[0]
        sy = pix_scale[1]
        degrees_per_pixel = np.sqrt(sx * sy)
    else:
        raise ValueError("Cannot show scalebar when WCS is not celestial")

    length = length / degrees_per_pixel

    corner = CORNERS[corner]

    scalebar = AnchoredSizeBar(
        ax.transData,
        length,
        label,
        corner,
        pad=pad,
        borderpad=borderpad,
        sep=5,
        frameon=frame,
        **kwargs,
    )

    ax.add_artist(scalebar)
