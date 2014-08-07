"""command-line script"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import argparse
from astropy import log
import astropy.io.fits as fits
import imageutils
import matplotlib.cm as cm
import matplotlib.image as mimg


def fits2bitmap(filename, exten=0, out_fn=None, scale='linear',
                power=1.0, noise_level=None, min_cut=None, max_cut=None,
                min_percent=None, max_percent=None, percent=None,
                cmap='Greys_r'):
    """
    Create a bitmap file from a FITS image, applying a
    scaling/stretching transform between minimum and maximum cut levels
    and a matplotlib colormap.

    Parameters
    ----------
    filename : str
        The filename of the FITS file.

    exten : int
        FITS extension number of the image to convert.  The default is 0.

    out_fn : str
        The filename of the output bitmap image.  The type of bitmap
        is determined by the filename extension (e.g. '.jpg', '.png').
        The default is a PNG file with the same name as the FITS file.

    scale : {{'linear', 'sqrt', 'power', log', 'asinh'}}
        The scaling/stretch function to apply to the image.  The default
        is 'linear'.

    power : float, optional
        The power index for the image scaling.  The default is 1.0.

    noise_level : float, optional
        The noise level of the image.  Pixel values less than
        ``noise_level`` will approximately be linearly scaled, while
        pixel values greater than ``noise_level`` will approximately be
        logarithmically scaled.  If ``noise_level`` is not input, an
        estimate of the 2-sigma level above the background will be used.

    min_cut : float, optional
        The pixel value of the minimum cut level.  Data values less than
        ``min_cut`` will set to ``min_cut`` before scaling the image.
        The default is the image minimum.  ``min_cut`` overrides
        ``min_percent``.

    max_cut : float, optional
        The pixel value of the maximum cut level.  Data values greater
        than ``min_cut`` will set to ``min_cut`` before scaling the
        image.  The default is the image maximum.  ``max_cut`` overrides
        ``max_percent``.

    min_percent : float, optional
        The percentile value used to determine the pixel value of
        minimum cut level.  The default is 0.0.  ``min_percent``
        overrides ``percent``.

    max_percent : float, optional
        The percentile value used to determine the pixel value of
        maximum cut level.  The default is 100.0.  ``max_percent``
        overrides ``percent``.

    percent : float, optional
        The percentage of the image values used to determine the pixel
        values of the minimum and maximum cut levels.  The lower cut
        level will set at the ``(100 - percent) / 2`` percentile, while
        the upper cut level will be set at the ``(100 + percent) / 2``
        percentile.  The default is 100.0.  ``percent`` is ignored if
        either ``min_percent`` or ``max_percent`` is input.

    cmap : str
        The matplotlib color map name.  The default is 'Greys_r'.
    """

    hdulist = fits.open(filename)
    image = hdulist[exten].data
    hdulist.close()
    if out_fn is None:
        out_fn = filename.replace('.fits', '.png')
    if cmap not in cm.datad.keys():
        log.critical('{0} is not a valid matplotlib colormap '
                     'name'.format(cmap))
        raise SystemExit()

    image_scaled = imageutils.scale_image(image, scale=scale, power=power,
                                          noise_level=noise_level,
                                          min_cut=min_cut, max_cut=max_cut,
                                          min_percent=min_percent,
                                          max_percent=max_percent,
                                          percent=percent)

    mimg.imsave(out_fn, image_scaled, cmap=cmap)
    log.info('Saved file to {0}'.format(out_fn))
    return


def main(args=None):
    parser = argparse.ArgumentParser(
        description='Create a bitmap file from a FITS image.')
    parser.add_argument('-e', '--exten', type=int, default=0,
                        help='FITS extension number of the image to convert')
    parser.add_argument('-o', metavar='filename', type=str, default=None,
                        help='Filename for the output image (Default is a '
                        'PNG file with the same name as the FITS file)')
    parser.add_argument('--scale', type=str, default='linear',
                        help='Type of image scaling ("linear", "sqrt", '
                        '"power", "log", or "asinh")')
    parser.add_argument('--power', type=float, default=1.0,
                        help='Power index for "power" scaling')
    parser.add_argument('--noise', type=float, default=None,
                        help=('The noise level of the image (used only for '
                              '"asinh" scaling)'))
    parser.add_argument('--min_cut', type=float, default=None,
                        help='The pixel value of the minimum cut level')
    parser.add_argument('--max_cut', type=float, default=None,
                        help='The pixel value of the maximum cut level')
    parser.add_argument('--min_percent', type=float, default=None,
                        help=('The percentile value used to determine the '
                              'minimum cut level'))
    parser.add_argument('--max_percent', type=float, default=None,
                        help=('The percentile value used to determine the '
                              'maximum cut level'))
    parser.add_argument('--percent', type=float, default=None,
                        help=('The percentage of the image values used to '
                              'determine the pixel values of the minimum and '
                              'maximum cut levels'))
    parser.add_argument('--cmap', metavar='colormap_name', type=str,
                        default='Greys_r', help='matplotlib color map name')
    parser.add_argument('filename', nargs='+',
                        help='Path to one or more FITS files to convert')
    args = parser.parse_args(args)

    for filename in args.filename:
        fits2bitmap(filename, exten=args.exten, out_fn=args.o,
                    scale=args.scale, min_cut=args.min_cut,
                    max_cut=args.max_cut, min_percent=args.min_percent,
                    max_percent=args.max_percent, percent=args.percent,
                    power=args.power, noise_level=args.noise, cmap=args.cmap)
