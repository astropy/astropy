# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
from distutils.version import LooseVersion

from ..mpl_normalize import simple_norm
from ... import log
from ...io.fits import getdata


def fits2bitmap(filename, ext=0, out_fn=None, stretch='linear',
                power=1.0, asinh_a=0.1, min_cut=None, max_cut=None,
                min_percent=None, max_percent=None, percent=None,
                cmap='Greys_r'):
    """
    Create a bitmap file from a FITS image, applying a stretching
    transform between minimum and maximum cut levels and a matplotlib
    colormap.

    Parameters
    ----------
    filename : str
        The filename of the FITS file.
    ext : int
        FITS extension name or number of the image to convert.  The
        default is 0.
    out_fn : str
        The filename of the output bitmap image.  The type of bitmap
        is determined by the filename extension (e.g. '.jpg', '.png').
        The default is a PNG file with the same name as the FITS file.
    stretch : {{'linear', 'sqrt', 'power', log', 'asinh'}}
        The stretching function to apply to the image.  The default is
        'linear'.
    power : float, optional
        The power index for ``stretch='power'``.  The default is 1.0.
    asinh_a : float, optional
        For ``stretch='asinh'``, the value where the asinh curve
        transitions from linear to logarithmic behavior, expressed as a
        fraction of the normalized image.  Must be in the range between
        0 and 1.  The default is 0.1.
    min_cut : float, optional
        The pixel value of the minimum cut level.  Data values less than
        ``min_cut`` will set to ``min_cut`` before stretching the image.
        The default is the image minimum.  ``min_cut`` overrides
        ``min_percent``.
    max_cut : float, optional
        The pixel value of the maximum cut level.  Data values greater
        than ``min_cut`` will set to ``min_cut`` before stretching the
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

    import matplotlib
    import matplotlib.cm as cm
    import matplotlib.image as mimg

    # __main__ gives ext as a string
    try:
        ext = int(ext)
    except ValueError:
        pass

    try:
        image = getdata(filename, ext)
    except Exception as e:
        log.critical(e)
        return 1

    if image.ndim != 2:
        log.critical('data in FITS extension {0} is not a 2D array'
                     .format(ext))

    if out_fn is None:
        out_fn = os.path.splitext(filename)[0]
        if out_fn.endswith('.fits'):
            out_fn = os.path.splitext(out_fn)[0]
        out_fn += '.png'

    # need to explicitly define the output format due to a bug in
    # matplotlib (<= 2.1), otherwise the format will always be PNG
    out_format = os.path.splitext(out_fn)[1][1:]

    # workaround for matplotlib 2.0.0 bug where png images are inverted
    # (mpl-#7656)
    if (out_format.lower() == 'png' and
            LooseVersion(matplotlib.__version__) == LooseVersion('2.0.0')):
        image = image[::-1]

    if cmap not in cm.datad:
        log.critical('{0} is not a valid matplotlib colormap name.'
                     .format(cmap))
        return 1

    norm = simple_norm(image, stretch=stretch, power=power, asinh_a=asinh_a,
                       min_cut=min_cut, max_cut=max_cut,
                       min_percent=min_percent, max_percent=max_percent,
                       percent=percent)

    mimg.imsave(out_fn, norm(image), cmap=cmap, origin='lower',
                format=out_format)
    log.info('Saved file to {0}.'.format(out_fn))


def main(args=None):

    import argparse

    parser = argparse.ArgumentParser(
        description='Create a bitmap file from a FITS image.')
    parser.add_argument('-e', '--ext', metavar='hdu', default=0,
                        help='Specify the HDU extension number or name '
                             '(Default is 0).')
    parser.add_argument('-o', metavar='filename', type=str, default=None,
                        help='Filename for the output image (Default is a '
                        'PNG file with the same name as the FITS file).')
    parser.add_argument('--stretch', type=str, default='linear',
                        help='Type of image stretching ("linear", "sqrt", '
                        '"power", "log", or "asinh") (Default is "linear").')
    parser.add_argument('--power', type=float, default=1.0,
                        help='Power index for "power" stretching (Default is '
                             '1.0).')
    parser.add_argument('--asinh_a', type=float, default=0.1,
                        help='The value in normalized image where the asinh '
                             'curve transitions from linear to logarithmic '
                             'behavior (used only for "asinh" stretch) '
                             '(Default is 0.1).')
    parser.add_argument('--min_cut', type=float, default=None,
                        help='The pixel value of the minimum cut level '
                             '(Default is the image minimum).')
    parser.add_argument('--max_cut', type=float, default=None,
                        help='The pixel value of the maximum cut level '
                             '(Default is the image maximum).')
    parser.add_argument('--min_percent', type=float, default=None,
                        help='The percentile value used to determine the '
                             'minimum cut level (Default is 0).')
    parser.add_argument('--max_percent', type=float, default=None,
                        help='The percentile value used to determine the '
                             'maximum cut level (Default is 100).')
    parser.add_argument('--percent', type=float, default=None,
                        help='The percentage of the image values used to '
                             'determine the pixel values of the minimum and '
                             'maximum cut levels (Default is 100).')
    parser.add_argument('--cmap', metavar='colormap_name', type=str,
                        default='Greys_r', help='matplotlib color map name '
                                                '(Default is "Greys_r").')
    parser.add_argument('filename', nargs='+',
                        help='Path to one or more FITS files to convert')
    args = parser.parse_args(args)

    for filename in args.filename:
        fits2bitmap(filename, ext=args.ext, out_fn=args.o,
                    stretch=args.stretch, min_cut=args.min_cut,
                    max_cut=args.max_cut, min_percent=args.min_percent,
                    max_percent=args.max_percent, percent=args.percent,
                    power=args.power, asinh_a=args.asinh_a, cmap=args.cmap)
