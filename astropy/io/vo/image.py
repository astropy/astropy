# Copyright (C) 2008-2010 Association of Universities for Research in Astronomy (AURA)

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

#     1. Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.

#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.

#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.

"""
Support VO simple image access capabilities.

Based on the `Simple Image Access Protocol 1.0 Recommendation
<http://www.ivoa.net/Documents/SIA/20091116/REC-SIA-1.0.html>`_.
"""

from __future__ import division, absolute_import

# STDLIB
import tempfile
import warnings

# THIRD-PARTY
import numpy as np

# LOCAL
from . import table
from . import util
from . import vos_catalog
from . import webquery

def query(catalog_db=None, pedantic=False, pos=None, size=None,
          intersect=None, naxis=None, cframe=None, equinox=None,
          crpix=None, crval=None, cdelt=None, rotang=None, proj=None,
          format=None, verb=None, **kwargs):
    """
    Perform a Simple Image Access data query on a specified
    catalog.

    %(catalog_db)s

    %(pedantic)s

    **Core parameters:**

    *pos*: The position of the region of interest, expressed as the
    right-ascension and declination of the field center, in decimal
    degrees using the ICRS coordinate system. May be a 2-tuple of
    floats or a string.  If a string, a comma should delimit the two
    values; embedded whitespace is not permitted. Example:
    ``pos="12.821,-33.4"``.

    *size*: The coordinate angular size of the region given in decimal
    degrees. The region may be specified using either one or two
    values:

      - If only one value is given, it applies to both coordinate
        axes.

      - If two values are given, the first value is the angular width
        in degrees of the right-ascension axis of the region, and the
        second value is the angular width in degrees of the declination
        axis. Note that the angular width along the right-ascension axis
        is not the coordinate width of the region in right-ascension,
        which is equal to the coordinate angular width multiplied by
        cos(DEC).

    A special case is ``size == 0``. For an atlas or pointed image
    archive, this tests whether the given point is in any image. For a
    cutout or mosaic service, this will cause a default-sized image to
    be returned. The default image size is service-defined and may be
    a value considered appropriate for the service, for the given
    image or data collection being accessed, or for the object (if
    any) at the given position.

    The *pos* and *size* parameters define a nonrotated, rectangular
    region on the sky, having the specified angular extent in right
    ascension and declination, using the cartesian (``CAR``)
    projection with the region center (*pos*) as the reference
    point. The cartesian projection is used here as it is simple and
    can scale to the whole sky, and works about as well as anything
    else for small regions. Note: the use of the ``CAR`` projection to
    define the ROI has nothing to do with what projection we choose
    for any actual generated images.

    *intersect*: A parameter that indicates how matched images should
    intersect the region of interest. The allowed values are the strings:

      - ``'COVERS'``: The candidate image covers or includes the
        entire ROI.

      - ``'ENCLOSED'``: The candidate image is entirely enclosed by
        the ROI.

      - ``'CENTER'``: The candidate image overlaps the center of the
        ROI.

      - ``'OVERLAPS'``: The candidate image overlaps some part of the
        ROI.

    If this parameter is not present, ``intersect == 'OVERLAPS'`` is
    assumed.

    **Image generation parameters:** (Not supported by all services)

    The default size and extent of the desired output image are
    determined from these parameters (*naxis*, *cframe*, *equinox*,
    *crpix*, *crval*, *cdelt*, *rotang* and *proj*) as follows: if a
    size or scale term is given explicitly this is the value used,
    otherwise a default is computed based on the angular extent of the
    query region (*size*). For example, if *naxis* and *cdelt* are
    given, then the angular size of the image is determined by these
    parameters. If only the *naxis* are given, then the default image
    scale is determined by *size*. If only *cdelt* is given, then the
    image size in pixels is determined by the angular extent of the
    query region (*size*). If only the region angular extent *size* is
    given then the size (in pixels) and scale of the image is
    determined by the service to best fit the data and the region of
    interest.

    In the simplest case, all image generation parameters may be
    omitted, and the generated output image defaults to the position
    and size of the query region. If any image generation parameters
    are specified these override the ROI-implied defaults.

    *naxis*: The size of the output image in pixels. May be a tuple
    of floats or a string.  If only one value is given it applies to
    both image axes. Default: determined from the ROI (see
    below).

    *cframe*: The coordinate system reference frame, a string selected
    from 'ICRS', 'FK5', 'FK4', 'ECL', 'GAL', and 'SGAL' (these
    abbreviations follow `CDS Aladin
    <http://aladin.u-strasbg.fr/java/doctech/cds.astro.Astroframe.html>`_).
    Default: 'ICRS'.

    *equinox*: Epoch of the mean equator and equinox for the specified
    coordinate system reference frame (*cframe*). Not required for
    ICRS. Default: B1950 for FK4, otherwise J2000.

    *crpix*: The coordinates of the reference pixel, expressed in the
    pixel coordinates of the output image, with [1,1] being the
    center of the first pixel of the first row of the image. This is
    a vector-valued quantity; if only one value is given it applies
    to both image axes. Default: the image center.

    *crval*: The world coordinates relative to *cframe* at the
    reference pixel. This is a vector-valued quantity; both array
    values are required. Default: the region center coordinates
    (*pos*) at the center of the image, transformed to the output
    coordinate system reference frame if other than ICRS. If *crpix*
    is specified to be other than the image center the corresponding
    *crval* can be computed, but should be specified explicitly by the
    client.

    *cdelt*: The scale of the output image in decimal degrees per
    pixel. A negative value implies an axis flip. Since the default
    image orientation is N up and E to the left, the default sign of
    CDELT is [-1,1]. This is a vector-valued quantity; if only one
    value is given it applies to both image axes, with the sign
    defaulting as specified above. Default: implied (see below),
    otherwise service-specific.

    *rotang*: The rotation angle of the image in degrees relative to
    *cframe* (an image which is unrotated in one reference frame may be
    rotated in another). This is the rotation of the WCS declination
    or latitude axis with respect to the second axis of the image,
    measured in the counterclockwise direction (as for FITS WCS,
    which is in turn based on the old AIPS convention). Default: 0
    (no rotation).

    *proj*: The celestial projection of the output image expressed as
    a three-character code as for FITS WCS, e.g., "TAN", "SIN", "ARC",
    and so forth. Default: TAN.

    **Image format:**

    *format*: A tuple or comma-delimited string where each element is
    any recognized MIME type.  Generally used MIME types include::

        image/fits
        image/png
        image/jpeg

    In addition these special values are defined:

        - ``'ALL'``: Denotes all formats supported by the service

        - ``'GRAPHIC'``: Denotes any of the following graphics
          formats: JPEG, PNG, GIF. These are types typically supported
          "in-line" by web-browsers.  There are more finegrained
          options for the ``'GRAPHIC'`` syntax defined in the Simple
          Image Access specification.

        - ``'METADATA'``: Denotes a Metadata Query: no images are
          requested; only metadata should be returned. This feature is
          described in more detail in `section 6.1 of the Simple Image
          Access specification
          <http://www.ivoa.net/Documents/SIA/20091116/REC-SIA-1.0.html#mdquery>`_.

    **Table verbosity:** (not supported by all services)

    *verb*: This parameter indicates the desired level of information
    to be returned in the output table, particularly the number of
    columns to be returned to describe each image.

        - 0: The output table should contain only the minimum columns
          required.

        - 1: In addition to level 0, the output table should contain
          columns sufficient for uniquely describing the image.

        - 2: In addition to level 1, the output table should contain,
          if possible, columns that contain values for all parameters
          supported as query constraints.

        - 3: The output table should return as much information about
          the images as possible. A table metadata query automatically
          implies the highest level of verbosity.
    """
    # Type-check and coerce arguments
    pos, pos_len = util.coerce_range_list_param(
        pos, frames=(), numeric=True)
    assert pos_len == 2

    size, size_len = util.coerce_range_list_param(
        size, frames=(), numeric=True)
    assert size_len in (0, 1, 2)

    if intersect is not None:
        assert isinstance(intersect, basestring)
        intersect = intersect.upper()
        assert intersect in ('COVERS', 'ENCLOSED', 'CENTER', 'OVERLAPS')

    naxis, naxis_len = util.coerce_range_list_param(
        naxis, frames=(), numeric=True)
    assert naxis_len in (0, 1, 2)

    if cframe is not None:
        assert isinstance(cframe, basestring)
        cframe = cframe.upper()
        assert cframe in ('ICRS', 'FK5', 'FK4', 'ECL', 'GAL', 'SGAL')

    if equinox is not None:
        assert isinstance(equinox, basestring)

    crpix, crpix_len = util.coerce_range_list_param(
        crpix, frames=(), numeric=True)
    assert crpix_len in (0, 1, 2)

    crval, crval_len = util.coerce_range_list_param(
        crval, frames=(), numeric=True)
    assert crval_len in (0, 2)

    cdelt, cdelt_len = util.coerce_range_list_param(
        cdelt, frames=(), numeric=True)
    assert cdelt_len in (0, 1, 2)

    if rotang is not None:
        rotang = float(rotang)

    if proj is not None:
        assert isinstance(proj, basestring)
        proj = proj.upper()
        assert len(proj) == 3

    if format is not None:
        format_correct = False
        if isinstance(format, basestring):
            format = format.upper()
            if format in ('ALL', 'METADATA'):
                format_correct = True
            elif format.startswith('GRAPHIC'):
                format_correct = True
            else:
                format = format.split(',')
        if not format_correct:
            for entry in format:
                assert xmlutil.check_mime_content_type(entry.lower())
            format = ','.join([x.lower() for x in format])

    args = {}
    util.dict_soft_update(args, kwargs)
    util.dict_soft_update(
        args, {
            'POS'       : pos,
            'SIZE'      : size,
            'INTERSECT' : intersect,
            'NAXIS'     : naxis,
            'CFRAME'    : cframe,
            'EQUINOX'   : equinox,
            'CRPIX'     : crpix,
            'CRVAL'     : crval,
            'CDELT'     : cdelt,
            'ROTANG'    : rotang,
            'PROJ'      : proj,
            'FORMAT'    : format,
            'VERB'      : verb
            })

    return vos_catalog.call_vo_service(
        'image', catalog_db=catalog_db, pedantic=pedantic, kwargs=args)

query.__doc__ = query.__doc__ % vos_catalog._doc_snippets


def get(url):
    """
    Retrieve an image defined in the results of a `query` call.

    *url* may either be a string containing a URL or a row from the
    `vo.table` array returned by `query`.  If the latter, the actual
    URL used is grabbed from the column named 'accessreference'.  If
    such a column does not exist, a `ValueError` is raised.

    If the resulting file is:

      - a fits file, a `pyfits.hdulist` object is returned

      - a VOTABLE file, a :class:`~vo.tree.VOTableFile` object is
        returned

      - otherwise, a read-only file-like object to the raw data is
        returned
    """
    if isinstance(url, np.flexible):
        url = url['accessreference']

    response = webquery.webget_open(url)
    try:
        if mimetype is None:
            if 'content-type' in response.info():
                mimetype = response.info()['content-type']
        if mimetype is not None:
            mimetype = mimetype[mimetype.find('/') + 1:]
    except:
        response.close()
        raise

    if mimetype == 'x-votable+xml':
        try:
            result = table.parse(response, filename=url, pedantic=pedantic)
        finally:
            response.close()
        for info in result.resources[0].infos:
            if info.name == 'QUERY_STATUS':
                if info.content is not None:
                    long_descr = ':\n%s' % info.content
                else:
                    long_descr = ''
                raise IOError(
                    "Image could not be retrieved from '%s'.\nStatus: '%s'%s" %
                    (url, info.value, long_descr))
            break
    elif mimetype == 'fits':
        try:
            import pyfits
        except ImportError:
            raise TypeError(
                "The response is a FITS file, but pyfits is not installed.")
        try:
            return pyfits.open(response)
        finally:
            response.close()
    else:
        return response

get.__doc__ = get.__doc__ % vos_catalog._doc_snippets


def list_catalogs():
    """
    Return the available Simple Image Access catalogs as a list of
    strings.  These can be used for the *catalog_db* argument to
    :func:`get_image`.
    """
    return vos_catalog.list_catalogs('image')
