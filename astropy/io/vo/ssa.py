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
Support VO simple spectral access capabilities.

Based on the `Simple Spectral Access Protocol 1.04 Recommendation
<http://www.ivoa.net/Documents/REC/DAL/SSA-20080201.html>`_.
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

# TODO: Support verification of named bands

def query_data(catalog_db=None, pedantic=False, pos=None, size=None, time=None,
               band=None, **kwargs):
    """
    Perform a Simple Spectral Access data query on a specified
    catalog.  To get the data itself, use `get_data`.

    %(catalog_db)s

    %(pedantic)s

    *pos*: The center of the region of interest.  The coordinate
    values may be specified either as a string or a 2-item sequence.
    If a string it is passed along to the service verbatim, and must
    be two numbers in list format (comma separated) with no embedded
    white space.

    *pos* defaults to right-ascension and declination in decimal
    degrees in the ICRS coordinate system. A coordinate system
    reference frame may optionally be specified to specify a
    coordinate system other than ICRS.  This may be the third element
    of a sequence, or if a string, separated from the two numbers by a
    semicolon.

    *size*: The diameter of the search region specified in decimal
    degrees.

    *band*: The spectral bandpass is given in `range-list-format`_,
    with each list element specified either numerically as a
    wavelength value or range, or textually as a spectral bandpass
    identifier, e.g., a filter or instrumental bandpass name.  May be
    either a single value or (for numerical ranges) an open or closed
    range.  Multiple element range-lists are supported by some
    services.  If a single numerical value is specified as a range
    element it matches any spectrum for which the spectral coverage
    includes the specified value.  If a two valued range is given, a
    dataset matches if any portion of it overlaps the given spectral
    region.

    For a numerical bandpass, the units are wavelength in vacuum in
    units of meters.  The spectral rest frame may be qualified as
    either ``source`` or ``observer``, using either the last element
    of the sequence or separated by a semi-colon if given as a
    string::

       band='1E-7/2E-6;source'

    For most queries the precision with which the spectral bandpass is
    specified in the query probably does not matter very much. A rough
    bandpass broad enough to find all the interesting data will
    generally suffice; the more precise spectral bandpass specified in
    the query response for each spectrum can then be used to refine
    the query. In some cases, for example a cutout service operating
    upon high resolution spectra, support at the service level for
    specifying the spectral rest frame could be important. If the
    service does not support specification of the spectral frame the
    syntax should be permitted but may be ignored.

    *time*: The time coverage (epoch) specified in
    `range-list-format`_, in `ISO 8601
    <http://www.iso.org/iso/date_and_time_format>`_.  If the time
    system used is not specified, UTC is assumed.  The value may be a
    single value or an open or closed range.  If a single value is
    specified, it matches any spectrum for which the time coverage
    includes the specified value.  If a two-valued range is given, a
    dataset matches if any portion of it overlaps the given temporal
    region.

    Additional kwargs may be provided to pass along to the server.
    These arguments are specific to the particular catalog being
    queried.  The standard set of optional arguments for SSA includes:

       APERTURE, SPECRP, SPATRES, TIMERES, SNR, REDSHIFT, VARAMPL,
       TARGETNAME, TARGETCLASS, FLUXCALIB, WAVECALIB, PUBDID,
       CREATORID, COLLECTION, TOP, MAXREC, MTIME, COMPRESS, RUNID
    """
    # Type-check and coerce arguments
    pos, pos_len = util.coerce_range_list_param(
        pos, frames=util.stc_reference_frames, numeric=True)
    assert pos_len == 2
    if size is not None:
        size = float(size)

    if band is not None:
        band, band_len = util.coerce_range_list_param(
            band, frames=('source', 'observer'), numeric=False)
    if time is not None:
        time, time_len = util.coerce_range_list_param(
            time, frames=(), numeric=True)

    args = {
        'REQUEST': 'queryData',

        # Passing format and version to some services seems to result
        # in 0 rows.  Not to spec???
        # 'FORMAT' : 'votable',
        }
    util.dict_soft_update(args, kwargs)
    util.dict_soft_update(
        args, {
            'POS': pos,
            'SIZE': size,
            'TIME': time,
            'BAND': band
            })

    return vos_catalog.call_vo_service(
        'ssa', catalog_db=catalog_db, pedantic=pedantic, kwargs=args)

query_data.__doc__ = query_data.__doc__ % vos_catalog._doc_snippets


def get_data(url, pedantic=False, mimetype=None):
    """
    Retrieve a set of spectral data defined in the results of a
    `query_data` call.

    *url* may either be a string containing a URL or a row from the
    `vo.table` array returned by `query_data`.  If the latter, the
    actual URL used is grabbed from the column named 'AcRef'.  If such
    a column does not exist, a `ValueError` is raised.

    %(pedantic)s

    *mimetype*, if provided, will override the type returned from the
    server.  Only use if you know what type of file you expect to be
    returned.

    If the resulting file is:

      - a fits file, a `pyfits.hdulist` object is returned

      - a VOTABLE file, a :class:`~vo.tree.VOTableFile` object is
        returned

      - otherwise, a read-only file-like object to the raw data is
        returned
    """
    if isinstance(url, np.flexible):
        url = url['AcRef']

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
            return table.parse(response, filename=url, pedantic=pedantic)
        finally:
            response.close()
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

get_data.__doc__ = get_data.__doc__ % vos_catalog._doc_snippets


def list_catalogs():
    """
    Return the available simple spectral access catalogs as a list of
    strings.  These can be used for the *catalog_db* argument to
    :func:`query_data`.
    """
    return vos_catalog.list_catalogs('ssa')
