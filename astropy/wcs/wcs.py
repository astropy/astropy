# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Under the hood, there are 3 separate classes that perform different
parts of the transformation:

   - `~astropy.wcs.Wcsprm`: Is a direct wrapper of the core WCS
     functionality in `wcslib`_.  (This includes TPV and TPD
     polynomial distortion, but not SIP distortion).

   - `~astropy.wcs.Sip`: Handles polynomial distortion as defined in the
     `SIP`_ convention.

   - `~astropy.wcs.DistortionLookupTable`: Handles `distortion paper`_
     lookup tables.

Additionally, the class `WCS` aggregates all of these transformations
together in a pipeline:

   - Detector to image plane correction (by a pair of
     `~astropy.wcs.DistortionLookupTable` objects).

   - `SIP`_ distortion correction (by an underlying `~astropy.wcs.Sip`
     object)

   - `distortion paper`_ table-lookup correction (by a pair of
     `~astropy.wcs.DistortionLookupTable` objects).

   - `wcslib`_ WCS transformation (by a `~astropy.wcs.Wcsprm` object)

"""
from __future__ import absolute_import, division, print_function, unicode_literals

# STDLIB
import copy
import io
import os
import textwrap
import warnings
import platform

# THIRD-PARTY
import numpy as np

# LOCAL
from ..extern import six
from ..io import fits
from . import _docutil as __
try:
    from . import _wcs
except ImportError:
    if not _ASTROPY_SETUP_:
        raise
    else:
        _wcs = None

if _wcs is not None:
    _parsed_version = _wcs.__version__.split('.')
    if int(_parsed_version[0]) == 5 and int(_parsed_version[1]) < 8:
        raise ImportError(
            "astropy.wcs is built with wcslib {0}, but only versions 5.8 and "
            "later on the 5.x series are known to work.  The version of wcslib "
            "that ships with astropy may be used.")

from ..utils import deprecated, deprecated_attribute
from ..utils.compat import possible_filename
from ..utils.exceptions import AstropyWarning, AstropyUserWarning, AstropyDeprecationWarning

if _wcs is not None:
    assert _wcs._sanity_check(), \
        "astropy.wcs did not pass its sanity check for your build " \
        "on your platform."


__all__ = ['FITSFixedWarning', 'WCS', 'find_all_wcs',
           'DistortionLookupTable', 'Sip', 'Tabprm', 'Wcsprm',
           'WCSBase', 'validate', 'WcsError', 'SingularMatrixError',
           'InconsistentAxisTypesError', 'InvalidTransformError',
           'InvalidCoordinateError', 'NoSolutionError',
           'InvalidSubimageSpecificationError',
           'NonseparableSubimageCoordinateSystemError',
           'NoWcsKeywordsFoundError', 'InvalidTabularParametersError']


if six.PY3 or platform.system() == 'Windows':
    __doctest_skip__ = ['WCS.all_world2pix']


if _wcs is not None:
    WCSBase = _wcs._Wcs
    DistortionLookupTable = _wcs.DistortionLookupTable
    Sip = _wcs.Sip
    Wcsprm = _wcs.Wcsprm
    Tabprm = _wcs.Tabprm
    WcsError = _wcs.WcsError
    SingularMatrixError = _wcs.SingularMatrixError
    InconsistentAxisTypesError = _wcs.InconsistentAxisTypesError
    InvalidTransformError = _wcs.InvalidTransformError
    InvalidCoordinateError = _wcs.InvalidCoordinateError
    NoSolutionError = _wcs.NoSolutionError
    InvalidSubimageSpecificationError = _wcs.InvalidSubimageSpecificationError
    NonseparableSubimageCoordinateSystemError = _wcs.NonseparableSubimageCoordinateSystemError
    NoWcsKeywordsFoundError = _wcs.NoWcsKeywordsFoundError
    InvalidTabularParametersError = _wcs.InvalidTabularParametersError

    # Copy all the constants from the C extension into this module's namespace
    for key, val in _wcs.__dict__.items():
        if (key.startswith('WCSSUB') or
            key.startswith('WCSHDR') or
            key.startswith('WCSHDO')):
            locals()[key] = val
            __all__.append(key)
else:
    WCSBase = object
    Wcsprm = object
    DistortionLookupTable = object
    Sip = object
    Tabprm = object
    WcsError = None
    SingularMatrixError = None
    InconsistentAxisTypesError = None
    InvalidTransformError = None
    InvalidCoordinateError = None
    NoSolutionError = None
    InvalidSubimageSpecificationError = None
    NonseparableSubimageCoordinateSystemError = None
    NoWcsKeywordsFoundError = None
    InvalidTabularParametersError = None


# Additional relax bit flags
WCSHDO_SIP = 0x10000


def _parse_keysel(keysel):
    keysel_flags = 0
    if keysel is not None:
        for element in keysel:
            if element.lower() == 'image':
                keysel_flags |= _wcs.WCSHDR_IMGHEAD
            elif element.lower() == 'binary':
                keysel_flags |= _wcs.WCSHDR_BIMGARR
            elif element.lower() == 'pixel':
                keysel_flags |= _wcs.WCSHDR_PIXLIST
            else:
                raise ValueError(
                    "keysel must be a list of 'image', 'binary' " +
                    "and/or 'pixel'")
    else:
        keysel_flags = -1

    return keysel_flags


class NoConvergence(Exception):
    """
    An error class used to report non-convergence and/or divergence
    of numerical methods. It is used to report errors in the
    iterative solution used by
    the :py:meth:`~astropy.wcs.WCS.all_world2pix`.

    Attributes
    ----------

    best_solution : numpy.ndarray
        Best solution achieved by the numerical method.

    accuracy : numpy.ndarray
        Accuracy of the :py:attr:`best_solution`.

    niter : int
        Number of iterations performed by the numerical method
        to compute :py:attr:`best_solution`.

    divergent : None, numpy.ndarray
        Indices of the points in :py:attr:`best_solution` array
        for which the solution appears to be divergent. If the
        solution does not diverge, `divergent` will be set to `None`.

    slow_conv : None, numpy.ndarray
        Indices of the solutions in :py:attr:`best_solution` array
        for which the solution failed to converge within the
        specified maximum number of iterations. If there are no
        non-converging solutions (i.e., if the required accuracy
        has been achieved for all input data points)
        then `slow_conv` will be set to `None`.

    """
    def __init__(self, *args, **kwargs):
        super(NoConvergence, self).__init__(*args)

        self.best_solution = kwargs.pop('best_solution', None)
        self.accuracy = kwargs.pop('accuracy', None)
        self.niter = kwargs.pop('niter', None)
        self.divergent = kwargs.pop('divergent', None)
        self.slow_conv = kwargs.pop('slow_conv', None)


class FITSFixedWarning(AstropyWarning):
    """
    The warning raised when the contents of the FITS header have been
    modified to be standards compliant.
    """
    pass


class WCS(WCSBase):
    """WCS objects perform standard WCS transformations, and correct for
    `SIP`_ and `distortion paper`_ table-lookup transformations, based
    on the WCS keywords and supplementary data read from a FITS file.

    Parameters
    ----------
    header : astropy.io.fits header object, string, dict-like, or None, optional
        If *header* is not provided or None, the object will be
        initialized to default values.

    fobj : An astropy.io.fits file (hdulist) object, optional
        It is needed when header keywords point to a `distortion
        paper`_ lookup table stored in a different extension.

    key : str, optional
        The name of a particular WCS transform to use.  This may be
        either ``' '`` or ``'A'``-``'Z'`` and corresponds to the
        ``\"a\"`` part of the ``CTYPEia`` cards.  *key* may only be
        provided if *header* is also provided.

    minerr : float, optional
        The minimum value a distortion correction must have in order
        to be applied. If the value of ``CQERRja`` is smaller than
        *minerr*, the corresponding distortion is not applied.

    relax : bool or int, optional
        Degree of permissiveness:

        - `True` (default): Admit all recognized informal extensions
          of the WCS standard.

        - `False`: Recognize only FITS keywords defined by the
          published WCS standard.

        - `int`: a bit field selecting specific extensions to accept.
          See :ref:`relaxread` for details.

    naxis : int or sequence, optional
        Extracts specific coordinate axes using
        :meth:`~astropy.wcs.Wcsprm.sub`.  If a header is provided, and
        *naxis* is not ``None``, *naxis* will be passed to
        :meth:`~astropy.wcs.Wcsprm.sub` in order to select specific
        axes from the header.  See :meth:`~astropy.wcs.Wcsprm.sub` for
        more details about this parameter.

    keysel : sequence of flags, optional
        A sequence of flags used to select the keyword types
        considered by wcslib.  When ``None``, only the standard image
        header keywords are considered (and the underlying wcspih() C
        function is called).  To use binary table image array or pixel
        list keywords, *keysel* must be set.

        Each element in the list should be one of the following
        strings:

        - 'image': Image header keywords

        - 'binary': Binary table image array keywords

        - 'pixel': Pixel list keywords

        Keywords such as ``EQUIna`` or ``RFRQna`` that are common to
        binary table image arrays and pixel lists (including
        ``WCSNna`` and ``TWCSna``) are selected by both 'binary' and
        'pixel'.

    colsel : sequence of int, optional
        A sequence of table column numbers used to restrict the WCS
        transformations considered to only those pertaining to the
        specified columns.  If `None`, there is no restriction.

    fix : bool, optional
        When `True` (default), call `~astropy.wcs.Wcsprm.fix` on
        the resulting object to fix any non-standard uses in the
        header.  `FITSFixedWarning` Warnings will be emitted if any
        changes were made.

    translate_units : str, optional
        Specify which potentially unsafe translations of non-standard
        unit strings to perform.  By default, performs none.  See
        `WCS.fix` for more information about this parameter.  Only
        effective when ``fix`` is `True`.

    Raises
    ------
    MemoryError
         Memory allocation failed.

    ValueError
         Invalid key.

    KeyError
         Key not found in FITS header.

    AssertionError
         Lookup table distortion present in the header but *fobj* was
         not provided.

    Notes
    -----

    1. astropy.wcs supports arbitrary *n* dimensions for the core WCS
       (the transformations handled by WCSLIB).  However, the
       `distortion paper`_ lookup table and `SIP`_ distortions must be
       two dimensional.  Therefore, if you try to create a WCS object
       where the core WCS has a different number of dimensions than 2
       and that object also contains a `distortion paper`_ lookup
       table or `SIP`_ distortion, a `~.exceptions.ValueError`
       exception will be raised.  To avoid this, consider using the
       *naxis* kwarg to select two dimensions from the core WCS.

    2. The number of coordinate axes in the transformation is not
       determined directly from the ``NAXIS`` keyword but instead from
       the highest of:

           - ``NAXIS`` keyword

           - ``WCSAXESa`` keyword

           - The highest axis number in any parameterized WCS keyword.
             The keyvalue, as well as the keyword, must be
             syntactically valid otherwise it will not be considered.

       If none of these keyword types is present, i.e. if the header
       only contains auxiliary WCS keywords for a particular
       coordinate representation, then no coordinate description is
       constructed for it.

       The number of axes, which is set as the ``naxis`` member, may
       differ for different coordinate representations of the same
       image.

    3. When the header includes duplicate keywords, in most cases the
       last encountered is used.

    4. `~astropy.wcs.Wcsprm.set` is called immediately after
       construction, so any invalid keywords or transformations will
       be raised by the constructor, not when subsequently calling a
       transformation method.

    """

    def __init__(self, header=None, fobj=None, key=' ', minerr=0.0,
                 relax=True, naxis=None, keysel=None, colsel=None,
                 fix=True, translate_units='', _do_set=True):
        close_fds = []

        if header is None:
            if naxis is None:
                naxis = 2
            wcsprm = _wcs.Wcsprm(header=None, key=key,
                                 relax=relax, naxis=naxis)
            self.naxis = wcsprm.naxis
            # Set some reasonable defaults.
            det2im = (None, None)
            cpdis = (None, None)
            sip = None
        else:
            keysel_flags = _parse_keysel(keysel)

            if isinstance(header, (six.text_type, six.binary_type)):
                try:
                    is_path = (possible_filename(header) and
                               os.path.exists(header))
                except (IOError, ValueError):
                    is_path = False

                if is_path:
                    if fobj is not None:
                        raise ValueError(
                            "Can not provide both a FITS filename to "
                            "argument 1 and a FITS file object to argument 2")
                    fobj = fits.open(header)
                    close_fds.append(fobj)
                    header = fobj[0].header
                    header_string = header.tostring()
                else:
                    header_string = header
            elif isinstance(header, fits.Header):
                header_string = header.tostring()
            else:
                try:
                    # Accept any dict-like object
                    new_header = fits.Header()
                    for dict_key in header.keys():
                        new_header[dict_key] = header[dict_key]
                    header_string = new_header.tostring()
                except TypeError:
                    raise TypeError(
                        "header must be a string, an astropy.io.fits.Header "
                        "object, or a dict-like object")

            # Importantly, header is a *copy* of the passed-in header
            # because we will be modifying it
            if isinstance(header_string, six.text_type):
                header_bytes = header_string.encode('ascii')
                header_string = header_string
            else:
                header_bytes = header_string
                header_string = header_string.decode('ascii')

            header = fits.Header.fromstring(header_string)

            if naxis is None:
                self.naxis = header.get('NAXIS', 2)
            else:
                self.naxis = naxis
            if self.naxis == 0:
                self.naxis = 2
            det2im = self._read_det2im_kw(header, fobj, err=minerr)
            cpdis = self._read_distortion_kw(
                header, fobj, dist='CPDIS', err=minerr)
            sip = self._read_sip_kw(header)

            header_string = header.tostring()
            header_string = header_string.replace('END' + ' ' * 77, '')

            if isinstance(header_string, six.text_type):
                header_bytes = header_string.encode('ascii')
                header_string = header_string
            else:
                header_bytes = header_string
                header_string = header_string.decode('ascii')

            try:
                wcsprm = _wcs.Wcsprm(header=header_bytes, key=key,
                                     relax=relax, keysel=keysel_flags,
                                     colsel=colsel)
            except _wcs.NoWcsKeywordsFoundError:
                # The header may have SIP or distortions, but no core
                # WCS.  That isn't an error -- we want a "default"
                # (identity) core Wcs transformation in that case.
                if colsel is None:
                    wcsprm = _wcs.Wcsprm(header=None, key=key,
                                         relax=relax, keysel=keysel_flags,
                                         colsel=colsel)
                else:
                    raise

            if naxis is not None:
                wcsprm = wcsprm.sub(naxis)
            self.naxis = wcsprm.naxis

            if (wcsprm.naxis != 2 and
                (det2im[0] or det2im[1] or cpdis[0] or cpdis[1] or sip)):
                raise ValueError(
                    """
FITS WCS distortion paper lookup tables and SIP distortions only work
in 2 dimensions.  However, WCSLIB has detected {0} dimensions in the
core WCS keywords.  To use core WCS in conjunction with FITS WCS
distortion paper lookup tables or SIP distortion, you must select or
reduce these to 2 dimensions using the naxis kwarg.
""".format(wcsprm.naxis))

            header_naxis = header.get('NAXIS', None)
            if header_naxis is not None and header_naxis < wcsprm.naxis:
                warnings.warn(
                    "The WCS transformation has more axes ({0:d}) than the "
                    "image it is associated with ({1:d})".format(
                        wcsprm.naxis, header_naxis), FITSFixedWarning)

        self._get_naxis(header)
        WCSBase.__init__(self, sip, cpdis, wcsprm, det2im)

        if fix:
            self.fix(translate_units=translate_units)

        if _do_set:
            self.wcs.set()

        for fd in close_fds:
            fd.close()

    def __copy__(self):
        new_copy = self.__class__()
        WCSBase.__init__(new_copy, self.sip,
                         (self.cpdis1, self.cpdis2),
                         self.wcs,
                         (self.det2im1, self.det2im2))
        new_copy.__dict__.update(self.__dict__)
        return new_copy

    def __deepcopy__(self, memo):
        new_copy = self.__class__()
        new_copy.naxis = copy.deepcopy(self.naxis, memo)
        WCSBase.__init__(new_copy, copy.deepcopy(self.sip, memo),
                         (copy.deepcopy(self.cpdis1, memo),
                          copy.deepcopy(self.cpdis2, memo)),
                         copy.deepcopy(self.wcs, memo),
                         (copy.deepcopy(self.det2im1, memo),
                          copy.deepcopy(self.det2im2, memo)))
        for key in self.__dict__:
            val = self.__dict__[key]
            new_copy.__dict__[key] = copy.deepcopy(val, memo)
        return new_copy

    def copy(self):
        """
        Return a shallow copy of the object.

        Convenience method so user doesn't have to import the
        :mod:`copy` stdlib module.
        """
        return copy.copy(self)

    def deepcopy(self):
        """
        Return a deep copy of the object.

        Convenience method so user doesn't have to import the
        :mod:`copy` stdlib module.
        """
        return copy.deepcopy(self)

    def sub(self, axes=None):
        copy = self.deepcopy()
        copy.wcs = self.wcs.sub(axes)
        copy.naxis = copy.wcs.naxis
        return copy
    if _wcs is not None:
        sub.__doc__ = _wcs.Wcsprm.sub.__doc__

    def _fix_scamp(self):
        """
        Remove SCAMP's PVi_m distortion parameters if SIP distortion parameters
        are also present. Some projects (e.g., Palomar Transient Factory)
        convert SCAMP's distortion parameters (which abuse the PVi_m cards) to
        SIP. However, wcslib gets confused by the presence of both SCAMP and
        SIP distortion parameters.

        See https://github.com/astropy/astropy/issues/299.
        """
        # Nothing to be done if no WCS attached
        if self.wcs is None:
            return

        # Nothing to be done if no PV parameters attached
        pv = self.wcs.get_pv()
        if not pv:
            return

        # Nothing to be done if axes don't use SIP distortion parameters
        if not all(ctype.endswith('-SIP') for ctype in self.wcs.ctype):
            return

        # Nothing to be done if any radial terms are present...
        # Loop over list to find any radial terms.
        # Certain values of the `j' index are used for storing
        # radial terms; refer to Equation (1) in
        # <http://web.ipac.caltech.edu/staff/shupe/reprints/SIP_to_PV_SPIE2012.pdf>.
        pv = np.asarray(pv)
        # Loop over distinct values of `i' index
        for i in set(pv[:, 0]):
            # Get all values of `j' index for this value of `i' index
            js = set(pv[:, 1][pv[:, 0] == i])
            # Find max value of `j' index
            max_j = max(js)
            for j in (3, 11, 23, 39):
                if j < max_j and j in js:
                    return

        self.wcs.set_pv([])
        warnings.warn("Removed redundant SCAMP distortion parameters " +
            "because SIP parameters are also present", FITSFixedWarning)

    def fix(self, translate_units='', naxis=None):
        """
        Perform the fix operations from wcslib, and warn about any
        changes it has made.

        Parameters
        ----------
        translate_units : str, optional
            Specify which potentially unsafe translations of
            non-standard unit strings to perform.  By default,
            performs none.

            Although ``"S"`` is commonly used to represent seconds,
            its translation to ``"s"`` is potentially unsafe since the
            standard recognizes ``"S"`` formally as Siemens, however
            rarely that may be used.  The same applies to ``"H"`` for
            hours (Henry), and ``"D"`` for days (Debye).

            This string controls what to do in such cases, and is
            case-insensitive.

            - If the string contains ``"s"``, translate ``"S"`` to
              ``"s"``.

            - If the string contains ``"h"``, translate ``"H"`` to
              ``"h"``.

            - If the string contains ``"d"``, translate ``"D"`` to
              ``"d"``.

            Thus ``''`` doesn't do any unsafe translations, whereas
            ``'shd'`` does all of them.

        naxis : int array[naxis], optional
            Image axis lengths.  If this array is set to zero or
            ``None``, then `~astropy.wcs.Wcsprm.cylfix` will not be
            invoked.
        """
        if self.wcs is not None:
            self._fix_scamp()
            fixes = self.wcs.fix(translate_units, naxis)
            for key, val in six.iteritems(fixes):
                if val != "No change":
                    warnings.warn(
                        ("'{0}' made the change '{1}'.").
                        format(key, val),
                        FITSFixedWarning)

    def calc_footprint(self, header=None, undistort=True, axes=None, center=True):
        """
        Calculates the footprint of the image on the sky.

        A footprint is defined as the positions of the corners of the
        image on the sky after all available distortions have been
        applied.

        Parameters
        ----------
        header : `~astropy.io.fits.Header` object, optional
            Used to get ``NAXIS1`` and ``NAXIS2``
            header and axes are mutually exclusive, alternative ways
            to provide the same information.

        undistort : bool, optional
            If `True`, take SIP and distortion lookup table into
            account

        axes : length 2 sequence ints, optional
            If provided, use the given sequence as the shape of the
            image.  Otherwise, use the ``NAXIS1`` and ``NAXIS2``
            keywords from the header that was used to create this
            `WCS` object.

        center : bool, optional
            If `True` use the center of the pixel, otherwise use the corner.

        Returns
        -------
        coord : (4, 2) array of (*x*, *y*) coordinates.
            The order is counter-clockwise starting with the bottom left corner.
        """
        if axes is not None:
            naxis1, naxis2 = axes
        else:
            if header is None:
                try:
                    # classes that inherit from WCS and define naxis1/2
                    # do not require a header parameter
                    naxis1 = self._naxis1
                    naxis2 = self._naxis2
                except AttributeError:
                    warnings.warn("Need a valid header in order to calculate footprint\n", AstropyUserWarning)
                    return None
            else:
                naxis1 = header.get('NAXIS1', None)
                naxis2 = header.get('NAXIS2', None)

        if naxis1 is None or naxis2 is None:
            raise ValueError(
                    "Image size could not be determined.")

        if center == True:
            corners = np.array([[1, 1],
                                [1, naxis2],
                                [naxis1, naxis2],
                                [naxis1, 1]], dtype = np.float64)
        else:
            corners = np.array([[0.5, 0.5],
                                [0.5, naxis2 + 0.5],
                                [naxis1 + 0.5, naxis2 + 0.5],
                                [naxis1 + 0.5, 0.5]], dtype = np.float64)

        if undistort:
            return self.all_pix2world(corners, 1)
        else:
            return self.wcs_pix2world(corners, 1)

    def _read_det2im_kw(self, header, fobj, err=0.0):
        """
        Create a `distortion paper`_ type lookup table for detector to
        image plane correction.
        """
        if fobj is None:
            return (None, None)

        if not isinstance(fobj, fits.HDUList):
            return (None, None)

        try:
            axiscorr = header[str('AXISCORR')]
            d2imdis = self._read_d2im_old_format(header, fobj, axiscorr)
            return d2imdis
        except KeyError:
            pass

        dist = 'D2IMDIS'
        d_kw = 'D2IM'
        err_kw = 'D2IMERR'
        tables = {}
        for i in range(1, self.naxis + 1):
            d_error = header.get(err_kw + str(i), 0.0)
            if d_error < err:
                tables[i] = None
                continue
            distortion = dist + str(i)
            if distortion in header:
                dis = header[distortion].lower()
                if dis == 'lookup':
                    assert isinstance(fobj, fits.HDUList), ('An astropy.io.fits.HDUList'
                                'is required for Lookup table distortion.')
                    dp = (d_kw + str(i)).strip()
                    d_extver = header.get(dp + '.EXTVER', 1)
                    if i == header[dp + '.AXIS.{0:d}'.format(i)]:
                        d_data = fobj[str('D2IMARR'), d_extver].data
                    else:
                        d_data = (fobj[str('D2IMARR'), d_extver].data).transpose()
                    d_header = fobj[str('D2IMARR'), d_extver].header
                    d_crpix = (d_header.get(str('CRPIX1'), 0.0), d_header.get(str('CRPIX2'), 0.0))
                    d_crval = (d_header.get(str('CRVAL1'), 0.0), d_header.get(str('CRVAL2'), 0.0))
                    d_cdelt = (d_header.get(str('CDELT1'), 1.0), d_header.get(str('CDELT2'), 1.0))
                    d_lookup = DistortionLookupTable(d_data, d_crpix,
                                                     d_crval, d_cdelt)
                    tables[i] = d_lookup
                else:
                    warnings.warn('Polynomial distortion is not implemented.\n', AstropyUserWarning)
            else:
                tables[i] = None
        if not tables:
            return (None, None)
        else:
            return (tables.get(1), tables.get(2))

    def _read_d2im_old_format(self, header, fobj, axiscorr):
        warnings.warn("The use of ``AXISCORR`` for D2IM correction has been deprecated."
                      "`~astropy.wcs` will read in files with ``AXISCORR`` but ``to_fits()`` will write "
                      "out files without it.",
                      AstropyDeprecationWarning)
        cpdis = [None, None]
        crpix = [0., 0.]
        crval = [0., 0.]
        cdelt = [1., 1.]
        try:
            d2im_data = fobj[(str('D2IMARR'), 1)].data
        except KeyError:
            return (None, None)
        except AttributeError:
            return (None, None)

        d2im_data = np.array([d2im_data])
        d2im_hdr = fobj[(str('D2IMARR'), 1)].header
        naxis = d2im_hdr[str('NAXIS')]

        for i in range(1, naxis + 1):
            crpix[i - 1] = d2im_hdr.get(str('CRPIX') + str(i), 0.0)
            crval[i - 1] = d2im_hdr.get(str('CRVAL') + str(i), 0.0)
            cdelt[i - 1] = d2im_hdr.get(str('CDELT') + str(i), 1.0)

        cpdis = DistortionLookupTable(d2im_data, crpix, crval, cdelt)

        if axiscorr == 1:
            return (cpdis, None)
        elif axiscorr == 2:
            return (None, cpdis)
        else:
            warnings.warn("Expected AXISCORR to be 1 or 2", AstropyUserWarning)
            return (None, None)

    def _write_det2im(self, hdulist):
        """
        Writes a `distortion paper`_ type lookup table to the given
        `astropy.io.fits.HDUList`.
        """

        if self.det2im1 is None and self.det2im2 is None:
            return
        dist = 'D2IMDIS'
        d_kw = 'D2IM'
        err_kw = 'D2IMERR'

        def write_d2i(num, det2im):
            if det2im is None:
                return
            str('{0}{1:d}').format(dist, num),
            hdulist[0].header[str('{0}{1:d}').format(dist, num)] = (
                'LOOKUP', 'Detector to image correction type')
            hdulist[0].header[str('{0}{1:d}.EXTVER').format(d_kw, num)] = (
                num, 'Version number of WCSDVARR extension')
            hdulist[0].header[str('{0}{1:d}.NAXES').format(d_kw, num)] = (
                len(det2im.data.shape), 'Number of independent variables in d2im function')
            for i in range(det2im.data.ndim):
                hdulist[0].header[str('{0}{1:d}.AXIS.{2:d}').format(d_kw, num, i + 1)] = (
                    i + 1, 'Axis number of the jth independent variable in a d2im function')

            image = fits.ImageHDU(det2im.data, name=str('D2IMARR'))
            header = image.header

            header[str('CRPIX1')] = (det2im.crpix[0],
                                     'Coordinate system reference pixel')
            header[str('CRPIX2')] = (det2im.crpix[1],
                                     'Coordinate system reference pixel')
            header[str('CRVAL1')] = (det2im.crval[0],
                                     'Coordinate system value at reference pixel')
            header[str('CRVAL2')] = (det2im.crval[1],
                                     'Coordinate system value at reference pixel')
            header[str('CDELT1')] = (det2im.cdelt[0],
                                     'Coordinate increment along axis')
            header[str('CDELT2')] = (det2im.cdelt[1],
                                     'Coordinate increment along axis')
            image.update_ext_version(
                int(hdulist[0].header[str('{0}{1:d}.EXTVER').format(d_kw, num)]))
            hdulist.append(image)
        write_d2i(1, self.det2im1)
        write_d2i(2, self.det2im2)

    def _read_distortion_kw(self, header, fobj, dist='CPDIS', err=0.0):
        """
        Reads `distortion paper`_ table-lookup keywords and data, and
        returns a 2-tuple of `~astropy.wcs.DistortionLookupTable`
        objects.

        If no `distortion paper`_ keywords are found, ``(None, None)``
        is returned.
        """
        if isinstance(header, (six.text_type, six.binary_type)):
            return (None, None)

        if dist == 'CPDIS':
            d_kw = str('DP')
            err_kw = str('CPERR')
        else:
            d_kw = str('DQ')
            err_kw = str('CQERR')

        tables = {}
        for i in range(1, self.naxis + 1):
            d_error = header.get(err_kw + str(i), 0.0)
            if d_error < err:
                tables[i] = None
                continue
            distortion = dist + str(i)
            if distortion in header:
                dis = header[distortion].lower()
                if dis == 'lookup':
                    assert isinstance(fobj, fits.HDUList), \
                        'An astropy.io.fits.HDUList is required for ' + \
                        'Lookup table distortion.'
                    dp = (d_kw + str(i)).strip()
                    d_extver = header.get(dp + str('.EXTVER'), 1)
                    if i == header[dp + str('.AXIS.') + str(i)]:
                        d_data = fobj[str('WCSDVARR'), d_extver].data
                    else:
                        d_data = (fobj[str('WCSDVARR'), d_extver].data).transpose()

                    d_header = fobj[str('WCSDVARR'), d_extver].header
                    d_crpix = (d_header.get(str('CRPIX1'), 0.0),
                               d_header.get(str('CRPIX2'), 0.0))
                    d_crval = (d_header.get(str('CRVAL1'), 0.0),
                               d_header.get(str('CRVAL2'), 0.0))
                    d_cdelt = (d_header.get(str('CDELT1'), 1.0),
                               d_header.get(str('CDELT2'), 1.0))
                    d_lookup = DistortionLookupTable(d_data, d_crpix, d_crval, d_cdelt)
                    tables[i] = d_lookup
                else:
                    warnings.warn('Polynomial distortion is not implemented.\n', AstropyUserWarning)
            else:
                tables[i] = None

        if not tables:
            return (None, None)
        else:
            return (tables.get(1), tables.get(2))

    def _write_distortion_kw(self, hdulist, dist='CPDIS'):
        """
        Write out `distortion paper`_ keywords to the given
        `fits.HDUList`.
        """
        if self.cpdis1 is None and self.cpdis2 is None:
            return

        if dist == 'CPDIS':
            d_kw = str('DP')
            err_kw = str('CPERR')
        else:
            d_kw = str('DQ')
            err_kw = str('CQERR')

        def write_dist(num, cpdis):
            if cpdis is None:
                return

            hdulist[0].header[str('{0}{1:d}').format(dist, num)] = (
                'LOOKUP', 'Prior distortion function type')
            hdulist[0].header[str('{0}{1:d}.EXTVER').format(d_kw, num)] = (
                num, 'Version number of WCSDVARR extension')
            hdulist[0].header[str('{0}{1:d}.NAXES').format(d_kw, num)] = (
                len(cpdis.data.shape), 'Number of independent variables in distortion function')

            for i in range(cpdis.data.ndim):
                hdulist[0].header[str('{0}{1:d}.AXIS.{2:d}').format(d_kw, num, i + 1)] = (
                    i + 1,
                    'Axis number of the jth independent variable in a distortion function')

            image = fits.ImageHDU(cpdis.data, name=str('WCSDVARR'))
            header = image.header

            header[str('CRPIX1')] = (cpdis.crpix[0], 'Coordinate system reference pixel')
            header[str('CRPIX2')] = (cpdis.crpix[1], 'Coordinate system reference pixel')
            header[str('CRVAL1')] = (cpdis.crval[0], 'Coordinate system value at reference pixel')
            header[str('CRVAL2')] = (cpdis.crval[1], 'Coordinate system value at reference pixel')
            header[str('CDELT1')] = (cpdis.cdelt[0], 'Coordinate increment along axis')
            header[str('CDELT2')] = (cpdis.cdelt[1], 'Coordinate increment along axis')
            image.update_ext_version(
                int(hdulist[0].header[str('{0}{1:d}.EXTVER').format(d_kw, num)]))
            hdulist.append(image)

        write_dist(1, self.cpdis1)
        write_dist(2, self.cpdis2)

    def _read_sip_kw(self, header):
        """
        Reads `SIP`_ header keywords and returns a `~astropy.wcs.Sip`
        object.

        If no `SIP`_ header keywords are found, ``None`` is returned.
        """
        if isinstance(header, (six.text_type, six.binary_type)):
            # TODO: Parse SIP from a string without pyfits around
            return None

        if str("A_ORDER") in header and header[str('A_ORDER')] > 1:
            if str("B_ORDER") not in header:
                raise ValueError(
                    "A_ORDER provided without corresponding B_ORDER "
                    "keyword for SIP distortion")

            m = int(header[str("A_ORDER")])
            a = np.zeros((m + 1, m + 1), np.double)
            for i in range(m + 1):
                for j in range(m - i + 1):
                    key = str("A_{0}_{1}").format(i, j)
                    if key in header:
                        a[i, j] = header[key]
                        del header[key]

            m = int(header[str("B_ORDER")])
            if m > 1:
                b = np.zeros((m + 1, m + 1), np.double)
                for i in range(m + 1):
                    for j in range(m - i + 1):
                        key = str("B_{0}_{1}").format(i, j)
                        if key in header:
                            b[i, j] = header[key]
                            del header[key]
            else:
                a = None
                b = None

            del header[str('A_ORDER')]
            del header[str('B_ORDER')]
        elif str("B_ORDER") in header and header[str('B_ORDER')] > 1:
            raise ValueError(
                "B_ORDER provided without corresponding A_ORDER " +
                "keyword for SIP distortion")
        else:
            a = None
            b = None

        if str("AP_ORDER") in header and header[str('AP_ORDER')] > 1:
            if str("BP_ORDER") not in header:
                raise ValueError(
                    "AP_ORDER provided without corresponding BP_ORDER "
                    "keyword for SIP distortion")

            m = int(header[str("AP_ORDER")])
            ap = np.zeros((m + 1, m + 1), np.double)
            for i in range(m + 1):
                for j in range(m - i + 1):
                    key = str("AP_{0}_{1}").format(i, j)
                    if key in header:
                        ap[i, j] = header[key]
                        del header[key]

            m = int(header[str("BP_ORDER")])
            if m > 1:
                bp = np.zeros((m + 1, m + 1), np.double)
                for i in range(m + 1):
                    for j in range(m - i + 1):
                        key = str("BP_{0}_{1}").format(i, j)
                        if key in header:
                            bp[i, j] = header[key]
                            del header[key]
            else:
                ap = None
                bp = None

            del header[str('AP_ORDER')]
            del header[str('BP_ORDER')]
        elif str("BP_ORDER") in header and header[str('BP_ORDER')] > 1:
            raise ValueError(
                "BP_ORDER provided without corresponding AP_ORDER "
                "keyword for SIP distortion")
        else:
            ap = None
            bp = None

        if a is None and b is None and ap is None and bp is None:
            return None

        if str("CRPIX1") not in header or str("CRPIX2") not in header:
            raise ValueError(
                "Header has SIP keywords without CRPIX keywords")

        crpix1 = header.get("CRPIX1")
        crpix2 = header.get("CRPIX2")

        return Sip(a, b, ap, bp, (crpix1, crpix2))

    def _write_sip_kw(self):
        """
        Write out SIP keywords.  Returns a dictionary of key-value
        pairs.
        """
        if self.sip is None:
            return {}

        keywords = {}

        def write_array(name, a):
            if a is None:
                return
            size = a.shape[0]
            keywords[str('{0}_ORDER').format(name)] = size - 1
            for i in range(size):
                for j in range(size - i):
                    if a[i, j] != 0.0:
                        keywords[
                            str('{0}_{1:d}_{2:d}').format(name, i, j)] = a[i, j]

        write_array(str('A'), self.sip.a)
        write_array(str('B'), self.sip.b)
        write_array(str('AP'), self.sip.ap)
        write_array(str('BP'), self.sip.bp)

        return keywords

    def _denormalize_sky(self, sky):
        if self.wcs.lngtyp != 'RA':
            raise ValueError(
                "WCS does not have longitude type of 'RA', therefore " +
                "(ra, dec) data can not be used as input")
        if self.wcs.lattyp != 'DEC':
            raise ValueError(
                "WCS does not have longitude type of 'DEC', therefore " +
                "(ra, dec) data can not be used as input")
        if self.wcs.naxis == 2:
            if self.wcs.lng == 0 and self.wcs.lat == 1:
                return sky
            elif self.wcs.lng == 1 and self.wcs.lat == 0:
                # Reverse the order of the columns
                return sky[:, ::-1]
            else:
                raise ValueError(
                    "WCS does not have longitude and latitude celestial " +
                    "axes, therefore (ra, dec) data can not be used as input")
        else:
            if self.wcs.lng < 0 or self.wcs.lat < 0:
                raise ValueError(
                    "WCS does not have both longitude and latitude "
                    "celestial axes, therefore (ra, dec) data can not be " +
                    "used as input")
            out = np.zeros((sky.shape[0], self.wcs.naxis))
            out[:, self.wcs.lng] = sky[:, 0]
            out[:, self.wcs.lat] = sky[:, 1]
            return out

    def _normalize_sky(self, sky):
        if self.wcs.lngtyp != 'RA':
            raise ValueError(
                "WCS does not have longitude type of 'RA', therefore " +
                "(ra, dec) data can not be returned")
        if self.wcs.lattyp != 'DEC':
            raise ValueError(
                "WCS does not have longitude type of 'DEC', therefore " +
                "(ra, dec) data can not be returned")
        if self.wcs.naxis == 2:
            if self.wcs.lng == 0 and self.wcs.lat == 1:
                return sky
            elif self.wcs.lng == 1 and self.wcs.lat == 0:
                # Reverse the order of the columns
                return sky[:, ::-1]
            else:
                raise ValueError(
                    "WCS does not have longitude and latitude celestial "
                    "axes, therefore (ra, dec) data can not be returned")
        else:
            if self.wcs.lng < 0 or self.wcs.lat < 0:
                raise ValueError(
                    "WCS does not have both longitude and latitude celestial "
                    "axes, therefore (ra, dec) data can not be returned")
            out = np.empty((sky.shape[0], 2))
            out[:, 0] = sky[:, self.wcs.lng]
            out[:, 1] = sky[:, self.wcs.lat]
            return out

    def _array_converter(self, func, sky, *args, **kwargs):
        """
        A helper function to support reading either a pair of arrays
        or a single Nx2 array.
        """
        ra_dec_order = kwargs.pop('ra_dec_order', False)
        if len(kwargs):
            raise TypeError("Unexpected keyword argument {0!r}".format(
                kwargs.keys()[0]))

        def _return_list_of_arrays(axes, origin):
            try:
                axes = np.broadcast_arrays(*axes)
            except ValueError:
                raise ValueError(
                    "Coordinate arrays are not broadcastable to each other")

            xy = np.hstack([x.reshape((x.size, 1)) for x in axes])

            if ra_dec_order and sky == 'input':
                xy = self._denormalize_sky(xy)
            output = func(xy, origin)
            if ra_dec_order and sky == 'output':
                output = self._normalize_sky(output)
                return (output[:, 0].reshape(axes[0].shape),
                        output[:, 1].reshape(axes[0].shape))
            return [output[:, i].reshape(axes[0].shape)
                    for i in range(output.shape[1])]

        def _return_single_array(xy, origin):
            if xy.shape[-1] != self.naxis:
                raise ValueError(
                    "When providing two arguments, the array must be "
                    "of shape (N, {0})".format(self.naxis))
            if ra_dec_order and sky == 'input':
                xy = self._denormalize_sky(xy)
            result = func(xy, origin)
            if ra_dec_order and sky == 'output':
                result = self._normalize_sky(result)
            return result

        if len(args) == 2:
            try:
                xy, origin = args
                xy = np.asarray(xy)
                origin = int(origin)
            except:
                raise TypeError(
                    "When providing two arguments, they must be "
                    "(coords[N][{0}], origin)".format(self.naxis))
            if self.naxis == 1 and len(xy.shape) == 1:
                return _return_list_of_arrays([xy], origin)
            return _return_single_array(xy, origin)

        elif len(args) == self.naxis + 1:
            axes = args[:-1]
            origin = args[-1]
            try:
                axes = [np.asarray(x) for x in axes]
                origin = int(origin)
            except:
                raise TypeError(
                    "When providing more than two arguments, they must be " +
                    "a 1-D array for each axis, followed by an origin.")

            return _return_list_of_arrays(axes, origin)

        raise TypeError(
            "WCS projection has {0} dimensions, so expected 2 (an Nx{0} array "
            "and the origin argument) or {1} arguments (the position in each "
            "dimension, and the origin argument). Instead, {2} arguments were "
            "given.".format(
                self.naxis, self.naxis + 1, len(args)))

    def all_pix2world(self, *args, **kwargs):
        return self._array_converter(
            self._all_pix2world, 'output', *args, **kwargs)
    all_pix2world.__doc__ = """
        Transforms pixel coordinates to world coordinates.

        Performs all of the following in series:

            - Detector to image plane correction (if present in the
              FITS file)

            - `SIP`_ distortion correction (if present in the FITS
              file)

            - `distortion paper`_ table-lookup correction (if present
              in the FITS file)

            - `wcslib`_ "core" WCS transformation

        Parameters
        ----------
        {0}

            For a transformation that is not two-dimensional, the
            two-argument form must be used.

        {1}

        Returns
        -------

        {2}

        Notes
        -----
        The order of the axes for the result is determined by the
        ``CTYPEia`` keywords in the FITS header, therefore it may not
        always be of the form (*ra*, *dec*).  The
        `~astropy.wcs.Wcsprm.lat`, `~astropy.wcs.Wcsprm.lng`,
        `~astropy.wcs.Wcsprm.lattyp` and `~astropy.wcs.Wcsprm.lngtyp`
        members can be used to determine the order of the axes.

        Raises
        ------
        MemoryError
            Memory allocation failed.

        SingularMatrixError
            Linear transformation matrix is singular.

        InconsistentAxisTypesError
            Inconsistent or unrecognized coordinate axis types.

        ValueError
            Invalid parameter value.

        ValueError
            Invalid coordinate transformation parameters.

        ValueError
            x- and y-coordinate arrays are not the same size.

        InvalidTransformError
            Invalid coordinate transformation parameters.

        InvalidTransformError
            Ill-conditioned coordinate transformation parameters.
        """.format(__.TWO_OR_MORE_ARGS('naxis', 8),
                   __.RA_DEC_ORDER(8),
                   __.RETURNS('sky coordinates, in degrees', 8))

    def wcs_pix2world(self, *args, **kwargs):
        if self.wcs is None:
            raise ValueError("No basic WCS settings were created.")
        return self._array_converter(
            lambda xy, o: self.wcs.p2s(xy, o)['world'],
            'output', *args, **kwargs)
    wcs_pix2world.__doc__ = """
        Transforms pixel coordinates to world coordinates by doing
        only the basic `wcslib`_ transformation.

        No `SIP`_ or `distortion paper`_ table lookup correction is
        applied.  To perform distortion correction, see
        `~astropy.wcs.WCS.all_pix2world`,
        `~astropy.wcs.WCS.sip_pix2foc`, `~astropy.wcs.WCS.p4_pix2foc`,
        or `~astropy.wcs.WCS.pix2foc`.

        Parameters
        ----------
        {0}

            For a transformation that is not two-dimensional, the
            two-argument form must be used.

        {1}

        Returns
        -------

        {2}

        Raises
        ------
        MemoryError
            Memory allocation failed.

        SingularMatrixError
            Linear transformation matrix is singular.

        InconsistentAxisTypesError
            Inconsistent or unrecognized coordinate axis types.

        ValueError
            Invalid parameter value.

        ValueError
            Invalid coordinate transformation parameters.

        ValueError
            x- and y-coordinate arrays are not the same size.

        InvalidTransformError
            Invalid coordinate transformation parameters.

        InvalidTransformError
            Ill-conditioned coordinate transformation parameters.

        Notes
        -----
        The order of the axes for the result is determined by the
        ``CTYPEia`` keywords in the FITS header, therefore it may not
        always be of the form (*ra*, *dec*).  The
        `~astropy.wcs.Wcsprm.lat`, `~astropy.wcs.Wcsprm.lng`,
        `~astropy.wcs.Wcsprm.lattyp` and `~astropy.wcs.Wcsprm.lngtyp`
        members can be used to determine the order of the axes.

        """.format(__.TWO_OR_MORE_ARGS('naxis', 8),
                   __.RA_DEC_ORDER(8),
                   __.RETURNS('world coordinates, in degrees', 8))


    def _all_world2pix(self, world, origin, tolerance, maxiter, adaptive,
                       detect_divergence, quiet):
        #############################################################
        ##          DESCRIPTION OF THE NUMERICAL METHOD            ##
        #############################################################
        # In this section I will outline the method of solving
        # the inverse problem of converting world coordinates to
        # pixel coordinates (*inverse* of the direct transformation
        # `all_pix2world`) and I will summarize some of the aspects
        # of the method proposed here and some of the issues of the
        # original `all_world2pix` (in relation to this method)
        # discussed in https://github.com/astropy/astropy/issues/1977
        # A more detailed discussion can be found here:
        # https://github.com/astropy/astropy/pull/2373
        #
        #
        #                  ### Background ###
        #
        #
        # I will refer here to the [SIP Paper]
        # (http://fits.gsfc.nasa.gov/registry/sip/SIP_distortion_v1_0.pdf).
        # According to this paper, the effect of distortions as
        # described in *their* equation (1) is:
        #
        # (1)   x = CD*(u+f(u)),
        #
        # where `x` is a *vector* of "intermediate spherical
        # coordinates" (equivalent to (x,y) in the paper) and `u`
        # is a *vector* of "pixel coordinates", and `f` is a vector
        # function describing geometrical distortions
        # (see equations 2 and 3 in SIP Paper.
        # However, I prefer to use `w` for "intermediate world
        # coordinates", `x` for pixel coordinates, and assume that
        # transformation `W` performs the **linear**
        # (CD matrix + projection onto celestial sphere) part of the
        # conversion from pixel coordinates to world coordinates.
        # Then we can re-write (1) as:
        #
        # (2)   w = W*(x+f(x)) = T(x)
        #
        # In `astropy.wcs.WCS` transformation `W` is represented by
        # the `wcs_pix2world` member, while the combined ("total")
        # transformation (linear part + distortions) is performed by
        # `all_pix2world`. Below I summarize the notations and their
        # equivalents in `astropy.wcs.WCS`:
        #
        #| Equation term | astropy.WCS/meaning          |
        #| ------------- | ---------------------------- |
        #| `x`           | pixel coordinates            |
        #| `w`           | world coordinates            |
        #| `W`           | `wcs_pix2world()`            |
        #| `W^{-1}`      | `wcs_world2pix()`            |
        #| `T`           | `all_pix2world()`            |
        #| `x+f(x)`      | `pix2foc()`                  |
        #
        #
        #      ### Direct Solving of Equation (2)  ###
        #
        #
        # In order to find the pixel coordinates that correspond to
        # given world coordinates `w`, it is necessary to invert
        # equation (2): `x=T^{-1}(w)`, or solve equation `w==T(x)`
        # for `x`. However, this approach has the following
        # disadvantages:
        #    1. It requires unnecessary transformations (see next
        #       section).
        #    2. It is prone to "RA wrapping" issues as described in
        # https://github.com/astropy/astropy/issues/1977
        # (essentially because `all_pix2world` may return points with
        # a different phase than user's input `w`).
        #
        #
        #      ### Description of the Method Used here ###
        #
        #
        # By applying inverse linear WCS transformation (`W^{-1}`)
        # to both sides of equation (2) and introducing notation `x'`
        # (prime) for the pixels coordinates obtained from the world
        # coordinates by applying inverse *linear* WCS transformation
        # ("focal plane coordinates"):
        #
        # (3)   x' = W^{-1}(w)
        #
        # we obtain the following equation:
        #
        # (4)   x' = x+f(x),
        #
        # or,
        #
        # (5)   x = x'-f(x)
        #
        # This equation is well suited for solving using the method
        # of fixed-point iterations
        # (http://en.wikipedia.org/wiki/Fixed-point_iteration):
        #
        # (6)   x_{i+1} = x'-f(x_i)
        #
        # As an initial value of the pixel coordinate `x_0` we take
        # "focal plane coordinate" `x'=W^{-1}(w)=wcs_world2pix(w)`.
        # We stop iterations when `|x_{i+1}-x_i|<tolerance`. We also
        # consider the process to be diverging if
        # `|x_{i+1}-x_i|>|x_i-x_{i-1}|`
        # **when** `|x_{i+1}-x_i|>=tolerance` (when current
        # approximation is close to the true solution,
        # `|x_{i+1}-x_i|>|x_i-x_{i-1}|` may be due to rounding errors
        # and we ignore such "divergences" when
        # `|x_{i+1}-x_i|<tolerance`). It may appear that checking for
        # `|x_{i+1}-x_i|<tolerance` in order to ignore divergence is
        # unnecessary since the iterative process should stop anyway,
        # however, the proposed implementation of this iterative
        # process is completely vectorized and, therefore, we may
        # continue iterating over *some* points even though they have
        # converged to within a specified tolerance (while iterating
        # over other points that have not yet converged to
        # a solution).
        #
        # In order to efficiently implement iterative process (6)
        # using available methods in `astropy.wcs.WCS`, we add and
        # subtract `x_i` from the right side of equation (6):
        #
        # (7)   x_{i+1} = x'-(x_i+f(x_i))+x_i = x'-pix2foc(x_i)+x_i,
        #
        # where `x'=wcs_world2pix(w)` and it is computed only *once*
        # before the beginning of the iterative process (and we also
        # set `x_0=x'`). By using `pix2foc` at each iteration instead
        # of `all_pix2world` we get about 25% increase in performance
        # (by not performing the linear `W` transformation at each
        # step) and we also avoid the "RA wrapping" issue described
        # above (by working in focal plane coordinates and avoiding
        # pix->world transformations).
        #
        # As an added benefit, the process converges to the correct
        # solution in just one iteration when distortions are not
        # present (compare to
        # https://github.com/astropy/astropy/issues/1977 and
        # https://github.com/astropy/astropy/pull/2294): in this case
        # `pix2foc` is the identical transformation
        # `x_i=pix2foc(x_i)` and from equation (7) we get:
        #
        # x' = x_0 = wcs_world2pix(w)
        # x_1 = x' - pix2foc(x_0) + x_0 = x' - pix2foc(x') + x' = x'
        #     = wcs_world2pix(w) = x_0
        # =>
        # |x_1-x_0| = 0 < tolerance (with tolerance > 0)
        #
        # However, for performance reasons, it is still better to
        # avoid iterations altogether and return the exact linear
        # solution (`wcs_world2pix`) right-away when non-linear
        # distortions are not present by checking that attributes
        # `sip`, `cpdis1`, `cpdis2`, `det2im1`, and `det2im2` are
        # *all* `None`.
        #
        #
        #         ### Outline of the Algorithm ###
        #
        #
        # While the proposed code is relatively long (considering
        # the simplicity of the algorithm), this is due to: 1)
        # checking if iterative solution is necessary at all; 2)
        # checking for divergence; 3) re-implementation of the
        # completely vectorized algorithm as an "adaptive" vectorized
        # algorithm (for cases when some points diverge for which we
        # want to stop iterations). In my tests, the adaptive version
        # of the algorithm is about 50% slower than non-adaptive
        # version for all HST images.
        #
        # The essential part of the vectorized non-adaptive algorithm
        # (without divergence and other checks) can be described
        # as follows:
        #
        #     pix0 = self.wcs_world2pix(world, origin)
        #     pix  = pix0.copy() # 0-order solution
        #
        #     for k in range(maxiter):
        #         # find correction to the previous solution:
        #         dpix = self.pix2foc(pix, origin) - pix0
        #
        #         # compute norm (L2) of the correction:
        #         dn = np.linalg.norm(dpix, axis=1)
        #
        #         # apply correction:
        #         pix -= dpix
        #
        #         # check convergence:
        #         if np.max(dn) < tolerance:
        #             break
        #
        #    return pix
        #
        # Here, the input parameter `world` can be a `MxN` array
        # where `M` is the number of coordinate axes in WCS and `N`
        # is the number of points to be converted simultaneously to
        # image coordinates.
        #
        #
        #                ###  IMPORTANT NOTE:  ###
        #
        # If, in the future releases of the `~astropy.wcs`,
        # `pix2foc` will not apply all the required distortion
        # corrections then in the code below, calls to `pix2foc` will
        # have to be replaced with
        # wcs_world2pix(all_pix2world(pix_list, origin), origin)
        #

        #############################################################
        ##            INITIALIZE ITERATIVE PROCESS:                ##
        #############################################################

        # initial approximation (linear WCS based only)
        pix0 = self.wcs_world2pix(world, origin)

        # Check that an iterative solution is required at all
        # (when any of the non-CD-matrix-based corrections are
        # present). If not required return the initial
        # approximation (pix0).
        if self.sip is None and \
           self.cpdis1 is None and self.cpdis2 is None and \
           self.det2im1 is None and self.det2im2 is None:
            # No non-WCS corrections detected so
            # simply return initial approximation:
            return pix0

        pix = pix0.copy()  # 0-order solution

        # initial correction:
        dpix = self.pix2foc(pix, origin) - pix0

        # Update initial solution:
        pix -= dpix

        # Norm (L2) squared of the correction:
        dn = np.sum(dpix*dpix, axis=1)
        dnprev = dn.copy()  # if adaptive else dn
        tol2 = tolerance**2

        # Prepare for iterative process
        k = 1
        ind = None
        inddiv = None

        # Turn off numpy runtime warnings for 'invalid' and 'over':
        old_invalid = np.geterr()['invalid']
        old_over = np.geterr()['over']
        np.seterr(invalid='ignore', over='ignore')

        #############################################################
        ##                NON-ADAPTIVE ITERATIONS:                 ##
        #############################################################
        if not adaptive:
            # Fixed-point iterations:
            while (np.nanmax(dn) >= tol2 and k < maxiter):
                # Find correction to the previous solution:
                dpix = self.pix2foc(pix, origin) - pix0

                # Compute norm (L2) squared of the correction:
                dn = np.sum(dpix*dpix, axis=1)

                # Check for divergence (we do this in two stages
                # to optimize performance for the most common
                # scenario when successive approximations converge):
                if detect_divergence:
                    divergent = (dn >= dnprev)
                    if np.any(divergent):
                        # Find solutions that have not yet converged:
                        slowconv = (dn >= tol2)
                        inddiv, = np.where(divergent & slowconv)

                        if inddiv.shape[0] > 0:
                            # Update indices of elements that
                            # still need correction:
                            conv = (dn < dnprev)
                            iconv = np.where(conv)

                            # Apply correction:
                            dpixgood = dpix[iconv]
                            pix[iconv] -= dpixgood
                            dpix[iconv] = dpixgood

                            # For the next iteration choose
                            # non-divergent points that have not yet
                            # converged to the requested accuracy:
                            ind, = np.where(slowconv & conv)
                            pix0 = pix0[ind]
                            dnprev[ind] = dn[ind]
                            k += 1

                            # Switch to adaptive iterations:
                            adaptive = True
                            break
                    # Save current correction magnitudes for later:
                    dnprev = dn

                # Apply correction:
                pix -= dpix
                k += 1

        #############################################################
        ##                  ADAPTIVE ITERATIONS:                   ##
        #############################################################
        if adaptive:
            if ind is None:
                ind, = np.where(np.isfinite(pix).all(axis=1))
                pix0 = pix0[ind]

            # "Adaptive" fixed-point iterations:
            while (ind.shape[0] > 0 and k < maxiter):
                # Find correction to the previous solution:
                dpixnew = self.pix2foc(pix[ind], origin) - pix0

                # Compute norm (L2) of the correction:
                dnnew = np.sum(np.square(dpixnew), axis=1)

                # Bookeeping of corrections:
                dnprev[ind] = dn[ind].copy()
                dn[ind] = dnnew

                if detect_divergence:
                    # Find indices of pixels that are converging:
                    conv = (dnnew < dnprev[ind])
                    iconv = np.where(conv)
                    iiconv = ind[iconv]

                    # Apply correction:
                    dpixgood = dpixnew[iconv]
                    pix[iiconv] -= dpixgood
                    dpix[iiconv] = dpixgood

                    # Find indices of solutions that have not yet
                    # converged to the requested accuracy
                    # AND that do not diverge:
                    subind, = np.where((dnnew >= tol2) & conv)

                else:
                    # Apply correction:
                    pix[ind] -= dpixnew
                    dpix[ind] = dpixnew

                    # Find indices of solutions that have not yet
                    # converged to the requested accuracy:
                    subind, = np.where(dnnew >= tol2)

                # Choose solutions that need more iterations:
                ind = ind[subind]
                pix0 = pix0[subind]

                k += 1

        #############################################################
        ##         FINAL DETECTION OF INVALID, DIVERGING,          ##
        ##         AND FAILED-TO-CONVERGE POINTS                   ##
        #############################################################
        # Identify diverging and/or invalid points:
        invalid = ((~np.all(np.isfinite(pix), axis=1)) &
                   (np.all(np.isfinite(world), axis=1)))

        # When detect_divergence==False, dnprev is outdated
        # (it is the norm of the very first correction).
        # Still better than nothing...
        inddiv, = np.where(((dn >= tol2) & (dn >= dnprev)) | invalid)
        if inddiv.shape[0] == 0:
            inddiv = None

        # Identify points that did not converge within 'maxiter'
        # iterations:
        if k >= maxiter:
            ind, = np.where((dn >= tol2) & (dn < dnprev) & (~invalid))
            if ind.shape[0] == 0:
                ind = None
        else:
            ind = None

        # Restore previous numpy error settings:
        np.seterr(invalid=old_invalid, over=old_over)

        #############################################################
        ##  RAISE EXCEPTION IF DIVERGING OR TOO SLOWLY CONVERGING  ##
        ##  DATA POINTS HAVE BEEN DETECTED:                        ##
        #############################################################
        if (ind is not None or inddiv is not None) and not quiet:
            if inddiv is None:
                raise NoConvergence(
                    "'WCS.all_world2pix' failed to "
                    "converge to the requested accuracy after {:d} "
                    "iterations.".format(k), best_solution=pix,
                    accuracy=np.abs(dpix), niter=k,
                    slow_conv=ind, divergent=None)
            else:
                raise NoConvergence(
                    "'WCS.all_world2pix' failed to "
                    "converge to the requested accuracy.\n"
                    "After {0:d} iterations, the solution is diverging "
                    "at least for one input point."
                    .format(k), best_solution=pix,
                    accuracy=np.abs(dpix), niter=k,
                    slow_conv=ind, divergent=inddiv)

        return pix

    def all_world2pix(self, *args, **kwargs):
        if self.wcs is None:
            raise ValueError("No basic WCS settings were created.")

        tolerance  = kwargs.pop('tolerance', 1e-4)
        maxiter    = kwargs.pop('maxiter', 20)
        adaptive   = kwargs.pop('adaptive', False)
        detect_div = kwargs.pop('detect_divergence', True)
        quiet      = kwargs.pop('quiet', False)

        return self._array_converter(
            lambda *args, **kwargs:
            self._all_world2pix(
                *args, tolerance=tolerance, maxiter=maxiter,
                adaptive=adaptive, detect_divergence=detect_div,
                quiet=quiet),
            'input', *args, **kwargs
        )

    all_world2pix.__doc__ = """
        all_world2pix(*arg, accuracy=1.0e-4, maxiter=20,
        adaptive=False, detect_divergence=True, quiet=False)

        Transforms world coordinates to pixel coordinates, using
        numerical iteration to invert the full forward transformation
        `~astropy.wcs.WCS.all_pix2world` with complete
        distortion model.


        Parameters
        ----------
        {0}

            For a transformation that is not two-dimensional, the
            two-argument form must be used.

        {1}

        tolerance : float, optional (Default = 1.0e-4)
            Tolerance of solution. Iteration terminates when the
            iterative solver estimates that the "true solution" is
            within this many pixels current estimate, more
            specifically, when the correction to the solution found
            during the previous iteration is smaller
            (in the sense of the L2 norm) than ``tolerance``.

        maxiter : int, optional (Default = 20)
            Maximum number of iterations allowed to reach a solution.

        quiet : bool, optional (Default = False)
            Do not throw :py:class:``NoConvergence`` exceptions when
            the method does not converge to a solution with the
            required accuracy within a specified number of maximum
            iterations set by ``maxiter`` parameter. Instead,
            simply return the found solution.

        Other Parameters
        ----------------
        adaptive : bool, optional (Default = False)
            Specifies whether to adaptively select only points that
            did not converge to a solution within the required
            accuracy for the next iteration. Default is recommended
            for HST as well as most other instruments.

            .. note::
               The :py:meth:`all_world2pix` uses a vectorized
               implementation of the method of consecutive
               approximations (see ``Notes`` section below) in which it
               iterates over *all* input points *regardless* until
               the required accuracy has been reached for *all* input
               points. In some cases it may be possible that
               *almost all* points have reached the required accuracy
               but there are only a few of input data points for
               which additional iterations may be needed (this
               depends mostly on the characteristics of the geometric
               distortions for a given instrument). In this situation
               it may be advantageous to set ``adaptive`` = `True` in
               which case :py:meth:`all_world2pix` will continue
               iterating *only* over the points that have not yet
               converged to the required accuracy. However, for the
               HST's ACS/WFC detector, which has the strongest
               distortions of all HST instruments, testing has
               shown that enabling this option would lead to a about
               50-100\% penalty in computational time (depending on
               specifics of the image, geometric distortions, and
               number of input points to be converted). Therefore,
               for HST and possibly instruments, it is recommended
               to set ``adaptive`` = `False`. The only danger in
               getting this setting wrong will be a performance
               penalty.

            .. note::
               When ``detect_divergence`` is `True`,
               :py:meth:`all_world2pix` will automatically switch
               to the adaptive algorithm once divergence has been
               detected.

        detect_divergence : bool, optional (Default = True)
            Specifies whether to perform a more detailed analysis
            of the convergence to a solution. Normally
            :py:meth:`all_world2pix` may not achieve the required
            accuracy if either the ``tolerance`` or ``maxiter`` arguments
            are too low. However, it may happen that for some
            geometric distortions the conditions of convergence for
            the the method of consecutive approximations used by
            :py:meth:`all_world2pix` may not be satisfied, in which
            case consecutive approximations to the solution will
            diverge regardless of the ``tolerance`` or ``maxiter``
            settings.

            When ``detect_divergence`` is `False`, these divergent
            points will be detected as not having achieved the
            required accuracy (without further details). In addition,
            if ``adaptive`` is `False` then the algorithm will not
            know that the solution (for specific points) is diverging
            and will continue iterating and trying to "improve"
            diverging solutions. This may result in ``NaN`` or
            ``Inf`` values in the return results (in addition to a
            performance penalties). Even when ``detect_divergence``
            is `False`, :py:meth:`all_world2pix`, at the end of the
            iterative process, will identify invalid results
            (``NaN`` or ``Inf``) as "diverging" solutions and will
            raise :py:class:``NoConvergence`` unless the ``quiet``
            parameter is set to `True`.

            When ``detect_divergence`` is `True`,
            :py:meth:`all_world2pix` will detect points for which
            current correction to the coordinates is larger than
            the correction applied during the previous iteration
            **if** the requested accuracy **has not yet been
            achieved**. In this case, if ``adaptive`` is `True`,
            these points will be excluded from further iterations and
            if ``adaptive`` is `False`, :py:meth:`all_world2pix` will
            automatically switch to the adaptive algorithm. Thus, the
            reported divergent solution will be the latest converging
            solution computed immediately *before* divergence
            has been detected.

            .. note::
               When accuracy has been achieved, small increases in
               current corrections may be possible due to rounding
               errors (when ``adaptive`` is `False`) and such
               increases will be ignored.

            .. note::
               Based on our testing using HST ACS/WFC images, setting
               ``detect_divergence`` to `True` will incur about 5-20\%
               performance penalty with the larger penalty
               corresponding to ``adaptive`` set to `True`.
               Because the benefits of enabling this
               feature outweigh the small performance penalty,
               especially when ``adaptive`` = `False`, it is
               recommended to set ``detect_divergence`` to `True`,
               unless extensive testing of the distortion models for
               images from specific instruments show a good stability
               of the numerical method for a wide range of
               coordinates (even outside the image itself).

            .. note::
               Indices of the diverging inverse solutions will be
               reported in the ``divergent`` attribute of the
               raised :py:class:``NoConvergence`` exception object.

        Returns
        -------

        {2}

        Notes
        -----
        The order of the axes for the input world array is determined by
        the ``CTYPEia`` keywords in the FITS header, therefore it may
        not always be of the form (*ra*, *dec*).  The
        `~astropy.wcs.Wcsprm.lat`, `~astropy.wcs.Wcsprm.lng`,
        `~astropy.wcs.Wcsprm.lattyp`, and
        `~astropy.wcs.Wcsprm.lngtyp`
        members can be used to determine the order of the axes.

        Using the method of fixed-point iterations approximations we
        iterate starting with the initial approximation, which is
        computed using the non-distortion-aware
        :py:meth:`wcs_world2pix` (or equivalent).

        The :py:meth:`all_world2pix` function uses a vectorized
        implementation of the method of consecutive approximations and
        therefore it is highly efficient (>30x) when *all* data points
        that need to be converted from sky coordinates to image
        coordinates are passed at *once*. Therefore, it is advisable,
        whenever possible, to pass as input a long array of all points
        that need to be converted to :py:meth:`all_world2pix` instead
        of calling :py:meth:`all_world2pix` for each data point. Also
        see the note to the ``adaptive`` parameter.

        Raises
        ------
        NoConvergence
            The method did not converge to a
            solution to the required accuracy within a specified
            number of maximum iterations set by the ``maxiter``
            parameter. To turn off this exception, set ``quiet`` to
            `True`. Indices of the points for which the requested
            accuracy was not achieved (if any) will be listed in the
            ``slow_conv`` attribute of the
            raised :py:class:``NoConvergence`` exception object.

            See :py:class:``NoConvergence`` documentation for
            more details.

        MemoryError
            Memory allocation failed.

        SingularMatrixError
            Linear transformation matrix is singular.

        InconsistentAxisTypesError
            Inconsistent or unrecognized coordinate axis types.

        ValueError
            Invalid parameter value.

        ValueError
            Invalid coordinate transformation parameters.

        ValueError
            x- and y-coordinate arrays are not the same size.

        InvalidTransformError
            Invalid coordinate transformation parameters.

        InvalidTransformError
            Ill-conditioned coordinate transformation parameters.

        Examples
        --------
        >>> import astropy.io.fits as fits
        >>> import astropy.wcs as wcs
        >>> import numpy as np
        >>> import os

        >>> filename = os.path.join(wcs.__path__[0], 'tests/data/j94f05bgq_flt.fits')
        >>> hdulist = fits.open(filename)
        >>> w = wcs.WCS(hdulist[('sci',1)].header, hdulist)
        >>> hdulist.close()

        >>> ra, dec = w.all_pix2world([1,2,3], [1,1,1], 1)
        >>> print(ra)
        [ 5.52645627  5.52649663  5.52653698]
        >>> print(dec)
        [-72.05171757 -72.05171276 -72.05170795]
        >>> radec = w.all_pix2world([[1,1], [2,1], [3,1]], 1)
        >>> print(radec)
        [[  5.52645627 -72.05171757]
         [  5.52649663 -72.05171276]
         [  5.52653698 -72.05170795]]
        >>> x, y = w.all_world2pix(ra, dec, 1)
        >>> print(x)
        [ 1.00000238  2.00000237  3.00000236]
        >>> print(y)
        [ 0.99999996  0.99999997  0.99999997]
        >>> xy = w.all_world2pix(radec, 1)
        >>> print(xy)
        [[ 1.00000238  0.99999996]
         [ 2.00000237  0.99999997]
         [ 3.00000236  0.99999997]]
        >>> xy = w.all_world2pix(radec, 1, maxiter=3,
        ...                      tolerance=1.0e-10, quiet=False)
        Traceback (most recent call last):
        ...
        NoConvergence: 'WCS.all_world2pix' failed to converge to the
        requested accuracy. After 3 iterations, the solution is
        diverging at least for one input point.

        >>> # Now try to use some diverging data:
        >>> divradec = w.all_pix2world([[1.0, 1.0],
        ...                             [10000.0, 50000.0],
        ...                             [3.0, 1.0]], 1)
        >>> print(divradec)
        [[  5.52645627 -72.05171757]
         [  7.15976932 -70.8140779 ]
         [  5.52653698 -72.05170795]]

        >>> # First, turn detect_divergence on:
        >>> try:
        ...   xy = w.all_world2pix(divradec, 1, maxiter=20,
        ...                        tolerance=1.0e-4, adaptive=False,
        ...                        detect_divergence=True,
        ...                        quiet=False)
        ... except wcs.wcs.NoConvergence as e:
        ...   print("Indices of diverging points: {{0}}"
        ...         .format(e.divergent))
        ...   print("Indices of poorly converging points: {{0}}"
        ...         .format(e.slow_conv))
        ...   print("Best solution:\\n{{0}}".format(e.best_solution))
        ...   print("Achieved accuracy:\\n{{0}}".format(e.accuracy))
        Indices of diverging points: [1]
        Indices of poorly converging points: None
        Best solution:
        [[  1.00000238e+00   9.99999965e-01]
         [ -1.99441636e+06   1.44309097e+06]
         [  3.00000236e+00   9.99999966e-01]]
        Achieved accuracy:
        [[  6.13968380e-05   8.59638593e-07]
         [  8.59526812e+11   6.61713548e+11]
         [  6.09398446e-05   8.38759724e-07]]
        >>> raise e
        Traceback (most recent call last):
        ...
        NoConvergence: 'WCS.all_world2pix' failed to converge to the
        requested accuracy.  After 5 iterations, the solution is
        diverging at least for one input point.

        >>> # This time turn detect_divergence off:
        >>> try:
        ...   xy = w.all_world2pix(divradec, 1, maxiter=20,
        ...                        tolerance=1.0e-4, adaptive=False,
        ...                        detect_divergence=False,
        ...                        quiet=False)
        ... except wcs.wcs.NoConvergence as e:
        ...   print("Indices of diverging points: {{0}}"
        ...         .format(e.divergent))
        ...   print("Indices of poorly converging points: {{0}}"
        ...         .format(e.slow_conv))
        ...   print("Best solution:\\n{{0}}".format(e.best_solution))
        ...   print("Achieved accuracy:\\n{{0}}".format(e.accuracy))
        Indices of diverging points: [1]
        Indices of poorly converging points: None
        Best solution:
        [[ 1.00000009  1.        ]
         [        nan         nan]
         [ 3.00000009  1.        ]]
        Achieved accuracy:
        [[  2.29417358e-06   3.21222995e-08]
         [             nan              nan]
         [  2.27407877e-06   3.13005639e-08]]
        >>> raise e
        Traceback (most recent call last):
        ...
        NoConvergence: 'WCS.all_world2pix' failed to converge to the
        requested accuracy.  After 6 iterations, the solution is
        diverging at least for one input point.

        """.format(__.TWO_OR_MORE_ARGS('naxis', 8),
                   __.RA_DEC_ORDER(8),
                   __.RETURNS('pixel coordinates', 8))


    def wcs_world2pix(self, *args, **kwargs):
        if self.wcs is None:
            raise ValueError("No basic WCS settings were created.")
        return self._array_converter(
            lambda xy, o: self.wcs.s2p(xy, o)['pixcrd'],
            'input', *args, **kwargs)
    wcs_world2pix.__doc__ = """
        Transforms world coordinates to pixel coordinates, using only
        the basic `wcslib`_ WCS transformation.  No `SIP`_ or
        `distortion paper`_ table lookup transformation is applied.

        Parameters
        ----------
        {0}

            For a transformation that is not two-dimensional, the
            two-argument form must be used.

        {1}

        Returns
        -------

        {2}

        Notes
        -----
        The order of the axes for the input world array is determined by
        the ``CTYPEia`` keywords in the FITS header, therefore it may
        not always be of the form (*ra*, *dec*).  The
        `~astropy.wcs.Wcsprm.lat`, `~astropy.wcs.Wcsprm.lng`,
        `~astropy.wcs.Wcsprm.lattyp` and `~astropy.wcs.Wcsprm.lngtyp`
        members can be used to determine the order of the axes.

        Raises
        ------
        MemoryError
            Memory allocation failed.

        SingularMatrixError
            Linear transformation matrix is singular.

        InconsistentAxisTypesError
            Inconsistent or unrecognized coordinate axis types.

        ValueError
            Invalid parameter value.

        ValueError
            Invalid coordinate transformation parameters.

        ValueError
            x- and y-coordinate arrays are not the same size.

        InvalidTransformError
            Invalid coordinate transformation parameters.

        InvalidTransformError
            Ill-conditioned coordinate transformation parameters.
        """.format(__.TWO_OR_MORE_ARGS('naxis', 8),
                   __.RA_DEC_ORDER(8),
                   __.RETURNS('pixel coordinates', 8))

    def pix2foc(self, *args):
        return self._array_converter(self._pix2foc, None, *args)
    pix2foc.__doc__ = """
        Convert pixel coordinates to focal plane coordinates using the
        `SIP`_ polynomial distortion convention and `distortion
        paper`_ table-lookup correction.

        The output is in absolute pixel coordinates, not relative to
        ``CRPIX``.

        Parameters
        ----------

        {0}

        Returns
        -------

        {1}

        Raises
        ------
        MemoryError
            Memory allocation failed.

        ValueError
            Invalid coordinate transformation parameters.
        """.format(__.TWO_OR_MORE_ARGS('2', 8),
                   __.RETURNS('focal coordinates', 8))

    def p4_pix2foc(self, *args):
        return self._array_converter(self._p4_pix2foc, None, *args)
    p4_pix2foc.__doc__ = """
        Convert pixel coordinates to focal plane coordinates using
        `distortion paper`_ table-lookup correction.

        The output is in absolute pixel coordinates, not relative to
        ``CRPIX``.

        Parameters
        ----------

        {0}

        Returns
        -------

        {1}

        Raises
        ------
        MemoryError
            Memory allocation failed.

        ValueError
            Invalid coordinate transformation parameters.
        """.format(__.TWO_OR_MORE_ARGS('2', 8),
                   __.RETURNS('focal coordinates', 8))

    def det2im(self, *args):
        return self._array_converter(self._det2im, None, *args)
    det2im.__doc__ = """
        Convert detector coordinates to image plane coordinates using
        `distortion paper`_ table-lookup correction.

        The output is in absolute pixel coordinates, not relative to
        ``CRPIX``.

        Parameters
        ----------

        {0}

        Returns
        -------

        {1}

        Raises
        ------
        MemoryError
            Memory allocation failed.

        ValueError
            Invalid coordinate transformation parameters.
        """.format(__.TWO_OR_MORE_ARGS('2', 8),
                   __.RETURNS('pixel coordinates', 8))

    def sip_pix2foc(self, *args):
        if self.sip is None:
            if len(args) == 2:
                return args[0]
            elif len(args) == 3:
                return args[:2]
            else:
                raise TypeError("Wrong number of arguments")
        return self._array_converter(self.sip.pix2foc, None, *args)
    sip_pix2foc.__doc__ = """
        Convert pixel coordinates to focal plane coordinates using the
        `SIP`_ polynomial distortion convention.

        The output is in pixel coordinates, relative to ``CRPIX``.

        FITS WCS `distortion paper`_ table lookup correction is not
        applied, even if that information existed in the FITS file
        that initialized this :class:`~astropy.wcs.WCS` object.  To
        correct for that, use `~astropy.wcs.WCS.pix2foc` or
        `~astropy.wcs.WCS.p4_pix2foc`.

        Parameters
        ----------

        {0}

        Returns
        -------

        {1}

        Raises
        ------
        MemoryError
            Memory allocation failed.

        ValueError
            Invalid coordinate transformation parameters.
        """.format(__.TWO_OR_MORE_ARGS('2', 8),
                   __.RETURNS('focal coordinates', 8))

    def sip_foc2pix(self, *args):
        if self.sip is None:
            if len(args) == 2:
                return args[0]
            elif len(args) == 3:
                return args[:2]
            else:
                raise TypeError("Wrong number of arguments")
        return self._array_converter(self.sip.foc2pix, None, *args)
    sip_foc2pix.__doc__ = """
        Convert focal plane coordinates to pixel coordinates using the
        `SIP`_ polynomial distortion convention.

        FITS WCS `distortion paper`_ table lookup distortion
        correction is not applied, even if that information existed in
        the FITS file that initialized this `~astropy.wcs.WCS` object.

        Parameters
        ----------

        {0}

        Returns
        -------

        {1}

        Raises
        ------
        MemoryError
            Memory allocation failed.

        ValueError
            Invalid coordinate transformation parameters.
        """.format(__.TWO_OR_MORE_ARGS('2', 8),
                   __.RETURNS('pixel coordinates', 8))

    def to_fits(self, relax=False, key=None):
        """
        Generate an `astropy.io.fits.HDUList` object with all of the
        information stored in this object.  This should be logically identical
        to the input FITS file, but it will be normalized in a number of ways.

        See `to_header` for some warnings about the output produced.

        Parameters
        ----------

        relax : bool or int, optional
            Degree of permissiveness:

            - `False` (default): Write all extensions that are
              considered to be safe and recommended.

            - `True`: Write all recognized informal extensions of the
              WCS standard.

            - `int`: a bit field selecting specific extensions to
              write.  See :ref:`relaxwrite` for details.

        key : str
            The name of a particular WCS transform to use.  This may be
            either ``' '`` or ``'A'``-``'Z'`` and corresponds to the ``"a"``
            part of the ``CTYPEia`` cards.

        Returns
        -------
        hdulist : `astropy.io.fits.HDUList`
        """

        header = self.to_header(relax=relax, key=key)

        hdu = fits.PrimaryHDU(header=header)
        hdulist = fits.HDUList(hdu)

        self._write_det2im(hdulist)
        self._write_distortion_kw(hdulist)

        return hdulist

    def to_header(self, relax=None, key=None):
        """Generate an `astropy.io.fits.Header` object with the basic WCS
        and SIP information stored in this object.  This should be
        logically identical to the input FITS file, but it will be
        normalized in a number of ways.

        .. warning::

          This function does not write out FITS WCS `distortion
          paper`_ information, since that requires multiple FITS
          header data units.  To get a full representation of
          everything in this object, use `to_fits`.

        Parameters
        ----------
        relax : bool or int, optional
            Degree of permissiveness:

            - `False` (default): Write all extensions that are
              considered to be safe and recommended.

            - `True`: Write all recognized informal extensions of the
              WCS standard.

            - `int`: a bit field selecting specific extensions to
              write.  See :ref:`relaxwrite` for details.

            If the ``relax`` keyword argument is not given and any
            keywords were omitted from the output, an
            `~astropy.utils.exceptions.AstropyWarning` is displayed.
            To override this, explicitly pass a value to ``relax``.

        key : str
            The name of a particular WCS transform to use.  This may be
            either ``' '`` or ``'A'``-``'Z'`` and corresponds to the ``"a"``
            part of the ``CTYPEia`` cards.

        Returns
        -------
        header : `astropy.io.fits.Header`

        Notes
        -----
        The output header will almost certainly differ from the input in a
        number of respects:

          1. The output header only contains WCS-related keywords.  In
             particular, it does not contain syntactically-required
             keywords such as ``SIMPLE``, ``NAXIS``, ``BITPIX``, or
             ``END``.

          2. Deprecated (e.g. ``CROTAn``) or non-standard usage will
             be translated to standard (this is partially dependent on
             whether ``fix`` was applied).

          3. Quantities will be converted to the units used internally,
             basically SI with the addition of degrees.

          4. Floating-point quantities may be given to a different decimal
             precision.

          5. Elements of the ``PCi_j`` matrix will be written if and
             only if they differ from the unit matrix.  Thus, if the
             matrix is unity then no elements will be written.

          6. Additional keywords such as ``WCSAXES``, ``CUNITia``,
             ``LONPOLEa`` and ``LATPOLEa`` may appear.

          7. The original keycomments will be lost, although
             `to_header` tries hard to write meaningful comments.

          8. Keyword order may be changed.

        """
        display_warning = False
        if relax is None:
            display_warning = True
            relax = False

        if key is not None:
            self.wcs.alt = key

        if relax not in (True, False):
            do_sip = relax & WCSHDO_SIP
            relax &= ~WCSHDO_SIP
        else:
            do_sip = relax

        if self.wcs is not None:
            header_string = self.wcs.to_header(relax)
            header = fits.Header.fromstring(header_string)
        else:
            header = fits.Header()

        if do_sip and self.sip is not None:
            for key, val in self._write_sip_kw().items():
                header[key] = val

        if display_warning:
            full_header = self.to_header(relax=True, key=key)
            missing_keys = []
            for key, val in full_header.items():
                if key not in header:
                    missing_keys.append(key)

            if len(missing_keys):
                warnings.warn(
                    "Some non-standard WCS keywords were excluded: {0} "
                    "Use the ``relax`` kwarg to control this.".format(
                        ', '.join(missing_keys)),
                    AstropyWarning)

        return header

    def to_header_string(self, relax=None):
        """
        Identical to `to_header`, but returns a string containing the
        header cards.
        """
        return str(self.to_header(relax))

    def footprint_to_file(self, filename=None, color='green', width=2):
        """
        Writes out a `ds9`_ style regions file. It can be loaded
        directly by `ds9`_.

        Parameters
        ----------
        filename : str, optional
            Output file name - default is ``'footprint.reg'``

        color : str, optional
            Color to use when plotting the line.

        width : int, optional
            Width of the region line.
        """
        if not filename:
            filename = 'footprint.reg'
        comments = '# Region file format: DS9 version 4.0 \n'
        comments += ('# global color=green font="helvetica 12 bold ' +
                     'select=1 highlite=1 edit=1 move=1 delete=1 ' +
                     'include=1 fixed=0 source\n')

        f = open(filename, 'a')
        f.write(comments)
        f.write('linear\n')
        f.write('polygon(')
        self.calc_footprint().tofile(f, sep=',')
        f.write(') # color={0}, width={1:d} \n'.format(color, width))
        f.close()

    def _get_naxis(self, header=None):
        self._naxis1 = 0
        self._naxis2 = 0
        if (header is not None and
            not isinstance(header, (six.text_type, six.binary_type))):
            self._naxis1 = header.get('NAXIS1', 0)
            self._naxis2 = header.get('NAXIS2', 0)

    def rotateCD(self, theta):
        _theta = np.deg2rad(theta)
        _mrot = np.zeros(shape=(2, 2), dtype=np.double)
        _mrot[0] = (np.cos(_theta), np.sin(_theta))
        _mrot[1] = (-np.sin(_theta), np.cos(_theta))
        new_cd = np.dot(self.wcs.cd, _mrot)
        self.wcs.cd = new_cd

    def printwcs(self):
        print(repr(self))

    def __repr__(self):
        '''
        Return a short description. Simply porting the behavior from
        the `printwcs()` method.
        '''
        description = ["WCS Keywords\n",
                       "Number of WCS axes: {0!r}".format(self.naxis)]
        sfmt = ' : ' +  "".join(["{"+"{0}".format(i)+"!r}  " for i in range(self.naxis)])

        keywords = ['CTYPE', 'CRVAL', 'CRPIX']
        values = [self.wcs.ctype, self.wcs.crval, self.wcs.crpix]
        for keyword, value in zip(keywords, values):
            description.append(keyword+sfmt.format(*value))

        if hasattr(self.wcs, 'pc'):
            for i in range(self.naxis):
                s = ''
                for j in range(self.naxis):
                    s += ''.join(['PC', str(i+1), '_', str(j+1), ' '])
                s += sfmt
                description.append(s.format(*self.wcs.pc[i]))
            s = 'CDELT' + sfmt
            description.append(s.format(*self.wcs.cdelt))
        elif hasattr(self.wcs, 'cd'):
            for i in range(self.naxis):
                s = ''
                for j in range(self.naxis):
                    s += "".join(['CD', str(i+1), '_', str(j+1), ' '])
                s += sfmt
                description.append(s.format(*self.wcs.cd[i]))

        description.append('NAXIS    : {0!r} {1!r}'.format(self._naxis1,
                           self._naxis2))
        return '\n'.join(description)

    def get_axis_types(self):
        """
        Similar to `self.wcsprm.axis_types <astropy.wcs.Wcsprm.axis_types>`
        but provides the information in a more Python-friendly format.

        Returns
        -------
        result : list of dicts

            Returns a list of dictionaries, one for each axis, each
            containing attributes about the type of that axis.

            Each dictionary has the following keys:

            - 'coordinate_type':

              - None: Non-specific coordinate type.

              - 'stokes': Stokes coordinate.

              - 'celestial': Celestial coordinate (including ``CUBEFACE``).

              - 'spectral': Spectral coordinate.

            - 'scale':

              - 'linear': Linear axis.

              - 'quantized': Quantized axis (``STOKES``, ``CUBEFACE``).

              - 'non-linear celestial': Non-linear celestial axis.

              - 'non-linear spectral': Non-linear spectral axis.

              - 'logarithmic': Logarithmic axis.

              - 'tabular': Tabular axis.

            - 'group'

              - Group number, e.g. lookup table number

            - 'number'

              - For celestial axes:

                - 0: Longitude coordinate.

                - 1: Latitude coordinate.

                - 2: ``CUBEFACE`` number.

              - For lookup tables:

                - the axis number in a multidimensional table.

            ``CTYPEia`` in ``"4-3"`` form with unrecognized algorithm code will
            generate an error.
        """
        if self.wcs is None:
            raise AttributeError(
                "This WCS object does not have a wcsprm object.")

        coordinate_type_map = {
            0: None,
            1: 'stokes',
            2: 'celestial',
            3: 'spectral'}

        scale_map = {
            0: 'linear',
            1: 'quantized',
            2: 'non-linear celestial',
            3: 'non-linear spectral',
            4: 'logarithmic',
            5: 'tabular'}

        result = []
        for axis_type in self.wcs.axis_types:
            subresult = {}

            coordinate_type = (axis_type // 1000) % 10
            subresult['coordinate_type'] = coordinate_type_map[coordinate_type]

            scale = (axis_type // 100) % 10
            subresult['scale'] = scale_map[scale]

            group = (axis_type // 10) % 10
            subresult['group'] = group

            number = axis_type % 10
            subresult['number'] = number

            result.append(subresult)

        return result

    def __reduce__(self):
        """
        Support pickling of WCS objects.  This is done by serializing
        to an in-memory FITS file and dumping that as a string.
        """

        hdulist = self.to_fits(relax=True)

        buffer = io.BytesIO()
        hdulist.writeto(buffer)

        return (__WCS_unpickle__,
                (self.__class__, self.__dict__, buffer.getvalue(),))

    def dropaxis(self, dropax):
        """
        Remove an axis from the WCS.

        Parameters
        ----------
        wcs : `~astropy.wcs.WCS`
            The WCS with naxis to be chopped to naxis-1
        dropax : int
            The index of the WCS to drop, counting from 0 (i.e., python convention,
            not FITS convention)

        Returns
        -------
        A new `~astropy.wcs.WCS` instance with one axis fewer
        """
        inds = list(range(self.wcs.naxis))
        inds.pop(dropax)

        # axis 0 has special meaning to sub
        # if wcs.wcs.ctype == ['RA','DEC','VLSR'], you want
        # wcs.sub([1,2]) to get 'RA','DEC' back
        return self.sub([i+1 for i in inds])

    def swapaxes(self, ax0, ax1):
        """
        Swap axes in a WCS.

        Parameters
        ----------
        wcs : `~astropy.wcs.WCS`
            The WCS to have its axes swapped
        ax0 : int
        ax1 : int
            The indices of the WCS to be swapped, counting from 0 (i.e., python
            convention, not FITS convention)

        Returns
        -------
        A new `~astropy.wcs.WCS` instance with the same number of axes, but two
        swapped
        """
        inds = list(range(self.wcs.naxis))
        inds[ax0],inds[ax1] = inds[ax1],inds[ax0]

        return self.sub([i+1 for i in inds])

    def reorient_celestial_first(self):
        """
        Reorient the WCS such that the celestial axes are first, followed by
        the spectral axis, followed by any others.
        Assumes at least celestial axes are present.
        """
        return self.sub([WCSSUB_CELESTIAL, WCSSUB_SPECTRAL, WCSSUB_STOKES])

    def slice(self, view, numpy_order=True):
        """
        Slice a WCS instance using a Numpy slice. The order of the slice should
        be reversed (as for the data) compared to the natural WCS order.

        Parameters
        ----------
        view : tuple
            A tuple containing the same number of slices as the WCS system.
            The ``step`` method, the third argument to a slice, is not
            presently supported.
        numpy_order : bool
            Use numpy order, i.e. slice the WCS so that an identical slice
            applied to a numpy array will slice the array and WCS in the same
            way. If set to `False`, the WCS will be sliced in FITS order,
            meaning the first slice will be applied to the *last* numpy index
            but the *first* WCS axis.

        Returns
        -------
        wcs_new : `~astropy.wcs.WCS`
            A new resampled WCS axis
        """
        if hasattr(view, '__len__') and len(view) > self.wcs.naxis:
            raise ValueError("Must have # of slices <= # of WCS axes")
        elif not hasattr(view, '__len__'): # view MUST be an iterable
            view = [view]

        if not all([isinstance(x, slice) for x in view]):
            raise ValueError("Cannot downsample a WCS with indexing.  Use "
                             "wcs.sub or wcs.dropaxis if you want to remove "
                             "axes.")

        wcs_new = self.deepcopy()
        for i, iview in enumerate(view):
            if iview.step is not None and iview.start is None:
                # Slice from "None" is equivalent to slice from 0 (but one
                # might want to downsample, so allow slices with
                # None,None,step or None,stop,step)
                iview = slice(0, iview.stop, iview.step)

            if iview.start is not None:
                if numpy_order:
                    wcs_index = self.wcs.naxis - 1 - i
                else:
                    wcs_index = i

                if iview.step not in (None, 1):
                    crpix = self.wcs.crpix[wcs_index]
                    cdelt = self.wcs.cdelt[wcs_index]
                    # equivalently (keep this comment so you can compare eqns):
                    # wcs_new.wcs.crpix[wcs_index] =
                    # (crpix - iview.start)*iview.step + 0.5 - iview.step/2.
                    crp = ((crpix - iview.start - 1.)/iview.step
                           + 0.5 + 1./iview.step/2.)
                    wcs_new.wcs.crpix[wcs_index] = crp
                    wcs_new.wcs.cdelt[wcs_index] = cdelt * iview.step
                else:
                    wcs_new.wcs.crpix[wcs_index] -= iview.start
        return wcs_new

    def __getitem__(self, item):
        # "getitem" is a shortcut for self.slice; it is very limited
        # there is no obvious and unambiguous interpretation of wcs[1,2,3]
        # We COULD allow wcs[1] to link to wcs.sub([2])
        # (wcs[i] -> wcs.sub([i+1])
        return self.slice(item)

    def __iter__(self):
        # Having __getitem__ makes Python think WCS is iterable. However,
        # Python first checks whether __iter__ is present, so we can raise an
        # exception here.
        raise TypeError("'{0}' object is not iterable".format(self.__class__.__name__))

    @property
    def axis_type_names(self):
        """
        World names for each coordinate axis

        Returns
        -------
        A list of names along each axis
        """
        names = list(self.wcs.cname)
        types = self.wcs.ctype
        for i in range(len(names)):
            if len(names[i]) > 0:
                continue
            names[i] = types[i].split('-')[0]
        return names

    @property
    def celestial(self):
        """
        A copy of the current WCS with only the celestial axes included
        """
        return self.sub([WCSSUB_CELESTIAL])

    @property
    def is_celestial(self):
        return self.has_celestial and self.naxis==2

    @property
    def has_celestial(self):
        try:
            return self.celestial.naxis == 2
        except InconsistentAxisTypesError:
            return False

    @property
    def pixel_scale_matrix(self):

        try:
            cdelt = np.matrix(np.diag(self.wcs.get_cdelt()))
            pc = np.matrix(self.wcs.get_pc())
        except InconsistentAxisTypesError:
            try:
                # for non-celestial axes, get_cdelt doesn't work
                cdelt = np.matrix(self.wcs.cd) * np.matrix(np.diag(self.wcs.cdelt))
            except AttributeError:
                cdelt = np.matrix(np.diag(self.wcs.cdelt))

            try:
                pc = np.matrix(self.wcs.pc)
            except AttributeError:
                pc = 1

        pccd = np.array(cdelt * pc)

        return pccd

    def _as_mpl_axes(self):
        """
        Compatibility hook for Matplotlib and WCSAxes.

        This functionality requires the WCSAxes package to work. The reason
        we include this here is that it allows users to use WCSAxes without
        having to explicitly import WCSAxes, which means that if in future we
        merge WCSAxes into the Astropy core package, the API will remain the
        same. With this method, one can do:

            from astropy.wcs import WCS
            import matplotlib.pyplot as plt

            wcs = WCS('filename.fits')

            fig = plt.figure()
            ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)
            ...

        and this will generate a plot with the correct WCS coordinates on the
        axes. See http://wcsaxes.readthedocs.org for more information.
        """

        try:
            from wcsaxes import WCSAxes
        except ImportError:
            raise ImportError("Using WCS instances as Matplotlib projections "
                              "requires the WCSAxes package to be installed. "
                              "See http://wcsaxes.readthedocs.org for more "
                              "details.")
        else:
            return WCSAxes, {'wcs': self}


def __WCS_unpickle__(cls, dct, fits_data):
    """
    Unpickles a WCS object from a serialized FITS string.
    """

    self = cls.__new__(cls)
    self.__dict__.update(dct)

    buffer = io.BytesIO(fits_data)
    hdulist = fits.open(buffer)

    WCS.__init__(self, hdulist[0].header, hdulist)

    return self


def find_all_wcs(header, relax=True, keysel=None, fix=True,
                 translate_units='',
                 _do_set=True):
    """
    Find all the WCS transformations in the given header.

    Parameters
    ----------
    header : str or astropy.io.fits header object.

    relax : bool or int, optional
        Degree of permissiveness:

        - `True` (default): Admit all recognized informal extensions of the
          WCS standard.

        - `False`: Recognize only FITS keywords defined by the
          published WCS standard.

        - `int`: a bit field selecting specific extensions to accept.
          See :ref:`relaxread` for details.

    keysel : sequence of flags, optional
        A list of flags used to select the keyword types considered by
        wcslib.  When ``None``, only the standard image header
        keywords are considered (and the underlying wcspih() C
        function is called).  To use binary table image array or pixel
        list keywords, *keysel* must be set.

        Each element in the list should be one of the following strings:

            - 'image': Image header keywords

            - 'binary': Binary table image array keywords

            - 'pixel': Pixel list keywords

        Keywords such as ``EQUIna`` or ``RFRQna`` that are common to
        binary table image arrays and pixel lists (including
        ``WCSNna`` and ``TWCSna``) are selected by both 'binary' and
        'pixel'.

    fix : bool, optional
        When `True` (default), call `~astropy.wcs.Wcsprm.fix` on
        the resulting objects to fix any non-standard uses in the
        header.  `FITSFixedWarning` warnings will be emitted if any
        changes were made.

    translate_units : str, optional
        Specify which potentially unsafe translations of non-standard
        unit strings to perform.  By default, performs none.  See
        `WCS.fix` for more information about this parameter.  Only
        effective when ``fix`` is `True`.

    Returns
    -------
    wcses : list of `WCS` objects
    """

    if isinstance(header, (six.text_type, six.binary_type)):
        header_string = header
    elif isinstance(header, fits.Header):
        header_string = header.tostring()
    else:
        raise TypeError(
            "header must be a string or astropy.io.fits.Header object")

    keysel_flags = _parse_keysel(keysel)

    if isinstance(header_string, six.text_type):
        header_bytes = header_string.encode('ascii')
    else:
        header_bytes = header_string

    wcsprms = _wcs.find_all_wcs(header_bytes, relax, keysel_flags)

    result = []
    for wcsprm in wcsprms:
        subresult = WCS(fix=False, _do_set=False)
        subresult.wcs = wcsprm
        result.append(subresult)

        if fix:
            subresult.fix(translate_units)

        if _do_set:
            subresult.wcs.set()

    return result


def validate(source):
    """
    Prints a WCS validation report for the given FITS file.

    Parameters
    ----------
    source : str path, readable file-like object or `astropy.io.fits.HDUList` object
        The FITS file to validate.

    Returns
    -------
    results : WcsValidateResults instance
        The result is returned as nested lists.  The first level
        corresponds to the HDUs in the given file.  The next level has
        an entry for each WCS found in that header.  The special
        subclass of list will pretty-print the results as a table when
        printed.
    """
    class _WcsValidateWcsResult(list):
        def __init__(self, key):
            self._key = key

        def __repr__(self):
            result = ["  WCS key '{0}':".format(self._key or ' ')]
            if len(self):
                for entry in self:
                    for i, line in enumerate(entry.splitlines()):
                        if i == 0:
                            initial_indent = '    - '
                        else:
                            initial_indent = '      '
                        result.extend(
                            textwrap.wrap(
                                line,
                                initial_indent=initial_indent,
                                subsequent_indent='      '))
            else:
                result.append("    No issues.")
            return '\n'.join(result)

    class _WcsValidateHduResult(list):
        def __init__(self, hdu_index, hdu_name):
            self._hdu_index = hdu_index
            self._hdu_name = hdu_name
            list.__init__(self)

        def __repr__(self):
            if len(self):
                if self._hdu_name:
                    hdu_name = ' ({0})'.format(self._hdu_name)
                else:
                    hdu_name = ''
                result = ['HDU {0}{1}:'.format(self._hdu_index, hdu_name)]
                for wcs in self:
                    result.append(repr(wcs))
                return '\n'.join(result)
            return ''

    class _WcsValidateResults(list):
        def __repr__(self):
            result = []
            for hdu in self:
                content = repr(hdu)
                if len(content):
                    result.append(content)
            return '\n\n'.join(result)

    global __warningregistry__

    if isinstance(source, fits.HDUList):
        hdulist = source
    else:
        hdulist = fits.open(source)

    results = _WcsValidateResults()

    for i, hdu in enumerate(hdulist):
        hdu_results = _WcsValidateHduResult(i, hdu.name)
        results.append(hdu_results)

        with warnings.catch_warnings(record=True) as warning_lines:
            wcses = find_all_wcs(
                hdu.header, relax=_wcs.WCSHDR_reject,
                fix=False, _do_set=False)

        for wcs in wcses:
            wcs_results = _WcsValidateWcsResult(wcs.wcs.alt)
            hdu_results.append(wcs_results)

            try:
                del __warningregistry__
            except NameError:
                pass

            with warnings.catch_warnings(record=True) as warning_lines:
                warnings.resetwarnings()
                warnings.simplefilter(
                    "always", FITSFixedWarning, append=True)

                try:
                    WCS(hdu.header,
                        key=wcs.wcs.alt or ' ',
                        relax=_wcs.WCSHDR_reject,
                        fix=True, _do_set=False)
                except WcsError as e:
                    wcs_results.append(str(e))

                wcs_results.extend([str(x.message) for x in warning_lines])

    return results
