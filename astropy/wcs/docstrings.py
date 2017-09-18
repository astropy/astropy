# Licensed under a 3-clause BSD style license - see LICENSE.rst

# It gets to be really tedious to type long docstrings in ANSI C
# syntax (since multi-line string literals are not valid).
# Therefore, the docstrings are written here in doc/docstrings.py,
# which are then converted by setup.py into docstrings.h, which is
# included by pywcs.c

from . import _docutil as __

a = """
``double array[a_order+1][a_order+1]`` Focal plane transformation
matrix.

The `SIP`_ ``A_i_j`` matrix used for pixel to focal plane
transformation.

Its values may be changed in place, but it may not be resized, without
creating a new `~astropy.wcs.Sip` object.
"""

a_order = """
``int`` (read-only) Order of the polynomial (``A_ORDER``).
"""

all_pix2world = """
all_pix2world(pixcrd, origin) -> ``double array[ncoord][nelem]``

Transforms pixel coordinates to world coordinates.

Does the following:

    - Detector to image plane correction (if present)

    - SIP distortion correction (if present)

    - FITS WCS distortion correction (if present)

    - wcslib "core" WCS transformation

The first three (the distortion corrections) are done in parallel.

Parameters
----------
pixcrd : double array[ncoord][nelem]
    Array of pixel coordinates.

{0}

Returns
-------
world : double array[ncoord][nelem]
    Returns an array of world coordinates.

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
    Invalid coordinate transformation.

InvalidTransformError
    Ill-conditioned coordinate transformation parameters.
""".format(__.ORIGIN())

alt = """
``str`` Character code for alternate coordinate descriptions.

For example, the ``"a"`` in keyword names such as ``CTYPEia``.  This
is a space character for the primary coordinate description, or one of
the 26 upper-case letters, A-Z.
"""

ap = """
``double array[ap_order+1][ap_order+1]`` Focal plane to pixel
transformation matrix.

The `SIP`_ ``AP_i_j`` matrix used for focal plane to pixel
transformation.  Its values may be changed in place, but it may not be
resized, without creating a new `~astropy.wcs.Sip` object.
"""

ap_order = """
``int`` (read-only) Order of the polynomial (``AP_ORDER``).
"""

axis_types = """
``int array[naxis]`` An array of four-digit type codes for each axis.

- First digit (i.e. 1000s):

  - 0: Non-specific coordinate type.

  - 1: Stokes coordinate.

  - 2: Celestial coordinate (including ``CUBEFACE``).

  - 3: Spectral coordinate.

- Second digit (i.e. 100s):

  - 0: Linear axis.

  - 1: Quantized axis (``STOKES``, ``CUBEFACE``).

  - 2: Non-linear celestial axis.

  - 3: Non-linear spectral axis.

  - 4: Logarithmic axis.

  - 5: Tabular axis.

- Third digit (i.e. 10s):

  - 0: Group number, e.g. lookup table number

- The fourth digit is used as a qualifier depending on the axis type.

  - For celestial axes:

    - 0: Longitude coordinate.

    - 1: Latitude coordinate.

    - 2: ``CUBEFACE`` number.

  - For lookup tables: the axis number in a multidimensional table.

``CTYPEia`` in ``"4-3"`` form with unrecognized algorithm code will
have its type set to -1 and generate an error.
"""

b = """
``double array[b_order+1][b_order+1]`` Pixel to focal plane
transformation matrix.

The `SIP`_ ``B_i_j`` matrix used for pixel to focal plane
transformation.  Its values may be changed in place, but it may not be
resized, without creating a new `~astropy.wcs.Sip` object.
"""

b_order = """
``int`` (read-only) Order of the polynomial (``B_ORDER``).
"""

bounds_check = """
bounds_check(pix2world, world2pix)

Enable/disable bounds checking.

Parameters
----------
pix2world : bool, optional
    When `True`, enable bounds checking for the pixel-to-world (p2x)
    transformations.  Default is `True`.

world2pix : bool, optional
    When `True`, enable bounds checking for the world-to-pixel (s2x)
    transformations.  Default is `True`.

Notes
-----
Note that by default (without calling `bounds_check`) strict bounds
checking is enabled.
"""

bp = """
``double array[bp_order+1][bp_order+1]`` Focal plane to pixel
transformation matrix.

The `SIP`_ ``BP_i_j`` matrix used for focal plane to pixel
transformation.  Its values may be changed in place, but it may not be
resized, without creating a new `~astropy.wcs.Sip` object.
"""

bp_order = """
``int`` (read-only) Order of the polynomial (``BP_ORDER``).
"""

cd = """
``double array[naxis][naxis]`` The ``CDi_ja`` linear transformation
matrix.

For historical compatibility, three alternate specifications of the
linear transformations are available in wcslib.  The canonical
``PCi_ja`` with ``CDELTia``, ``CDi_ja``, and the deprecated
``CROTAia`` keywords.  Although the latter may not formally co-exist
with ``PCi_ja``, the approach here is simply to ignore them if given
in conjunction with ``PCi_ja``.

`~astropy.wcs.Wcsprm.has_pc`, `~astropy.wcs.Wcsprm.has_cd` and
`~astropy.wcs.Wcsprm.has_crota` can be used to determine which of
these alternatives are present in the header.

These alternate specifications of the linear transformation matrix are
translated immediately to ``PCi_ja`` by `~astropy.wcs.Wcsprm.set` and
are nowhere visible to the lower-level routines.  In particular,
`~astropy.wcs.Wcsprm.set` resets `~astropy.wcs.Wcsprm.cdelt` to unity
if ``CDi_ja`` is present (and no ``PCi_ja``).  If no ``CROTAia`` is
associated with the latitude axis, `~astropy.wcs.Wcsprm.set` reverts
to a unity ``PCi_ja`` matrix.
"""

cdelt = """
``double array[naxis]`` Coordinate increments (``CDELTia``) for each
coord axis.

If a ``CDi_ja`` linear transformation matrix is present, a warning is
raised and `~astropy.wcs.Wcsprm.cdelt` is ignored.  The ``CDi_ja``
matrix may be deleted by::

  del wcs.wcs.cd

An undefined value is represented by NaN.
"""

cdfix = """
cdfix()

Fix erroneously omitted ``CDi_ja`` keywords.

Sets the diagonal element of the ``CDi_ja`` matrix to unity if all
``CDi_ja`` keywords associated with a given axis were omitted.
According to Paper I, if any ``CDi_ja`` keywords at all are given in a
FITS header then those not given default to zero.  This results in a
singular matrix with an intersecting row and column of zeros.

Returns
-------
success : int
    Returns ``0`` for success; ``-1`` if no change required.
"""

cel_offset = """
``boolean`` Is there an offset?

If `True`, an offset will be applied to ``(x, y)`` to force ``(x, y) =
(0, 0)`` at the fiducial point, (phi_0, theta_0).  Default is `False`.
"""

celfix = """
Translates AIPS-convention celestial projection types, ``-NCP`` and
``-GLS``.

Returns
-------
success : int
    Returns ``0`` for success; ``-1`` if no change required.
"""

cname = """
``list of strings`` A list of the coordinate axis names, from
``CNAMEia``.
"""

colax = """
``int array[naxis]`` An array recording the column numbers for each
axis in a pixel list.
"""

colnum = """
``int`` Column of FITS binary table associated with this WCS.

Where the coordinate representation is associated with an image-array
column in a FITS binary table, this property may be used to record the
relevant column number.

It should be set to zero for an image header or pixel list.
"""

compare = """
compare(other, cmp=0, tolerance=0.0)

Compare two Wcsprm objects for equality.

Parameters
----------

other : Wcsprm
    The other Wcsprm object to compare to.

cmp : int, optional
    A bit field controlling the strictness of the comparison.  When 0,
    (the default), all fields must be identical.

    The following constants may be or'ed together to loosen the
    comparison.

    - ``WCSCOMPARE_ANCILLARY``: Ignores ancillary keywords that don't
      change the WCS transformation, such as ``DATE-OBS`` or
      ``EQUINOX``.

    - ``WCSCOMPARE_TILING``: Ignore integral differences in
      ``CRPIXja``.  This is the 'tiling' condition, where two WCSes
      cover different regions of the same map projection and align on
      the same map grid.

    - ``WCSCOMPARE_CRPIX``: Ignore any differences at all in
      ``CRPIXja``.  The two WCSes cover different regions of the same
      map projection but may not align on the same grid map.
      Overrides ``WCSCOMPARE_TILING``.

tolerance : float, optional
    The amount of tolerance required.  For example, for a value of
    1e-6, all floating-point values in the objects must be equal to
    the first 6 decimal places.  The default value of 0.0 implies
    exact equality.

Returns
-------
equal : bool
"""

convert = """
convert(array)

Perform the unit conversion on the elements of the given *array*,
returning an array of the same shape.
"""

coord = """
``double array[K_M]...[K_2][K_1][M]`` The tabular coordinate array.

Has the dimensions::

    (K_M, ... K_2, K_1, M)

(see `~astropy.wcs.Tabprm.K`) i.e. with the `M` dimension
varying fastest so that the `M` elements of a coordinate vector are
stored contiguously in memory.
"""

copy = """
Creates a deep copy of the WCS object.
"""

cpdis1 = """
`~astropy.wcs.DistortionLookupTable`

The pre-linear transformation distortion lookup table, ``CPDIS1``.
"""

cpdis2 = """
`~astropy.wcs.DistortionLookupTable`

The pre-linear transformation distortion lookup table, ``CPDIS2``.
"""

crder = """
``double array[naxis]`` The random error in each coordinate axis,
``CRDERia``.

An undefined value is represented by NaN.
"""

crota = """
``double array[naxis]`` ``CROTAia`` keyvalues for each coordinate
axis.

For historical compatibility, three alternate specifications of the
linear transformations are available in wcslib.  The canonical
``PCi_ja`` with ``CDELTia``, ``CDi_ja``, and the deprecated
``CROTAia`` keywords.  Although the latter may not formally co-exist
with ``PCi_ja``, the approach here is simply to ignore them if given
in conjunction with ``PCi_ja``.

`~astropy.wcs.Wcsprm.has_pc`, `~astropy.wcs.Wcsprm.has_cd` and
`~astropy.wcs.Wcsprm.has_crota` can be used to determine which of
these alternatives are present in the header.

These alternate specifications of the linear transformation matrix are
translated immediately to ``PCi_ja`` by `~astropy.wcs.Wcsprm.set` and
are nowhere visible to the lower-level routines.  In particular,
`~astropy.wcs.Wcsprm.set` resets `~astropy.wcs.Wcsprm.cdelt` to unity
if ``CDi_ja`` is present (and no ``PCi_ja``).  If no ``CROTAia`` is
associated with the latitude axis, `~astropy.wcs.Wcsprm.set` reverts
to a unity ``PCi_ja`` matrix.
"""

crpix = """
``double array[naxis]`` Coordinate reference pixels (``CRPIXja``) for
each pixel axis.
"""

crval = """
``double array[naxis]`` Coordinate reference values (``CRVALia``) for
each coordinate axis.
"""

crval_tabprm = """
``double array[M]`` Index values for the reference pixel for each of
the tabular coord axes.
"""

csyer = """
``double array[naxis]`` The systematic error in the coordinate value
axes, ``CSYERia``.

An undefined value is represented by NaN.
"""

ctype = """
``list of strings[naxis]`` List of ``CTYPEia`` keyvalues.

The `~astropy.wcs.Wcsprm.ctype` keyword values must be in upper case
and there must be zero or one pair of matched celestial axis types,
and zero or one spectral axis.
"""

cubeface = """
``int`` Index into the ``pixcrd`` (pixel coordinate) array for the
``CUBEFACE`` axis.

This is used for quadcube projections where the cube faces are stored
on a separate axis.

The quadcube projections (``TSC``, ``CSC``, ``QSC``) may be
represented in FITS in either of two ways:

    - The six faces may be laid out in one plane and numbered as
      follows::


                                       0

                              4  3  2  1  4  3  2

                                       5

      Faces 2, 3 and 4 may appear on one side or the other (or both).
      The world-to-pixel routines map faces 2, 3 and 4 to the left but
      the pixel-to-world routines accept them on either side.

    - The ``COBE`` convention in which the six faces are stored in a
      three-dimensional structure using a ``CUBEFACE`` axis indexed
      from 0 to 5 as above.

These routines support both methods; `~astropy.wcs.Wcsprm.set`
determines which is being used by the presence or absence of a
``CUBEFACE`` axis in `~astropy.wcs.Wcsprm.ctype`.
`~astropy.wcs.Wcsprm.p2s` and `~astropy.wcs.Wcsprm.s2p` translate the
``CUBEFACE`` axis representation to the single plane representation
understood by the lower-level projection routines.
"""

cunit = """
``list of astropy.UnitBase[naxis]`` List of ``CUNITia`` keyvalues as
`astropy.units.UnitBase` instances.

These define the units of measurement of the ``CRVALia``, ``CDELTia``
and ``CDi_ja`` keywords.

As ``CUNITia`` is an optional header keyword,
`~astropy.wcs.Wcsprm.cunit` may be left blank but otherwise is
expected to contain a standard units specification as defined by WCS
Paper I.  `~astropy.wcs.Wcsprm.unitfix` is available to translate
commonly used non-standard units specifications but this must be done
as a separate step before invoking `~astropy.wcs.Wcsprm.set`.

For celestial axes, if `~astropy.wcs.Wcsprm.cunit` is not blank,
`~astropy.wcs.Wcsprm.set` uses ``wcsunits`` to parse it and scale
`~astropy.wcs.Wcsprm.cdelt`, `~astropy.wcs.Wcsprm.crval`, and
`~astropy.wcs.Wcsprm.cd` to decimal degrees.  It then resets
`~astropy.wcs.Wcsprm.cunit` to ``"deg"``.

For spectral axes, if `~astropy.wcs.Wcsprm.cunit` is not blank,
`~astropy.wcs.Wcsprm.set` uses ``wcsunits`` to parse it and scale
`~astropy.wcs.Wcsprm.cdelt`, `~astropy.wcs.Wcsprm.crval`, and
`~astropy.wcs.Wcsprm.cd` to SI units.  It then resets
`~astropy.wcs.Wcsprm.cunit` accordingly.

`~astropy.wcs.Wcsprm.set` ignores `~astropy.wcs.Wcsprm.cunit` for
other coordinate types; `~astropy.wcs.Wcsprm.cunit` may be used to
label coordinate values.
"""

cylfix = """
cylfix()

Fixes WCS keyvalues for malformed cylindrical projections.

Returns
-------
success : int
    Returns ``0`` for success; ``-1`` if no change required.
"""

data = """
``float array`` The array data for the
`~astropy.wcs.DistortionLookupTable`.
"""

data_wtbarr = """
``double array``

The array data for the BINTABLE.
"""

dateavg = """
``string`` Representative mid-point of the date of observation.

In ISO format, ``yyyy-mm-ddThh:mm:ss``.

See also
--------
astropy.wcs.Wcsprm.dateobs
"""

dateobs = """
``string`` Start of the date of observation.

In ISO format, ``yyyy-mm-ddThh:mm:ss``.

See also
--------
astropy.wcs.Wcsprm.dateavg
"""

datfix = """
datfix()

Translates the old ``DATE-OBS`` date format to year-2000 standard form
``(yyyy-mm-ddThh:mm:ss)`` and derives ``MJD-OBS`` from it if not
already set.

Alternatively, if `~astropy.wcs.Wcsprm.mjdobs` is set and
`~astropy.wcs.Wcsprm.dateobs` isn't, then `~astropy.wcs.Wcsprm.datfix`
derives `~astropy.wcs.Wcsprm.dateobs` from it.  If both are set but
disagree by more than half a day then `ValueError` is raised.

Returns
-------
success : int
    Returns ``0`` for success; ``-1`` if no change required.
"""

delta = """
``double array[M]`` (read-only) Interpolated indices into the coord
array.

Array of interpolated indices into the coordinate array such that
Upsilon_m, as defined in Paper III, is equal to
(`~astropy.wcs.Tabprm.p0` [m] + 1) + delta[m].
"""

det2im = """
Convert detector coordinates to image plane coordinates.
"""

det2im1 = """
A `~astropy.wcs.DistortionLookupTable` object for detector to image plane
correction in the *x*-axis.
"""

det2im2 = """
A `~astropy.wcs.DistortionLookupTable` object for detector to image plane
correction in the *y*-axis.
"""

dims = """
``int array[ndim]`` (read-only)

The dimensions of the tabular array
`~astropy.wcs.Wtbarr.data`.
"""

DistortionLookupTable = """
DistortionLookupTable(*table*, *crpix*, *crval*, *cdelt*)

Represents a single lookup table for a `distortion paper`_
transformation.

Parameters
----------
table : 2-dimensional array
    The distortion lookup table.

crpix : 2-tuple
    The distortion array reference pixel

crval : 2-tuple
    The image array pixel coordinate

cdelt : 2-tuple
    The grid step size
"""

equinox = """
``double`` The equinox associated with dynamical equatorial or
ecliptic coordinate systems.

``EQUINOXa`` (or ``EPOCH`` in older headers).  Not applicable to ICRS
equatorial or ecliptic coordinates.

An undefined value is represented by NaN.
"""

extlev = """
``int`` (read-only)

``EXTLEV`` identifying the binary table extension.
"""

extnam = """
``str`` (read-only)

``EXTNAME`` identifying the binary table extension.
"""

extrema = """
``double array[K_M]...[K_2][2][M]`` (read-only)

An array recording the minimum and maximum value of each element of
the coordinate vector in each row of the coordinate array, with the
dimensions::

    (K_M, ... K_2, 2, M)

(see `~astropy.wcs.Tabprm.K`).  The minimum is recorded
in the first element of the compressed K_1 dimension, then the
maximum.  This array is used by the inverse table lookup function to
speed up table searches.
"""

extver = """
``int`` (read-only)

``EXTVER`` identifying the binary table extension.
"""

find_all_wcs = """
find_all_wcs(relax=0, keysel=0)

Find all WCS transformations in the header.

Parameters
----------

header : str
    The raw FITS header data.

relax : bool or int
    Degree of permissiveness:

    - `False`: Recognize only FITS keywords defined by the published
      WCS standard.

    - `True`: Admit all recognized informal extensions of the WCS
      standard.

    - `int`: a bit field selecting specific extensions to accept.  See
      :ref:`relaxread` for details.

keysel : sequence of flags
    Used to restrict the keyword types considered:

    - ``WCSHDR_IMGHEAD``: Image header keywords.

    - ``WCSHDR_BIMGARR``: Binary table image array.

    - ``WCSHDR_PIXLIST``: Pixel list keywords.

    If zero, there is no restriction.  If -1, `wcspih` is called,
    rather than `wcstbh`.

Returns
-------
wcs_list : list of `~astropy.wcs.Wcsprm` objects
"""

fix = """
fix(translate_units='', naxis=0)

Applies all of the corrections handled separately by
`~astropy.wcs.Wcsprm.datfix`, `~astropy.wcs.Wcsprm.unitfix`,
`~astropy.wcs.Wcsprm.celfix`, `~astropy.wcs.Wcsprm.spcfix`,
`~astropy.wcs.Wcsprm.cylfix` and `~astropy.wcs.Wcsprm.cdfix`.

Parameters
----------

translate_units : str, optional
    Specify which potentially unsafe translations of non-standard unit
    strings to perform.  By default, performs all.

    Although ``"S"`` is commonly used to represent seconds, its
    translation to ``"s"`` is potentially unsafe since the standard
    recognizes ``"S"`` formally as Siemens, however rarely that may be
    used.  The same applies to ``"H"`` for hours (Henry), and ``"D"``
    for days (Debye).

    This string controls what to do in such cases, and is
    case-insensitive.

    - If the string contains ``"s"``, translate ``"S"`` to ``"s"``.

    - If the string contains ``"h"``, translate ``"H"`` to ``"h"``.

    - If the string contains ``"d"``, translate ``"D"`` to ``"d"``.

    Thus ``''`` doesn't do any unsafe translations, whereas ``'shd'``
    does all of them.

naxis : int array[naxis], optional
    Image axis lengths.  If this array is set to zero or ``None``,
    then `~astropy.wcs.Wcsprm.cylfix` will not be invoked.

Returns
-------
status : dict

    Returns a dictionary containing the following keys, each referring
    to a status string for each of the sub-fix functions that were
    called:

    - `~astropy.wcs.Wcsprm.cdfix`

    - `~astropy.wcs.Wcsprm.datfix`

    - `~astropy.wcs.Wcsprm.unitfix`

    - `~astropy.wcs.Wcsprm.celfix`

    - `~astropy.wcs.Wcsprm.spcfix`

    - `~astropy.wcs.Wcsprm.cylfix`
"""

get_offset = """
get_offset(x, y) -> (x, y)

Returns the offset as defined in the distortion lookup table.

Returns
-------
coordinate : coordinate pair
    The offset from the distortion table for pixel point (*x*, *y*).
"""

get_cdelt = """
get_cdelt() -> double array[naxis]

Coordinate increments (``CDELTia``) for each coord axis.

Returns the ``CDELT`` offsets in read-only form.  Unlike the
`~astropy.wcs.Wcsprm.cdelt` property, this works even when the header
specifies the linear transformation matrix in one of the alternative
``CDi_ja`` or ``CROTAia`` forms.  This is useful when you want access
to the linear transformation matrix, but don't care how it was
specified in the header.
"""

get_pc = """
get_pc() -> double array[naxis][naxis]

Returns the ``PC`` matrix in read-only form.  Unlike the
`~astropy.wcs.Wcsprm.pc` property, this works even when the header
specifies the linear transformation matrix in one of the alternative
``CDi_ja`` or ``CROTAia`` forms.  This is useful when you want access
to the linear transformation matrix, but don't care how it was
specified in the header.
"""

get_ps = """
get_ps() -> list of tuples

Returns ``PSi_ma`` keywords for each *i* and *m*.

Returns
-------
ps : list of tuples

    Returned as a list of tuples of the form (*i*, *m*, *value*):

    - *i*: int.  Axis number, as in ``PSi_ma``, (i.e. 1-relative)

    - *m*: int.  Parameter number, as in ``PSi_ma``, (i.e. 0-relative)

    - *value*: string.  Parameter value.

See also
--------
astropy.wcs.Wcsprm.set_ps : Set ``PSi_ma`` values
"""

get_pv = """
get_pv() -> list of tuples

Returns ``PVi_ma`` keywords for each *i* and *m*.

Returns
-------

    Returned as a list of tuples of the form (*i*, *m*, *value*):

    - *i*: int.  Axis number, as in ``PVi_ma``, (i.e. 1-relative)

    - *m*: int.  Parameter number, as in ``PVi_ma``, (i.e. 0-relative)

    - *value*: string. Parameter value.

See also
--------
astropy.wcs.Wcsprm.set_pv : Set ``PVi_ma`` values

Notes
-----

Note that, if they were not given, `~astropy.wcs.Wcsprm.set` resets
the entries for ``PVi_1a``, ``PVi_2a``, ``PVi_3a``, and ``PVi_4a`` for
longitude axis *i* to match (``phi_0``, ``theta_0``), the native
longitude and latitude of the reference point given by ``LONPOLEa``
and ``LATPOLEa``.
"""

has_cd = """
has_cd() -> bool

Returns `True` if ``CDi_ja`` is present.

``CDi_ja`` is an alternate specification of the linear transformation
matrix, maintained for historical compatibility.

Matrix elements in the IRAF convention are equivalent to the product
``CDi_ja = CDELTia * PCi_ja``, but the defaults differ from that of
the ``PCi_ja`` matrix.  If one or more ``CDi_ja`` keywords are present
then all unspecified ``CDi_ja`` default to zero.  If no ``CDi_ja`` (or
``CROTAia``) keywords are present, then the header is assumed to be in
``PCi_ja`` form whether or not any ``PCi_ja`` keywords are present
since this results in an interpretation of ``CDELTia`` consistent with
the original FITS specification.

While ``CDi_ja`` may not formally co-exist with ``PCi_ja``, it may
co-exist with ``CDELTia`` and ``CROTAia`` which are to be ignored.

See also
--------
astropy.wcs.Wcsprm.cd : Get the raw ``CDi_ja`` values.
"""

has_cdi_ja = """
has_cdi_ja() -> bool

Alias for `~astropy.wcs.Wcsprm.has_cd`.  Maintained for backward
compatibility.
"""

has_crota = """
has_crota() -> bool

Returns `True` if ``CROTAia`` is present.

``CROTAia`` is an alternate specification of the linear transformation
matrix, maintained for historical compatibility.

In the AIPS convention, ``CROTAia`` may only be associated with the
latitude axis of a celestial axis pair.  It specifies a rotation in
the image plane that is applied *after* the ``CDELTia``; any other
``CROTAia`` keywords are ignored.

``CROTAia`` may not formally co-exist with ``PCi_ja``.  ``CROTAia`` and
``CDELTia`` may formally co-exist with ``CDi_ja`` but if so are to be
ignored.

See also
--------
astropy.wcs.Wcsprm.crota : Get the raw ``CROTAia`` values
"""

has_crotaia = """
has_crotaia() -> bool

Alias for `~astropy.wcs.Wcsprm.has_crota`.  Maintained for backward
compatibility.
"""

has_pc = """
has_pc() -> bool

Returns `True` if ``PCi_ja`` is present.  ``PCi_ja`` is the
recommended way to specify the linear transformation matrix.

See also
--------
astropy.wcs.Wcsprm.pc : Get the raw ``PCi_ja`` values
"""

has_pci_ja = """
has_pci_ja() -> bool

Alias for `~astropy.wcs.Wcsprm.has_pc`.  Maintained for backward
compatibility.
"""

i = """
``int`` (read-only)

Image axis number.
"""

imgpix_matrix = """
``double array[2][2]`` (read-only) Inverse of the ``CDELT`` or ``PC``
matrix.

Inverse containing the product of the ``CDELTia`` diagonal matrix and
the ``PCi_ja`` matrix.
"""

is_unity = """
is_unity() -> bool

Returns `True` if the linear transformation matrix
(`~astropy.wcs.Wcsprm.cd`) is unity.
"""

K = """
``int array[M]`` (read-only) The lengths of the axes of the coordinate
array.

An array of length `M` whose elements record the lengths of the axes of
the coordinate array and of each indexing vector.
"""

kind = """
``str`` (read-only)

Character identifying the wcstab array type:

    - ``'c'``: coordinate array,
    - ``'i'``: index vector.
"""

lat = """
``int`` (read-only) The index into the world coord array containing
latitude values.
"""

latpole = """
``double`` The native latitude of the celestial pole, ``LATPOLEa`` (deg).
"""

lattyp = """
``string`` (read-only) Celestial axis type for latitude.

For example, "RA", "DEC", "GLON", "GLAT", etc. extracted from "RA--",
"DEC-", "GLON", "GLAT", etc. in the first four characters of
``CTYPEia`` but with trailing dashes removed.
"""

lng = """
``int`` (read-only) The index into the world coord array containing
longitude values.
"""

lngtyp = """
``string`` (read-only) Celestial axis type for longitude.

For example, "RA", "DEC", "GLON", "GLAT", etc. extracted from "RA--",
"DEC-", "GLON", "GLAT", etc. in the first four characters of
``CTYPEia`` but with trailing dashes removed.
"""

lonpole = """
``double`` The native longitude of the celestial pole.

``LONPOLEa`` (deg).
"""

M = """
``int`` (read-only) Number of tabular coordinate axes.
"""

m = """
``int`` (read-only)

Array axis number for index vectors.
"""

map = """
``int array[M]`` Association between axes.

A vector of length `~astropy.wcs.Tabprm.M` that defines
the association between axis *m* in the *M*-dimensional coordinate
array (1 <= *m* <= *M*) and the indices of the intermediate world
coordinate and world coordinate arrays.

When the intermediate and world coordinate arrays contain the full
complement of coordinate elements in image-order, as will usually be
the case, then ``map[m-1] == i-1`` for axis *i* in the *N*-dimensional
image (1 <= *i* <= *N*).  In terms of the FITS keywords::

    map[PVi_3a - 1] == i - 1.

However, a different association may result if the intermediate
coordinates, for example, only contains a (relevant) subset of
intermediate world coordinate elements.  For example, if *M* == 1 for
an image with *N* > 1, it is possible to fill the intermediate
coordinates with the relevant coordinate element with ``nelem`` set to
1.  In this case ``map[0] = 0`` regardless of the value of *i*.
"""

mix = """
mix(mixpix, mixcel, vspan, vstep, viter, world, pixcrd, origin)

Given either the celestial longitude or latitude plus an element of
the pixel coordinate, solves for the remaining elements by iterating
on the unknown celestial coordinate element using
`~astropy.wcs.Wcsprm.s2p`.

Parameters
----------
mixpix : int
    Which element on the pixel coordinate is given.

mixcel : int
    Which element of the celestial coordinate is given. If *mixcel* =
    ``1``, celestial longitude is given in ``world[self.lng]``,
    latitude returned in ``world[self.lat]``.  If *mixcel* = ``2``,
    celestial latitude is given in ``world[self.lat]``, longitude
    returned in ``world[self.lng]``.

vspan : pair of floats
    Solution interval for the celestial coordinate, in degrees.  The
    ordering of the two limits is irrelevant.  Longitude ranges may be
    specified with any convenient normalization, for example
    ``(-120,+120)`` is the same as ``(240,480)``, except that the
    solution will be returned with the same normalization, i.e. lie
    within the interval specified.

vstep : float
    Step size for solution search, in degrees.  If ``0``, a sensible,
    although perhaps non-optimal default will be used.

viter : int
    If a solution is not found then the step size will be halved and
    the search recommenced.  *viter* controls how many times the step
    size is halved.  The allowed range is 5 - 10.

world : double array[naxis]
    World coordinate elements.  ``world[self.lng]`` and
    ``world[self.lat]`` are the celestial longitude and latitude, in
    degrees.  Which is given and which returned depends on the value
    of *mixcel*.  All other elements are given.  The results will be
    written to this array in-place.

pixcrd : double array[naxis].
    Pixel coordinates.  The element indicated by *mixpix* is given and
    the remaining elements will be written in-place.

{0}

Returns
-------
result : dict

    Returns a dictionary with the following keys:

    - *phi* (double array[naxis])

    - *theta* (double array[naxis])

        - Longitude and latitude in the native coordinate system of
          the projection, in degrees.

    - *imgcrd* (double array[naxis])

        - Image coordinate elements.  ``imgcrd[self.lng]`` and
          ``imgcrd[self.lat]`` are the projected *x*- and
          *y*-coordinates, in decimal degrees.

    - *world* (double array[naxis])

        - Another reference to the *world* argument passed in.

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

InvalidTransformError
    Invalid coordinate transformation parameters.

InvalidTransformError
    Ill-conditioned coordinate transformation parameters.

InvalidCoordinateError
    Invalid world coordinate.

NoSolutionError
    No solution found in the specified interval.

See also
--------
astropy.wcs.Wcsprm.lat, astropy.wcs.Wcsprm.lng
    Get the axes numbers for latitude and longitude

Notes
-----

Initially, the specified solution interval is checked to see if it's a
\"crossing\" interval.  If it isn't, a search is made for a crossing
solution by iterating on the unknown celestial coordinate starting at
the upper limit of the solution interval and decrementing by the
specified step size.  A crossing is indicated if the trial value of
the pixel coordinate steps through the value specified.  If a crossing
interval is found then the solution is determined by a modified form
of \"regula falsi\" division of the crossing interval.  If no crossing
interval was found within the specified solution interval then a
search is made for a \"non-crossing\" solution as may arise from a
point of tangency.  The process is complicated by having to make
allowance for the discontinuities that occur in all map projections.

Once one solution has been determined others may be found by
subsequent invocations of `~astropy.wcs.Wcsprm.mix` with suitably
restricted solution intervals.

Note the circumstance that arises when the solution point lies at a
native pole of a projection in which the pole is represented as a
finite curve, for example the zenithals and conics.  In such cases two
or more valid solutions may exist but `~astropy.wcs.Wcsprm.mix` only
ever returns one.

Because of its generality, `~astropy.wcs.Wcsprm.mix` is very
compute-intensive.  For compute-limited applications, more efficient
special-case solvers could be written for simple projections, for
example non-oblique cylindrical projections.
""".format(__.ORIGIN())

mjdavg = """
``double`` Modified Julian Date corresponding to ``DATE-AVG``.

``(MJD = JD - 2400000.5)``.

An undefined value is represented by NaN.

See also
--------
astropy.wcs.Wcsprm.mjdobs
"""

mjdobs = """
``double`` Modified Julian Date corresponding to ``DATE-OBS``.

``(MJD = JD - 2400000.5)``.

An undefined value is represented by NaN.

See also
--------
astropy.wcs.Wcsprm.mjdavg
"""

name = """
``string`` The name given to the coordinate representation
``WCSNAMEa``.
"""

naxis = """
``int`` (read-only) The number of axes (pixel and coordinate).

Given by the ``NAXIS`` or ``WCSAXESa`` keyvalues.

The number of coordinate axes is determined at parsing time, and can
not be subsequently changed.

It is determined from the highest of the following:

  1. ``NAXIS``

  2. ``WCSAXESa``

  3. The highest axis number in any parameterized WCS keyword.  The
     keyvalue, as well as the keyword, must be syntactically valid
     otherwise it will not be considered.

If none of these keyword types is present, i.e. if the header only
contains auxiliary WCS keywords for a particular coordinate
representation, then no coordinate description is constructed for it.

This value may differ for different coordinate representations of the
same image.
"""

nc = """
``int`` (read-only) Total number of coord vectors in the coord array.

Total number of coordinate vectors in the coordinate array being the
product K_1 * K_2 * ... * K_M.
"""

ndim = """
``int`` (read-only)

Expected dimensionality of the wcstab array.
"""

obsgeo = """
``double array[3]`` Location of the observer in a standard terrestrial
reference frame.

``OBSGEO-X``, ``OBSGEO-Y``, ``OBSGEO-Z`` (in meters).

An undefined value is represented by NaN.
"""

p0 = """
``int array[M]`` Interpolated indices into the coordinate array.

Vector of length `~astropy.wcs.Tabprm.M` of interpolated
indices into the coordinate array such that Upsilon_m, as defined in
Paper III, is equal to ``(p0[m] + 1) + delta[m]``.
"""

p2s = """
p2s(pixcrd, origin)

Converts pixel to world coordinates.

Parameters
----------

pixcrd : double array[ncoord][nelem]
    Array of pixel coordinates.

{0}

Returns
-------
result : dict
    Returns a dictionary with the following keys:

    - *imgcrd*: double array[ncoord][nelem]

      - Array of intermediate world coordinates.  For celestial axes,
        ``imgcrd[][self.lng]`` and ``imgcrd[][self.lat]`` are the
        projected *x*-, and *y*-coordinates, in pseudo degrees.  For
        spectral axes, ``imgcrd[][self.spec]`` is the intermediate
        spectral coordinate, in SI units.

    - *phi*: double array[ncoord]

    - *theta*: double array[ncoord]

      - Longitude and latitude in the native coordinate system of the
        projection, in degrees.

    - *world*: double array[ncoord][nelem]

      - Array of world coordinates.  For celestial axes,
        ``world[][self.lng]`` and ``world[][self.lat]`` are the
        celestial longitude and latitude, in degrees.  For spectral
        axes, ``world[][self.spec]`` is the intermediate spectral
        coordinate, in SI units.

    - *stat*: int array[ncoord]

      - Status return value for each coordinate. ``0`` for success,
        ``1+`` for invalid pixel coordinate.

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
    *x*- and *y*-coordinate arrays are not the same size.

InvalidTransformError
    Invalid coordinate transformation parameters.

InvalidTransformError
    Ill-conditioned coordinate transformation parameters.

See also
--------
astropy.wcs.Wcsprm.lat, astropy.wcs.Wcsprm.lng
    Definition of the latitude and longitude axes
""".format(__.ORIGIN())

p4_pix2foc = """
p4_pix2foc(*pixcrd, origin*) -> double array[ncoord][nelem]

Convert pixel coordinates to focal plane coordinates using `distortion
paper`_ lookup-table correction.

Parameters
----------
pixcrd : double array[ncoord][nelem].
    Array of pixel coordinates.

{0}

Returns
-------
foccrd : double array[ncoord][nelem]
    Returns an array of focal plane coordinates.

Raises
------
MemoryError
    Memory allocation failed.

ValueError
    Invalid coordinate transformation parameters.
""".format(__.ORIGIN())

pc = """
``double array[naxis][naxis]`` The ``PCi_ja`` (pixel coordinate)
transformation matrix.

The order is::

  [[PC1_1, PC1_2],
   [PC2_1, PC2_2]]

For historical compatibility, three alternate specifications of the
linear transformations are available in wcslib.  The canonical
``PCi_ja`` with ``CDELTia``, ``CDi_ja``, and the deprecated
``CROTAia`` keywords.  Although the latter may not formally co-exist
with ``PCi_ja``, the approach here is simply to ignore them if given
in conjunction with ``PCi_ja``.

`~astropy.wcs.Wcsprm.has_pc`, `~astropy.wcs.Wcsprm.has_cd` and
`~astropy.wcs.Wcsprm.has_crota` can be used to determine which of
these alternatives are present in the header.

These alternate specifications of the linear transformation matrix are
translated immediately to ``PCi_ja`` by `~astropy.wcs.Wcsprm.set` and
are nowhere visible to the lower-level routines.  In particular,
`~astropy.wcs.Wcsprm.set` resets `~astropy.wcs.Wcsprm.cdelt` to unity
if ``CDi_ja`` is present (and no ``PCi_ja``).  If no ``CROTAia`` is
associated with the latitude axis, `~astropy.wcs.Wcsprm.set` reverts
to a unity ``PCi_ja`` matrix.
"""

phi0 = """
``double`` The native latitude of the fiducial point.

The point whose celestial coordinates are given in ``ref[1:2]``.  If
undefined (NaN) the initialization routine, `~astropy.wcs.Wcsprm.set`,
will set this to a projection-specific default.

See also
--------
astropy.wcs.Wcsprm.theta0
"""

pix2foc = """
pix2foc(*pixcrd, origin*) -> double array[ncoord][nelem]

Perform both `SIP`_ polynomial and `distortion paper`_ lookup-table
correction in parallel.

Parameters
----------
pixcrd : double array[ncoord][nelem]
    Array of pixel coordinates.

{0}

Returns
-------
foccrd : double array[ncoord][nelem]
    Returns an array of focal plane coordinates.

Raises
------
MemoryError
    Memory allocation failed.

ValueError
    Invalid coordinate transformation parameters.
""".format(__.ORIGIN())

piximg_matrix = """
``double array[2][2]`` (read-only) Matrix containing the product of
the ``CDELTia`` diagonal matrix and the ``PCi_ja`` matrix.
"""

print_contents = """
print_contents()

Print the contents of the `~astropy.wcs.Wcsprm` object to stdout.
Probably only useful for debugging purposes, and may be removed in the
future.

To get a string of the contents, use `repr`.
"""

print_contents_tabprm = """
print_contents()

Print the contents of the `~astropy.wcs.Tabprm` object to
stdout.  Probably only useful for debugging purposes, and may be
removed in the future.

To get a string of the contents, use `repr`.
"""

radesys = """
``string`` The equatorial or ecliptic coordinate system type,
``RADESYSa``.
"""

restfrq = """
``double`` Rest frequency (Hz) from ``RESTFRQa``.

An undefined value is represented by NaN.
"""

restwav = """
``double`` Rest wavelength (m) from ``RESTWAVa``.

An undefined value is represented by NaN.
"""

row = """
``int`` (read-only)

Table row number.
"""

s2p = """
s2p(world, origin)

Transforms world coordinates to pixel coordinates.

Parameters
----------
world : double array[ncoord][nelem]
    Array of world coordinates, in decimal degrees.

{0}

Returns
-------
result : dict
    Returns a dictionary with the following keys:

    - *phi*: double array[ncoord]

    - *theta*: double array[ncoord]

        - Longitude and latitude in the native coordinate system of
          the projection, in degrees.

    - *imgcrd*: double array[ncoord][nelem]

       - Array of intermediate world coordinates.  For celestial axes,
         ``imgcrd[][self.lng]`` and ``imgcrd[][self.lat]`` are the
         projected *x*-, and *y*-coordinates, in pseudo \"degrees\".
         For quadcube projections with a ``CUBEFACE`` axis, the face
         number is also returned in ``imgcrd[][self.cubeface]``.  For
         spectral axes, ``imgcrd[][self.spec]`` is the intermediate
         spectral coordinate, in SI units.

    - *pixcrd*: double array[ncoord][nelem]

        - Array of pixel coordinates.  Pixel coordinates are
          zero-based.

    - *stat*: int array[ncoord]

        - Status return value for each coordinate. ``0`` for success,
          ``1+`` for invalid pixel coordinate.

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

InvalidTransformError
   Invalid coordinate transformation parameters.

InvalidTransformError
    Ill-conditioned coordinate transformation parameters.

See also
--------
astropy.wcs.Wcsprm.lat, astropy.wcs.Wcsprm.lng
    Definition of the latitude and longitude axes
""".format(__.ORIGIN())

sense = """
``int array[M]`` +1 if monotonically increasing, -1 if decreasing.

A vector of length `~astropy.wcs.Tabprm.M` whose elements
indicate whether the corresponding indexing vector is monotonically
increasing (+1), or decreasing (-1).
"""

set = """
set()

Sets up a WCS object for use according to information supplied within
it.

Note that this routine need not be called directly; it will be invoked
by `~astropy.wcs.Wcsprm.p2s` and `~astropy.wcs.Wcsprm.s2p` if
necessary.

Some attributes that are based on other attributes (such as
`~astropy.wcs.Wcsprm.lattyp` on `~astropy.wcs.Wcsprm.ctype`) may not
be correct until after `~astropy.wcs.Wcsprm.set` is called.

`~astropy.wcs.Wcsprm.set` strips off trailing blanks in all string
members.

`~astropy.wcs.Wcsprm.set` recognizes the ``NCP`` projection and
converts it to the equivalent ``SIN`` projection and it also
recognizes ``GLS`` as a synonym for ``SFL``.  It does alias
translation for the AIPS spectral types (``FREQ-LSR``, ``FELO-HEL``,
etc.) but without changing the input header keywords.

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

InvalidTransformError
    Invalid coordinate transformation parameters.

InvalidTransformError
    Ill-conditioned coordinate transformation parameters.
"""

set_tabprm = """
set()

Allocates memory for work arrays.

Also sets up the class according to information supplied within it.

Note that this routine need not be called directly; it will be invoked
by functions that need it.

Raises
------
MemoryError
    Memory allocation failed.

InvalidTabularParameters
    Invalid tabular parameters.
"""

set_ps = """
set_ps(ps)

Sets ``PSi_ma`` keywords for each *i* and *m*.

Parameters
----------
ps : sequence of tuples

    The input must be a sequence of tuples of the form (*i*, *m*,
    *value*):

    - *i*: int.  Axis number, as in ``PSi_ma``, (i.e. 1-relative)

    - *m*: int.  Parameter number, as in ``PSi_ma``, (i.e. 0-relative)

    - *value*: string.  Parameter value.

See also
--------
astropy.wcs.Wcsprm.get_ps
"""

set_pv = """
set_pv(pv)

Sets ``PVi_ma`` keywords for each *i* and *m*.

Parameters
----------
pv : list of tuples

    The input must be a sequence of tuples of the form (*i*, *m*,
    *value*):

    - *i*: int.  Axis number, as in ``PVi_ma``, (i.e. 1-relative)

    - *m*: int.  Parameter number, as in ``PVi_ma``, (i.e. 0-relative)

    - *value*: float.  Parameter value.

See also
--------
astropy.wcs.Wcsprm.get_pv
"""

sip = """
Get/set the `~astropy.wcs.Sip` object for performing `SIP`_ distortion
correction.
"""

Sip = """
Sip(*a, b, ap, bp, crpix*)

The `~astropy.wcs.Sip` class performs polynomial distortion correction
using the `SIP`_ convention in both directions.

Parameters
----------
a : double array[m+1][m+1]
    The ``A_i_j`` polynomial for pixel to focal plane transformation.
    Its size must be (*m* + 1, *m* + 1) where *m* = ``A_ORDER``.

b : double array[m+1][m+1]
    The ``B_i_j`` polynomial for pixel to focal plane transformation.
    Its size must be (*m* + 1, *m* + 1) where *m* = ``B_ORDER``.

ap : double array[m+1][m+1]
    The ``AP_i_j`` polynomial for pixel to focal plane transformation.
    Its size must be (*m* + 1, *m* + 1) where *m* = ``AP_ORDER``.

bp : double array[m+1][m+1]
    The ``BP_i_j`` polynomial for pixel to focal plane transformation.
    Its size must be (*m* + 1, *m* + 1) where *m* = ``BP_ORDER``.

crpix : double array[2]
    The reference pixel.

Notes
-----
Shupe, D. L., M. Moshir, J. Li, D. Makovoz and R. Narron.  2005.
"The SIP Convention for Representing Distortion in FITS Image
Headers."  ADASS XIV.
"""

sip_foc2pix = """
sip_foc2pix(*foccrd, origin*) -> double array[ncoord][nelem]

Convert focal plane coordinates to pixel coordinates using the `SIP`_
polynomial distortion convention.

Parameters
----------
foccrd : double array[ncoord][nelem]
    Array of focal plane coordinates.

{0}

Returns
-------
pixcrd : double array[ncoord][nelem]
    Returns an array of pixel coordinates.

Raises
------
MemoryError
    Memory allocation failed.

ValueError
    Invalid coordinate transformation parameters.
""".format(__.ORIGIN())

sip_pix2foc = """
sip_pix2foc(*pixcrd, origin*) -> double array[ncoord][nelem]

Convert pixel coordinates to focal plane coordinates using the `SIP`_
polynomial distortion convention.

Parameters
----------
pixcrd : double array[ncoord][nelem]
    Array of pixel coordinates.

{0}

Returns
-------
foccrd : double array[ncoord][nelem]
    Returns an array of focal plane coordinates.

Raises
------
MemoryError
    Memory allocation failed.

ValueError
    Invalid coordinate transformation parameters.
""".format(__.ORIGIN())

spcfix = """
spcfix() -> int

Translates AIPS-convention spectral coordinate types.  {``FREQ``,
``VELO``, ``FELO``}-{``OBS``, ``HEL``, ``LSR``} (e.g. ``FREQ-LSR``,
``VELO-OBS``, ``FELO-HEL``)

Returns
-------
success : int
    Returns ``0`` for success; ``-1`` if no change required.
"""

spec = """
``int`` (read-only) The index containing the spectral axis values.
"""

specsys = """
``string`` Spectral reference frame (standard of rest), ``SPECSYSa``.

See also
--------
astropy.wcs.Wcsprm.ssysobs, astropy.wcs.Wcsprm.velosys
"""

sptr = """
sptr(ctype, i=-1)

Translates the spectral axis in a WCS object.

For example, a ``FREQ`` axis may be translated into ``ZOPT-F2W`` and
vice versa.

Parameters
----------
ctype : str
    Required spectral ``CTYPEia``, maximum of 8 characters.  The first
    four characters are required to be given and are never modified.
    The remaining four, the algorithm code, are completely determined
    by, and must be consistent with, the first four characters.
    Wildcarding may be used, i.e.  if the final three characters are
    specified as ``\"???\"``, or if just the eighth character is
    specified as ``\"?\"``, the correct algorithm code will be
    substituted and returned.

i : int
    Index of the spectral axis (0-relative).  If ``i < 0`` (or not
    provided), it will be set to the first spectral axis identified
    from the ``CTYPE`` keyvalues in the FITS header.

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

InvalidTransformError
    Invalid coordinate transformation parameters.

InvalidTransformError
    Ill-conditioned coordinate transformation parameters.

InvalidSubimageSpecificationError
    Invalid subimage specification (no spectral axis).
"""

ssysobs = """
``string`` Spectral reference frame.

The spectral reference frame in which there is no differential
variation in the spectral coordinate across the field-of-view,
``SSYSOBSa``.

See also
--------
astropy.wcs.Wcsprm.specsys, astropy.wcs.Wcsprm.velosys
"""

ssyssrc = """
``string`` Spectral reference frame for redshift.

The spectral reference frame (standard of rest) in which the redshift
was measured, ``SSYSSRCa``.
"""

sub = """
sub(axes)

Extracts the coordinate description for a subimage from a
`~astropy.wcs.WCS` object.

The world coordinate system of the subimage must be separable in the
sense that the world coordinates at any point in the subimage must
depend only on the pixel coordinates of the axes extracted.  In
practice, this means that the ``PCi_ja`` matrix of the original image
must not contain non-zero off-diagonal terms that associate any of the
subimage axes with any of the non-subimage axes.

`sub` can also add axes to a wcsprm object.  The new axes will be
created using the defaults set by the Wcsprm constructor which produce
a simple, unnamed, linear axis with world coordinates equal to the
pixel coordinate.  These default values can be changed before
invoking `set`.

Parameters
----------
axes : int or a sequence.

    - If an int, include the first *N* axes in their original order.

    - If a sequence, may contain a combination of image axis numbers
      (1-relative) or special axis identifiers (see below).  Order is
      significant; ``axes[0]`` is the axis number of the input image
      that corresponds to the first axis in the subimage, etc.  Use an
      axis number of 0 to create a new axis using the defaults.

    - If ``0``, ``[]`` or ``None``, do a deep copy.

    Coordinate axes types may be specified using either strings or
    special integer constants.  The available types are:

    - ``'longitude'`` / ``WCSSUB_LONGITUDE``: Celestial longitude

    - ``'latitude'`` / ``WCSSUB_LATITUDE``: Celestial latitude

    - ``'cubeface'`` / ``WCSSUB_CUBEFACE``: Quadcube ``CUBEFACE`` axis

    - ``'spectral'`` / ``WCSSUB_SPECTRAL``: Spectral axis

    - ``'stokes'`` / ``WCSSUB_STOKES``: Stokes axis

    - ``'celestial'`` / ``WCSSUB_CELESTIAL``: An alias for the
      combination of ``'longitude'``, ``'latitude'`` and ``'cubeface'``.

Returns
-------
new_wcs : `~astropy.wcs.WCS` object

Raises
------
MemoryError
    Memory allocation failed.

InvalidSubimageSpecificationError
    Invalid subimage specification (no spectral axis).

NonseparableSubimageCoordinateSystem
    Non-separable subimage coordinate system.

Notes
-----
Combinations of subimage axes of particular types may be extracted in
the same order as they occur in the input image by combining the
integer constants with the 'binary or' (``|``) operator.  For
example::

    wcs.sub([WCSSUB_LONGITUDE | WCSSUB_LATITUDE | WCSSUB_SPECTRAL])

would extract the longitude, latitude, and spectral axes in the same
order as the input image.  If one of each were present, the resulting
object would have three dimensions.

For convenience, ``WCSSUB_CELESTIAL`` is defined as the combination
``WCSSUB_LONGITUDE | WCSSUB_LATITUDE | WCSSUB_CUBEFACE``.

The codes may also be negated to extract all but the types specified,
for example::

    wcs.sub([
      WCSSUB_LONGITUDE,
      WCSSUB_LATITUDE,
      WCSSUB_CUBEFACE,
      -(WCSSUB_SPECTRAL | WCSSUB_STOKES)])

The last of these specifies all axis types other than spectral or
Stokes.  Extraction is done in the order specified by ``axes``, i.e. a
longitude axis (if present) would be extracted first (via ``axes[0]``)
and not subsequently (via ``axes[3]``).  Likewise for the latitude and
cubeface axes in this example.

The number of dimensions in the returned object may be less than or
greater than the length of ``axes``.  However, it will never exceed the
number of axes in the input image.
"""

tab = """
``list of Tabprm`` Tabular coordinate objects.

A list of tabular coordinate objects associated with this WCS.
"""

Tabprm = """
A class to store the information related to tabular coordinates,
i.e., coordinates that are defined via a lookup table.

This class can not be constructed directly from Python, but instead is
returned from `~astropy.wcs.Wcsprm.tab`.
"""

theta0 = """
``double``  The native longitude of the fiducial point.

The point whose celestial coordinates are given in ``ref[1:2]``.  If
undefined (NaN) the initialization routine, `~astropy.wcs.Wcsprm.set`,
will set this to a projection-specific default.

See also
--------
astropy.wcs.Wcsprm.phi0
"""

to_header = """
to_header(relax=False)

`to_header` translates a WCS object into a FITS header.

The details of the header depends on context:

    - If the `~astropy.wcs.Wcsprm.colnum` member is non-zero then a
      binary table image array header will be produced.

    - Otherwise, if the `~astropy.wcs.Wcsprm.colax` member is set
      non-zero then a pixel list header will be produced.

    - Otherwise, a primary image or image extension header will be
      produced.

The output header will almost certainly differ from the input in a
number of respects:

    1. The output header only contains WCS-related keywords.  In
       particular, it does not contain syntactically-required keywords
       such as ``SIMPLE``, ``NAXIS``, ``BITPIX``, or ``END``.

    2. Deprecated (e.g. ``CROTAn``) or non-standard usage will be
       translated to standard (this is partially dependent on whether
       ``fix`` was applied).

    3. Quantities will be converted to the units used internally,
       basically SI with the addition of degrees.

    4. Floating-point quantities may be given to a different decimal
       precision.

    5. Elements of the ``PCi_j`` matrix will be written if and only if
       they differ from the unit matrix.  Thus, if the matrix is unity
       then no elements will be written.

    6. Additional keywords such as ``WCSAXES``, ``CUNITia``,
       ``LONPOLEa`` and ``LATPOLEa`` may appear.

    7. The original keycomments will be lost, although
       `~astropy.wcs.Wcsprm.to_header` tries hard to write meaningful
       comments.

    8. Keyword order may be changed.

Keywords can be translated between the image array, binary table, and
pixel lists forms by manipulating the `~astropy.wcs.Wcsprm.colnum` or
`~astropy.wcs.Wcsprm.colax` members of the `~astropy.wcs.WCS`
object.

Parameters
----------

relax : bool or int
    Degree of permissiveness:

    - `False`: Recognize only FITS keywords defined by the published
      WCS standard.

    - `True`: Admit all recognized informal extensions of the WCS
      standard.

    - `int`: a bit field selecting specific extensions to write.
      See :ref:`relaxwrite` for details.

Returns
-------
header : str
    Raw FITS header as a string.
"""

ttype = """
``str`` (read-only)

``TTYPEn`` identifying the column of the binary table that contains
the wcstab array.
"""

unitfix = """
unitfix(translate_units='')

Translates non-standard ``CUNITia`` keyvalues.

For example, ``DEG`` -> ``deg``, also stripping off unnecessary
whitespace.

Parameters
----------
translate_units : str, optional
    Do potentially unsafe translations of non-standard unit strings.

    Although ``\"S\"`` is commonly used to represent seconds, its
    recognizes ``\"S\"`` formally as Siemens, however rarely that may
    be translation to ``\"s\"`` is potentially unsafe since the
    standard used.  The same applies to ``\"H\"`` for hours (Henry),
    and ``\"D\"`` for days (Debye).

    This string controls what to do in such cases, and is
    case-insensitive.

    - If the string contains ``\"s\"``, translate ``\"S\"`` to ``\"s\"``.

    - If the string contains ``\"h\"``, translate ``\"H\"`` to ``\"h\"``.

    - If the string contains ``\"d\"``, translate ``\"D\"`` to ``\"d\"``.

    Thus ``''`` doesn't do any unsafe translations, whereas ``'shd'``
    does all of them.

Returns
-------
success : int
    Returns ``0`` for success; ``-1`` if no change required.
"""

velangl = """
``double`` Velocity angle.

The angle in degrees that should be used to decompose an observed
velocity into radial and transverse components.

An undefined value is represented by NaN.
"""

velosys = """
``double`` Relative radial velocity.

The relative radial velocity (m/s) between the observer and the
selected standard of rest in the direction of the celestial reference
coordinate, ``VELOSYSa``.

An undefined value is represented by NaN.

See also
--------
astropy.wcs.Wcsprm.specsys, astropy.wcs.Wcsprm.ssysobs
"""

velref = """
``int`` AIPS velocity code.

From ``VELREF`` keyword.
"""

wcs = """
A `~astropy.wcs.Wcsprm` object to perform the basic `wcslib`_ WCS
transformation.
"""

Wcs = """
Wcs(*sip, cpdis, wcsprm, det2im*)

Wcs objects amalgamate basic WCS (as provided by `wcslib`_), with
`SIP`_ and `distortion paper`_ operations.

To perform all distortion corrections and WCS transformation, use
``all_pix2world``.

Parameters
----------
sip : `~astropy.wcs.Sip` object or `None`

cpdis : A pair of `~astropy.wcs.DistortionLookupTable` objects, or
  ``(None, None)``.

wcsprm : `~astropy.wcs.Wcsprm` object

det2im : A pair of `~astropy.wcs.DistortionLookupTable` objects, or
   ``(None, None)``.
"""

Wcsprm = """
Wcsprm(header=None, key=' ', relax=False, naxis=2, keysel=0, colsel=None)

`~astropy.wcs.Wcsprm` performs the core WCS transformations.

.. note::
    The members of this object correspond roughly to the key/value
    pairs in the FITS header.  However, they are adjusted and
    normalized in a number of ways that make performing the WCS
    transformation easier.  Therefore, they can not be relied upon to
    get the original values in the header.  For that, use
    `astropy.io.fits.Header` directly.

The FITS header parsing enforces correct FITS "keyword = value" syntax
with regard to the equals sign occurring in columns 9 and 10.
However, it does recognize free-format character (NOST 100-2.0,
Sect. 5.2.1), integer (Sect. 5.2.3), and floating-point values
(Sect. 5.2.4) for all keywords.

Parameters
----------
header : An `astropy.io.fits.Header`, string, or `None`.
  If ``None``, the object will be initialized to default values.

key : str, optional
    The key referring to a particular WCS transform in the header.
    This may be either ``' '`` or ``'A'``-``'Z'`` and corresponds to
    the ``\"a\"`` part of ``\"CTYPEia\"``.  (*key* may only be
    provided if *header* is also provided.)

relax : bool or int, optional

    Degree of permissiveness:

    - `False`: Recognize only FITS keywords defined by the published
      WCS standard.

    - `True`: Admit all recognized informal extensions of the WCS
      standard.

    - `int`: a bit field selecting specific extensions to accept.  See
      :ref:`relaxread` for details.

naxis : int, optional
    The number of world coordinates axes for the object.  (*naxis* may
    only be provided if *header* is `None`.)

keysel : sequence of flag bits, optional
    Vector of flag bits that may be used to restrict the keyword types
    considered:

        - ``WCSHDR_IMGHEAD``: Image header keywords.

        - ``WCSHDR_BIMGARR``: Binary table image array.

        - ``WCSHDR_PIXLIST``: Pixel list keywords.

    If zero, there is no restriction.  If -1, the underlying wcslib
    function ``wcspih()`` is called, rather than ``wcstbh()``.

colsel : sequence of int
    A sequence of table column numbers used to restrict the keywords
    considered.  `None` indicates no restriction.

Raises
------
MemoryError
     Memory allocation failed.

ValueError
     Invalid key.

KeyError
     Key not found in FITS header.
"""

Wtbarr = """
Classes to construct coordinate lookup tables from a binary table
extension (BINTABLE).

This class can not be constructed directly from Python, but instead is
returned from `~astropy.wcs.Wcsprm.wtb`.
"""

zsource = """
``double`` The redshift, ``ZSOURCEa``, of the source.

An undefined value is represented by NaN.
"""

WcsError = """
Base class of all invalid WCS errors.
"""

SingularMatrix = """
SingularMatrixError()

The linear transformation matrix is singular.
"""

InconsistentAxisTypes = """
InconsistentAxisTypesError()

The WCS header inconsistent or unrecognized coordinate axis type(s).
"""

InvalidTransform = """
InvalidTransformError()

The WCS transformation is invalid, or the transformation parameters
are invalid.
"""

InvalidCoordinate = """
InvalidCoordinateError()

One or more of the world coordinates is invalid.
"""

NoSolution = """
NoSolutionError()

No solution can be found in the given interval.
"""

InvalidSubimageSpecification = """
InvalidSubimageSpecificationError()

The subimage specification is invalid.
"""

NonseparableSubimageCoordinateSystem = """
NonseparableSubimageCoordinateSystemError()

Non-separable subimage coordinate system.
"""

NoWcsKeywordsFound = """
NoWcsKeywordsFoundError()

No WCS keywords were found in the given header.
"""

InvalidTabularParameters = """
InvalidTabularParametersError()

The given tabular parameters are invalid.
"""
