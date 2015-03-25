1.0.2 (unreleased)
------------------

New Features
^^^^^^^^^^^^

- ``astropy.config``

- ``astropy.constants``

- ``astropy.convolution``

- ``astropy.coordinates``

- ``astropy.cosmology``

- ``astropy.io.ascii``

- ``astropy.io.fits``

- ``astropy.io.misc``

- ``astropy.io.registry``

- ``astropy.io.votable``

- ``astropy.modeling``

  - Added support for polynomials with degree 0 or degree greater than 15.
    [#3574, 3589]

- ``astropy.nddata``

- ``astropy.stats``

- ``astropy.sphinx``

- ``astropy.table``

- ``astropy.time``

- ``astropy.units``

- ``astropy.utils``

- ``astropy.vo``

- ``astropy.wcs``

API Changes
^^^^^^^^^^^

- ``astropy.config``

- ``astropy.constants``

- ``astropy.convolution``

- ``astropy.coordinates``

- ``astropy.cosmology``

- ``astropy.io.ascii``

- ``astropy.io.fits``

- ``astropy.io.misc``

- ``astropy.io.registry``

- ``astropy.io.votable``

- ``astropy.modeling``

- ``astropy.nddata``

- ``astropy.stats``

- ``astropy.table``

- ``astropy.time``

- ``astropy.units``

- ``astropy.utils``

- ``astropy.vo``

- ``astropy.wcs``

Bug Fixes
^^^^^^^^^

- ``astropy.config``

- ``astropy.constants``

- ``astropy.convolution``

- ``astropy.coordinates``

- ``astropy.cosmology``

- ``astropy.io.ascii``

- ``astropy.io.fits``

  - Fixed handling of BINTABLE with TDIMn of size 1. [#3580]

- ``astropy.io.misc``

- ``astropy.io.registry``

- ``astropy.io.votable``

- ``astropy.modeling``

- ``astropy.nddata``

- ``astropy.stats``

- ``astropy.table``

  - Ensure ``QTable`` can be pickled [#3590]

- ``astropy.time``

  - Ensure a ``Column`` without units is treated as an ``array``, not as an
    dimensionless ``Quantity``. [#3648]

- ``astropy.units``

  - Ensure equivalencies that do more than just scale a ``Quantity`` are
    properly handled also in ``ufunc`` evaluations. [#2496, #3586]

  - The LaTeX representation of the Angstrom unit has changed from
    ``\overset{\circ}{A}`` to ``\mathring{A}``, which should have
    better support across regular LaTeX, MathJax and matplotlib (as of
    version 1.5) [#3617]

- ``astropy.utils``

- ``astropy.vo``

  - Using HTTPS/SSL for communication between SAMP hubs now works
    correctly on all supported versions of Python [#3613]

- ``astropy.wcs``

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Nothing changed yet.


1.0.1 (2015-03-06)
------------------

Bug Fixes
^^^^^^^^^

- ``astropy.constants``

  - Ensure constants can be turned into ``Quantity`` safely. [#3537, #3538]

- ``astropy.io.ascii``

  - Fix a segfault in the fast C parser when one of the column headers
    is empty [#3545].

  - Fixed support for reading inf and nan values with the fast reader in
    Windows.  Also fixed in the case of using ``use_fast_converter=True``
    with the fast reader. [#3525]

  - Fixed use of mmap in the fast reader on Windows. [#3525]

  - Fixed issue where commented header would treat comments defining the table
    (i.e. column headers) as purely information comments, leading to problems
    when trying to round-trip the table. [#3562]

- ``astropy.modeling``

  - Fixed propagation of parameter constraints ('fixed', 'bounds', 'tied')
    between compound models and their components.  There is may still be some
    difficulty defining 'tied' constraints properly for use with compound
    models, however. [#3481]

- ``astropy.nddata``

  - Restore several properties to the compatibility class ``NDDataArray`` that
    were inadvertently omitted [#3466].

- ``astropy.time``

  - Time objects now always evalutate to ``True``, except when empty. [#3530]

Miscellaneous
^^^^^^^^^^^^^

- ``astropy._erfa``

  - The ERFA wrappers are now written directly in the Python/C API
    rather than using Cython, for greater performance. [#3521]
- Miscellaneous

  - Improve import time of astropy [#3488].

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Updated bundled astropy-helpers version to v1.0.1 to address installation
  issues with some packages that depend on Astropy. [#3541]


1.0 (2015-02-18)
----------------

General
^^^^^^^

Astropy now requires a Numpy 1.6.0 or later.

New Features
^^^^^^^^^^^^

- ``astropy.analytic_functions``

  - The ``astropy.analytic_functions`` was added to contain analytic functions
    useful for astronomy [#3077].

- ``astropy.coordinates``

  - ``astropy.coordinates`` now has a full stack of frames allowing
    transformations from ICRS or other celestial systems down to Alt/Az
    coordinates. [#3217]

  - ``astropy.coordinates`` now has a ``get_sun`` function that gives
    the coordinates  of the Sun at a specified time. [#3217]

  - ``SkyCoord`` now has ``to_pixel`` and ``from_pixel`` methods that convert
    between celestial coordinates as ``SkyCoord`` objects and pixel coordinates
    given an ``astropy.wcs.WCS`` object. [#3002]

  - ``SkyCoord`` now has ``search_around_sky`` and ``search_around_3d``
    convenience methods that allow searching for all coordinates within
    a certain distance of another ``SkyCoord``. [#2953]

  - ``SkyCoord`` can now accept a frame instance for the ``frame=`` keyword
    argument. [#3063]

  - ``SkyCoord`` now has a ``guess_from_table`` method that can be used to
    quickly create ``SkyCoord`` objects from an ``astropy.table.Table``
    object. [#2951]

  - ``astropy.coordinates`` now has a ``Galactocentric`` frame, a coordinate
    frame centered on a (user specified) center of the Milky Way. [#2761, #3286]

  - ``SkyCoord`` now accepts more formats of the coordinate string when the
    representation has ``ra`` and ``dec`` attributes. [#2920]

  - ``SkyCoord`` can now accept lists of ``SkyCoord`` objects, frame objects,
    or representation objects and will combine them into a single object.
    [#3285]

  - Frames and ``SkyCoord`` instances now have a method ``is_equivalent_frame``
    that can be used to check that two frames are equivalent (ignoring the
    data).  [#3330]

  - The ``__repr__`` of coordinate objects now shows scalar coordinates in the
    same format as vector coordinates. [#3350, 3448]

- ``astropy.cosmology``

  - Added ``lookback_distance``, which is ``c * lookback_time``. [#3145]

  - Add baryonic matter density and dark matter only density parameters
    to cosmology objects [#2757].

  - Add a ``clone`` method to cosmology objects to allow copies
    of cosmological objects to be created with the specified variables
    modified [#2592].

  - Increase default numerical precision of ``z_at_value`` following
    the accurate by default, fast by explicit request model [#3074].

  - Cosmology functions that take a single (redshift) input now
    broadcast like numpy ufuncs.  So, passing an arbitrarily shaped
    array of inputs will produce an output of the same shape. [#3178, #3194]

- ``astropy.io.ascii``

  - Simplify the way new Reader classes are defined, allowing custom behavior
    entirely by overriding inherited class attributes instead of setting
    instance attributes in the Reader ``__init__`` method. [#2812]

  - There is now a faster C/Cython engine available for reading and writing
    simple ASCII formats like CSV. Both are enabled by default, and fast
    reading will fall back on an ordinary reader in case of a parsing
    failure. Their behavior can be altered with the parameter ``fast_reader``
    in ``read`` and ``fast_writer`` in ``write``. [#2716]

  - Make Latex/AASTex tables use unit attribute of Column for output. [#3064]

  - Store comment lines encountered during reading in metadata of the
    output table via ``meta['comment_lines']``. [#3222]

  - Write comment lines in Table metadata during output for all basic formats,
    IPAC, and fast writers. This functionality can be disabled with
    ``comment=False``. [#3255]

  - Add reader / writer for the Enhanced CSV format which stores table and
    column meta data, in particular data type and unit. [#2319]

  - Automatically use ``guess=False`` when reading if the file ``format`` is
    provided and the format parameters are uniquely specified.  This update
    also removes duplicate format guesses to improve performance. [#3418]

- ``astropy.io.fits``

  - The ``fitsdiff`` script ignores some things by default when comparing fits
    files (e.g. empty header lines). This adds a ``--exact`` option where
    nothing is ignored. [#2782, #3110]

  - The ``fitsheader`` script now takes a ``--keyword`` option to extract a
    specific keyword from the header of a FITS file, and a ``--table`` option
    to export headers into any of the data formats supported by
    ``astropy.table``. [#2555, #2588]

  - ``Section`` now supports all advanced indexing features ``ndarray`` does
    (slices with any steps, integer arrays, boolean arrays, None, Ellipsis).
    It also properly returns scalars when this is appropriate. [#3148]

- ``astropy.io.votable``

  - ``astropy.io.votable.parse`` now takes a ``datatype_mapping``
    keyword argument to map invalid datatype names to valid ones in
    order to support non-compliant files. [#2675]

- ``astropy.modeling``

  - Added the capability of creating new "compound" models by combining
    existing models using arithmetic operators.  See the "What's New in 1.0"
    page in the Astropy documentation for more details. [#3231]

  - A new ``custom_model`` decorator/factory function has been added for
    converting normal functions to ``Model`` classes that can work within
    the Astropy modeling framework.  This replaces the old ``custom_model_1d``
    function which is now deprecated.  The new function works the same as
    the old one but is less limited in the types of models it can be used to
    created.  [#1763]

  - The ``Model`` and ``Fitter`` classes have ``.registry`` attributes which
    provide sets of all loaded ``Model`` and ``Fitter`` classes (this is
    useful for building UIs for models and fitting). [#2725]

  - A dict-like ``meta`` member was added to ``Model``. it is to be used to
    store any optional information which is relevant to a project and is not
    in the standard ``Model`` class. [#2189]

  - Added ``Ellipse2D`` model. [#3124]

- ``astropy.nddata``

  - New array-related utility functions in ``astropy.nddata.utils`` for adding
    and removing arrays from other arrays with different sizes/shapes. [#3201]

  - New metaclass ``NDDataBase`` for enforcing the nddata interface in
    subclasses without restricting implementation of the data storage. [#2905]

  - New mixin classes ``NDSlicingMixin`` for slicing, ``NDArithmeticMixin``
    for arithmetic operations, and ``NDIOMixin`` for input/ouput in NDData. [#2905]

  - Added a decorator ``support_nddata`` that can be used to write functions
    that can either take separate arguments or NDData objects. [#2855]

- ``astropy.stats``

  - Added ``mad_std()`` function. [#3208]

  - Added ``gaussian_fwhm_to_sigma`` and ``gaussian_sigma_to_fwhm``
    constants. [#3208]

  - New function ``sigma_clipped_stats`` which can be used to quickly get
    common statistics for an array, using sigma clipping at the same time.
    [#3201]

- ``astropy.table``

  - Changed the internal implementation of the ``Table`` class changed so that
    it no longer uses numpy structured arrays as the core table data container.
    [#2790, #3179]

  - Tables can now be written to an html file that includes interactive
    browsing capabilities. To write out to this format, use
    ``Table.write('filename.html', format='jsviewer')``. [#2875]

  - A ``quantity`` property and ``to`` method were added to ``Table``
    columns that allow the column values to be easily converted to
    ``astropy.units.Quantity`` objects. [#2950]

  - Add ``unique`` convenience method to table. [#3185]

- ``astropy.tests``

  - Added a new Quantity-aware ``assert_quantity_allclose``. [#3273]

- ``astropy.time``

  - ``Time`` can now handle arbitrary array dimensions, with operations
    following standard numpy broadcasting rules. [#3138]

- ``astropy.units``

  - Support for VOUnit has been updated to be compliant with version
    1.0 of the standard. [#2901]

  - Added an ``insert`` method to insert values into a ``Quantity`` object.
    This is similar to the ``numpy.insert`` function. [#3049]

  - When viewed in IPython, ``Quantity`` objects with array values now render
    using LaTeX and scientific notation. [#2271]

  - Added ``units.quantity_input`` decorator to validate quantity inputs to a
    function for unit compatibility. [#3072]

  - Added ``units.astronomical_unit`` as a long form for ``units.au``. [#3303]

- ``astropy.utils``

  - Added a new decorator ``astropy.utils.wraps`` which acts as a replacement
    for the standard library's ``functools.wraps``, the only difference being
    that the decorated function also preserves the wrapped function's call
    signature. [#2849]

  - ``astropy.utils.compat.numpy`` has been revised such that it can include
    patched versions of routines from newer ``numpy`` versions.  The first
    addition is a version of ``broadcast_arrays`` that can be used with
    ``Quantity`` and other ``ndarray`` subclasses (using the ``subok=True``
    flag). [#2327]

  - Added ``astropy.utils.resolve_name`` which returns a member of a module
    or class given the fully qualified dotted name of that object as a
    string. [#3389]

  - Added ``astropy.utils.minversion`` which can be used to check minimum
    version requirements of Python modules (to test for specific features and/
    or bugs and the like). [#3389]

- ``astropy.visualization``

  - Created ``astropy.visualization`` module and added functionality relating
    to image normalization (i.e. stretching and scaling) as well as a new
    script ``fits2bitmap`` that can produce a bitmap image from a FITS file.
    [#3201]

  - Added dictionary ``astropy.visualization.mpl_style.astropy_mpl_style``
    which can be used to set a uniform plotstyle specifically for tutorials
    that is improved compared to matplotlib defaults. [#2719, #2787, #3200]

- ``astropy.wcs``

  - ``wcslib`` has been upgraded to version 4.25.  This brings a
    single new feature:

    - ``equinox`` and ``radesys`` will now be given default values
      conforming with the WCS specification if ``EQUINOXa`` and
      ``RADESYSa``, respectively, are not present in the header.

  - The minimum required version of ``wcslib`` is now 4.24. [#2503]

  - Added a new function ``wcs_to_celestial_frame`` that can be used to find
    the astropy.coordinates celestial frame corresponding to a particular WCS.
    [#2730]

  - ``astropy.wcs.WCS.compare`` now supports a ``tolerance`` keyword argument
    to allow for approximate comparison of floating-point values. [#2503]

  - added ``pixel_scale_matrix``, ``celestial``, ``is_celestial``, and
    ``has_celestial`` convenience attributes. Added
    ``proj_plane_pixel_scales``, ``proj_plane_pixel_area``, and
    ``non_celestial_pixel_scales`` utility functions for retrieving WCS pixel
    scale and area information [#2832, #3304]

  - Added two functions ``pixel_to_skycoord`` and
    ``skycoord_to_pixel`` that make it easy to convert between
    SkyCoord objects and pixel coordinates. [#2885]

  - ``all_world2pix`` now uses a much more sophisticated and complete
    algorithm to iteratively compute the inverse WCS transform. [#2816]

  - Add ability to use ``WCS`` object to define projections in Matplotlib,
    using the ``WCSAxes`` package. [#3183]

  - Added ``is_proj_plane_distorted`` for testing if pixels are
    distorted. [#3329]

- Misc

  - ``astropy._erfa`` was added as a new subpackage wrapping the functionality
    of the ERFA library in python.  This is primarily of use for other astropy
    subpackages, but the API may be made more public in the future. [#2992]


API Changes
^^^^^^^^^^^

- ``astropy.coordinates``

  - Subclasses of ``BaseCoordinateFrame`` which define a custom ``repr`` should
    be aware of the format expected in ``SkyCoord.__repr__()``, which changed in
    this release. [#2704, #2882]

  - The ``CartesianPoints`` class (deprecated in v0.4) has now been removed.
    [#2990]

  - The previous ``astropy.coordinates.builtin_frames`` module is now a
    subpackage.  Everything that was in the
    ``astropy.coordinates.builtin_frames`` module is still accessible from the
    new package, but the classes are now in separate modules.  This should have
    no direct impact at the user level. [#3120]

  - Support for passing a frame as a positional argument in the ``SkyCoord``
    class has now been deprecated, except in the case where a frame with data
    is passed as the sole positional argument. [#3152]

  - Improved ``__repr__`` of coordinate objects representing a single
    coordinate point for the sake of easier copy/pasting. [#3350]

- ``astropy.cosmology``

  - The functional interface to the cosmological routines as well as
    ``set_current`` and ``get_current`` (deprecated in v0.4) have now been
    removed. [#2990]

- ``astropy.io.ascii``

  - Added a new argument to ``htmldict`` in the HTML reader named
    ``parser``, which allows the user to specify which parser
    BeautifulSoup should use as a backend. [#2815]

  - Add ``FixedWidthTwoLine`` reader to guessing. This will allows to read
    tables that a copied from screen output like ``print my_table`` to be read
    automatically. Discussed in #3025 and #3099 [#3109]

- ``astropy.io.fits``

  - A new optional argument ``cache`` has been added to
    ``astropy.io.fits.open()``.  When opening a FITS file from a URL,
    ``cache`` is a boolean value specifying whether or not to save the
    file locally in Astropy's download cache (``True`` by default). [#3041]

- ``astropy.modeling``

  - Model classes should now specify ``inputs`` and ``outputs`` class
    attributes instead of the old ``n_inputs`` and ``n_outputs``.  These
    should be tuples providing human-readable *labels* for all inputs and
    outputs of the model.  The length of the tuple indicates the numbers
    of inputs and outputs.  See "What's New in Astropy 1.0" for more
    details. [#2835]

  - It is no longer necessary to include ``__init__`` or ``__call__``
    definitions in ``Model`` subclasses if all they do is wrap the
    super-method in order to provide a nice call signature to the docs.
    The ``inputs`` class attribute is now used to generate a nice call
    signature, so these methods should only be overridden by ``Model``
    subclasses in order to provide new functionality. [#2835]

  - Most models included in Astropy now have sensible default values for most
    or all of their parameters.  Call ``help(ModelClass)`` on any model to
    check what those defaults are.  Most of them time they should be
    overridden, but some of them are useful (for example spatial offsets are
    always set at the origin by default). Another rule of thumb is that, where
    possible, default parameters are set so that the model is a no-op, or
    close to it, by default. [#2932]

  - The ``Model.inverse`` method has been changed to a *property*, so that
    now accessing ``model.inverse`` on a model returns a new model that
    implements that model's inverse, and *calling* ``model.inverse(...)``` on
    some independent variable computes the value of the inverse (similar to what
    the old ``Model.invert()`` method was meant to do).  [#3024]

  - The ``Model.invert()`` method has been removed entirely (it was never
    implemented and there should not be any existing code that relies on it).
    [#3024]

  - ``custom_model_1d`` is deprecated in favor of the new ``custom_model``
    (see "New Features" above).  [#1763]

  - The ``Model.param_dim`` property (deprecated in v0.4) has now been removed.
    [#2990]

  - The ``Beta1D`` and ``Beta2D`` models have been renamed to ``Moffat1D`` and
    ``Moffat2D``. [#3029]

- ``astropy.nddata``

  - ``flags``, ``shape``, ``size``, ``dtype`` and ``ndim`` properties removed
    from ``astropy.nddata.NDData``. [#2905]

  - Arithmetic operations, uncertainty propagation, slicing and automatic
    conversion to a numpy array removed from ``astropy.nddata.NDData``. The
    class ``astropy.nddata.NDDataArray`` is functionally equivalent to the
    old ``NDData``.  [#2905]

- ``astropy.table``

  - The ``Column.units`` property (deprecated in v0.3) has now been removed.
    [#2990]

  - The ``Row.data`` and ``Table._data`` attributes have been deprecated
    related to the change in Table implementation.  They are replaced by
    ``Row.as_void()`` and ``Table.as_array()`` methods, respectively. [#2790]

  - The ``Table.create_mask`` method has been removed.  This undocumented
    method was a development orphan and would cause corruption of the
    table if called. [#2790]

  - The return type for integer item access to a Column (e.g. col[12] or
    t['a'][12]) is now always a numpy scalar, numpy ``ndarray``, or numpy
    ``MaskedArray``.  Previously if the column was multidimensional then a
    Column object would be returned. [#3095]

  - The representation of Table and Column objects has been changed to
    be formatted similar to the print output. [#3239]

- ``astropy.time``

  - The ``Time.val`` and ``Time.vals`` properties (deprecated in v0.3) and the
    ``Time.lon``, and ``Time.lat`` properties (deprecated in v0.4) have now
    been removed. [#2990]

  - Add ``decimalyear`` format that represents time as a decimal year. [#3265]

- ``astropy.units``

  - Support for VOUnit has been updated to be compliant with version
    1.0 of the standard. This means that some VOUnit strings that were
    rejected before are now acceptable. [#2901] Notably:

      - SI prefixes are supported on most units
      - Binary prefixes are supported on "bits" and "bytes"
      - Custom units can be defined "inline" by placing them between single
        quotes.

  - ``Unit.get_converter`` has been deprecated.  It is not strictly
    necessary for end users, and it was confusing due to lack of
    support for ``Quantity`` objects. [#3456]

- ``astropy.utils``

  - Some members of ``astropy.utils.misc`` were moved into new submodules.
    Specifically:

    - ``deprecated``, ``deprecated_attribute``, and ``lazyproperty`` ->
      ``astropy.utils.decorators``

    - ``find_current_module``, ``find_mod_objs`` ->
      ``astropy.utils.introspection``

    All of these functions can be imported directly from ``astropy.utils``
    which should be preferred over referencing individual submodules of
    ``astropy.utils``.  [#2857]

  - The ProgressBar.iterate class method (deprecated in v0.3) has now been
    removed. [#2990]

  - Updated ``astropy/utils/console.py`` ProgressBar() module to
    display output to IPython notebook with the addition of an
    ``interactive`` kwarg. [#2658] [#2789]

- ``astropy.wcs``

  - The ``WCS.calcFootprint`` method (deprecated in v0.4) has now been removed.
    [#2990]

  - An invalid unit in a ``CUNITn`` keyword now displays a warning and
    returns a ``UnrecognizedUnit`` instance rather than raising an
    exception [#3190]

Bug Fixes
^^^^^^^^^

- ``astropy.convolution``

  - ``astropy.convolution.discretize_model`` now handles arbitrary callables
    correctly [#2274].

- ``astropy.coordinates``

  - ``Angle.to_string`` now outputs unicode arrays instead of object arrays.
    [#2981]

  - ``SkyCoord.to_string`` no longer gives an error when used with an array
    coordinate with more than one dimension. [#3340]

  - Fixed support for subclasses of ``UnitSphericalRepresentation`` and
    ``SphericalRepresentation`` [#3354, #3366]

  - Fixed latex display of array angles in IPython notebook. [#3480]

- ``astropy.io.ascii``

  - In the ``CommentedHeader`` the ``data_start`` parameter now defaults to
    ``0``, which is the first uncommented line. Discussed in #2692. [#3054]

  - Position lines in ``FixedWidthTwoLine`` reader could consist of many characters.
    Now, only one character in addition to the delimiter is allowed. This bug was
    discovered as part of [#3109]

  - The IPAC table writer now consistently uses the ``fill_values`` keyword to
    specify the output null values.  Previously the behavior was inconsistent
    or incorrect. [#3259]

  - The IPAC table reader now correctly interprets abbreviated column types.
    [#3279]

  - Tables that look almost, but not quite like DAOPhot tables could cause
    guessing to fail. [#3342]

- ``astropy.io.fits``

  - Fixed the problem in ``fits.open`` of some filenames with colon (``:``) in
    the name being recognized as URLs instead of file names. [#3122]

  - Setting ``memmap=True`` in ``fits.open`` and related functions now raises
    a ValueError if opening a file in memory-mapped mode is impossible. [#2298]

  - CONTINUE cards no longer end the value of the final card in the series with
    an ampersand, per the specification of the CONTINUE card convention. [#3282]

  - Fixed a crash that occurred when reading an ASCII table containing
    zero-precision floating point fields. [#3422]

  - When a float field for an ASCII table has zero-precision a decimal point
    (with no digits following it) is still written to the field as long as
    there is space for it, as recommended by the FITS standard.  This makes it
    less ambiguous that these columns should be interpreted as floats. [#3422]

- ``astropy.logger``

  - Fix a bug that occurred when displaying warnings that produced an error
    message ``dictionary changed size during iteration``. [#3353]

- ``astropy.modeling``

  - Fixed a bug in ``SLSQPLSQFitter`` where the ``maxiter`` argument was not
    passed correctly to the optimizer. [#3339]

- ``astropy.table``

  - Fix a problem where ``table.hstack`` fails to stack multiple references to
    the same table, e.g. ``table.hstack([t, t])``. [#2995]

  - Fixed a problem where ``table.vstack`` and ``table.hstack`` failed to stack
    a single table, e.g. ``table.vstack([t])``. [#3313]

  - Fix a problem when doing nested iterators on a single table. [#3358]

  - Fix an error when an empty list, tuple, or ndarray is used for item access
    within a table.  This now returns the table with no rows. [#3442]

- ``astropy.time``

  - When creating a Time object from a datetime object the time zone
    info is now correctly used. [#3160]

  - For Time objects, it is now checked that numerical input is finite. [#3396]

- ``astropy.units``

  - Added a ``latex_inline`` unit format that returns the units in LaTeX math
    notation with negative exponents instead of fractions [#2622].

  - When using a unit that is deprecated in a given unit format,
    non-deprecated alternatives will be suggested. [#2806] For
    example::

      >>> import astropy.units as u
      >>> u.Unit('Angstrom', format='fits')
      WARNING: UnitsWarning: The unit 'Angstrom' has been deprecated
      in the FITS standard. Suggested: nm (with data multiplied by
      0.1).  [astropy.units.format.utils]

- ``astropy.utils``

  - ``treat_deprecations_as_exceptions`` has been fixed to recognize Astropy
    deprecation warnings. [#3015]

- ``astropy.wcs``

  - ``astropy.wcs.WCS.sub`` now accepts unicode strings as input on
    Python 2.x [#3356]

- Misc

  - Some modules and tests that would crash upon import when using a non-final
    release of Numpy (e.g. 1.9.0rc1). [#3471]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- The bundled copy of astropy-helpers has been updated to v1.0. [#3515]
- The bundled copy of astropy-helpers has been updated to v1.0. [#3513]

- Updated ``astropy.extern.configobj`` to Version 5. Version 5 uses ``six``
  and the same code covers both Python 2 and Python 3. [#3149]

- ``astropy.coordinates``

  - The ``repr`` of ``SkyCoord`` and coordinate frame classes now separate
    frame attributes and coordinate information.  [#2704, #2882]

- ``astropy.io.fits``

  - Overwriting an existing file using the ``clobber=True`` option no longer
    displays a warning message. [#1963]

  - ``fits.open`` no longer catches ``OSError`` exceptions on missing or
    unreadable files-- instead it raises the standard Python exceptions in such
    cases. [#2756, #2785]

- ``astropy.table``

  - Sped up setting of ``Column`` slices by an order of magnitude. [#2994, #3020]

- Updated the bundled ``six`` module to version 1.7.3 and made 1.7.3 the
  minimum acceptable version of ``six``. [#2814]

- The version of ERFA included with Astropy is now v1.1.1 [#2971]

- The code base is now fully Python 2 and 3 compatible and no longer requires
  2to3. [#2033]

- `funcsigs <https://pypi.python.org/pypi/funcsigs>`_ is included in
  utils.compat, but defaults to the inspect module components where available
  (3.3+) [#3151].

- The list of modules displayed in the pytest header can now be customized.
  [#3157]

- `jinja2 <http://jinja.pocoo.org/docs/dev/>`_>=2.7 is now required to build the
  source code from the git repository, in order to allow the ERFA wrappers to
  be generated. [#3166]


0.4.4 (2015-01-21)
------------------

Bug Fixes
^^^^^^^^^

- ``astropy.vo.samp``

  - ``astropy.vo.samp`` is now usable on Python builds that do not
    support the SSLv3 protocol (which depends both on the version of
    Python and the version of OpenSSL or LibreSSL that it is built
    against.) [#3308]

API Changes
^^^^^^^^^^^

- ``astropy.vo.samp``

  - The default SSL protocol used is now determined from the default
    used in the Python ``ssl`` standard library.  This default may be
    different depending on the exact version of Python you are using.
    [#3308]

- ``astropy.wcs``

  - WCS allows slices of the form slice(None, x, y), which previously resulted
    in an unsliced copy being returned (note: this was previously incorrectly
    reported as fixed in v0.4.3) [#2909]


0.4.3 (2015-01-15)
------------------

Bug Fixes
^^^^^^^^^

- ``astropy.coordinates``

  - The ``Distance`` class has been fixed to no longer rely on the deprecated
    cosmology functions. [#2991]

  - Ensure ``float32`` values can be used in coordinate representations. [#2983]

  - Fix frame attribute inheritance in ``SkyCoord.transform_to()`` method so
    that the default attribute value (e.g. equinox) for the destination frame
    gets used if no corresponding value was explicitly specified. [#3106]

  - ``Angle`` accepts hours:mins or deg:mins initializers (without
     seconds). In these cases float minutes are also accepted. [#2843]

  - ``astropy.coordinates.SkyCoord`` objects are now copyable. [#2888]

  - ``astropy.coordinates.SkyCoord`` object attributes are now
    immutable.  It is still technically possible to change the
    internal data for an array-valued coordinate object but this leads
    to inconsistencies [#2889] and should not be done. [#2888]

- ``astropy.cosmology``

  - The ``ztol`` keyword argument to z_at_value now works correctly [#2993].

- ``astropy.io.ascii``

  - Fix a bug in Python 3 when guessing file format using a file object as
    input.  Also improve performance in same situation for Python 2. [#3132]

  - Fix a problem where URL was being downloaded for each guess. [#2001]

- ``astropy.io.fits``

  - The ``in`` operator now works correctly for checking if an extension
    is in an ``HDUList`` (as given via EXTNAME, (EXTNAME, EXTVER) tuples,
    etc.) [#3060]

  - Added workaround for bug in MacOS X <= 10.8 that caused np.fromfile to
    fail. [#3078]

  - Added support for the ``RICE_ONE`` compression type synonym. [#3115]

- ``astropy.modeling``

  - Fixed a test failure on Debian/PowerPC and Debian/s390x. [#2708]

  - Fixed crash in evaluating models that have more outputs than inputs--this
    case may not be handled as desired for all conceivable models of this
    format (some may have to implement custom ``prepare_inputs`` and
    ``prepare_outputs`` methods).  But as long as all outputs can be assumed
    to have a shape determined from the broadcast of all inputs with all
    parameters then this can be used safely. [#3250]

- ``astropy.table``

  - Fix a bug that caused join to fail for multi-dimensional columns. [#2984]

  - Fix a bug where MaskedColumn attributes which had been changed since
    the object was created were not being carried through when slicing. [#3023]

  - Fix a bug that prevented initializing a table from a structured array
    with multi-dimensional columns with copy=True. [#3034]

  - Fixed unnecessarily large unicode columns when instantiating a table from
    row data on Python 3. [#3052]

  - Improved the warning message when unable to aggregate non-numeric
    columns. [#2700]

- ``astropy.units``

  - Operations on quantities with incompatible types now raises a much
    more informative ``TypeError``. [#2934]

  - ``Quantity.tolist`` now overrides the ``ndarray`` method to give a
    ``NotImplementedError`` (by renaming the previous ``list`` method). [#3050]

  - ``Quantity.round`` now always returns a ``Quantity`` (previously it
    returned an ``ndarray`` for ``decimals>0``). [#3062]

  - Ensured ``np.squeeze`` always returns a ``Quantity`` (it only worked if
    no dimensions were removed). [#3045]

  - Input to ``Quantity`` with a ``unit`` attribute no longer can get mangled
    with ``copy=False``. [#3051]

  - Remove trailing space in ``__format__`` calls for dimensionless quantities.
    [#3097]

  - Comparisons between units and non-unit-like objects now works
    correctly. [#3108]

  - Units with fractional powers are now correctly multiplied together
    by using rational arithmetic.  [#3121]

  - Removed a few entries from spectral density equivalencies which did not
    make sense. [#3153]

- ``astropy.utils``

  - Fixed an issue with the ``deprecated`` decorator on classes that invoke
    ``super()`` in their ``__init__`` method. [#3004]

  - Fixed a bug which caused the ``metadata_conflicts`` parameter to be
    ignored in the ``astropy.utils.metadata.merge`` function. [#3294]

- ``astropy.vo``

  - Fixed an issue with reconnecting to a SAMP Hub. [#2674]

- ``astropy.wcs``

  - Invalid or out of range values passed to ``wcs_world2pix`` will
    now be correctly identified and returned as ``nan``
    values. [#2965]

  - Fixed an issue which meant that Python thought ``WCS`` objects were
    iterable. [#3066]

- Misc

  - Astropy will now work if your Python interpreter does not have the
    ``bz2`` module installed. [#3104]

  - Fixed ``ResourceWarning`` for ``astropy/extern/bundled/six.py`` that could
    occur sometimes after using Astropy in Python 3.4. [#3156]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``astropy.coordinates``

  - Improved the agreement of the FK5 <-> Galactic conversion with other
    codes, and with the FK5 <-> FK4 <-> Galactic route. [#3107]


0.4.2 (2014-09-23)
------------------

Bug Fixes
^^^^^^^^^

- ``astropy.coordinates``

  - ``Angle`` accepts hours:mins or deg:mins initializers (without
     seconds). In these cases float minutes are also accepted.

  - The ``repr`` for coordinate frames now displays the frame attributes
    (ex: ra, dec) in a consistent order.  It should be noted that as part of
    this fix, the ``BaseCoordinateFrame.get_frame_attr_names()`` method now
    returns an ``OrderedDict`` instead of just a ``dict``. [#2845]

- ``astropy.io.fits``

  - Fixed a crash when reading scaled float data out of a FITS file that was
    loaded from a string (using ``HDUList.fromfile``) rather than from a file.
    [#2710]

  - Fixed a crash when reading data from an HDU whose header contained in
    invalid value for the BLANK keyword (e.g., a string value instead of an
    integer as required by the FITS Standard). Invalid BLANK keywords are now
    warned about, but are otherwise ignored. [#2711]

  - Fixed a crash when reading the header of a tile-compressed HDU if that
    header contained invalid duplicate keywords resulting in a ``KeyError``
    [#2750]

  - Fixed crash when reading gzip-compressed FITS tables through the Astropy
    ``Table`` interface. [#2783]

  - Fixed corruption when writing new FITS files through to gzipped files.
    [#2794]

  - Fixed crash when writing HDUs made with non-contiguous data arrays to
    file-like objects. [#2794]

  - It is now possible to create ``astropy.io.fits.BinTableHDU``
    objects with a table with zero rows. [#2916]

- ``astropy.io.misc``

  - Fixed a bug that prevented h5py ``Dataset`` objects from being
    automatically recognized by ``Table.read``. [#2831]

- ``astropy.modeling``

  - Make ``LevMarLSQFitter`` work with ``weights`` keyword. [#2900]

- ``astropy.table``

  - Fixed reference cycle in tables that could prevent ``Table`` objects
    from being freed from memory. [#2879]

  - Fixed an issue where ``Table.pprint()`` did not print the header to
    ``stdout`` when ``stdout`` is redirected (say, to a file). [#2878]

  - Fixed printing of masked values when a format is specified. [#1026]

  - Ensured that numpy ufuncs that return booleans return plain ``ndarray``
    instances, just like the comparison operators. [#2963]

- ``astropy.time``

  - Ensure bigendian input to Time works on a little-endian machine
    (and vice versa).  [#2942]

- ``astropy.units``

  - Ensure unit is kept when adding 0 to quantities. [#2968]

- ``astropy.utils``

  - Fixed color printing on Windows with IPython 2.0. [#2878]

- ``astropy.vo``

  - Improved error message on Cone Search time out. [#2687]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed a couple issues with files being inappropriately included and/or
  excluded from the source archive distributions of Astropy. [#2843, #2854]

- As part of fixing the fact that masked elements of table columns could not be
  printed when a format was specified, the column format string options were
  expanded to allow simple specifiers such as ``'5.2f'``. [#2898]

- Ensure numpy 1.9 is supported. [#2917]

- Ensure numpy master is supported, by making ``np.cbrt`` work with quantities.
  [#2937]

0.4.1 (2014-08-08)
------------------

Bug Fixes
^^^^^^^^^

- ``astropy.config``

  - Fixed a bug where an unedited configuration file from astropy
    0.3.2 would not be correctly identified as unedited. [#2772] This
    resulted in the warning::

      WARNING: ConfigurationChangedWarning: The configuration options
      in astropy 0.4 may have changed, your configuration file was not
      updated in order to preserve local changes.  A new configuration
      template has been saved to
      '~/.astropy/config/astropy.0.4.cfg'. [astropy.config.configuration]

  - Fixed the error message that is displayed when an old
    configuration item has moved.  Before, the destination
    section was wrong.  [#2772]

  - Added configuration settings for ``io.fits``, ``io.votable`` and
    ``table.jsviewer`` that were missing from the configuration file
    template. [#2772]

  - The configuration template is no longer rewritten on every import
    of astropy, causing race conditions. [#2805]

- ``astropy.convolution``

  - Fixed the multiplication of ``Kernel`` with numpy floats. [#2174]

- ``astropy.coordinates``

  - ``Distance`` can now take a list of quantities. [#2261]

  - For in-place operations for ``Angle`` instances in which the result unit
    is not an angle, an exception is raised before the instance is corrupted.
    [#2718]

  - ``CartesianPoints`` are now deprecated in favor of
    ``CartesianRepresentation``. [#2727]

- ``astropy.io.misc``

  - An existing table within an HDF5 file can be overwritten without affecting
    other datasets in the same HDF5 file by simultaneously using
    ``overwrite=True`` and ``append=True`` arguments to the ``Table.write``
    method. [#2624]

- ``astropy.logger``

  - Fixed a crash that could occur in rare cases when (such as in bundled
    apps) where submodules of the ``email`` package are not importable. [#2671]

- ``astropy.nddata``

  - ``astropy.nddata.NDData()`` no longer raises a ``ValueError`` when passed
    a numpy masked array which has no masked entries. [#2784]

- ``astropy.table``

  - When saving a table to a FITS file containing a unit that is not
    supported by the FITS standard, a warning rather than an exception
    is raised. [#2797]

- ``astropy.units``

  - By default, ``Quantity`` and its subclasses will now convert to float also
    numerical types such as ``decimal.Decimal``, which are stored as objects
    by numpy. [#1419]

  - The units ``count``, ``pixel``, ``voxel`` and ``dbyte`` now output
    to FITS, OGIP and VOUnit formats correctly. [#2798]

- ``astropy.utils``

  - Restored missing information from deprecation warning messages
    from the ``deprecated`` decorator. [#2811]

  - Fixed support for ``staticmethod`` deprecation in the ``deprecated``
    decorator. [#2811]

- ``astropy.wcs``

  - Fixed a memory leak when ``astropy.wcs.WCS`` objects are copied
    [#2754]

  - Fixed a crash when passing ``ra_dec_order=True`` to any of the
    ``*2world`` methods. [#2791]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Bundled copy of astropy-helpers upgraded to v0.4.1. [#2825]

- General improvements to documentation and docstrings [#2722, #2728, #2742]

- Made it easier for third-party packagers to have Astropy use their own
  version of the ``six`` module (so long as it meets the minimum version
  requirement) and remove the copy bundled with Astropy.  See the
  astropy/extern/README file in the source tree.  [#2623]


0.4 (2014-07-16)
----------------

New Features
^^^^^^^^^^^^

- ``astropy.constants``

  - Added ``b_wien`` to represent Wien wavelength displacement law constant.
    [#2194]

- ``astropy.convolution``

  - Changed the input parameter in ``Gaussian1DKernel`` and
    ``Gaussian2DKernel`` from ``width`` to ``stddev`` [#2085].

- ``astropy.coordinates``

  - The coordinates package has undergone major changes to implement
    `APE5 <https://github.com/astropy/astropy-APEs/blob/master/APE5.rst>`_ .
    These include backwards-incompatible changes, as the underlying framework
    has changed substantially. See the APE5 text and the package documentation
    for more details. [#2422]

  - A ``position_angle`` method has been added to the new ``SkyCoord``. [#2487]

  - Updated ``Angle.dms`` and ``Angle.hms`` to return ``namedtuple`` -s instead
    of regular tuples, and added ``Angle.signed_dms`` attribute that gives the
    absolute value of the ``d``, ``m``, and ``s`` along with the sign.  [#1988]

  - By default, ``Distance`` objects are now required to be positive. To
    allow negative values, set ``allow_negative=True`` in the ``Distance``
    constructor when creating a ``Distance`` instance.

  - ``Longitude`` (resp. ``Latitude``) objects cannot be used any more to
    initialize or set ``Latitude`` (resp. ``Longitude``) objects. An explicit
    conversion to ``Angle`` is now required. [#2461]

  - The deprecated functions for pre-0.3 coordinate object names like
    ``ICRSCoordinates`` have been removed. [#2422]

  - The ``rotation_matrix`` and ``angle_axis`` functions in
    ``astropy.coordinates.angles`` were made more numerically consistent and
    are now tested explicitly [#2619]

- ``astropy.cosmology``

  - Added ``z_at_value`` function to find the redshift at which a cosmology
    function matches a desired value. [#1909]

  - Added ``FLRW.differential_comoving_volume`` method to give the differential
    comoving volume at redshift z. [#2103]

  - The functional interface is now deprecated in favor of the more-explicit
    use of methods on cosmology objects. [#2343]

  - Updated documentation to reflect the removal of the functional
    interface. [#2507]

- ``astropy.io.ascii``

  - The ``astropy.io.ascii`` output formats ``latex`` and ``aastex`` accept a
    dictionary called ``latex_dict`` to specify options for LaTeX output.  It is
    now possible to specify the table alignment within the text via the
    ``tablealign`` keyword. [#1838]

  - If ``header_start`` is specified in a call to ``ascii.get_reader`` or any
    method that calls ``get_reader`` (e.g. ``ascii.read``) but ``data_start``
    is not specified at the same time, then ``data_start`` is calculated so
    that the data starts after the header. Before this, the default was
    that the header line was read again as the first data line
    [#855 and #1844].

  - A new ``csv`` format was added as a convenience for handling CSV (comma-
    separated values) data. [#1935]
    This format also recognises rows with an inconsistent number of elements.
    [#1562]

  - An option was added to guess the start of data for CDS format files when
    they do not strictly conform to the format standard. [#2241]

  - Added an HTML reader and writer to the ``astropy.io.ascii`` package.
    Parsing requires the installation of BeautifulSoup and is therefore
    an optional feature. [#2160]

  - Added support for inputting column descriptions and column units
    with the ``io.ascii.SExtractor`` reader. [#2372]

  - Allow the use of non-local ReadMe files in the CDS reader. [#2329]

  - Provide a mechanism to select how masked values are printed. [#2424]

  - Added support for reading multi-aperture daophot file. [#2656]

- ``astropy.io.fits``

  - Included a new command-line script called ``fitsheader`` to display the
    header(s) of a FITS file from the command line. [#2092]

  - Added new verification options ``fix+ignore``, ``fix+warn``,
    ``fix+exception``, ``silentfix+ignore``, ``silentfix+warn``, and
    ``silentfix+exception`` which give more control over how to report fixable
    errors as opposed to unfixable errors.

- ``astropy.modeling``

  - Prototype implementation of fitters that treat optimization algorithms
    separately from fit statistics, allowing new fitters to be created by
    mixing and matching optimizers and statistic functions. [#1914]

  - Slight overhaul to how inputs to and outputs from models are handled with
    respect to array-valued parameters and variables, as well as sets of
    multiple models.  See the associated PR and the modeling section of the
    v0.4 documentation for more details. [#2634]

  - Added a new ``SimplexLSQFitter`` which uses a downhill simplex optimizer
    with a least squares statistic. [#1914]

  - Changed ``Gaussian2D`` model such that ``theta`` now increases
    counterclockwise. [#2199]

  - Replaced the ``MatrixRotation2D`` model with a new model called simply
    ``Rotation2D`` which requires only an angle to specify the rotation.
    The new ``Rotation2D`` rotates in a counter-clockwise sense whereas
    the old ``MatrixRotation2D`` increased the angle clockwise.
    [#2266, #2269]

  - Added a new ``AffineTransformation2D`` model which serves as a
    replacement for the capability of ``MatrixRotation2D`` to accept an
    arbitrary matrix, while also adding a translation capability. [#2269]

  - Added ``GaussianAbsorption1D`` model. [#2215]

  - New ``Redshift`` model [#2176].

- ``astropy.nddata``

  - Allow initialization ``NDData`` or ``StdDevUncertainty`` with a
    ``Quantity``. [#2380]

- ``astropy.stats``

  - Added flat prior to binom_conf_interval and binned_binom_proportion

  - Change default in ``sigma_clip`` from ``np.median`` to ``np.ma.median``.
    [#2582]

- ``astropy.sphinx``

  - Note, the following new features are included in astropy-helpers as well:

  - The ``automodapi`` and ``automodsumm`` extensions now include sphinx
    configuration options to write out what ``automodapi`` and ``automodsumm``
    generate, mainly for debugging purposes. [#1975, #2022]

  - Reference documentation now shows functions/class docstrings at the
    inteded user-facing API location rather than the actual file where
    the implementation is found. [#1826]

  - The ``automodsumm`` extension configuration was changed to generate
    documentation of class ``__call__`` member functions. [#1817, #2135]

  - ``automodapi`` and ``automodsumm`` now have an ``:allowed-package-names:``
    option that make it possible to document functions and classes that
    are in a different namespace.  [#2370]

- ``astropy.table``

  - Improved grouped table aggregation by using the numpy ``reduceat()`` method
    when possible. This can speed up the operation by a factor of at least 10
    to 100 for large unmasked tables and columns with relatively small
    group sizes.  [#2625]

  - Allow row-oriented data input using a new ``rows`` keyword argument.
    [#850]

  - Allow subclassing of ``Table`` and the component classes ``Row``, ``Column``,
    ``MaskedColumn``, ``TableColumns``, and ``TableFormatter``. [#2287]

  - Fix to allow numpy integer types as valid indices into tables in
    Python 3.x [#2477]

  - Remove transition code related to the order change in ``Column`` and
    ``MaskedColumn`` arguments ``name`` and ``data`` from Astropy 0.2
    to 0.3. [#2511]

  - Change HTML table representation in IPython notebook to show all
    table columns instead of restricting to 80 column width.  [#2651]

- ``astropy.time``

  - Mean and apparent sidereal time can now be calculated using the
    ``sidereal_time`` method [#1418].

  - The time scale now defaults to UTC if no scale is provided. [#2091]

  - ``TimeDelta`` objects can have all scales but UTC, as well as, for
    consistency with time-like quantities, undefined scale (where the
    scale is taken from the object one adds to or subtracts from).
    This allows, e.g., to work consistently in TDB.  [#1932]

  - ``Time`` now supports ISO format strings that end in "Z". [#2211, #2203]

- ``astropy.units``

  - Support for the unit format `Office of Guest Investigator Programs (OGIP)
    FITS files
    <http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/>`__
    has been added. [#377]

  - The ``spectral`` equivalency can now handle angular wave number. [#1306 and
    #1899]

  - Added ``one`` as a shorthand for ``dimensionless_unscaled``. [#1980]

  - Added ``dex`` and ``dB`` units. [#1628]

  - Added ``temperature()`` equivalencies to support conversion between
    Kelvin, Celsius, and Fahrenheit. [#2209]

  - Added ``temperature_energy()`` equivalencies to support conversion
    between electron-volt and Kelvin. [#2637]

  - The runtime of ``astropy.units.Unit.compose`` is greatly improved
    (by a factor of 2 in most cases) [#2544]

  - Added ``electron`` unit. [#2599]

- ``astropy.utils``

  - ``timer.RunTimePredictor`` now uses ``astropy.modeling`` in its
    ``do_fit()`` method. [#1896]

- ``astropy.vo``

  - A new sub-package, ``astropy.vo.samp``, is now available (this was
    previously the SAMPy package, which has been refactored for use in
    Astropy). [#1907]

  - Enhanced functionalities for ``VOSCatalog`` and ``VOSDatabase``. [#1206]

- ``astropy.wcs``

  - astropy now requires wcslib version 4.23.  The version of wcslib
    included with astropy has been updated to version 4.23.

  - Bounds checking is now performed on native spherical
    coordinates.  Any out-of-bounds values will be returned as
    ``NaN``, and marked in the ``stat`` array, if using the
    low-level ``wcslib`` interface such as
    ``astropy.wcs.Wcsprm.p2s``. [#2107]

  - A new method, ``astropy.wcs.WCS.compare()``, compares two wcsprm
    structs for equality with varying degrees of strictness. [#2361]

  - New ``astropy.wcs.utils`` module, with a handful of tools for manipulating
    WCS objects, including dropping, swapping, and adding axes.

- Misc

  - Includes the new astropy-helpers package which separates some of Astropy's
    build, installation, and documentation infrastructure out into an
    independent package, making it easier for Affiliated Packages to depend on
    these features.  astropy-helpers replaces/deprecates some of the submodules
    in the ``astropy`` package (see API Changes below).  See also
    `APE 4 <https://github.com/astropy/astropy-APEs/blob/master/APE4.rst>`_
    for more details on the motivation behind and implementation of
    astropy-helpers.  [#1563]


API Changes
^^^^^^^^^^^

- ``astropy.config``

  - The configuration system received a major overhaul, as part of APE3.  It is
    no longer possible to save configuration items from Python, but instead
    users must edit the configuration file directly.  The locations of
    configuration items have moved, and some have been changed to science state
    values.  The old locations should continue to work until astropy 0.5, but
    deprecation warnings will be displayed.  See the `Configuration transition
    <http://astropy.readthedocs.org/en/v0.4/config/config_0_4_transition.html>`_
    docs for a detailed description of the changes and how to update existing
    code. [#2094]

- ``astropy.io.fits``

  - The ``astropy.io.fits.new_table`` function is now fully deprecated (though
    will not be removed for a long time, considering how widely it is used).

    Instead please use the more explicit ``BinTableHDU.from_columns`` to create
    a new binary table HDU, and the similar ``TableHDU.from_columns`` to create
    a new ASCII table.  These otherwise accept the same arguments as
    ``new_table`` which is now just a wrapper for these.

  - The ``.fromstring`` classmethod of each HDU type has been simplified such
    that, true to its namesake, it only initializes an HDU from a string
    containing its header *and* data.

  - Fixed an issue where header wildcard matching (for example
    ``header['DATE*']``) can be used to match *any* characters that might
    appear in a keyword.  Previously this only matched keywords containing
    characters in the set ``[0-9A-Za-z_]``.  Now this can also match a hyphen
    ``-`` and any other characters, as some conventions like ``HIERARCH`` and
    record-valued keyword cards allow a wider range of valid characters than
    standard FITS keywords.

  - This will be the *last* release to support the following APIs that have
    been marked deprecated since Astropy v0.1/PyFITS v3.1:

    - The ``CardList`` class, which was part of the old header implementation.

    - The ``Card.key`` attribute.  Use ``Card.keyword`` instead.

    - The ``Card.cardimage`` and ``Card.ascardimage`` attributes.  Use simply
      ``Card.image`` or ``str(card)`` instead.

    - The ``create_card`` factory function.  Simply use the normal ``Card``
      constructor instead.

    - The ``create_card_from_string`` factory function.  Use ``Card.fromstring``
      instead.

    - The ``upper_key`` function.  Use ``Card.normalize_keyword`` method
      instead (this is not unlikely to be used outside of PyFITS itself, but it
      was technically public API).

    - The usage of ``Header.update`` with ``Header.update(keyword, value,
      comment)`` arguments.  ``Header.update`` should only be used analogously
      to ``dict.update``.  Use ``Header.set`` instead.

    - The ``Header.ascard`` attribute.  Use ``Header.cards`` instead for a list
      of all the ``Card`` objects in the header.

    - The ``Header.rename_key`` method.  Use ``Header.rename_keyword`` instead.

    - The ``Header.get_history`` method.  Use ``header['HISTORY']`` instead
      (normal keyword lookup).

    - The ``Header.get_comment`` method.  Use ``header['COMMENT']`` instead.

    - The ``Header.toTxtFile`` method.  Use ``header.totextfile`` instead.

    - The ``Header.fromTxtFile`` method.  Use ``Header.fromtextfile`` instead.

    - The ``tdump`` and ``tcreate`` functions.  Use ``tabledump`` and
      ``tableload`` respectively.

    - The ``BinTableHDU.tdump`` and ``tcreate`` methods.  Use
      ``BinTableHDU.dump`` and ``BinTableHDU.load`` respectively.

    - The ``txtfile`` argument to the ``Header`` constructor.  Use
      ``Header.fromfile`` instead.

    - The ``startColumn`` and ``endColumn`` arguments to the ``FITS_record``
      constructor.  These are unlikely to be used by any user code.

    These deprecated interfaces will be removed from the development version of
    Astropy following the v0.4 release (they will still be available in any
    v0.4.x bugfix releases, however).

- ``astropy.modeling``

  - The method computing the derivative of the model with respect
    to parameters was renamed from ``deriv`` to ``fit_deriv``. [#1739]

  - ``ParametricModel`` and the associated ``Parametric1DModel`` and
    ``Parametric2DModel`` classes have been renamed ``FittableModel``,
    ``Fittable1DModel``, and ``Fittable2DModel`` respectively.  The base
    ``Model`` class has subsumed the functionality of the old

    ``ParametricModel`` class so that all models support parameter constraints.
    The only distinction of ``FittableModel`` is that anything which subclasses
    it is assumed "safe" to use with Astropy fitters. [#2276]

  - ``NonLinearLSQFitter`` has been renamed ``LevMarLSQFitter`` to emphasise
    that it uses the Levenberg-Marquardt optimization algorithm with a
    least squares statistic function. [#1914]

  - The ``SLSQPFitter`` class has been renamed ``SLSQPLSQFitter`` to emphasize
    that it uses the Sequential Least Squares Programming optimization
    algorithm with a least squares statistic function. [#1914]

  - The ``Fitter.errorfunc`` method has been renamed to the more general
    ``Fitter.objective_function``. [#1914]

- ``astropy.nddata``

  - Issue warning if unit is changed from a non-trivial value by directly
    setting ``NDData.unit``. [#2411]

  - The ``mask`` and ``flag`` attributes of ``astropy.nddata.NDData`` can now
    be set with any array-like object instead of requiring that they be set
    with a ``numpy.ndarray``. [#2419]

- ``astropy.sphinx``

  - Use of the ``astropy.sphinx`` module is deprecated; all new development of
    this module is in ``astropy_helpers.sphinx`` which should be used instead
    (therefore documentation builds that made use of any of the utilities in
    ``astropy.sphinx`` now have ``astropy_helpers`` as a documentation
    dependency).

- ``astropy.table``

  - The default table printing function now shows a table header row for units
    if any columns have the unit attribute set.  [#1282]

  - Before, an unmasked ``Table`` was automatically converted to a masked
    table if generated from a masked Table or a ``MaskedColumn``.
    Now, this conversion is only done if explicitly requested or if any
    of the input values is actually masked. [#1185]

  - The repr() function of ``astropy.table.Table`` now shows the units
    if any columns have the unit attribute set.  [#2180]

  - The semantics of the config options ``table.max_lines`` and
    ``table.max_width`` has changed slightly.  If these values are not
    set in the config file, astropy will try to determine the size
    automatically from the terminal. [#2683]

- ``astropy.time``

  - Correct use of UT in TDB calculation [#1938, #1939].

  - ``TimeDelta`` objects can have scales other than TAI [#1932].

  - Location information should now be passed on via an ``EarthLocation``
    instance or anything that initialises it, e.g., a tuple containing
    either geocentric or geodetic coordinates. [#1928]

- ``astropy.units``

  - ``Quantity`` now converts input to float by default, as this is physically
    most sensible for nearly all units [#1776].

  - ``Quantity`` comparisons with ``==`` or ``!=`` now always return ``True``
    or ``False``, even if units do not match (for which case a ``UnitsError``
    used to be raised).  [#2328]

  - Applying ``float`` or ``int`` to a ``Quantity`` now works for all
    dimensionless quantities; they are automatically converted to unscaled
    dimensionless. [#2249]

  - The exception ``astropy.units.UnitException``, which was
    deprecated in astropy 0.2, has been removed.  Use
    ``astropy.units.UnitError`` instead [#2386]

  - Initializing a ``Quantity`` with a valid number/array with a ``unit``
    attribute now interprets that attribute as the units of the input value.
    This makes it possible to initialize a ``Quantity`` from an Astropy
    ``Table`` column and have it correctly pick up the units from the column.
    [#2486]

- ``astropy.wcs``

  - ``calcFootprint`` was deprecated. It is replaced by
    ``calc_footprint``.  An optional boolean keyword ``center`` was
    added to ``calc_footprint``.  It controls whether the centers or
    the corners of the pixels are used in the computation. [#2384]

  - ``astropy.wcs.WCS.sip_pix2foc`` and
    ``astropy.wcs.WCS.sip_foc2pix`` formerly did not conform to the
    ``SIP`` standard: ``CRPIX`` was added to the ``foc`` result so
    that it could be used as input to "core FITS WCS".  As of astropy
    0.4, ``CRPIX`` is no longer added to the result, so the ``foc``
    space is correct as defined in the `SIP convention
    <http://adsabs.harvard.edu/abs/2005ASPC..347..491S>`__. [#2360]

  - ``astropy.wcs.UnitConverter``, which was deprecated in astropy
    0.2, has been removed.  Use the ``astropy.units`` module
    instead. [#2386]

  - The following methods on ``astropy.wcs.WCS``, which were
    deprecated in astropy 0.1, have been removed [#2386]:

    - ``all_pix2sky`` -> ``all_pix2world``
    - ``wcs_pix2sky`` -> ``wcs_pix2world``
    - ``wcs_sky2pix`` -> ``wcs_world2pix``

  - The ``naxis1`` and ``naxis2`` attributes and the ``get_naxis``
    method of ``astropy.wcs.WCS``, which were deprecated in astropy
    0.2, have been removed.  Use the shape of the underlying FITS data
    array instead.  [#2386]

- Misc

  - The ``astropy.setup_helpers`` and ``astropy.version_helpers`` modules are
    deprecated; any non-critical fixes and development to those modules should
    be in ``astropy_helpers`` instead.  Packages that use these modules in
    their ``setup.py`` should depend on ``astropy_helpers`` following the same
    pattern as in the Astropy package template.


Bug Fixes
^^^^^^^^^

- ``astropy.constants``

  - ``astropy.constants.Contant`` objects can now be deep
    copied. [#2601]

- ``astropy.cosmology``

  - The distance modulus function in ``astropy.cosmology`` can now handle
    negative distances, which can occur in certain closed cosmologies. [#2008]

  - Removed accidental imports of some extraneous variables in
    ``astropy.cosmology`` [#2025]

- ``astropy.io.ascii``

  - ``astropy.io.ascii.read`` would fail to read lists of strings where some of
    the strings consisted of just a newline ("\n"). [#2648]

- ``astropy.io.fits``

  - Use NaN for missing values in FITS when using Table.write for float
    columns. Earlier the default fill value was close to 1e20.[#2186]

  - Fixes for checksums on 32-bit platforms.  Results may be different
    if writing or checking checksums in "nonstandard" mode.  [#2484]

  - Additional minor bug fixes ported from PyFITS.  [#2575]

- ``astropy.io.votable``

  - It is now possible to save an ``astropy.table.Table`` object as a
    VOTable with any of the supported data formats, ``tabledata``,
    ``binary`` and ``binary2``, by using the ``tabledata_format``
    kwarg. [#2138]

  - Fixed a crash writing out variable length arrays. [#2577]

- ``astropy.nddata``

  - Indexing ``NDData`` in a way that results in a single element returns that
    element. [#2170]

  - Change construction of result of arithmetic and unit conversion to allow
    subclasses to require the presence of attribute like unit. [#2300]

  - Scale uncertainties to correct units in arithmetic operations and unit
    conversion. [#2393]

  - Ensure uncertainty and mask members are copied in arithmetic and
    convert_unit_to. [#2394]

  - Mask result of arithmetic if either of the operands is masked. [#2403]

  - Copy all attributes of input object if ``astropy.nddata.NDData`` is
    initialized with an ``NDData`` object. [#2406]

  - Copy ``flags`` to new object in ``convert_unit_to``. [#2409]

  - Result of ``NDData`` arithmetic makes a copy of any WCS instead of using
    a reference. [#2410]

  - Fix unit handling for multiplication/division and use
    ``astropy.units.Quantity`` for units arithmetic. [#2413]

  - A masked ``NDData`` is now converted to a masked array when used in an
    operation or ufunc with a numpy array. [#2414]

  - An unmasked ``NDData`` now uses an internal representation of its mask
    state that ``numpy.ma`` expects so that an ``NDData`` behaves as an
    unmasked array. [#2417]

- ``astropy.sphinx``

  - Fix crash in smart resolver when the resolution doesn't work. [#2591]

- ``astropy.table``

  - The ``astropy.table.Column`` object can now use both functions and callable
    objects as formats. [#2313]

  - Fixed a problem on 64 bit windows that caused errors
    "expected 'DTYPE_t' but got 'long long'" [#2490]

  - Fix initialisation of ``TableColumns`` with lists or tuples.  [#2647]

  - Fix removal of single column using ``remove_columns``. [#2699]

  - Fix a problem that setting a row element within a masked table did not
    update the corresponding table element. [#2734]

- ``astropy.time``

  - Correct UT1->UTC->UT1 round-trip being off by 1 second if UT1 is
    on a leap second. [#2077]

- ``astropy.units``

  - ``Quantity.copy`` now behaves identically to ``ndarray.copy``, and thus
    supports the ``order`` argument (for numpy >=1.6). [#2284]

  - Composing base units into identical composite units now works. [#2382]

  - Creating and composing/decomposing units is now substantially faster [#2544]

  - ``Quantity`` objects now are able to be assigned NaN [#2695]

- ``astropy.wcs``

  - Astropy now requires wcslib version 4.23.  The version of wcslib
    included with astropy has been updated to version 4.23.

  - Bug fixes in the projection routines: in ``hpxx2s`` [the
    cartesian-to-spherical operation of the ``HPX`` projection]
    relating to bounds checking, bug introduced at wcslib 4.20; in
    ``parx2s`` and molx2s`` [the cartesion-to-spherical operation of
    the ``PAR`` and ``MOL`` projections respectively] relating to
    setting the stat vector; in ``hpxx2s`` relating to implementation
    of the vector API; and in ``xphx2s`` relating to setting an
    out-of-bounds value of *phi*.

  - In the ``PCO`` projection, use alternative projection equations
    for greater numerical precision near theta == 0.  In the ``COP``
    projection, return an exact result for theta at the poles.
    Relaxed the tolerance for bounds checking a little in ``SFL``
    projection.

  - Fix a bug allocating insufficient memory in
    ``astropy.wcs.WCS.sub`` [#2468]

  - A new method, ``Wcsprm.bounds_check`` (corresponding to wcslib's
    ``wcsbchk``) has been added to control what bounds checking is performed by
    wcslib.

  - ``WCS.to_header`` will now raise a more meaningful exception when the WCS
    information is invalid or inconsistent in some way. [#1854]

  - In ``WCS.to_header``, ``RESTFRQ`` and ``RESTWAV`` are no longer
    rewritten if zero. [#2468]

  - In ``WCS.to_header``, floating point values will now always be written
    with an exponent or fractional part, i.e. ``.0`` being appended if necessary
    to acheive this. [#2468]

  - If the C extension for ``astropy.wcs`` was not built or fails to import for
    any reason, ``import astropy.wcs`` will result in an ``ImportError``,
    rather than getting obscure errors once the ``astropy.wcs`` is used.
    [#2061]

  - When the C extension for ``astropy.wcs`` is built using a version of
    ``wscslib`` already present in the system, the package does not try
    to install ``wcslib`` headers under ``astropy/wcs/include``. [#2536]

  - Fixes an unresolved external symbol error in the
    `astropy.wcs._wcs` C extension on Microsoft Windows when built
    with a Microsoft compiler. [#2478]

- Misc

  - Running the test suite with ``python setup.py test`` now works if
    the path to the source contains spaces. [#2488]

  - The version of ERFA included with Astropy is now v1.1.0 [#2497]

  - Removed deprecated option from travis configuration and force use of
    wheels rather than allowing build from source. [#2576]

  - The short option ``-n`` to run tests in parallel was broken
    (conflicts with the distutils built-in option of "dry-run").
    Changed to ``-j``. [#2566]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``python setup.py test --coverage`` will now give more accurate
  results, because the coverage analysis will include early imports of
  astropy.  There doesn't seem to be a way to get this to work when
  doing ``import astropy; astropy.test()``, so the ``coverage``
  keyword to ``astropy.test`` has been removed.  Coverage testing now
  depends only on `coverage.py
  <http://nedbatchelder.com/code/coverage/>`__, not
  ``pytest-cov``. [#2112]

- The included version of py.test has been upgraded to 2.5.1. [#1970]

- The included version of six.py has been upgraded to 1.5.2. [#2006]

- Where appropriate, tests are now run both with and without the
  ``unicode_literals`` option to ensure that we support both cases. [#1962]

- Running the Astropy test suite from within the IPython REPL is disabled for
  now due to bad interaction between the test runner and IPython's logging
  and I/O handler.  For now, run the Astropy tests should be run in the basic
  Python interpreter. [#2684]

- Added support for numerical comparison of floating point values appearing in
  the output of doctests using a ``+FLOAT_CMP`` doctest flag. [#2087]

- A monkey patch is performed to fix a bug in Numpy version 1.7 and
  earlier where unicode fill values on masked arrays are not
  supported.  This may cause unintended side effects if your
  application also monkey patches ``numpy.ma`` or relies on the broken
  behavior.  If unicode support of masked arrays is important to your
  application, upgrade to Numpy 1.8 or later for best results. [#2059]

- The developer documentation has been extensively rearranged and
  rewritten. [#1712]

- The ``human_time`` function in ``astropy.utils`` now returns strings
  without zero padding. [#2420]

- The ``bdist_dmg`` command for ``setup.py`` has now been removed. [#2553]

- Many broken API links have been fixed in the documentation, and the
  ``nitpick`` Sphinx option is now used to avoid broken links in future.
  [#1221, #2019, #2109, #2161, #2162, #2192, #2200, #2296, #2448, #2456,
  #2460, #2467, #2476, #2508, #2509]


0.3.2 (2014-05-13)
------------------

Bug Fixes
^^^^^^^^^

- ``astropy.coordinates``

  - if ``sep`` argument is specified to be a single character in
    ``sexagisimal_to_string``, it now includes seperators only between
    items [#2183]

  - Ensure comparisons involving ``Distance`` objects do not raise exceptions;
    also ensure operations that lead to units other than length return
    ``Quantity``. [#2206, #2250]

  - Multiplication and division of ``Angle`` objects is now
    supported. [#2273]

  - Fixed ``Angle.to_string`` functionality so that negative angles have the
    correct amount of padding when ``pad=True``. [#2337]

  - Mixing strings and quantities in the ``Angle`` constructor now
    works.  For example: ``Angle(['1d', 1. * u.d])``.  [#2398]

  - If ``Longitude`` is given a ``Longitude`` as input, use its ``wrap_angle``
    by default [#2705]

- ``astropy.cosmology``

  - Fixed ``format()`` compatibility with Python 2.6. [#2129]

  - Be more careful about converting to floating point internally [#1815, #1818]

- ``astropy.io.ascii``

  - The CDS reader in ``astropy.io.ascii`` can now handle multiple
    description lines in ReadMe files. [#2225]

  - When reading a table with values that generate an overflow error during
    type conversion (e.g. overflowing the native C long type), fall through to
    using string. Previously this generated an exception [#2234].

  - Some CDS files mark missing values with ``"---"``, others with ``"--"``.
    Recognize any string with one to four dashes as null value. [#1335]

- ``astropy.io.fits``

  - Allow pickling of ``FITS_rec`` objects. [#1597]

  - Improved behavior when writing large compressed images on OSX by removing
    an unnecessary check for platform architecture. [#2345]

  - Fixed an issue where Astropy ``Table`` objects containing boolean columns
    were not correctly written out to FITS files. [#1953]

  - Several other bug fixes ported from PyFITS v3.2.3 [#2368]

  - Fixed a crash on Python 2.x when writing a FITS file directly to a
    ``StringIO.StringIO`` object. [#2463]

- ``astropy.io.registry``

  - Allow readers/writers with the same name to be attached to different
    classes. [#2312]

- ``astropy.io.votable``

  - By default, floating point values are now written out using
    ``repr`` rather than ``str`` to preserve precision [#2137]

- ``astropy.modeling``

  - Fixed the ``SIP`` and ``InverseSIP`` models both so that they work in the
    first place, and so that they return results consistent with the SIP
    functions in ``astropy.wcs``. [#2177]

- ``astropy.stats``

  - Ensure the ``axis`` keyword in ``astropy.stats.funcs`` can now be used for
    all axes. [#2173]

- ``astropy.table``

  - Ensure nameless columns can be printed, using 'None' for the header. [#2213]

- ``astropy.time``

  - Fixed pickling of ``Time`` objects. [#2123]

- ``astropy.units``

  - ``Quantity._repr_latex_()`` returns ``NotImplementedError`` for quantity
    arrays instead of an uninformative formatting exception. [#2258]

  - Ensure ``Quantity.flat`` always returns ``Quantity``. [#2251]

  - Angstrom unit renders better in MathJax [#2286]

- ``astropy.utils``

  - Progress bars will now be displayed inside the IPython
    qtconsole. [#2230]

  - ``data.download_file()`` now evaluates ``REMOTE_TIMEOUT()`` at runtime
    rather than import time. Previously, setting ``REMOTE_TIMEOUT`` after
    import had no effect on the function's behavior. [#2302]

  - Progressbar will be limited to 100% so that the bar does not exceed the
    terminal width.  The numerical display can still exceed 100%, however.

  - Converted representation of progress bar units without suffix
    from float to int in console.human_file_size. [#2201,#2202,#2721,#3299]

- ``astropy.vo``

  - Fixed ``format()`` compatibility with Python 2.6. [#2129]

  - Cone Search validation no longer raises ``ConeSearchError`` for positive RA.
    [#2240, #2242]

- ``astropy.wcs``

  - Fixed a bug where calling ``astropy.wcs.Wcsprm.sub`` with
    ``WCSSUB_CELESTIAL`` may cause memory corruption due to
    underallocation of a temporary buffer. [#2350]

  - Fixed a memory allocation bug in ``astropy.wcs.Wcsprm.sub`` and
    ``astropy.wcs.Wcsprm.copy``.  [#2439]

- Misc

  - Fixes for compatibility with Python 3.4. [#1945]

  - ``import astropy; astropy.test()`` now correctly uses the same test
    configuration as ``python setup.py test`` [#1811]


0.3.1 (2014-03-04)
------------------

Bug Fixes
^^^^^^^^^

- ``astropy.config``

  - Fixed a bug where ``ConfigurationItem.set_temp()`` does not reset to
    default value when exception is raised within ``with`` block. [#2117]

- ``astropy.convolution``

  - Fixed a bug where ``_truncation`` was left undefined for ``CustomKernel``.
    [#2016]

  - Fixed a bug with ``_normalization`` when ``CustomKernel`` input array
    sums to zero. [#2016]

- ``astropy.coordinates``

  - Fixed a bug where using ``==`` on two array coordinates wouldn't
    work. [#1832]

  - Fixed bug which caused ``len()`` not to work for coordinate objects and
    added a ``.shape`` property to get appropriately array-like behavior.
    [#1761, #2014]

  - Fixed a bug where sexagesimal notation would sometimes include
    exponential notation in the last field. [#1908, #1913]

  - ``CompositeStaticMatrixTransform`` no longer attempts to reference the
    undefined variable ``self.matrix`` during instantiation. [#1944]

  - Fixed pickling of ``Longitude``, ensuring ``wrap_angle`` is preserved
    [#1961]

  - Allow ``sep`` argument in ``Angle.to_string`` to be empty (resulting in no
    separators) [#1989]

- ``astropy.io.ascii``

  - Allow passing unicode delimiters when reading or writing tables.  The
    delimiter must be convertible to pure ASCII.  [#1949]

  - Fix a problem when reading a table and renaming the columns to names that
    already exist. [#1991]

- ``astropy.io.fits``

  - Ported all bug fixes from PyFITS 3.2.1.  See the PyFITS changelog at
    http://pyfits.readthedocs.org/en/v3.2.1/ [#2056]

- ``astropy.io.misc``

  - Fixed issues in the HDF5 Table reader/writer functions that occurred on
    Windows. [#2099]

- ``astropy.io.votable``

  - The ``write_null_values`` kwarg to ``VOTable.to_xml``, when set to `False`
    (the default) would produce non-standard VOTable files.  Therefore, this
    functionality has been replaced by a better understanding that knows which
    fields in a VOTable may be left empty (only ``char``, ``float`` and
    ``double`` in VOTable 1.1 and 1.2, and all fields in VOTable 1.3).  The
    kwarg is still accepted but it will be ignored, and a warning is emitted.
    [#1809]

  - Printing out a ``astropy.io.votable.tree.Table`` object using `repr` or
    `str` now uses the pretty formatting in ``astropy.table``, so it's possible
    to easily preview the contents of a ``VOTable``. [#1766]

- ``astropy.modeling``

  - Fixed bug in computation of model derivatives in ``LinearLSQFitter``.
    [#1903]

  - Raise a ``NotImplementedError`` when fitting composite models. [#1915]

  - Fixed bug in the computation of the ``Gaussian2D`` model. [#2038]

  - Fixed bug in the computation of the ``AiryDisk2D`` model. [#2093]

- ``astropy.sphinx``

  - Added slightly more useful debug info for AstropyAutosummary. [#2024]

- ``astropy.table``

  - The column string representation for n-dimensional cells with only
    one element has been fixed. [#1522]

  - Fix a problem that caused ``MaskedColumn.__getitem__`` to not preserve
    column metadata. [#1471, #1872]

  - With Numpy prior to version 1.6.2, tables with Unicode columns now
    sort correctly. [#1867]

  - ``astropy.table`` can now print out tables with Unicode columns containing
    non-ascii characters. [#1864]

  - Columns can now be named with Unicode strings, as long as they contain only
    ascii characters.  This makes using ``astropy.table`` easier on Python 2
    when ``from __future__ import unicode_literals`` is used. [#1864]

  - Allow pickling of ``Table``, ``Column``, and ``MaskedColumn`` objects. [#792]

  - Fix a problem where it was not possible to rename columns after sorting or
    adding a row. [#2039]

- ``astropy.time``

  - Fix a problem where scale conversion problem in TimeFromEpoch
    was not showing a useful error [#2046]

  - Fix a problem when converting to one of the formats ``unix``, ``cxcsec``,
    ``gps`` or ``plot_date`` when the time scale is ``UT1``, ``TDB`` or ``TCB``
    [#1732]

  - Ensure that ``delta_ut1_utc`` gets calculated when accessed directly,
    instead of failing and giving a rather obscure error message [#1925]

  - Fix a bug when computing the TDB to TT offset.  The transform routine was
    using meters instead of kilometers for the Earth vector.  [#1929]

  - Increase ``__array_priority__`` so that ``TimeDelta`` can convert itself
    to a ``Quantity`` also in reverse operations [#1940]

  - Correct hop list from TCG to TDB to ensure that conversion is
    possible [#2074]

- ``astropy.units``

  - ``Quantity`` initialisation rewritten for speed [#1775]

  - Fixed minor string formatting issue for dimensionless quantities. [#1772]

  - Fix error for inplace operations on non-contiguous quantities [#1834].

  - The definition of the unit ``bar`` has been corrected to "1e5
    Pascal" from "100 Pascal" [#1910]

  - For units that are close to known units, but not quite, for
    example due to differences in case, the exception will now include
    recommendations. [#1870]

  - The generic and FITS unit parsers now accept multiple slashes in
    the unit string.  There are multiple ways to interpret them, but
    the approach taken here is to convert "m/s/kg" to "m s-1 kg-1".
    Multiple slashes are accepted, but discouraged, by the FITS
    standard, due to the ambiguity of parsing, so a warning is raised
    when it is encountered. [#1911]

  - The use of "angstrom" (with a lower case "a") is now accepted in FITS unit
    strings, since it is in common usage.  However, since it is not officially
    part of the FITS standard, a warning will be issued when it is encountered.
    [#1911]

  - Pickling unrecognized units will not raise a ``AttributeError``. [#2047]

  - ``astropy.units`` now correctly preserves the precision of
    fractional powers. [#2070]

  - If a ``Unit`` or ``Quantity`` is raised to a floating point power
    that is very close to a rational number with a denominator less
    than or equal to 10, it is converted to a ``Fraction`` object to
    preserve its precision through complex unit conversion operations.
    [#2070]

- ``astropy.utils``

  - Fixed crash in ``timer.RunTimePredictor.do_fit``. [#1905]

  - Fixed ``astropy.utils.compat.argparse`` for Python 3.1. [#2017]

- ``astropy.wcs``

  - ``astropy.wcs.WCS``, ``astropy.wcs.WCS.fix`` and
    ``astropy.wcs.find_all_wcs`` now have a ``translate_units`` keyword
    argument that is passed down to ``astropy.wcs.Wcsprm.fix``.  This can be
    used to specify any unsafe translations of units from rarely used ones to
    more commonly used ones.

    Although ``"S"`` is commonly used to represent seconds, its translation to
    ``"s"`` is potentially unsafe since the standard recognizes ``"S"``
    formally as Siemens, however rarely that may be used.  The same applies to
    ``"H"`` for hours (Henry), and ``"D"`` for days (Debye).

    When these sorts of changes are performed, a warning is emitted.
    [#1854]

  - When a unit is "fixed" by ``astropy.wcs.WCS.fix`` or
    ``astropy.wcs.Wcsprm.unitfix``, it now correctly reports the ``CUNIT``
    field that was changed. [#1854]

  - ``astropy.wcs.Wcs.printwcs`` will no longer warn that ``cdelt`` is being
    ignored when none was present in the FITS file. [#1845]

  - ``astropy.wcs.Wcsprm.set`` is called from within the ``astropy.wcs.WCS``
    constructor, therefore any invalid information in the keywords will be
    raised from the constructor, rather than on a subsequent call to a
    transformation method. [#1918]

  - Fix a memory corruption bug when using ``astropy.wcs.Wcs.sub`` with
    ``astropy.wcs.WCSSUB_CELESTIAL``. [#1960]

  - Fixed the ``AttributeError`` exception that was raised when using
    ``astropy.wcs.WCS.footprint_to_file``. [#1912]

  - Fixed a ``NameError`` exception that was raised when using
    ``astropy.wcs.validate`` or the ``wcslint`` script. [#2053]

  - Fixed a bug where named WCSes may be erroneously reported as ``' '`` when
    using ``astropy.wcs.validate`` or the ``wcslint`` script. [#2053]

  - Fixed a bug where error messages about incorrect header keywords
    may not be propagated correctly, resulting in a "NULL error object
    in wcslib" message. [#2106]

- Misc

  - There are a number of improvements to make Astropy work better on big
    endian platforms, such as MIPS, PPC, s390x and SPARC. [#1849]

  - The test suite will now raise exceptions when a deprecated feature of
    Python or Numpy is used.  [#1948]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- A new function, ``astropy.wcs.get_include``, has been added to get the
  location of the ``astropy.wcs`` C header files. [#1755]

- The doctests in the ``.rst`` files in the ``docs`` folder are now
  tested along with the other unit tests.  This is in addition to the
  testing of doctests in docstrings that was already being performed.
  See ``docs/development/testguide.rst`` for more information. [#1771]

- Fix a problem where import fails on Python 3 if setup.py exists
  in current directory. [#1877]


0.3 (2013-11-20)
----------------

New Features
^^^^^^^^^^^^

- General

  - A top-level configuration item, ``unicode_output`` has been added to
    control whether the Unicode string representation of certain
    objects will contain Unicode characters.  For example, when
    ``use_unicode`` is `False` (default)::

        >>> from astropy import units as u
        >>> print(unicode(u.degree))
        deg

    When ``use_unicode`` is `True`::

        >>> from astropy import units as u
        >>> print(unicode(u.degree))
        °

    See `handling-unicode
    <http://docs.astropy.org/en/v0.3/development/codeguide.html#unicode-guidelines>`_
    for more information. [#1441]

    - ``astropy.utils.misc.find_api_page`` is now imported into the top-level.
      This allows usage like ``astropy.find_api_page(astropy.units.Quantity)``.
      [#1779]

- ``astropy.convolution``

  - New class-based system for generating kernels, replacing ``make_kernel``.
    [#1255] The ``astropy.nddata.convolution`` sub-package has now been moved
    to ``astropy.convolution``. [#1451]

- ``astropy.coordinates``

  - Two classes ``astropy.coordinates.Longitude`` and
    ``astropy.coordinates.Latitude`` have been added.  These are derived from
    the new ``Angle`` class and used for all longitude-like (RA, azimuth,
    galactic L) and latitude-like coordinates (Dec, elevation, galactic B)
    respectively.  The ``Longitude`` class provides auto-wrapping capability
    and ``Latitude`` performs bounds checking.

  - ``astropy.coordinates.Distance`` supports conversion to and from distance
    modulii. [#1472]

  - ``astropy.coordinates.SphericalCoordinateBase`` and derived classes now
    support arrays of coordinates, enabling large speed-ups for some operations
    on multiple coordinates at the same time. These coordinates can also be
    indexed using standard slicing or any Numpy-compatible indexing. [#1535,
    #1615]

  - Array coordinates can be matched to other array coordinates, finding the
    closest matches between the two sets of coordinates (see the
    ``astropy.coordinates.matching.match_coordinates_3d`` and
    ``astropy.coordinates.matching.match_coordinates_sky`` functions). [#1535]

- ``astropy.cosmology``

  - Added support for including massive Neutrinos in the cosmology classes. The
    Planck (2013) cosmology has been updated to use this. [#1364]

  - Calculations now use and return ``Quantity`` objects where appropriate.
    [#1237]

- ``astropy.io.ascii``

  - Added support for writing IPAC format tables [#1152].

- ``astropy.io.fits``

  - Added initial support for table columns containing pseudo-unsigned
    integers.  This is currently enabled by using the ``uint=True`` option when
    opening files; any table columns with the correct BZERO value will be
    interpreted and returned as arrays of unsigned integers. [#906]

  - Upgraded vendored copy of CFITSIO to v3.35, though backwards compatibility
    back to version v3.28 is maintained.

  - Added support for reading and writing tables using the Q format for columns.
    The Q format is identical to the P format (variable-length arrays) except
    that it uses 64-bit integers for the data descriptors, allowing more than
    4 GB of variable-length array data in a single table.

  - Some refactoring of the table and ``FITS_rec`` modules in order to better
    separate the details of the FITS binary and ASCII table data structures from
    the HDU data structures that encapsulate them.  Most of these changes should
    not be apparent to users (but see API Changes below).

- ``astropy.io.votable``

  - Updated to support the VOTable 1.3 draft. [#433]

  - Added the ability to look up and group elements by their utype attribute.
    [#622]

  - The format of the units of a VOTable file can be specified using the
    ``unit_format`` parameter.  Note that units are still always written out
    using the CDS format, to ensure compatibility with the standard.

- ``astropy.modeling``

  - Added a new framework for representing and evaluating mathematical models
    and for fitting data to models.  See "What's New in Astropy 0.3" in the
    documentation for further details. [#493]

- ``astropy.stats``

  - Added robust statistics functions
    ``astropy.stats.funcs.median_absolute_deviation``,
    ``astropy.stats.funcs.biweight_location``, and
    ``astropy.stats.funcs.biweight_midvariance``. [#621]

  - Added ``astropy.stats.funcs.signal_to_noise_oir_ccd`` for computing the
    signal to noise ratio for source being observed in the optical/IR using a
    CCD. [#870]

  - Add ``axis=int`` option to ``stropy.stats.funcs.sigma_clip`` to allow
    clipping along a given axis for multidimensional data. [#1083]

- ``astropy.table``

  - New columns can be added to a table via assignment to a non-existing
    column by name. [#726]

  - Added ``join`` function to perform a database-like join on two tables. This
    includes support for inner, left, right, and outer joins as well as
    metadata merging.  [#903]

  - Added ``hstack`` and ``vstack`` functions to stack two or more tables.
    [#937]

  - Tables now have a ``.copy`` method and include support for ``copy`` and
    ``deepcopy``. [#1208]

  - Added support for selecting and manipulating groups within a table with
    a database style ``group_by`` method. [#1424]

  - Table ``read`` and ``write`` functions now include rudimentary support
    reading and writing of FITS tables via the unified reading/writing
    interface. [#591]

  - The ``units`` and ``dtypes`` attributes and keyword arguments in Column,
    MaskedColumn, Row, and Table are now deprecated in favor of the
    single-tense ``unit`` and ``dtype``. [#1174]

  - Setting a column from a Quantity now correctly sets the unit on the Column
    object. [#732]

  - Add ``remove_row`` and ``remove_rows`` to remove table rows. [#1230]

  - Added a new ``Table.show_in_browser`` method that opens a web browser
    and displays the table rendered as HTML. [#1342]

  - New tables can now be instantiated using a single row from an existing
    table. [#1417]

- ``astropy.time``

  - New ``Time`` objects can be instantiated from existing ``Time`` objects
    (but with different format, scale, etc.) [#889]

  - Added a ``Time.now`` classmethod that returns the current UTC time,
    similarly to Python's ``datetime.now``. [#1061]

  - Update internal time manipulations so that arithmetic with Time and
    TimeDelta objects maintains sub-nanosecond precision over a time span
    longer than the age of the universe. [#1189]

  - Use ``astropy.utils.iers`` to provide ``delta_ut1_utc``, so that
    automatic calculation of UT1 becomes possible. [#1145]

  - Add ``datetime`` format which allows converting to and from standard
    library ``datetime.datetime`` objects. [#860]

  - Add ``plot_date`` format which allows converting to and from the date
    representation used when plotting dates with matplotlib via the
    ``matplotlib.pyplot.plot_date`` function. [#860]

  - Add ``gps`` format (seconds since 1980-01-01 00:00:00 UTC,
    including leap seconds) [#1164]

  - Add array indexing to Time objects [#1132]

  - Allow for arithmetic of multi-element and single-element Time and TimeDelta
    objects. [#1081]

  - Allow multiplication and division of TimeDelta objects by
    constants and arrays, as well as changing sign (negation) and
    taking the absolute value of TimeDelta objects. [#1082]

  - Allow comparisons of Time and TimeDelta objects. [#1171]

  - Support interaction of Time and Quantity objects that represent a time
    interval. [#1431]

- ``astropy.units``

  - Added parallax equivalency for length-angle. [#985]

  - Added mass-energy equivalency. [#1333]

  - Added a new-style format method which will use format specifiers
    (like ``0.03f``) in new-style format strings for the Quantity's value.
    Specifiers which can't be applied to the value will fall back to the
    entire string representation of the quantity. [#1383]

  - Added support for complex number values in quantities. [#1384]

  - Added new spectroscopic equivalencies for velocity conversions
    (relativistic, optical, and radio conventions are supported) [#1200]

  - The ``spectral`` equivalency now also handles wave number.

  - The ``spectral_density`` equivalency now also accepts a Quantity for the
    frequency or wavelength. It also handles additional flux units.

  - Added Brightness Temperature (antenna gain) equivalency for conversion
    between :math:`T_B` and flux density. [#1327]

  - Added percent unit, and allowed any string containing just a number to be
    interpreted as a scaled dimensionless unit. [#1409]

  - New-style format strings can be used to set the unit output format.  For
    example, ``"{0:latex}".format(u.km)`` will print with the latex formatter.
    [#1462]

  - The ``Unit.is_equivalent`` method can now take a tuple. In this case, the
    method returns ``True`` if the unit is equivalent to any of the units
    listed in the tuple. [#1521]

  - ``def_unit`` can now take a 2-tuple of names of the form (short, long),
    where each entry is a list.  This allows for handling strange units that
    might have multiple short names. [#1543]

  - Added ``dimensionless_angles`` equivalency, which allows conversion of any
    power of radian to dimensionless. [#1161]

  - Added the ability to enable set of units, or equivalencies that are used by
    default.  Also provided context managers for these cases. [#1268]

  - Imperial units are disabled by default. [#1593, #1662]

  - Added an ``astropy.units.add_enabled_units`` context manager, which allows
    creating a temporary context with additional units temporarily enabled in
    the global units namespace. [#1662]

  - ``Unit`` instances now have ``.si`` and ``.cgs`` properties a la
    ``Quantity``.  These serve as shortcuts for ``Unit.to_system(cgs)[0]``
    etc. [#1610]

- ``astropy.vo``

  - New package added to support Virtual Observatory Simple Cone Search query
    and service validation. [#552]

- ``astropy.wcs``

  - Fixed attribute error in ``astropy.wcs.Wcsprm`` (lattype->lattyp) [#1463]

  - Included a new command-line script called ``wcslint`` and accompanying API
    for validating the WCS in a given FITS file or header. [#580]

  - Upgraded included version of WCSLIB to 4.19.

- ``astropy.utils``

  - Added a new set of utilities in ``astropy.utils.timer`` for analyzing the
    runtime of functions and making runtime predections for larger inputs.
    [#743]

  - ``ProgressBar`` and ``Spinner`` classes can now be used directly to return
    generator expressions. [#771]

  - Added ``astropy.utils.iers`` which allows reading in of IERS A or IERS B
    bulletins and interpolation in UT1-UTC.

  - Added a function ``astropy.utils.find_api_page``--given a class or object
    from the ``astropy`` package, this will open that class's API documentation
    in a web browser. [#663]

  - Data download functions such as ``download_file`` now accept a
    ``show_progress`` argument to suppress console output, and a ``timeout``
    argument. [#865, #1258]

- ``astropy.extern.six``

  - Added `six <https://pypi.python.org/pypi/six/>`_ for python2/python3
    compatibility

- Astropy now uses the ERFA library instead of the IAU SOFA library for
  fundamental time transformation routines.  The ERFA library is derived, with
  permission, from the IAU SOFA library but is distributed under a BSD license.
  See ``license/ERFA.rst`` for details. [#1293]

- ``astropy.logger``

  - The Astropy logger now no longer catches exceptions by default, and also
    only captures warnings emitted by Astropy itself (prior to this change,
    following an import of Astropy, any warning got re-directed through the
    Astropy logger). Logging to the Astropy log file has also been disabled by
    default. However, users of Astropy 0.2 will likely still see the previous
    behavior with Astropy 0.3 for exceptions and logging to file since the
    default configuration file installed by 0.2 set the exception logging to be
    on by default. To get the new behavior, set the ``log_exceptions`` and
    ``log_to_file`` configuration items to ``False`` in the ``astropy.cfg``
    file. [#1331]

API Changes
^^^^^^^^^^^

- General

  - The configuration option ``utils.console.use_unicode`` has been
    moved to the top level and renamed to ``unicode_output``.  It now
    not only affects console widgets, such as progress bars, but also
    controls whether calling `unicode` on certain classes will return a
    string containing unicode characters.

- ``astropy.coordinates``

  - The ``astropy.coordinates.Angle`` class is now a subclass of
    ``astropy.units.Quantity``. This means it has all of the methods of a
    `numpy.ndarray`. [#1006]

  - The ``astropy.coordinates.Distance`` class is now a subclass of
    ``astropy.units.Quantity``. This means it has all of the methods of a
    `numpy.ndarray`. [#1472]

    - All angular units are now supported, not just ``radian``, ``degree`` and
      ``hour``, but now ``arcsecond`` and ``arcminute`` as well.  The object
      will retain its native unit, so when printing out a value initially
      provided in hours, its ``to_string()`` will, by default, also be
      expressed in hours.

    - The ``Angle`` class now supports arrays of angles.

    - To be consistent with ``units.Unit``, ``Angle.format`` has been
      deprecated and renamed to ``Angle.to_string``.

    - To be consistent with ``astropy.units``, all plural forms of unit names
      have been removed.  Therefore, the following properties of
      ``astropy.coordinates.Angle`` should be renamed:

      - ``radians`` -> ``radian``
      - ``degrees`` -> ``degree``
      - ``hours`` -> ``hour``

    - Multiplication and division of two ``Angle`` objects used to raise
      ``NotImplementedError``.  Now they raise ``TypeError``.

  - The ``astropy.coordinates.Angle`` class no longer has a ``bounds``
    attribute so there is no bounds-checking or auto-wrapping at this level.
    This allows ``Angle`` objects to be used in arbitrary arithmetic
    expressions (e.g. coordinate distance computation).

  - The ``astropy.coordinates.RA`` and ``astropy.coordinates.Dec`` classes have
    been removed and replaced with ``astropy.coordinates.Longitude`` and
    ``astropy.coordinates.Latitude`` respectively.  These are now used for the
    components of Galactic and Horizontal (Alt-Az) coordinates as well instead
    of plain ``Angle`` objects.

  - ``astropy.coordinates.angles.rotation_matrix`` and
    ``astropy.coordinates.angles.angle_axis`` now take a ``unit`` kwarg instead
    of ``degrees`` kwarg to specify the units of the angles.
    ``rotation_matrix`` will also take the unit from the given ``Angle`` object
    if no unit is provided.

  - The ``AngularSeparation`` class has been removed.  The output of the
    coordinates ``separation()`` method is now an
    ``astropy.coordinates.Angle``.  [#1007]

  - The coordinate classes have been renamed in a way that remove the
    ``Coordinates`` at the end of the class names.  E.g., ``ICRSCoordinates``
    from previous versions is now called ``ICRS``. [#1614]

  - ``HorizontalCoordinates`` are now named ``AltAz``, to reflect more common
    terminology.

- ``astropy.cosmology``

  - The Planck (2013) cosmology will likely give slightly different (and more
    accurate) results due to the inclusion of Neutrino masses. [#1364]

  - Cosmology class properties now return ``Quantity`` objects instead of
    simple floating-point values. [#1237]

  - The names of cosmology instances are now truly optional, and are set to
    ``None`` rather than the name of the class if the user does not provide
    them.  [#1705]

- ``astropy.io.ascii``

  - In the ``read`` method of ``astropy.io.ascii``, empty column values in an
    ASCII table are now treated as missing values instead of the previous
    treatment as a zero-length string "".  This now corresponds to the behavior
    of other table readers like ``numpy.genfromtxt``.  To restore the previous
    behavior set ``fill_values=None`` in the call to ``ascii.read()``. [#919]

  - The ``read`` and ``write`` methods of ``astropy.io.ascii`` now have a
    ``format`` argument for specifying the file format.  This is the preferred
    way to choose the format instead of the ``Reader`` and ``Writer``
    arguments. [#961]

  - The ``include_names`` and ``exclude_names`` arguments were removed from
    the ``BaseHeader`` initializer, and now instead handled by the reader and
    writer classes directly. [#1350]

  - Allow numeric and otherwise unusual column names when reading a table
    where the ``format`` argument is specified, but other format details such
    as the delimiter or quote character are being guessed. [#1692]

  - When reading an ASCII table using the ``Table.read()`` method, the default
    has changed from ``guess=False`` to ``guess=True`` to allow auto-detection
    of file format.  This matches the default behavior of ``ascii.read()``.

- ``astropy.io.fits``

  - The ``astropy.io.fits.new_table`` function is marked "pending deprecation".
    This does not mean it will be removed outright or that its functionality
    has changed.  It will likely be replaced in the future for a function with
    similar, if not subtly different functionality.  A better, if not slightly
    more verbose approach is to use ``pyfits.FITS_rec.from_columns`` to create
    a new ``FITS_rec`` table--this has the same interface as
    ``pyfits.new_table``.  The difference is that it returns a plan
    ``FITS_rec`` array, and not an HDU instance.  This ``FITS_rec`` object can
    then be used as the data argument in the constructors for ``BinTableHDU``
    (for binary tables) or ``TableHDU`` (for ASCII tables).  This is analogous
    to creating an ``ImageHDU`` by passing in an image array.
    ``pyfits.FITS_rec.from_columns`` is just a simpler way of creating a
    FITS-compatible recarray from a FITS column specification.

  - The ``updateHeader``, ``updateHeaderData``, and ``updateCompressedData``
    methods of the ``CompDataHDU`` class are pending deprecation and moved to
    internal methods.  The operation of these methods depended too much on
    internal state to be used safely by users; instead they are invoked
    automatically in the appropriate places when reading/writing compressed
    image HDUs.

  - The ``CompDataHDU.compData`` attribute is pending deprecation in favor of
    the clearer and more PEP-8 compatible ``CompDataHDU.compressed_data``.

  - The constructor for ``CompDataHDU`` has been changed to accept new keyword
    arguments.  The new keyword arguments are essentially the same, but are in
    underscore_separated format rather than camelCase format.  The old
    arguments are still pending deprecation.

  - The internal attributes of HDU classes ``_hdrLoc``, ``_datLoc``, and
    ``_datSpan`` have been replaced with ``_header_offset``, ``_data_offset``,
    and ``_data_size`` respectively.  The old attribute names are still pending
    deprecation.  This should only be of interest to advanced users who have
    created their own HDU subclasses.

  - The following previously deprecated functions and methods have been removed
    entirely: ``createCard``, ``createCardFromString``, ``upperKey``,
    ``ColDefs.data``, ``setExtensionNameCaseSensitive``, ``_File.getfile``,
    ``_TableBaseHDU.get_coldefs``, ``Header.has_key``, ``Header.ascardlist``.

  - Interfaces that were pending deprecation are now fully deprecated.  These
    include: ``create_card``, ``create_card_from_string``, ``upper_key``,
    ``Header.get_history``, and ``Header.get_comment``.

  - The ``.name`` attribute on HDUs is now directly tied to the HDU's header, so
    that if ``.header['EXTNAME']`` changes so does ``.name`` and vice-versa.

- ``astropy.io.registry``

  - Identifier functions for reading/writing Table and NDData objects should
    now accept ``(origin, *args, **kwargs)`` instead of ``(origin, args,
    kwargs)``. [#591]

  - Added a new ``astropy.io.registry.get_formats`` function for listing
    registered I/O formats and details about the their readers/writers. [#1669]

- ``astropy.io.votable``

  - Added a new option ``use_names_over_ids`` option to use when converting
    from VOTable objects to Astropy Tables. This can prevent a situation where
    column names are not preserved when converting from a VOTable. [#609]

- ``astropy.nddata``

  - The ``astropy.nddata.convolution`` sub-package has now been moved to
    ``astropy.convolution``, and the ``make_kernel`` function has been removed.
    (the kernel classes should be used instead) [#1451]

- ``astropy.stats.funcs``

  - For ``sigma_clip``, the ``maout`` optional parameter has been removed, and
    the function now always returns a masked array.  A new boolean parameter
    ``copy`` can be used to indicated whether the input data should be copied
    (``copy=True``, default) or used by reference (``copy=False``) in the
    output masked array. [#1083]

- ``astropy.table``

  - The first argument to the ``Column`` and ``MaskedColumn`` classes is now
    the data array--the ``name`` argument has been changed to an optional
    keyword argument. [#840]

  - Added support for instantiating a ``Table`` from a list of dict, each one
    representing a single row with the keys mapping to column names. [#901]

  - The plural 'units' and 'dtypes' have been switched to 'unit' and 'dtype'
    where appropriate. The original attributes are still present in this
    version as deprecated attributes, but will be removed in the next version.
    [#1174]

  - The ``copy`` methods of ``Column`` and ``MaskedColumn`` were changed so
    that the first argument is now ``order='C'``.  This is required for
    compatibility with Numpy 1.8 which is currently in development. [#1250]

  - Comparing a column (with == or !=) to a scalar, an array, or another column
    now always returns a boolean Numpy array (which is a masked array if either
    of the arguments in the comparison was masked). This is in contrast to the
    previous behavior, which in some cases returned a boolean Numpy array, and
    in some cases returned a boolean Column object. [#1446]

- ``astropy.time``

  - For consistency with ``Quantity``, the attributes ``val`` and
    ``is_scalar`` have been renamed to ``value`` and ``isscalar``,
    respectively, and the attribute ``vals`` has been dropped. [#767]

  - The double-float64 internal representation of time is used more
    efficiently to enable better accuracy. [#366]

  - Format and scale arguments are now allowed to be case-insensitive. [#1128]

- ``astropy.units``

  - The ``Quantity`` class now inherits from the Numpy array class, and
    includes the following API changes [#929]:

    - Using ``float(...)``, ``int(...)``, and ``long(...)`` on a quantity will
      now only work if the quantity is dimensionless and unscaled.

    - All Numpy ufuncs should now treat units correctly (or raise an exception
      if not supported), rather than extract the value of quantities and
      operate on this, emitting a warning about the implicit loss of units.

    - When using relevant Numpy ufuncs on dimensionless quantities (e.g.
      ``np.exp(h * nu / (k_B * T))``), or combining dimensionless quantities
      with Python scalars or plain Numpy arrays ``1 + v / c``, the
      dimensionless Quantity will automatically be converted to an unscaled
      dimensionless Quantity.

    - When initializing a quantity from a value with no unit, it is now set to
      be dimensionless and unscaled by default. When initializing a Quantity
      from another Quantity and with no unit specified in the initializer, the
      unit is now taken from the unit of the Quantity being initialized from.

  - Strings are no longer allowed as the values for Quantities. [#1005]

  - Quantities are always comparable with zero regardless of their units.
    [#1254]

  - The exception ``astropy.units.UnitsException`` has been renamed to
    ``astropy.units.UnitsError`` to be more consistent with the naming
    of built-in Python exceptions. [#1406]

  - Multiplication with and division by a string now always returns a Unit
    (rather than a Quantity when the string was first) [#1408]

  - Imperial units are disabled by default.

- ``astropy.wcs``

  - For those including the ``astropy.wcs`` C headers in their project, they
    should now include it as:

       #include "astropy_wcs/astropy_wcs_api.h"

    instead of:

       #include "astropy_wcs_api.h"

    [#1631]

- The ``--enable-legacy`` option for ``setup.py`` has been removed. [#1493]

Bug Fixes
^^^^^^^^^

- ``astropy.io.ascii``

  - The ``write()`` function was ignoring the ``fill_values`` argument. [#910]

  - Fixed an issue in ``DefaultSplitter.join`` where the delimiter attribute
    was ignored when writing the CSV. [#1020]

  - Fixed writing of IPAC tables containing null values. [#1366]

  - When a table with no header row was read without specifying the format and
    using the ``names`` argument, then the first row could be dropped. [#1692]

- ``astropy.io.fits``

  - Binary tables containing compressed images may, optionally, contain other
    columns unrelated to the tile compression convention. Although this is an
    uncommon use case, it is permitted by the standard.

  - Reworked some of the file I/O routines to allow simpler, more consistent
    mapping between OS-level file modes ('rb', 'wb', 'ab', etc.) and the more
    "PyFITS-specific" modes used by PyFITS like "readonly" and "update".  That
    is, if reading a FITS file from an open file object, it doesn't matter as
    much what "mode" it was opened in so long as it has the right capabilities
    (read/write/etc.)  Also works around bugs in the Python io module in 2.6+
    with regard to file modes.

  - Fixed a long-standing issue where writing binary tables did not correctly
    write the TFORMn keywords for variable-length array columns (they omitted
    the max array length parameter of the format).  This was thought fixed in
    an earlier version, but it was only fixed for compressed image HDUs and
    not for binary tables in general.

- ``astropy.nddata``

  - Fixed crash when trying to multiple or divide ``NDData`` objects with
    uncertainties. [#1547]

- ``astropy.table``

  - Using a list of strings to index a table now correctly returns a new table
    with the columns named in the list. [#1454]

  - Inequality operators now work properly with ``Column`` objects. [#1685]

- ``astropy.time``

  - ``Time`` scale and format attributes are now shown when calling ``dir()``
    on a ``Time`` object. [#1130]

- ``astropy.wcs``

  - Fixed assignment to string-like WCS attributes on Python 3. [#956]

- ``astropy.units``

  - Fixed a bug that caused the order of multiplication/division of plain
    Numpy arrays with Quantities to matter (i.e. if the plain array comes
    first the units were not preserved in the output). [#899]

  - Directly instantiated ``CompositeUnits`` were made printable without
    crashing. [#1576]

- Misc

  - Fixed various modules that hard-coded ``sys.stdout`` as default arguments
    to functions at import time, rather than using the runtime value of
    ``sys.stdout``. [#1648]

  - Minor documentation fixes and enhancements [#922, #1034, #1210, #1217,
    #1491, #1492, #1498, #1582, #1608, #1621, #1646, #1670, #1756]

  - Fixed a crash that could sometimes occur when running the test suite on
    systems with platform names containing non-ASCII characters. [#1698]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- General

  - Astropy now follows the PSF Code of Conduct. [#1216]

  - Astropy's test suite now tests all doctests in inline docstrings.  Support
    for running doctests in the reST documentation is planned to follow in
    v0.3.1.

  - Astropy's test suite can be run on multiple CPUs in parallel, often
    greatly improving runtime, using the ``--parallel`` option. [#1040]

  - A warning is now issued when using Astropy with Numpy < 1.5--much of
    Astropy may still work in this case but it shouldn't be expected to
    either. [#1479]

  - Added automatic download/build/installation of Numpy during Astropy
    installation if not already found. [#1483]

  - Handling of metadata for the ``NDData`` and ``Table`` classes has been
    unified by way of a common ``MetaData`` descriptor--it allows instantiating
    an object with metadata of any mapping type, and subsequently prevents
    replacing the mapping stored in the ``.meta`` attribute (only direct
    updates to that object are allowed). [#1686]

- ``astropy.coordinates``

  - Angles containing out of bounds minutes or seconds (e.g. 60) can be
    parsed--the value modulo 60 is used with carry to the hours/minutes, and a
    warning is issued rather than raising an exception. [#990]

- ``astropy.io.fits``

  - The new compression code also adds support for the ZQUANTIZ and ZDITHER0
    keywords added in more recent versions of this FITS Tile Compression spec.
    This includes support for lossless compression with GZIP. (#198) By default
    no dithering is used, but the ``SUBTRACTIVE_DITHER_1`` and
    ``SUBTRACTIVE_DITHER_2`` methods can be enabled by passing the correct
    constants to the ``quantize_method`` argument to the ``CompImageHDU``
    constructor.  A seed can be manually specified, or automatically generated
    using either the system clock or checksum-based methods via the
    ``dither_seed`` argument.  See the documentation for ``CompImageHDU`` for
    more details.

  - Images compressed with the Tile Compression standard can now be larger than
    4 GB through support of the Q format.

  - All HDUs now have a ``.ver`` ``.level`` attribute that returns the value of
    the EXTVAL and EXTLEVEL keywords from that HDU's header, if the exist.
    This was added for consistency with the ``.name`` attribute which returns
    the EXTNAME value from the header.

  - Then ``Column`` and ``ColDefs`` classes have new ``.dtype`` attributes
    which give the Numpy dtype for the column data in the first case, and the
    full Numpy compound dtype for each table row in the latter case.

  - There was an issue where new tables created defaulted the values in all
    string columns to '0.0'.  Now string columns are filled with empty strings
    by default--this seems a less surprising default, but it may cause
    differences with tables created with older versions of PyFITS or Astropy.

- ``astropy.io.misc``

  - The HDF5 reader can now refer to groups in the path as well as datasets;
    if given a group, the first dataset in that group is read. [#1159]

- ``astropy.nddata``

  - ``NDData`` objects have more helpful, though still rudimentary ``__str__`
    and ``__repr__`` displays. [#1313]

- ``astropy.units``

  - Added 'cycle' unit. [#1160]

  - Extended units supported by the CDS formatter/parser. [#1468]

  - Added unicode an LaTeX symbols for liter. [#1618]

- ``astropy.wcs``

  - Redundant SCAMP distortion parameters are removed with SIP distortions are
    also present. [#1278]

  - Added iterative implementation of ``all_world2pix`` that can be reliably
    inverted. [#1281]


0.2.5 (2013-10-25)
------------------

Bug Fixes
^^^^^^^^^

- ``astropy.coordinates``

  - Fixed incorrect string formatting of Angles using ``precision=0``. [#1319]

  - Fixed string formatting of Angles using ``decimal=True`` which ignored the
    ``precision`` argument. [#1323]

  - Fixed parsing of format strings using appropriate unicode characters
    instead of the ASCII ``-`` for minus signs. [#1429]

- ``astropy.io.ascii``

  - Fixed a crash in the IPAC table reader when the ``include/exclude_names``
    option is set. [#1348]

  - Fixed writing AASTex tables to honor the ``tabletype`` option. [#1372]

- ``astropy.io.fits``

  - Improved round-tripping and preservation of manually assigned column
    attributes (``TNULLn``, ``TSCALn``, etc.) in table HDU headers. (Note: This
    issue was previously reported as fixed in Astropy v0.2.2 by mistake; it is
    not fixed until v0.3.) [#996]

  - Fixed a bug that could cause a segfault when trying to decompress an
    compressed HDU whose contents are truncated (due to a corrupt file, for
    example). This still causes a Python traceback but better that than a
    segfault. [#1332]

  - Newly created ``CompImageHDU`` HDUs use the correct value of the
    ``DEFAULT_COMPRESSION_TYPE`` module-level constant instead of hard-coding
    "RICE_1" in the header.

  - Fixed a corner case where when extra memory is allocated to compress an
    image, it could lead to unnecessary in-memory copying of the compressed
    image data and a possible memory leak through Numpy.

  - Fixed a bug where assigning from an mmap'd array in one FITS file over
    the old (also mmap'd) array in another FITS file failed to update the
    destination file. Corresponds to PyFITS issue 25.

  - Some miscellaneous documentation fixes.

- ``astropy.io.votable``

  - Added a warning for when a VOTable 1.2 file contains no ``RESOURCES``
    elements (at least one should be present). [#1337]

  - Fixed a test failure specific to MIPS architecture caused by an errant
    floating point warning. [#1179]

- ``astropy.nddata.convolution``

  - Prevented in-place modification of the input arrays to ``convolve()``.
    [#1153]

- ``astropy.table``

  - Added HTML escaping for string values in tables when outputting the table
    as HTML. [#1347]

  - Added a workaround in a bug in Numpy that could cause a crash when
    accessing a table row in a masked table containing ``dtype=object``
    columns. [#1229]

  - Fixed an issue similar to the one in #1229, but specific to unmasked
    tables. [#1403]

- ``astropy.units``

  - Improved error handling for unparseable units and fixed parsing CDS units
    without mantissas in the exponent. [#1288]

  - Added a physical type for spectral flux density. [#1410]

  - Normalized conversions that should result in a scale of exactly 1.0 to
    round off slight floating point imprecisions. [#1407]

  - Added support in the CDS unit parser/formatter for unusual unit prefixes
    that are nonetheless required to be supported by that convention. [#1426]

  - Fixed the parsing of ``sqrt()`` in unit format strings which was returning
    ``unit ** 2`` instead of ``unit ** 0.5``. [#1458]

- ``astropy.wcs``

  - When passing a single array to the wcs transformation functions,
    (``astropy.wcs.Wcs.all_pix2world``, etc.), its second dimension must now
    exactly match the number of dimensions in the transformation. [#1395]

  - Improved error message when incorrect arguments are passed to
    ``WCS.wcs_world2pix``. [#1394]

  - Fixed a crash when trying to read WCS from FITS headers on Python 3.3
    in Windows. [#1363]

  - Only headers that are required as part of the WCSLIB C API are installed
    by the package, per request of system packagers. [#1666]

- Misc

  - Fixed crash when the ``COLUMNS`` environment variable is set to a
    non-integer value. [#1291]

  - Fixed a bug in ``ProgressBar.map`` where ``multiprocess=True`` could cause
    it to hang on waiting for the process pool to be destroyed. [#1381]

  - Fixed a crash on Python 3.2 when affiliated packages try to use the
    ``astropy.utils.data.get_pkg_data_*`` functions. [#1256]

  - Fixed a minor path normalization issue that could occur on Windows in
    ``astropy.utils.data.get_pkg_data_filename``. [#1444]

  - Fixed an annoyance where configuration items intended only for testing
    showed up in users' astropy.cfg files. [#1477]

  - Prevented crashes in exception logging in unusual cases where no traceback
    is associated with the exception. [#1518]

  - Fixed a crash when running the tests in unusual environments where
    ``sys.stdout.encoding`` is ``None``. [#1530]

  - Miscellaneous documentation fixes and improvements [#1308, #1317, #1377,
    #1393, #1362, #1516]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Astropy installation now requests setuptools >= 0.7 during build/installation
  if neither distribute or setuptools >= 0.7 is already installed.  In other
  words, if ``import setuptools`` fails, ``ez_setup.py`` is used to bootstrap
  the latest setuptools (rather than using ``distribute_setup.py`` to bootstrap
  the now obsolete distribute package). [#1197]

- When importing Astropy from a source checkout without having built the
  extension modules first an ``ImportError`` is raised rather than a
  ``SystemExit`` exception. [#1269]


0.2.4 (2013-07-24)
------------------

Bug Fixes
^^^^^^^^^

- ``astropy.coordinates``

  - Fixed the angle parser to support parsing the string "1 degree". [#1168]

- ``astropy.cosmology``

  - Fixed a crash in the ``comoving_volume`` method on non-flat cosmologies
    when passing it an array of redshifts.

- ``astropy.io.ascii``

  - Fixed a bug that prevented saving changes to the comment symbol when
    writing changes to a table. [#1167]

- ``astropy.io.fits``

  - Added a workaround for a bug in 64-bit OSX that could cause truncation when
    writing files greater than 2^32 bytes in size. [#839]

- ``astropy.io.votable``

  - Fixed incorrect reading of tables containing multiple ``<RESOURCE>``
    elements. [#1223]

- ``astropy.table``

  - Fixed a bug where ``Table.remove_column`` and ``Table.rename_column``
    could cause a masked table to lose its masking. [#1120]

  - Fixed bugs where subclasses of ``Table`` did not preserver their class in
    certain operations. [#1142]

  - Fixed a bug where slicing a masked table did not preserve the mask. [#1187]

- ``astropy.units``

  - Fixed a bug where the ``.si`` and ``.cgs`` properties of dimensionless
    ``Quantity`` objects raised a ``ZeroDivisionError``. [#1150]

  - Fixed a bug where multiple subsequent calls to the ``.decompose()`` method
    on array quantities applied a scale factor each time. [#1163]

- Misc

  - Fixed an installation crash that could occur sometimes on Debian/Ubuntu
    and other \*NIX systems where ``pkg_resources`` can be installed without
    installing ``setuptools``. [#1150]

  - Updated the ``distribute_setup.py`` bootstrapper to use setuptools >= 0.7
    when installing on systems that don't already have an up to date version
    of distribute/setuptools. [#1180]

  - Changed the ``version.py`` template so that Astropy affiliated packages can
    (and they should) use their own ``cython_version.py`` and
    ``utils._compiler`` modules where appropriate. This issue only pertains to
    affiliated package maintainers. [#1198]

  - Fixed a corner case where the default config file generation could crash
    if building with matplotlib but *not* Sphinx installed in a virtualenv.
    [#1225]

  - Fixed a crash that could occur in the logging module on systems that
    don't have a default preferred encoding (in particular this happened
    in some versions of PyCharm). [#1244]

  - The Astropy log now supports passing non-string objects (and calling
    ``str()`` on them by default) to the logging methods, in line with Python's
    standard logging API. [#1267]

  - Minor documentation fixes [#582, #696, #1154, #1194, #1212, #1213, #1246,
    #1252]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``astropy.cosmology``

  - Added a new ``Plank13`` object representing the Plank 2013 results. [#895]

- ``astropy.units``

  - Performance improvements in initialization of ``Quantity`` objects with
    a large number of elements. [#1231]


0.2.3 (2013-05-30)
------------------

Bug Fixes
^^^^^^^^^

- ``astropy.time``

  - Fixed inaccurate handling of leap seconds when converting from UTC to UNIX
    timestamps. [#1118]

  - Tightened required accuracy in many of the time conversion tests. [#1121]

- Misc

  - Fixed a regression that was introduced in v0.2.2 by the fix to issue #992
    that was preventing installation of Astropy affiliated packages that use
    Astropy's setup framework. [#1124]


0.2.2 (2013-05-21)
------------------

Bug Fixes
^^^^^^^^^

- ``astropy.io``

  - Fixed issues in both the ``fits`` and ``votable`` sub-packages where array
    byte order was not being handled consistently, leading to possible crashes
    especially on big-endian systems. [#1003]

- ``astropy.io.fits``

  - When an error occurs opening a file in fitsdiff the exception message will
    now at least mention which file had the error.

  - Fixed a couple cases where creating a new table using TDIMn in some of the
    columns could cause a crash.

  - Slightly refactored how tables containing variable-length array columns are
    handled to add two improvements: Fixes an issue where accessing the data
    after a call to the ``astropy.io.fits.getdata`` convenience function caused
    an exception, and allows the VLA data to be read from an existing mmap of
    the FITS file.

  - Fixed a bug on Python 3 where attempting to open a non-existent file on
    Python 3 caused a seemingly unrelated traceback.

  - Fixed an issue in the tests that caused some tests to fail if Astropy is
    installed with read-only permissions.

  - Fixed a bug where instantiating a ``BinTableHDU`` from a numpy array
    containing boolean fields converted all the values to ``False``.

  - Fixed an issue where passing an array of integers into the constructor of
    ``Column()`` when the column type is floats of the same byte width caused
    the column array to become garbled.

  - Fixed inconsistent behavior in creating CONTINUE cards from byte strings
    versus unicode strings in Python 2--CONTINUE cards can now be created
    properly from unicode strings (so long as they are convertable to ASCII).

  - Fixed a bug in parsing HIERARCH keywords that do not have a space after the
    first equals sign (before the value).

  - Prevented extra leading whitespace on HIERARCH keywords from being treated
    as part of the keyword.

  - Fixed a bug where HIERARCH keywords containing lower-case letters was
    mistakenly marked as invalid during header validation along with an
    ancillary issue where the ``Header.index()`` method id not work correctly
    with HIERARCH keywords containing lower-case letters.

  - Disallowed assigning NaN and Inf floating point values as header values,
    since the FITS standard does not define a way to represent them in. Because
    this is undefined, the previous behavior did not make sense and produced
    invalid FITS files. [#954]

  - Fixed an obscure issue that can occur on systems that don't have flush to
    memory-mapped files implemented (namely GNU Hurd). [#968]

- ``astropy.io.votable``

  - Stopped deprecation warnings from the ``astropy.io.votable`` package that
    could occur during setup. [#970]

  - Fixed an issue where INFO elements were being incorrectly dropped when
    occurring inside a TABLE element. [#1000]

  - Fixed obscure test failures on MIPS platforms. [#1010]

- ``astropy.nddata.convolution``

  - Fixed an issue in ``make_kernel()`` when using an Airy function kernel.
    Also removed the superfluous 'brickwall' option. [#939]

- ``astropy.table``

  - Fixed a crash that could occur when adding a row to an empty (rowless)
    table with masked columns. [#973]

  - Made it possible to assign to one table row from the value of another row,
    effectively making it easier to copy rows, for example. [#1019]

- ``astropy.time``

  - Added appropriate ``__copy__`` and ``__deepcopy__`` behavior; this
    omission caused a seemingly unrelated error in FK5 coordinate separation.
    [#891]

- ``astropy.units``

  - Fixed an issue where the ``isiterable()`` utility returned ``True`` for
    quantities with scalar values.  Added an ``__iter__`` method for the
    ``Quantity`` class and fixed ``isiterable()`` to catch false positives.
    [#878]

  - Fixed previously undefined behavior when multiplying a unit by a string.
    [#949]

  - Added 'time' as a physical type--this was a simple omission. [#959]

  - Fixed issues with pickling unit objects so as to play nicer with the
    multiprocessing module. [#974]

  - Made it more difficult to accidentally override existing units with a new
    unit of the same name. [#1070]

  - Added several more physical types and units that were previously omitted,
    including 'mass density', 'specific volume', 'molar volume', 'momentum',
    'angular momentum', 'angular speed', 'angular acceleration', 'electric
    current', 'electric current density', 'electric field strength', 'electric
    flux density', 'electric charge density', 'permittivity', 'electromagnetic
    field strength', 'radiant intensity', 'data quantity', 'bandwidth'; and
    'knots', 'nautical miles', 'becquerels', and 'curies' respectively. [#1072]

- Misc

  - Fixed a permission error that could occur when running ``astropy.test()``
    on Python 3 when Astropy is installed as root. [#811]

  - Made it easier to filter warnings from the ``convolve()`` function and
    from ``Quantity`` objects. [#853]

  - Fixed a crash that could occur in Python 3 when generation of the default
    config file fails during setup. [#952]

  - Fixed an unrelated error message that could occur when trying to import
    astropy from a source checkout without having build the extension modules
    first. This issue was claimed to be fixed in v0.2.1, but the fix itself had
    a bug. [#971]

  - Fixed a crash that could occur when running the ``build_sphinx`` setup
    command in Python 3. [#977]

  - Added a more helpful error message when trying to run the
    ``setup.py build_sphinx`` command when Sphinx is not installed. [#1027]

  - Minor documentation fixes and restructuring.
    [#935, #967, #978, #1004, #1028, #1047]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Some performance improvements to the ``astropy.units`` package, in particular
  improving the time it takes to import the sub-package. [#1015]


0.2.1 (2013-04-03)
------------------

Bug Fixes
^^^^^^^^^

- ``astropy.coordinates``

  - Fixed encoding errors that could occur when formatting coordinate objects
    in code using ``from __future__ import unicode_literals``. [#817]

  - Fixed a bug where the minus sign was dropped when string formatting dms
    coordinates with -0 degrees. [#875]

- ``astropy.io.fits``

  - Properly supports the ZQUANTIZ keyword used to support quantization
    level--this includes working support for lossless GZIP compression of
    images.

  - Fixed support for opening gzipped FITS files in a writeable mode. [#256]

  - Added a more helpful exception message when trying to read invalid values
    from a table when the required ``TNULLn`` keyword is missing. [#309]

  - More refactoring of the tile compression handling to work around a
    potential memory access violation that was particularly prevalent on
    Windows. [#507]

  - Fixed an integer size mismatch in the compression module that could affect
    32-bit systems. [#786]

  - Fixed malformatting of the ``TFORMn`` keywords when writing compressed
    image tables (they omitted the max array length parameter from the
    variable-length array format).

  - Fixed a crash that could occur when writing a table containing multi-
    dimensional array columns from an existing file into a new file.

  - Fixed a bug in fitsdiff that reported two header keywords containing NaN
    as having different values.

- ``astropy.io.votable``

  - Fixed links to the ``astropy.io.votable`` documentation in the VOTable
    validator output. [#806]

  - When reading VOTables containing integers that are out of range for their
    column type, display a warning rather than raising an exception. [#825]

  - Changed the default string format for floating point values for better
    round-tripping. [#856]

  - Fixed opening VOTables through the ``Table.read()`` interface for tables
    that have no names. [#927]

  - Fixed creation of VOTables from an Astropy table that does not have a data
    mask. [#928]

  - Minor documentation fixes. [#932]

- ``astropy.nddata.convolution``

  - Added better handling of ``inf`` values to the ``convolve_fft`` family of
    functions. [#893]

- ``astropy.table``

  - Fixed silent failure to assign values to a row on multiple columns. [#764]

  - Fixed various buggy behavior when viewing a table after sorting by one of
    its columns. [#829]

  - Fixed using ``numpy.where()`` with table indexing. [#838]

  - Fixed a bug where opening a remote table with ``Table.read()`` could cause
    the entire table to be downloaded twice. [#845]

  - Fixed a bug where ``MaskedColumn`` no longer worked if the column being
    masked is renamed. [#916]

- ``astropy.units``

  - Added missing capability for array ``Quantity``\s to be initializable by
    a list of ``Quantity``\s. [#835]

  - Fixed the definition of year and lightyear to be in terms of Julian year
    per the IAU definition. [#861]

  - "degree" was removed from the list of SI base units. [#863]

- ``astropy.wcs``

  - Fixed ``TypeError`` when calling ``WCS.to_header_string()``. [#822]

  - Added new method ``WCS.all_world2pix`` for converting from world
    coordinates to pixel space, including inversion of the astrometric
    distortion correction. [#1066, #1281]


- Misc

  - Fixed a minor issue when installing with ``./setup.py develop`` on a fresh
    git clone.  This is likely only of interest to developers on Astropy.
    [#725]

  - Fixes a crash with ``ImportError: No module named 'astropy.version'`` when
    running setup.py from a source checkout for the first time on OSX with
    Python 3.3. [#820]

  - Fixed an installation issue where running ``./setup.py install`` or when
    installing with pip the ``.astropy`` directory gets created in the home
    directory of the user running the command.  The user's ``.astropy``
    directory should only be created when they use Astropy, not when they
    install it. [#867]

  - Fixed an exception when creating a ``ProgressBar`` with a "total" of 0.
    [#752]

  - Added better documentation of behavior that can occur when trying to import
    the astropy package from within a source checkout without first building
    the extension modules. [#795, #864]

  - Added link to the installation instructions in the README. [#797]

  - Catches segfaults in xmllint which can occur sometimes and is otherwise out
    of our control. [#803]

  - Minor changes to the documentation template. [#805]

  - Fixed a minor exception handling bug in ``download_file()``. [#808]

  - Added cleanup of any temporary files if an error occurs in
    ``download_file()``. [#857]

  - Filesystem free space is checked for before attempting to download a file
    with ``download_file()``. [#858]

  - Fixed package data locating to work across symlinks--required to work with
    some OS packaging layouts. [#827]

  - Fixed a bug when building Cython extensions where hidden files containing
    ``.pyx`` extensions could cause the build to crash. This can be an issue
    with software and filesystems that autogenerate hidden files. [#834]

  - Fixed bug that could cause a "script" called README.rst to be installed
    in a bin directory. [#852]

  - Fixed some miscellaneous and mostly rare reference leaks caught by
    cpychecker. [#914]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added logo and branding for Windows binary installers. [#741]

- Upgraded included version libexpat to 2.1.0. [#781]

- ~25% performance improvement in unit composition/decomposition. [#836]

- Added previously missing LaTeX formatting for ``L_sun`` and ``R_sun``. [#841]

- ``ConfigurationItem``\s now have a more useful and informative ``__repr__``
  and improved documentation for how to use them. [#855]

- Added a friendlier error message when trying to import astropy from a source
  checkout without first building the extension modules inplace. [#864]

- ``py.test`` now outputs more system information for help in debugging issues
  from users. [#869]

- Added unit definitions "mas" and "uas" for "milliarcsecond" and
  "microarcsecond" respectively. [#892]


0.2 (2013-02-19)
----------------

New Features
^^^^^^^^^^^^

This is a brief overview of the new features included in Astropy 0.2--please
see the "What's New" section of the documentation for more details.

- ``astropy.coordinates``

  - This new subpackage contains a representation of celestial coordinates,
    and provides a wide range of related functionality.  While
    fully-functional, it is a work in progress and parts of the API may
    change in subsequent releases.

- ``astropy.cosmology``

  - Update to include cosmologies with variable dark energy equations of state.
    (This introduces some API incompatibilities with the older Cosmology
    objects).

  - Added parameters for relativistic species (photons, neutrinos) to the
    astropy.cosmology classes. The current treatment assumes that neutrinos are
    massless. [#365]

  - Add a WMAP9 object using the final (9-year) WMAP parameters from
    Hinshaw et al. 2013. It has also been made the default cosmology.
    [#629, #724]

- ``astropy.table`` I/O infrastructure for custom readers/writers
  implemented. [#305]

  - Added support for reading/writing HDF5 files [#461]

  - Added support for masked tables with missing or invalid data [#451]

- New ``astropy.time`` sub-package. [#332]

- New ``astropy.units`` sub-package that includes a class for units
  (``astropy.units.Unit``) and scalar quantities that have units
  (``astropy.units.Quantity``). [#370, #445]

  This has the following effects on other sub-packages:

  - In ``astropy.wcs``, the ``wcs.cunit`` list now takes and returns
    ``astropy.units.Unit`` objects. [#379]

  - In ``astropy.nddata``, units are now stored as ``astropy.units.Unit``
    objects. [#382]

  - In ``astropy.table``, units on columns are now stored as
    ``astropy.units.Unit`` objects. [#380]

  - In ``astropy.constants``, constants are now stored as
    ``astropy.units.Quantity`` objects. [#529]

- ``astropy.io.ascii``

  - Improved integration with the ``astropy.table`` Table class so that
    table and column metadata (e.g. keywords, units, description,
    formatting) are directly available in the output table object.  The
    CDS, DAOphot, and IPAC format readers now provide this type of
    integrated metadata.

  - Changed to using `astropy.table` masked tables instead of NumPy
    masked arrays for tables with missing values.

  - Added SExtractor table reader to ``astropy.io.ascii`` [#420]

  - Removed the Memory reader class which was used to convert data input
    passed to the ``write`` function into an internal table.  Instead
    ``write`` instantiates an astropy Table object using the data
    input to ``write``.

  - Removed the NumpyOutputter as the output of reading a table is now
    always a ``Table`` object.

  - Removed the option of supplying a function as a column output
    formatter.

  - Added a new ``strip_whitespace`` keyword argument to the ``write``
    function.  This controls whether whitespace is stripped from
    the left and right sides of table elements before writing.
    Default is True.

  - Fixed a bug in reading IPAC tables with null values.

- Generalized I/O infrastructure so that ``astropy.nddata`` can also have
  custom readers/writers [#659]

- ``astropy.wcs``

  - From updating the underlying wcslib 4.16:

    - When ``astropy.wcs.WCS`` constructs a default coordinate representation
      it will give it the special name "DEFAULTS", and will not report "Found
      one coordinate representation".

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- A configuration file with all options set to their defaults is now generated
  when astropy is installed.  This file will be pulled in as the users'
  astropy configuration file the first time they ``import astropy``.  [#498]

- Astropy doc themes moved into ``astropy.sphinx`` to allow affiliated packages
  to access them.

- Added expanded documentation for the ``astropy.cosmology`` sub-package.
  [#272]

- Added option to disable building of "legacy" packages (pyfits, vo, etc.).

- The value of the astronomical unit (au) has been updated to that adopted by
  IAU 2012 Resolution B2, and the values of the pc and kpc constants have been
  updated to reflect this. [#368]

- Added links to the documentation pages to directly edit the documentation on
  GitHub. [#347]

- Several updates merged from ``pywcs`` into ``astropy.wcs`` [#384]:

  - Improved the reading of distortion images.

  - Added a new option to choose whether or not to write SIP coefficients.

  - Uses the ``relax`` option by default so that non-standard keywords are
    allowed. [#585]


- Added HTML representation of tables in IPython notebook [#409]

- Rewrote CFITSIO-based backend for handling tile compression of FITS files.
  It now uses a standard CFITSIO instead of heavily modified pieces of CFITSIO
  as before.  Astropy ships with its own copy of CFITSIO v3.30, but system
  packagers may choose instead to strip this out in favor of a
  system-installed version of CFITSIO.  This corresponds to PyFITS ticket 169.
  [#318]

- Moved ``astropy.config.data`` to ``astropy.utils.data`` and re-factored the
  I/O routines to separate out the generic I/O code that can be used to open
  any file or resource from the code used to access Astropy-related data. The
  'core' I/O routine is now ``get_readable_fileobj``, which can be used to
  access any local as well as remote data, supports caching, and can decompress
  gzip and bzip2 files on-the-fly. [#425]

- Added a classmethod to
  ``astropy.coordinates.coordsystems.SphericalCoordinatesBase`` that performs a
  name resolve query using Sesame to retrieve coordinates for the requested
  object. This works for any subclass of ``SphericalCoordinatesBase``, but
  requires an internet connection. [#556]

- ``astropy.nddata.convolution`` removed requirement of PyFFTW3; uses Numpy's
  FFT by default instead with the added ability to specify an FFT
  implementation to use. [#660]


Bug Fixes
^^^^^^^^^

- ``astropy.io.ascii``

  - Fixed crash when pprinting a row with INDEF values. [#511]

  - Fixed failure when reading DAOphot files with empty keyword values. [#666]

- ``astropy.io.fits``

  - Improved handling of scaled images and pseudo-unsigned integer images in
    compressed image HDUs.  They now work more transparently like normal image
    HDUs with support for the ``do_not_scale_image_data`` and ``uint`` options,
    as well as ``scale_back`` and ``save_backup``.  The ``.scale()`` method
    works better too. Corresponds to PyFITS ticket 88.

  - Permits non-string values for the EXTNAME keyword when reading in a file,
    rather than throwing an exception due to the malformatting.  Added
    verification for the format of the EXTNAME keyword when writing.
    Corresponds to PyFITS ticket 96.

  - Added support for EXTNAME and EXTVER in PRIMARY HDUs.  That is, if EXTNAME
    is specified in the header, it will also be reflected in the ``.name``
    attribute and in ``fits.info()``.  These keywords used to be verboten in
    PRIMARY HDUs, but the latest version of the FITS standard allows them.
    Corresponds to PyFITS ticket 151.

  - HCOMPRESS can again be used to compress data cubes (and higher-dimensional
    arrays) so long as the tile size is effectively 2-dimensional. In fact,
    compatible tile sizes will automatically be used even if they're not
    explicitly specified. Corresponds to PyFITS ticket 171.

  - Fixed a bug that could cause a deadlock in the filesystem on OSX when
    reading the data from certain types of FITS files. This only occurred
    when used in conjunction with Numpy 1.7. [#369]

  - Added support for the optional ``endcard`` parameter in the
    ``Header.fromtextfile()`` and ``Header.totextfile()`` methods.  Although
    ``endcard=False`` was a reasonable default assumption, there are still text
    dumps of FITS headers that include the END card, so this should have been
    more flexible. Corresponds to PyFITS ticket 176.

  - Fixed a crash when running fitsdiff on two empty (that is, zero row) tables.
    Corresponds to PyFITS ticket 178.

  - Fixed an issue where opening a FITS file containing a random group HDU in
    update mode could result in an unnecessary rewriting of the file even if
    no changes were made. This corresponds to PyFITS ticket 179.

  - Fixed a crash when generating diff reports from diffs using the
    ``ignore_comments`` options. Corresponds to PyFITS ticket 181.

  - Fixed some bugs with WCS Paper IV record-valued keyword cards:

    - Cards that looked kind of like RVKCs but were not intended to be were
      over-permissively treated as such--commentary keywords like COMMENT and
      HISTORY were particularly affected. Corresponds to PyFITS ticket 183.

    - Looking up a card in a header by its standard FITS keyword only should
      always return the raw value of that card.  That way cards containing
      values that happen to valid RVKCs but were not intended to be will still
      be treated like normal cards. Corresponds to PyFITS ticket 184.

    - Looking up a RVKC in a header with only part of the field-specifier (for
      example "DP1.AXIS" instead of "DP1.AXIS.1") was implicitly treated as a
      wildcard lookup. Corresponds to PyFITS ticket 184.

  - Fixed a crash when diffing two FITS files where at least one contains a
    compressed image HDU which was not recognized as an image instead of a
    table. Corresponds to PyFITS ticket 187.

  - Fixed a bug where opening a file containing compressed image HDUs in
    'update' mode and then immediately closing it without making any changes
    caused the file to be rewritten unnecessarily.

  - Fixed two memory leaks that could occur when writing compressed image data,
    or in some cases when opening files containing compressed image HDUs in
    'update' mode.

  - Fixed a bug where ``ImageHDU.scale(option='old')`` wasn't working at
    all--it was not restoring the image to its original BSCALE and BZERO
    values.

  - Fixed a bug when writing out files containing zero-width table columns,
    where the TFIELDS keyword would be updated incorrectly, leaving the table
    largely unreadable.

  - Fixed a minor string formatting issue.

  - Fixed bugs in the backwards compatibility layer for the ``CardList.index``
    and ``CardList.count`` methods. Corresponds to PyFITS ticket 190.

  - Improved ``__repr__`` and text file representation of cards with long
    values that are split into CONTINUE cards. Corresponds to PyFITS ticket
    193.

  - Fixed a crash when trying to assign a long (> 72 character) value to blank
    ('') keywords. This also changed how blank keywords are represented--there
    are still exactly 8 spaces before any commentary content can begin; this
    *may* affect the exact display of header cards that assumed there could be
    fewer spaces in a blank keyword card before the content begins. However,
    the current approach is more in line with the requirements of the FITS
    standard. Corresponds to PyFITS ticket 194.

- ``astropy.io.votable``

  - The ``Table`` class now maintains a single array object which is a
    Numpy masked array.  For variable-length columns, the object that
    is stored there is also a Numpy masked array.

  - Changed the ``pedantic`` configuration option to be ``False`` by default
    due to the vast proliferation of non-compliant VO Tables. [#296]

  - Renamed ``astropy.io.vo`` to ``astropy.io.votable``.

- ``astropy.table``

  - Added a workaround for an upstream bug in Numpy 1.6.2 that could cause
    a maximum recursion depth RuntimeError when printing table rows. [#341]

- ``astropy.wcs``

  - Updated to wcslib 4.15 [#418]

  - Fixed a problem with handling FITS headers on locales that do not use
    dot as a decimal separator. This required an upstream fix to wcslib which
    is included in wcslib 4.14. [#313]

- Fixed some tests that could fail due to missing/incorrect logging
  configuration--ensures that tests don't have any impact on the default log
  location or contents. [#291]

- Various minor documentation fixes [#293 and others]

- Fixed a bug where running the tests with the ``py.test`` command still tried
  to replace the system-installed pytest with the one bundled with Astropy.
  [#454]

- Improved multiprocessing compatibility for file downloads. [#615]

- Fixed handling of Cython modules when building from a source checkout of a
  tagged release version. [#594]

- Added a workaround for a bug in Sphinx that could occur when using the
  ``:tocdepth:`` directive. [#595]

- Minor VOTable fixes [#596]

- Fixed how ``setup.py`` uses ``distribute_setup.py`` to prevent possible
  ``VersionConflict`` errors when an older version of distribute is already
  installed on the user's system. [#616][#640]

- Changed use of ``log.warn`` in the logging module to ``log.warning`` since
  the former is deprecated. [#624]


0.1 (2012-06-19)
----------------

- Initial release.
