Version 6.1.0 (2024-05-03)
==========================

New Features
------------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- ``BaseCoordinateFrame`` now has a ``position_angle()`` method, which is the
  same as the ``position_angle`` method of ``SkyCoord`` instances. [#15737]

- By default the ``SkyCoord`` and ``BaseCoordinateFrame`` ``separation()``
  methods now emit a warning if they have to perform a coordinate transformation
  that is not a pure rotation to inform the user that the angular separation can
  depend on the direction of the transformation.
  It is possible to modify this behaviour with the new optional keyword-only
  ``origin_mismatch`` argument.
  Specifying ``origin_mismatch="ignore"`` allows any transformation to
  succeed without warning, which has been the behaviour so far.
  ``origin_mismatch="error"`` forbids all transformations that are not
  pure rotations. [#16246]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Clearer error message in reading ASCII tables when there is
  a mismatch between converter type and column type. [#15991]

astropy.io.registry
^^^^^^^^^^^^^^^^^^^

- The module ``astropy.io.typing`` has been added to provide type annotations for
  I/O-related functionality. [#15916]

astropy.samp
^^^^^^^^^^^^

- SAMP web profile CORS HTTP server implements `Private Network Access proposal <https://wicg.github.io/private-network-access>`_. [#16193]

astropy.table
^^^^^^^^^^^^^

- ``Table`` now has a ``setdefault()`` method, analogous to
  ``dict.setdefault()``. [#16188]

astropy.units
^^^^^^^^^^^^^

- Added a new module ``astropy.units.typing`` that provides support for type annotations related to
  ``astropy.units``. [#15860]

- Added a new CGS unit Oersted. [#15962]

- Added "surface brightness", "surface brightness wav", "photon surface brightness", and "photon surface brightness wav" to recognized physical types. [#16032]

- Added magnetic helicity as a physical type. [#16101]

astropy.utils
^^^^^^^^^^^^^

- For gufuncs on ``Masked`` instances, add support for the ``axes`` argument. [#16121]

- ``Masked`` instances now support the various numpy array set operations, such
  as ``np.unique`` and ``np.isin``. [#16224]

astropy.wcs
^^^^^^^^^^^

- Added support for slicing WCS objects containing ``cpdis`` or ``det2im`` distortions, which previously were ignored. [#16163]


API Changes
-----------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- The ``astropy.coordinates.transformations`` module has been refactored into a module.
  There should be no user-visible changes, but if you notice any, please open an
  Issue. [#15895]

- Changed the default value of the ``copy`` argument in
  ``astropy.coordinates.representation.CylindricalDifferential.__init__`` from
  ``False`` to ``True``, which is the intended behaviour for all subclasses of
  ``astropy.coordinates.representation.BaseDifferential``. [#16198]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- ``Cosmology`` and its subclasses are now frozen ``dataclass`` objects. [#15484]

- The argument ``verbose`` in the function ``z_at_value`` is now keyword-only. [#15855]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- The ``io.ascii`` Python and C table readers were updated to use a 64-bit integer field by
  default when reading a column of integer numeric data. This changes the default behavior
  on Windows and potentially 32-bit architectures. Previously on those platforms, table
  columns with any long integers which overflowed the 32-bit integer would be returned
  as string columns. The new default behavior is consistent with ``numpy`` v2 and ``pandas``. [#16005]

- The parallel fast-reader parser for reading ASCII files has been removed.
  Since astropy v4.0.4 requesting this option has issued a warning that
  this option is broken and that the serial parser will be used.
  The ``parallel`` key in the ``fast_reader`` argument for reading
  ASCII tables is no longer available. [#16103]

astropy.table
^^^^^^^^^^^^^

- ``show_in_notebook`` is deprecated and it is recommended to use dedicated
  tools in the Jupyter ecosystem to create interactive plots in notebooks. [#15905]

- A warning is now emitted when ``Quantity`` values are inserted into empty ``Column`` objects
  via ``Table.insert_row`` or ``Table.add_row``. [#16038]

- ``show_in_browser`` is deprecated (pending feedback from the community).
  Please use https://github.com/astropy/astropy/issues/16067 if you are
  actively using the function. [#16068]

- ``TableColumns.setdefault()``  and ``TableColumns.update()`` methods (which
  would typically be called as ``Table.columns.setdefault()`` and
  ``Table.columns.update()``) have been deprecated because they can easily
  corrupt the ``Table`` instance the ``TableColumns`` instance is attached to.
  The ``Table.setdefault()`` and ``Table.update()`` methods are safe. [#16154]

astropy.time
^^^^^^^^^^^^

- ``TIME_FORMATS`` and ``TIME_DELTA_FORMATS`` in ``astropy.time.formats``
  are changed from ``OrderedDict`` to Python ``dict``. [#15491]

- A ``FutureWarning`` is now emitted when mutating ``Time.location`` post-initialization. [#16063]

- Following the removal of ``np.ndarray.ptp`` in Numpy v2, ``Time.ptp`` is now
  deprecated in favor of ``np.ptp``. [#16212]

astropy.units
^^^^^^^^^^^^^

- If any iterable such as a list of tuple was input to ``Quantity``, a check was
  done to see if they contained only quantities, and, if so, the quantities were
  concatenated.  This makes sense for list and tuple, but is not necessarily
  logical for all iterables and indeed was broken for those that do not have a
  length (such as ``array_api`` array instances). Hence, the check will now be
  done only for values where it makes sense, i.e., instances of list and tuple. [#15752]

- Units now exposes ``get_converter`` which returns a function that
  will convert a scalar or array from one unit to another. This can be
  useful to speed up code that converts many quantities with the same
  unit to another one, especially if the quantity has not many elements,
  so that the overhead of creating a conversion function is relatively large. [#16139]

astropy.utils
^^^^^^^^^^^^^

- Deprecate importing ``ErfaError`` and ``ErfaWarning`` from ``astropy.utils.exceptions``.
  They should be imported directly from ``erfa`` instead. [#15777]

- ``introspection.isinstancemethod()`` and ``introspection.find_mod_objs()`` are
  deprecated. [#15934]

- ``astropy.utils.console.terminal_size`` is now deprecated in favour of
  ``shutil.get_terminal_size`` from the standard library. [#16045]

- ``indent()`` is deprecated.
  Use ``textwrap.indent()`` from Python standard library instead. [#16223]

- Unmasked ``Masked`` scalar instances are now considered hashable, to match the
  implicit behaviour of regular arrays, where if an operation leads to a scalar,
  a hashable array scalar is returned. [#16224]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Renamed the ``min_cut`` and ``max_cut`` keywords in ``simple_norm`` and
  ``fits2bitmap`` to ``vmin`` and ``vmax``. The old names are deprecated. [#15621]

- If ``vmin == vmax``, the ``ImageNormalize`` class now maps the input
  data to 0. If ``vmin > vmax``, the ``ImageNormalize`` class now raises a
  ``ValueError``. [#15622]


Bug Fixes
---------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Avoid a segfault when calling ``astropy.convolution.convolve`` on an empty array.
  An exception is now raised instead. [#15840]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Previously passing a ``SkyCoord`` instance to the ``BaseCoordinateFrame``
  ``separation()`` or ``separation_3d()`` methods could produce wrong results,
  depending on what additional frame attributes were defined on the ``SkyCoord``,
  but now ``SkyCoord`` input can be used safely. [#15659]

- ``Distance`` now accepts as ``parallax`` any angle-like value.
  This includes types like ``Column`` which have a unit but are not ``Quantity`` subclasses. [#15712]

- The new default for the class method ``SkyCoord.from_name()``
  is to look for coordinates first in SIMBAD, then in NED, and then in VizieR,
  instead of having no specific order. [#16046]

- Fix ``Angle.to_string()`` for angles in degrees represented in 'hms' and angles in hours represented in 'dms'. [#16085]

- Fix a bug where ``SkyCoord.spherical_offsets_by`` would crash when a wrap
  was needed. [#16241]

- ``search_around_3d()`` now always raises a ``UnitConversionError`` if the units
  of the distances in ``coord1`` and ``coord2`` and the unit of ``distlimit`` do
  not agree.
  Previously the error was not raised if at least one of the coordinates was
  empty. [#16280]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Fixed a bug where the attribute ``ParametersAttribute.attr_name`` could be None
  instead of a string. [#15882]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Reading of CDS header files with multi-line descriptions where the continued line started with a number was broken. This is now fixed. [#15617]

- Ensure that the names of mixin columns are properly propagated as
  labels for the MRT format. [#15848]

- Fixed reading IPAC tables for ``long`` column type on some platforms, e.g., Windows. [#16005]

astropy.io.fits
^^^^^^^^^^^^^^^

- Avoid ``WinError 1455`` in opening some large files with memory
  mapping on windows. [#15388]

- Fix TDISP parsing for floating numbers. [#16007]

- Fix a crash when calling FITS ``writeto`` methods with stdout as the output stream. [#16008]

- Fix TDISP parsing for floating numbers in formats ES / EN. [#16015]

- Fix conversion of ``Table`` to ``BinTableHDU`` with ``character_as_bytes=True``. [#16358]

- Improved error message when instantiating a fits table with an ill-formed array. [#16363]

astropy.io.misc
^^^^^^^^^^^^^^^

- Reading an empty table stored in parquet format now creates an empty
  table instead of raising an unexpected error. [#16237]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- When reading a VOTable, if some user-requested columns were not present then the
  resulting error message previously listed all the requested column names.
  Now only columns that are actually missing are shown. [#15956]

astropy.stats
^^^^^^^^^^^^^

- Fix a spurious warning when calling ``sigma_clipped_stats`` on a ``MaskedColumn``. [#15844]

astropy.table
^^^^^^^^^^^^^

- Fix a Table bug when setting items (via slice or index list) in a ``bytes`` type
  ``MaskedColumn`` would cause the column mask to be set to all ``False``. A common way to
  trigger this bug was reading a FITS file with masked string data and then sorting the
  table. [#15669]

- Fix slicing logic for Row.
  Previously, slicing a ``astropy.table.row.Row`` object would incorrectly return a column,
  now it correctly returns a list of values from that row. [#15733]

- Fix a ``ValueError`` raised by ``table.join`` when fed with large tables.
  This would typically happen in situations when the result joined table would be
  too large to fit in memory. In those situations, the error message is now much more
  clearly about the necessary memory size. [#15734]

- Fix an unintended exception being raised when attempting to compare two unequal ``Table`` instances. [#15845]

- Ensure that if a ``Column`` is initialized with a ``Quantity`` it will use by
  default a possible name defined on the quantity's ``.info``. [#15848]

- Fix a bug where columns with ``dtype=object`` wouldn't be properly deep-copied using ``copy.deepcopy``. [#15871]

- Fix ``hasattr(Table, "iloc")`` raising an exception, preventing use of tables e.g. with scikit-learn. [#15913]

- Calling ``Table.group_by`` on an empty table no longer raises an exception. [#16093]

- The unit conversion ``convert_unit_to`` with MaskedColumn was
  broken as it was storing the old unit in a dictionary attached
  to underlying np.ma.MaskedArray. This fixes it by overwriting
  the old unit after unit conversion. [#16118]

- ``astropy.table.vstack`` will no longer modify the input list even when it
  contains non-Table objects like ``astropy.table.Row``. [#16130]

- Update old dataTables.js version.
  This should not affect the end user. [#16315]

astropy.time
^^^^^^^^^^^^

- Fix comparing NaN ``Quantity`` with ``TimeDelta`` object. [#15830]

- Scalar ``Time`` instances are now hashable if they are not masked, also if one
  uses ``Masked`` internally, matching the behaviour prior to astropy 6.0 (and
  the current behaviour when masking using ``np.ma.MaskedArray``). [#16224]

astropy.units
^^^^^^^^^^^^^

- Fix rare signature incompatibilities between helper and helped array functions.
  Most involve cases where the corresponding numpy function has had its
  arguments renamed between numpy versions. Since all those generally changed
  the first arguments, which are typically passed as positional arguments,
  this should not affect user code.
  Affected functions:
  - ``numpy.array_str``
  - ``numpy.choose``
  - ``numpy.convolve``
  - ``numpy.correlate``
  - ``numpy.histogram``
  - ``numpy.histogramdd``
  - ``numpy.histogram2d``
  - ``numpy.isin``
  - ``numpy.inner``
  - ``numpy.nanmedian``
  - ``numpy.unique``
  - ``numpy.matrix_rank``
  - ``numpy.unwrap``
  - ``numpy.vdot``
  - ``numpy.lib.recfunctions.unstructured_to_structured`` [#15710]

- Fix an issue with unicode string representations of units shown as
  superscripts (like degree) when raised to some power. Like for
  LaTeX representations, now the superscript unicode character is
  replaced by the literal short name before adding the power. [#15755]

- Fix a missing ``Sun`` unit in the list of VOUnits simple_units. [#15832]

- Fix an unhelpful ``TypeError`` when attempting truediv, ``lshift`` (``<<``) or ``mul`` (``*``) or ``truediv`` (``/``) with a ``Unit`` for right operand and a numpy array with non-numerical dtype for left operand. [#15883]

- Fix write/read roundtrips with empty ``Table`` dumped to ECSV. [#15885]

- Fix a bug where LaTeX formatter would return empty strings for unity (1) input. [#15923]

- Fix extraneous space in LaTeX repr for ``Quantity`` objects with superscript
  units (e.g. angles or temperatures in degree Celsius). [#16043]

- Ensure powers of units are consistently as simple as possible. So, an
  integer if possible, otherwise a float, or a fraction if the float is
  really close to that. This also ensures the hash of a unit is unique
  for any given unit (previously, the same power could be represented as
  float, int or fraction, which made the hash different). [#16058]

- Ensure that ``find_equivalent_units`` only returns actual units, not units
  that raised to some power match the requested one.  With this fix,
  ``(u.m**-3).find_equivalent_units()`` properly finds nothing, rather than all
  units of length. [#16127]

- Using a dimensionless ``Quantity`` as an exponent works anew.
  In astropy 6.0.1 an exception was erroneously raised. [#16261]

astropy.utils
^^^^^^^^^^^^^

- Fix rare signature incompatibilities between helper and helped array functions.
  These typically cover corner cases and should not affect user code.
  Some arguments weren't being re-exposed correctly or at all, depending on
  numpy's version.
  Affected functions:
  - ``numpy.broadcast_arrays``
  - ``numpy.median``
  - ``numpy.quantile``
  - ``numpy.empty_like``
  - ``numpy.ones_like``
  - ``numpy.zeros_like``
  - ``numpy.full_like`` [#16025]

- Fix a bug where ``astropy.utils.console.Spinner`` would leak newlines for
  messages longer than terminal width. [#16040]

- Update ``report_diff_values`` so the diff no longer depends on the
  console terminal size. [#16065]

- Fix support in ``Masked`` for generalized ufuncs with more than a
  single core dimension (such as ``erfa.rxp``). [#16120]

- ``Masked`` array instances now deal more properly with structured dtypes,
  combining field masks to get element masks for generalized ufuncs, and
  allowing ``.view()`` any time the mask can be viewed as well. This allows a
  larger number of ``erfa`` routines to work with masked data. [#16125]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- ``WCSAxes`` will correctly set certain defaults when ``wcs.world_axis_physical_types`` contains ``custom:`` prefixes. [#15626]

- Fix an edge case where ``quantity_support`` would produce duplicate tick labels for small data ranges. [#15841]

- Fix a bug where ``AngleFormatterLocator`` and ``ScalarFormatterLocator`` wouldn't respect matplotlib.rc's ``axes.unicode_minus`` parameter. [#15902]

- Fixed a bug in ``CoordinateHelper.grid`` method to properly handle ``draw_grid=False`` and ``draw_grid=None``,
  ensuring grid lines are controlled correctly even when not explicitly called. [#15985]

astropy.wcs
^^^^^^^^^^^

- Updated bundled WCSLIB version to 8.2.2. This update fixes character buffer
  overflows in the comment string for the longitude and latitude axes triggered
  by some projections in ``wcshdo()``, and also the formatting for generic
  coordinate systems. For a full list of changes - see
  http://www.atnf.csiro.au/people/mcalabre/WCS/CHANGES or
  ``astropy/cextern/wcslib/CHANGES`` [#15795]

- Fixed a bug in ``fit_wcs_from_points`` that does not set the default value of the ``cdelt`` of the returned WCS object. [#16027]

- Fixed a bug in ``DistortionLookupTable`` (which implements ``cpdis`` and ``det2im`` projection corrections to a WCS) in which image pixels received an incorrect distortion value, from a location in the lookup table incorrectly offset by about 1 table pixel. [#16163]


Other Changes and Additions
---------------------------

- Update minimum supported Python version to 3.10 [#15603]

- The minimum required NumPy version is now 1.23 and the minimum required SciPy version is 1.8. [#15706]

- Fix loading parser tabs on pyc-only installations.

  Fix a bug in the wrappers for the lex and yacc wrappers that are
  used for parsing Astropy units so that they work on pyc-only
  installations.

  According to the Python module loading
  `flow chart <https://peps.python.org/pep-3147/#flow-chart>`_, when evaluating
  ``import foo`` and ``foo.py`` is not found, Python then reads ``foo.pyc``.

  One can take advantage of this fact to strip source files and leave only Python
  bytecode files for deployment inspace-constrained execution environments such
  as AWS Lambda. Astropy is now compatible with pyc-only deployments. [#16159]

- Change the default value of ``copy`` arguments in public APIs from ``False`` to
  ``None`` if Numpy 2.0 or newer is installed.
  For details, see the "Copy semantics" section on the What's New page for Astropy 6.1 . [#16181]

- astropy is now compiled against NumPy 2.0, enabling runtime compatibility
  with this new major release. Compatibility with NumPy 1.23 and newer
  versions of NumPy 1.x is preserved through this change. [#16252]

Version 6.0.1 (2024-03-25)
==========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Previously passing a ``SkyCoord`` instance to the ``BaseCoordinateFrame``
  ``separation()`` or ``separation_3d()`` methods could produce wrong results,
  depending on what additional frame attributes were defined on the ``SkyCoord``,
  but now ``SkyCoord`` input can be used safely. [#15659]

- ``Distance`` now accepts as ``parallax`` any angle-like value.
  This includes types like ``Column`` which have a unit but are not ``Quantity`` subclasses. [#15712]

- The new default for the class method ``SkyCoord.from_name()``
  is to look for coordinates first in SIMBAD, then in NED, and then in VizieR,
  instead of having no specific order. [#16046]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Reading of CDS header files with multi-line descriptions where the continued line started with a number was broken. This is now fixed. [#15617]

- Ensure that the names of mixin columns are properly propagated as
  labels for the MRT format. [#15848]

- Fixed reading IPAC tables for ``long`` column type on some platforms, e.g., Windows. [#15992]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fix TDISP parsing for floating numbers. [#16007]

- Fix a crash when calling FITS ``writeto`` methods with stdout as the output stream. [#16008]

- Fix TDISP parsing for floating numbers in formats ES / EN. [#16015]

astropy.stats
^^^^^^^^^^^^^

- Fix a spurious warning when calling ``sigma_clipped_stats`` on a ``MaskedColumn``. [#15844]

astropy.table
^^^^^^^^^^^^^

- Fix a Table bug when setting items (via slice or index list) in a ``bytes`` type
  ``MaskedColumn`` would cause the column mask to be set to all ``False``. A common way to
  trigger this bug was reading a FITS file with masked string data and then sorting the
  table. [#15669]

- Fix slicing logic for Row.
  Previously, slicing a ``astropy.table.row.Row`` object would incorrectly return a column,
  now it correctly returns a list of values from that row. [#15733]

- Fix a ``ValueError`` raised by ``table.join`` when fed with large tables.
  This would typically happen in situations when the result joined table would be
  too large to fit in memory. In those situations, the error message is now much more
  clearly about the necessary memory size. [#15734]

- Fix an unintended exception being raised when attempting to compare two unequal ``Table`` instances. [#15845]

- Ensure that if a ``Column`` is initialized with a ``Quantity`` it will use by
  default a possible name defined on the quantity's ``.info``. [#15848]

- The unit conversion ``convert_unit_to`` with MaskedColumn was
  broken as it was storing the old unit in a dictionary attached
  to underlying np.ma.MaskedArray. This fixes it by overwriting
  the old unit after unit conversion. [#16118]

- ``astropy.table.vstack`` will no longer modify the input list even when it
  contains non-Table objects like ``astropy.table.Row``. [#16130]

astropy.units
^^^^^^^^^^^^^

- Fix an issue with unicode string representations of units shown as
  superscripts (like degree) when raised to some power. Like for
  LaTeX representations, now the superscript unicode character is
  replaced by the literal short name before adding the power. [#15755]

- Fix a missing ``Sun`` unit in the list of VOUnits simple_units. [#15832]

- Fix write/read roundtrips with empty ``Table`` dumped to ECSV. [#15885]

- Fix a bug where LaTeX formatter would return empty strings for unity (1) input. [#15923]

- Ensure powers of units are consistently as simple as possible. So, an
  integer if possible, otherwise a float, or a fraction if the float is
  really close to that. This also ensures the hash of a unit is unique
  for any given unit (previously, the same power could be represented as
  float, int or fraction, which made the hash different). [#16058]

- Ensure that ``find_equivalent_units`` only returns actual units, not units
  that raised to some power match the requested one.  With this fix,
  ``(u.m**-3).find_equivalent_units()`` properly finds nothing, rather than all
  units of length. [#16127]

astropy.utils
^^^^^^^^^^^^^

- Fix a bug where ``astropy.utils.console.Spinner`` would leak newlines for
  messages longer than terminal width. [#16040]

- Update ``report_diff_values`` so the diff no longer depends on the
  console terminal size. [#16065]

- Fix support in ``Masked`` for generalized ufuncs with more than a
  single core dimension (such as ``erfa.rxp``). [#16120]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Fix an edge case where ``quantity_support`` would produce duplicate tick labels for small data ranges. [#15841]

astropy.wcs
^^^^^^^^^^^

- Updated bundled WCSLIB version to 8.2.2. This update fixes character buffer
  overflows in the comment string for the longitude and latitude axes triggered
  by some projections in ``wcshdo()``, and also the formatting for generic
  coordinate systems. For a full list of changes - see
  http://www.atnf.csiro.au/people/mcalabre/WCS/CHANGES or
  ``astropy/cextern/wcslib/CHANGES`` [#15795]

- Fixed a bug in ``fit_wcs_from_points`` that does not set the default value of the ``cdelt`` of the returned WCS object. [#16027]

Other Changes and Additions
---------------------------

- Given the potential breaking changes with the upcoming Numpy 2.0 release,
  this release pins Numpy<2.0 and support for Numpy 2.0 will be added in the
  v6.1.0 release.

Version 6.0.0 (2023-11-25)
==========================

New Features
------------

astropy.config
^^^^^^^^^^^^^^

- The new ``ConfigNamespace.help()`` method provides a convenient way to get
  information about configuration items. [#13499]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Support has been added to create geodetic representations not just for existing ellipsoids
  from ERFA, but also with explicitly provided values, by defining a subclass of
  ``BaseGeodeticRepresentation`` with the equatorial radius and flattening assigned to
  ``_equatorial_radius`` and ``_flattening`` attributes. [#14763]

- Add ``BaseBodycentricRepresentation``, a new spheroidal representation for bodycentric
  latitudes and longitudes. [#14851]

- Support Numpy broadcasting over frame data and attributes. [#15121]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Registered a ``latex`` writer for exporting a Cosmology object to a LaTex table. [#14701]

- Added argument ``rename`` to Cosmology's I/O, allowing for input and output symbols to
  be renamed. [#14780]

- All non-abstract Cosmology subclasses are now automatically registered to work with
  Astropy's YAML serialization. [#14979]

- Cosmology I/O now auto-identifies the '.tex' suffix with the 'ascii.latex' format. [#15088]

- The ``Cosmology`` class now has a new property to access the parameters of the
  cosmology: ``.parameters``. This property return a read-only dictionary of all the
  non-derived parameter values on the cosmology object. When accessed from the class (not
  an instance) the dictionary contains ``Parameter`` instances, not the values. [#15168]

- The field ``default`` has been added to ``Parameter``. This can be used to introspect
  the default value of a parameter on a cosmology class e.g. ``LambdaCDM.H0.default``. [#15400]

astropy.io.fits
^^^^^^^^^^^^^^^

- Add new option ``decompress_in_memory`` to ``fits.open``, to decompress the
  whole file in memory at once, instead of decompressing the file progressively
  as data is needed.  Default behavior is better for memory usage but sometimes
  slow, especially for files with many small HDUs. [#15501]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Add support for Parquet serialization of VOTables. Writing of this
  serialization is available with using the new ``'votable.parquet'`` format. [#15281]

- Added MIVOT feature through the ``MivotBlock`` class
  that allows model annotations reading and writing in VOTable. [#15390]

astropy.modeling
^^^^^^^^^^^^^^^^

- Added a ``GeneralSersic2D`` model that can have "boxy" or "disky"
  isophotes. [#15545]

astropy.nddata
^^^^^^^^^^^^^^

- A more flexible and/or compact string representation is now available for
  ``NDData`` objects which visually indicates masked entries, and provides for
  better for dask array support. [#14438]

astropy.table
^^^^^^^^^^^^^

- The new ``Row.get()`` method, analogous to ``dict.get()``, returns the value of
  the specified column from the row if the column present, otherwise it returns a
  fallback value, which by default is ``None``. [#14878]

astropy.time
^^^^^^^^^^^^

- Masked ``Time`` instances now use astropy's own ``Masked`` class internally.
  This means that ``Masked`` input is now properly recognized, and that masks
  get propagated also to ``Quantity`` output (such as from a ``TimeDelta``
  converted to a unit of time), creating ``MaskedQuantity`` instances. [#15231]

- Added a ``TimeDelta`` format ``quantity_str`` that represents the time delta as a string
  with one or more ``Quantity`` components. This format provides a human-readable
  multi-scale string representation of a time delta. The default output sub-format is not
  considered stable in this release, please see https://github.com/astropy/astropy/issues/15485
  for more information. [#15264]

astropy.uncertainty
^^^^^^^^^^^^^^^^^^^

- Uncertainty ``Distribution`` now support structured data types, and as
  a result it now works also with ``EarthLocation``. [#15304]

- Uncertainty ``Distribution`` can now be used inside representations, which
  also allows basic support in ``SkyCoord``. While most calculations work, there
  are remaining issues.  For instance, the ``repr`` does not show that the
  coordinates are distributions. [#15395]

astropy.units
^^^^^^^^^^^^^

- Add support for gc2gde and gd2gce erfa functions to allow geodetic representations
  using equatorial radius and flattening. [#14729]

astropy.utils
^^^^^^^^^^^^^

- The ``astropy.utils.metadata.MetaData`` default dictionary can now be
  set with the ``default_factory`` keyword argument. [#15265]

- ``astropy.utils.decorators.deprecated`` now adds the ``__deprecated__`` attribute to
  the objects it wraps, following the practice in https://peps.python.org/pep-0702/. [#15310]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Add ``WCSAxes.text_coord`` method to print text using ``SkyCoord`` objects
  parallel to plotting data points with ``WCSAxes.plot_coord``. [#14661]

astropy.wcs
^^^^^^^^^^^

- Support WCS descriptions of basic planetary coordinate frames. [#14820]

- Updated bundled WCSLIB version to 8.1. This update adds support planetary keywords ``A_RADIUS``, ``B_RADIUS``, ``C_RADIUS``, ``BLON_OBS``, ``BLAT_OBS``, and ``BDIS_OBS`` in ``auxprm`` and adds ``wcsprm::time`` to the ``wcsprm`` struct to record the ``TIME`` axis. This update also includes several bug fixes. For a full list of changes - see http://www.atnf.csiro.au/people/mcalabre/WCS/CHANGES [#15035]


API Changes
-----------

astropy.config
^^^^^^^^^^^^^^

- Removed deprecated ``ConfigurationMissingWarning`` class and ``update_default_config`` function;
  There are no replacements as they should no be used anymore. [#15466]

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Invalid kernel arithmetic operations now raise a ``KernelArithmeticError`` instead of a
  bare ``Exception``. [#14728]

- Added base ``KernelError`` error class and removed ``DiscretizationError`` error class (a ``ValueError`` will be raised instead). [#14732]

- ``discretize_model`` will now raise a ``ValueError`` if
  ``mode='oversample'`` and ``factor`` does not have an integer value. [#14794]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Removed deprecated angle parsing and formatting utilities from ``angle_utilities``.
  Use the functions from ``angle_formats`` instead. [#14675]

- The deprecated functionality of initializing ``Angle`` or ``Longitude`` from a
  ``tuple`` is no longer supported. [#15205]

- Angle-related classes and functions have been moved within ``astropy.coordinates``.
  There is no change to public API as everything moved should still be imported from
  ``astropy.coordinates``, not a sub-module. If you are using private API, try importing
  from ``astropy.coordinates`` instead. If you need something that has been moved and is
  not available in ``astropy.coordinates``, please open an issue on the Astropy issue
  tracker. [#15220]

- It is no longer possible to pass frame classes to the ``transform_to()`` method
  of a low-level coordinate-frame class. It is still possible to pass frame
  instances. The ``transform_to()`` method of the high-level ``SkyCoord`` class
  is unaffected. [#15500]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Removed support of importing private constants and functions from ``astropy.cosmology.flrw``. [#14672]

- Removed deprecated Cosmology Parameter argument ``fmt``. [#14673]

- Removed deprecated ``vectorize_if_needed`` and ``inf_like`` from ``cosmology.utils``. [#14677]

- Removed deprecated import paths from ``astropy.cosmology.core``. [#14782]

- Cosmology ``Parameter`` is now a ``dataclass``, and can work with all of Python's dataclasses
  machinery, like field introspection and type conversion. [#14874]

- A new property -- ``scale_factor0`` -- has been added to Cosmology objects.
  This is the scale factor at redshift 0, and is defined to be 1.0. [#14931]

- Added registration label ``ascii.latex`` to Cosmology IO. [#14938]

- The private module ``astropy.cosmology.utils`` has been deprecated. [#14980]

- Removed deprecated ``get_cosmology_from_string`` class method in ``default_cosmology``; use ``get`` instead. [#15467]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Several arguments in functions within ``astropy.io.ascii`` have been deprecated and
  are either renamed or scheduled to be removed.

  ``read()``:
  - ``Reader`` will be removed. Instead supply the equivalent ``format`` argument.
  - ``Inputter`` has been renamed to ``inputter_cls``.
  - ``Outputter`` has been renamed to ``outputter_cls``.

  ``get_reader()``:
  - ``Reader`` has been renamed to ``reader_cls``.
  - ``Inputter`` has been renamed to ``inputter_cls``.
  - ``Outputter`` has been renamed to ``outputter_cls``.

  ``write()``:
  - ``Writer`` will be removed. Instead supply the equivalent ``format`` argument.

  ``get_writer()``:
  - ``Writer`` has been renamed to ``writer_cls``. [#14914]

- Removed deprecated ``astropy.io.ascii.tests.common.raises`` test helper; use ``pytest.raises`` instead. [#15470]

astropy.io.fits
^^^^^^^^^^^^^^^

- Deprecate ``_ExtensionHDU`` and ``_NonstandardExtHDU`` (use ``ExtensionHDU`` or
  ``NonstandardExtHDU`` instead). [#15396]

- Remove special handling of TCTYP TCUNI TCRPX TCRVL TCDLT TRPOS (#7157). [#15396]

- Rename and deprecate ``TableHDU.update`` to ``TableHDU.update_header``, for
  consistency with ``ImageHDU``. [#15396]

astropy.io.misc
^^^^^^^^^^^^^^^

- Removed deprecated ``astropy.io.misc.asdf`` subpackage. Use ``asdf-astropy`` package instead. [#14668]

- ``fnunpickle`` and ``fnpickle`` are deprecated because they are not used anywhere within ``astropy``.
  If you must, use the module from Python standard library but be advised that pickle is insecure
  so you should only unpickle data that you trust. [#15418]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Removed deprecated ``pedantic`` option from the
  ``astropy.io.votable.table.parse()`` function and the corresponding configuration
  setting. Use the ``verify`` option instead. [#14669]

- Class ``astropy.io.votable.tree.Table`` has been renamed to ``TableElement``
  to avoid sharing the name with ``astropy.table.Table``. [#15372]

- Fully removed support for version = '1.0' on ``VOTableFile__init__()`` and changed its tests to check correctly.
  It was raising a ``DeprecationWarning`` and now is raising a ``ValueError``. [#15490]

astropy.modeling
^^^^^^^^^^^^^^^^

- Removed the ``AliasDict`` class from ``modeling.utils``. [#12943]

- Creating a model instance with parameters that have incompatible shapes will
  now raise a ``ValueError`` rather than an ``IncompatibleShapeError``. [#15209]

- Removal of deprecated code ``_model_to_fit_params`` and ``_fitter_to_model_params`` from ``fitting.py``. [#15461]

astropy.stats
^^^^^^^^^^^^^

- The ``BoxLeastSquares``, ``BoxLeastSquaresResults`` and ``LombScargle`` classes
  are not available from ``astropy.stats`` anymore, they are now available only
  from ``astropy.timeseries``. [#15530]

astropy.tests
^^^^^^^^^^^^^

- Removed deprecated deprecation, warning, and exception handling functionality provided by ``astropy.tests.helper``. [#14670]

- ``astropy.tests.command.FixRemoteDataOption`` and ``astropy.tests.command.AstropyTest`` are deprecated.
  They are no longer necessary after sunsetting ``astropy-helpers``. [#15204]

astropy.time
^^^^^^^^^^^^

- ``Time`` has switched to use ``Masked`` arrays internally, instead of
  indicating masked values using NaN in the internal ``jd2`` attribute.  As a
  result, any output from instances, such as one gets with, say, the ``.isot``
  format, will also use ``Masked`` by default.

  For backwards compatibility, a new configuration item,
  ``astropy.time.conf.masked_array_type`` is introduced which is set to
  "astropy" by default (which indicates one wants to use ``Masked``), but can
  also be set to "numpy", in which case ``numpy.ma.MaskedArray`` will be used
  where possible (essentially, for all but ``Quantity``). [#15231]

- Changed the ``TimeDelta`` init signature to be consistent with that of ``Time``.
  Previously the argument order was ``val, val2, format, scale, copy``. Now the order is
  ``val, val2, format, scale, *, precision, in_subfmt, out_subfmt, copy``, where the
  arguments after the ``*`` must be specified by keyword. [#15264]

astropy.timeseries
^^^^^^^^^^^^^^^^^^

- Removed deprecated ``midpoint_epoch`` in ``fold`` function; use ``epoch_time`` instead. [#15462]

astropy.uncertainty
^^^^^^^^^^^^^^^^^^^

- The ``.dtype`` attribute exposed by ``Distribution`` is now that of
  the samples, rather than one that has a "samples" entry.  This makes
  quantities with structured data types and units easier to support, and
  generally makes the ``Distribution`` appear more similar to regular
  arrays.  It should have little effect on code.  For instance,
  ``distribution["samples"]`` still will return the actual distribution.

  As a consequence of this refactoring, most arrays that are not
  C-contiguous can now be viewed and will thus not be copied on input
  any more.  The only exceptions are arrays for which the strides are
  negative.

  Note that the true data type is considered an implementation detail.
  But for reference, it now is a structured data type with a single
  field, "samples", which itself is an array of "sample" fields, which
  contain the actual data. [#15304]

astropy.units
^^^^^^^^^^^^^

- Like ``np.ndarray``, under numpy 2.0 ``Quantity`` and all its subclasses
  (``Angle``, ``Masked``, etc.) will no longer support the ``.ptp()`` method.
  Use ``np.ptp(...)`` instead.

  Similarly, support for the much less frequently used ``.newbyteorder()`` and
  ``.itemset()`` methods has been removed. [#15378]

- The following deprecated functionality has been removed:

    * ``littleh`` unit and ``with_H0`` equivalency. They are still available from
      ``cosmology.units``.
    * ``brightness_temperature`` equivalency no longer automatically swaps the
      order of its arguments if it does not match the expectation.
    * ``PhysicalType`` no longer supports ``str`` methods and attributes. [#15514]

astropy.utils
^^^^^^^^^^^^^

- Removed deprecated ``OrderedDescriptor``, ``OrderedDescriptorContainer``, and ``set_locale`` in ``astropy.utils.misc``. [#14679]

- ``is_path_hidden()`` and ``walk_skip_hidden()`` are deprecated. [#14759]

- The structure of ``utils.metadata`` has been refactored, but all the available
  functions and classes are still present and should be imported as before. [#15166]

- The ``astropy.utils.metadata.MetaData`` class, which is used throughout astropy
  to carry metadata on tables, columns, etc., can now also be used on dataclasses.

  When accessing the meta attribute on a class ``astropy.utils.metadata.MetaData``
  now returns None instead of itself. [#15237]

- The ``astropy.utils.metadata.MetaData`` class, which is used throughout astropy
  to carry metadata on tables, columns, etc., can now also be used on frozen dataclasses. [#15404]

- Removed deprecated ``version_path`` in ``minversion`` function; it is no longer used. [#15468]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- The ``bboxes``, ``ticklabels_bbox``, and ``tick_out_size`` arguments to ``astropy.visualization.wcaxes.ticklabels.TickLabels.draw()`` now have no effect and are deprecated.
  This is to allow rasterized ticks to be drawn correctly on WCSAxes. [#14760]

- It is now not possible to pass any keyword arguments to ``astropy.visualization.wcsaxes.WCSAxes.draw()``.
  Previously passing any keyword arguments would have errored anyway, as ``matplotlib.axes.Axes.draw()`` does not accept keyword arguments. [#14772]

- Deprecated the ``exp`` attribute in the ``LogStretch``,
  ``InvertedLogStretch``, ``PowerDistStretch``, and
  ``InvertedPowerDistStretch`` stretch classes, and the ``power``
  attribute in the ``PowerStretch``. Instead, use the ``a`` attribute,
  which matches the input keyword. [#15538]

- Removed the maximum value of the ``a`` parameter in the ``AsinhStretch``
  and ``SinhStretch`` stretch classes. [#15539]

astropy.wcs
^^^^^^^^^^^

- Removed deprecated ``accuracy`` from ``all_world2pix`` method in ``WCS``; use ``tolerance`` instead. [#15464]

- ``NoConvergence`` no longer accepts arbitrary keyword arguments. [#15504]


Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Fixed minor bug when getting solar system positions of objects from Type 3 SPICE kernel files. [#15612]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- The exponent in ``w0wzCDM.de_density_scale`` has been corrected to 3, from -3.
  This correction has also been made to the scalar ``inv_efunc`` cpython functions. [#14991]

- ``pandas.Series`` are now uniformly converted to their underlying data type when given
  as an argument to a Cosmology method. [#15600]

astropy.io.fits
^^^^^^^^^^^^^^^

- Reading a table from FITS now respects the TNULL property of a column, passing
  it into the column's ``fill_value``. [#14723]

- Fix crash when a PrimaryHDU has a GROUPS keyword with a non-boolean value (i.e.
  not a random-groups HDU). [#14998]

- Fixed a bug that caused ``Cutout2D`` to not work correctly with ``CompImageHDU.section`` [#14999]

- Fixed a bug that caused compressed images with TFORM missing the optional '1' prefix to not be readable. [#15001]

- Ensure that tables written to FITS with both masked and unmasked columns
  roundtrip properly (previously, all integer columns would become masked
  if any column was masked). [#15473]

- Fix segfault with error report in tile decompression. [#15489]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Output of ``repr`` for VOTable instance now clearly shows it is a VOTable and not generic astropy Table. [#14702]

astropy.modeling
^^^^^^^^^^^^^^^^

- All models can be pickled now. [#14902]

astropy.nddata
^^^^^^^^^^^^^^

- Restore bitmask propagation behavior in ``NDData.mask``, plus a fix
  for arithmetic between masked and unmasked ``NDData`` objects. [#14995]

astropy.table
^^^^^^^^^^^^^

- ``Table.as_array`` now respects the ``fill_value`` property of masked columns. [#14723]

- Fix a bug where table indexes were not using a stable sort order. This was causing the
  order of rows within groups to not match the original table order when an indexed table
  was grouped. [#14907]

- Fixed issue #14964 that when grouping a Table on a mixin column such as ``Quantity`` or
  ``Time``, the grouped table keys did not reflect the original column values. For
  ``Quantity`` this meant that the key values were pure float values without the unit,
  while for ``Time`` the key values were the pair of ``jd1`` and ``jd2`` float values. [#14966]

astropy.time
^^^^^^^^^^^^

- Ensure that the ``Time`` caches of formats and scales do not get out
  of sync with the actual data, even if another instance, holding a view
  of the data is written to.  E.g., if one does ``t01 = t[:2]``, and
  sets ``t[0]`` after, it is now guaranteed that ``t01.value`` will
  correctly reflect that change in value. [#15453]

astropy.units
^^^^^^^^^^^^^

- In VOunits, "pix", "au", "a", and "ct" are removed from the list of deprecated units. [#14885]

astropy.utils
^^^^^^^^^^^^^

- Ufuncs with more than 2 operands (such as ``erfa.dtf2d``) now work
  also if all inputs are scalars and more than two inputs have masks. [#15450]

- Ensured that ``str(masked_array)`` looks like ``str(unmasked_array)`` also for
  array scalars. Thus, like regular array scalars, the precision is ignored for
  float, and strings do not include extra quoting. [#15451]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- The location of ticklabels on a WCSAxes is now correctly calculated when the figure is rasterized. [#14760]

- Fixed a bug where a ``ValueError`` would be raised in the
  ``AsinhStretch`` and ``SinhStretch`` classes for valid ``a`` parameter
  values. [#15539]

astropy.wcs
^^^^^^^^^^^

- ``wcs.validate(filename)`` now properly closes the file handler. [#15054]

- Fix a regression in custom WCS mapping due to the recent introduction of
  Solar System frames. [#15630]


Other Changes and Additions
---------------------------

- The minimum supported version of NumPy is now 1.22. [#15006]

- Moved International Earth Rotation and Reference Systems (IERS) and Leap Second
  files out into standalone astropy-iers-data package, maintaining full
  backward-compatibility in the ``astropy.utils.iers`` API. Deprecation
  warnings may be issued when certain files are accessed directly. [#14819]

- Switch from using ``setup.cfg`` for project configuration to using ``pyproject.toml``. [#15247]

- Update bundled expat to 2.5.0. [#15585]

Version 5.3.4 (2023-10-03)
==========================

Bug Fixes
---------

astropy.io.misc
^^^^^^^^^^^^^^^

- Updated ``astropy.io.misc.yaml`` so ``dump()` with a numpy object array or
  ``load()`` with YAML representing a Numpy object array both raise
  ``TypeError``. This prevents problems like a segmentation fault. [#15373]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Fixed a bug in ``convert_to_writable_filelike`` where ``GzipFile`` was not
  closed properly. [#15359]

astropy.units
^^^^^^^^^^^^^

- In VOUnit, the spaces around the slash were removed in the formatting of
  fractions, and fractional powers now also use the "**" operator. [#15282]

- We now ensure that the unit ``u.cgs.cm`` is just an alias of ``u.si.cm``,
  instead of a redefinition.  This ensures that ``u.Unit("cm") / u.cm``
  will reliably cancel to dimensionless (instead of some "cm / cm"). [#15368]

astropy.utils
^^^^^^^^^^^^^

- For ``Masked``, ``np.ptp`` and the ``.ptp()`` method now properly account for
  the mask, ensuring the result is identical to subtracting the maximum and
  minimum (with the same arguments). [#15380]

Other Changes and Additions
---------------------------

- Compatibility with Python 3.12. [#14784]

- Replaced the URL of ``IETF_LEAP_SECOND_URL`` because the original is now
  defunct and IETF now defers to IANA for such look-up. [#15421]


Version v5.3.3 (2023-09-07)
===========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- ``TransformGraph.to_dot_graph()`` now throws an exception for invalid ``savelayout``.

astropy.cosmology
^^^^^^^^^^^^^^^^^

- The exponent of ``w0wzCDM`` functions in ``inv_efunc`` has been corrected to 3, from -3. [#15224]

astropy.modeling
^^^^^^^^^^^^^^^^

- Astropy modeling can filter non-finite data values using the ``filter_non_finite``
  keyword argument in a fitter call. Now when ``filter_non_finite`` is True,
  non-finite *weights* will also be filtered to prevent crashes in ``LevMarLSQFitter``. [#15215]

astropy.units
^^^^^^^^^^^^^

- Fixed ``astropy.units.Quantity``'s implementation of ``numpy.nanmedian()``,
  where for Numpy >= 1.25 an exception was raised for some array shapes and axis
  combinations. [#15228]


Other Changes and Additions
---------------------------

- v5.3.x will not support NumPy 2.0 or later. [#15234]


Version 5.3.2 (2023-08-11)
==========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Fixed import when called with Python ``-OO`` flag. [#15037]

astropy.nddata
^^^^^^^^^^^^^^

- Fix for collapse operations on ``NDData`` without masks or units. [#15082]

astropy.units
^^^^^^^^^^^^^

- Modified the implementation of ``np.power()`` for instances of ``Quantity`` to
  allow any array as the second operand if all its elements have the same value. [#15101]

Version 5.3.1 (2023-07-06)
==========================

Bug Fixes
---------

astropy.cosmology
^^^^^^^^^^^^^^^^^

- The exponent in ``wowzCDM.de_density_scale`` has been corrected to 3, from -3. [#14991]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fix crash when a PrimaryHDU has a GROUPS keyword with a non-boolean value (i.e.
  not a random-groups HDU). [#14998]

- Fixed a bug that caused ``Cutout2D`` to not work correctly with ``CompImageHDU.section`` [#14999]

- Fixed a bug that caused compressed images with TFORM missing the optional '1' prefix to not be readable. [#15001]

astropy.modeling
^^^^^^^^^^^^^^^^

- All models can be pickled now. [#14902]

astropy.nddata
^^^^^^^^^^^^^^

- Restore bitmask propagation behavior in ``NDData.mask``, plus a fix
  for arithmetic between masked and unmasked ``NDData`` objects. [#14995]

astropy.table
^^^^^^^^^^^^^

- Fix a bug where table indexes were not using a stable sort order. This was causing the
  order of rows within groups to not match the original table order when an indexed table
  was grouped. [#14907]

astropy.units
^^^^^^^^^^^^^

- In VOunits, "pix", "au", "a", and "ct" are removed from the list of deprecated units. [#14885]

Version 5.3 (2023-05-22)
========================

New Features
------------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Add optional parameter ``refresh_cache`` to ``EarthLocation.of_site()`` and
  ``EarthLocation.get_site_names()`` to force the download of the latest site
  registry. [#13993]

- Added ``atol`` argument to function ``is_O3`` and ``is_rotation`` in matrix utilities. [#14371]

- A new class ``astropy.coordinates.StokesCoord`` has been added to represent world coordinates describing polarization state.
  This change introduces a breaking change to the return value of ``astropy.wcs.WCS.pixel_to_world`` where before a ``u.Quantity`` object would be returned containing numerical values representing a Stokes profile now a ``StokesCoord`` object is returned. The previous numerical values can be accessed with ``StokesCoord.value``. [#14482]

- Add an optional parameter ``location`` to ``EarthLocation.get_itrs()``
  to allow the generation of topocentric ITRS coordinates with respect
  to a specific location. [#14628]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Two new cosmologies have been added, ``FlatwpwaCDM`` and ``Flatw0wzCDM``, which are the
  flat variants of ``wpwaCDM`` and ``w0wzCDM``, respectively. [#12353]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Add ability to read and write an RST (reStructuredText) ASCII table that
  includes additional header rows specifying any or all of the column dtype, unit,
  format, and description. This is available via the new ``header_rows`` keyword
  argument. [#14182]

astropy.io.fits
^^^^^^^^^^^^^^^

- Added support for >3D data in CompImageHDU [#14252]

- Added a ``CompImageHDU.section`` property which can be used to
  efficiently access subsets of the data, similarly to ``ImageHDU.section``.
  When using this, only the tiles required to cover the section are
  read from disk and decompressed. [#14353]

- Added support for ``'NOCOMPRESS'`` for the ``compression_type`` option in ``CompImageHDU``. [#14408]

- Added new properties ``compression_type`` and ``tile_shape`` on
  ``CompImageHDU``, giving the name of the compression algorithm
  and the shape of the tiles in the tiled compression respectively. [#14428]

- Do not call ``gc.collect()`` when closing a ``CompImageHDU`` object as it has a
  large performance penalty. [#14576]

- VLA tables can now be written with the unified I/O interface.
  When object types are present or the VLA contains different types a `TypeError`
  is thrown. [#14578]

astropy.io.misc
^^^^^^^^^^^^^^^

- Add support for writing/reading fixed-size and variable-length array columns to the parquet formatter. [#14237]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Added a method ``get_infos_by_name`` to make it easier to implement
  DALI-compliant protocols [#14212]

- Updating the built-in UCD list to upstream 1.5 (which requires a minor
  update to the parser) [#14554]

astropy.modeling
^^^^^^^^^^^^^^^^

- Enable check for poorly conditioned fits in ``LinearLSQFitter`` for polynomial
  models with fixed inputs. [#14037]

astropy.nddata
^^^^^^^^^^^^^^

- ``astropy.nddata.NDDataArray`` now has collapsing methods like ``sum``,
  ``mean``, ``min``, and ``max`` which operate along any axes, and better
  support for ``astropy.utils.Masked`` objects. [#14175]

astropy.stats
^^^^^^^^^^^^^

- ``vonmisesmle`` has now functioning "weights" and "axis" parameters that work equivalently
  to the rest of the functions in the ``circstats`` module (``circmean``, ``rayleightest``, etc.) [#14533]

astropy.table
^^^^^^^^^^^^^

- ``Table`` and ``QTable`` can now use the ``|`` and ``|=`` operators for
  dictionary-style merge and update. [#14187]

astropy.time
^^^^^^^^^^^^

- Add a ``leap_second_strict`` argument to the ``Time.to_datetime()`` method. This
  controls the behavior when converting a time within a leap second to the ``datetime``
  format and can take the values ``raise`` (the default), ``warn``, or ``silent``. [#14606]

astropy.timeseries
^^^^^^^^^^^^^^^^^^

- Adds the ``astropy.timeseries.LombScargleMultiband`` class, which is an
  extension of the ``astropy.timeseries.LombScargle`` class. It enables the
  generation of periodograms for datasets with measurements taken in more than
  one photometric band. [#14016]

- Add ``unit_parse_strict`` parameter to the Kepler reader to control the warnings
  emitted when reading files. [#14294]

astropy.units
^^^^^^^^^^^^^

- Add support for degrees Celsius for FITS. Parsing "Celsius" and "deg C" is now
  supported and astropy will output "Celsius" into FITS.

  Note that "deg C" is only provided for compatibility with existing FITS files,
  as it does not conform to the normal unit standard, where this should be read
  as "degree * Coulomb". Indeed, compound units like "deg C kg-1" will still be
  parsed as "Coulomb degree per kilogram". [#14042]

- Enabled the ``equal_nan`` keyword argument for ``np.array_equal()`` when the
  arguments are ``astropy.units.Quantity`` instances. [#14135]

- Allow "console" and "unicode" formats for conversion to string of
  function units. [#14407]

- Add a "fraction" options to all the unit ``format`` classes, which determine
  whether, if a unit has bases raised to a negative power, a string
  representation should just show the negative powers (``fraction=False``) or
  use a fraction, and, in the latter case, whether to use a single-line
  representation using a solidus (``fraction='inline'`` or ``fraction=True``)
  or, if the format supports it, a multi-line presentation with the numerator
  and denominator separated by a horizontal line (``fraction='multiline'``). [#14449]

astropy.utils
^^^^^^^^^^^^^

- The ``mean`` method on ``NDDataArray`` now avoids a division by zero
  warning when taking the mean of a fully-masked slice (and still
  returns ``np.nan``). [#14341]

- Ensure we can read the newer ``IERS_B`` files produced by the International
  Earth Rotation and Reference Systems Service, and point
  ``astropy.utils.iers.IERS_B_URL`` to the new location. [#14382]


API Changes
-----------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- ``get_moon()`` is deprecated and may be removed in a future version of
  ``astropy``. Calling ``get_moon(...)`` should be replaced with
  ``get_body("moon", ...)``. [#14354]

astropy.io.fits
^^^^^^^^^^^^^^^

- Deprecate the auto-fixing of tile sizes for HCOMPRESS_1 tiled
  image compression when the tile size could be changed by +1
  to make it acceptable. [#14410]

- The ``tile_size=`` argument to ``CompImageHDU`` has been deprecated
  as it was confusing that it was required to be in the opposite
  order to the data shape (it was in header rather than Numpy order).
  Instead, users should make use of the ``tile_shape=`` argument which
  is in Numpy shape order. [#14428]

astropy.modeling
^^^^^^^^^^^^^^^^

- Deprecate the ``humlicek2`` method for `~astropy.modeling.functional_models.Voigt1D` in favor
  of using the ``wofz`` method using the `scipy.special.wofz` implementation of the
  Fadeeva function whenever `scipy` is installed. [#14013]

- Deprecated ``astropy.modeling.utils.comb()`` function in favor of ``comb()``
  from ``math`` standard library. [#14038]

- Propagate measurement uncertainties via the ``weights`` keyword argument into the
  parameter covariances. [#14519]

astropy.units
^^^^^^^^^^^^^

- The conversion of ``astropy.units.Quantity`` to ``bool``
  that was deprecated since astropy 3.0 now raises a ``ValueError``.
  This affects statements like ``if quantity``.
  Use explicit comparisons like ``if quantity.value != 0``
  or ``if quantity is not None`` instead. [#14124]

- Operations on ``Quantity`` in tables are sped up by only copying ``info`` when
  it makes sense (i.e., when the object can still logically be thought of as the
  same, such as in unit changes or slicing). ``info`` is no longer copied if a
  ``Quantity`` is part of an operation. [#14253]

- The ``Quantity.nansum`` method has been deprecated. It was always weird that it
  was present, since ``ndarray`` does not have a similar method, and the other
  ``nan*`` functions such as ``nanmean`` did not have a corresponding method.
  Use ``np.nansum(quantity)`` instead. [#14267]

- The unused ``units.format.Unscaled`` format class has been deprecated. [#14417]

- The order in which unit bases are displayed has been changed to match the
  order bases are stored in internally, which is by descending power to which
  the base is raised, and alphabetical after. This helps avoid monstrosities
  like ``beam^-1 Jy`` for ``format='fits'``.

  Note that this may affect doctests that use quantities with complicated units. [#14439]

astropy.utils
^^^^^^^^^^^^^

- For ``Masked`` instances, the ``where`` argument for any ufunc can now
  also be masked (with any masked elements masked in the output as well).
  This is not very useful in itself, but avoids problems in conditional
  functions (like ``np.add(ma, 1, where=ma>10)``). [#14590]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- The pixel attribute of ``astropy.visualization.wcsaxes.frame.Spine`` is deprecated
  and will be removed in a future astropy version.
  Because it is (in general) not possible to correctly calculate pixel
  coordinates before Matplotlib is drawing a figure, instead set the world or data
  coordinates of the ``Spine`` using the appropriate setters. [#13989]

- Passing a bare number as the ``coord_wrap`` argument to ``CoordinateHelper.set_coord_type`` is deprecated.
  Pass a ``Quantity`` with units equivalent to angular degrees instead.

  The ``.coord_wrap`` attribute of ``CoordinateHelper`` is now a ``Quantity`` instead of a bare number. [#14050]


Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- ``Angle.to_string()`` was changed to ensure it matches the behaviour of
  ``Quantity.to_string()`` in having a space between the value and the unit
  for display with non-degree and hourangle units (i.e., the case in which
  units are displayed by their name; the sexagesimal case for degrees or
  hourangle that uses symbols is not changed). [#14379]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fix an issue in the ``io.ascii`` QDP format reader to allow lower-case commands in the
  table data file. Previously it required all upper case in order to parse QDP files. [#14365]

astropy.io.fits
^^^^^^^^^^^^^^^

- Compressing/decompressing a floating point dataset containing NaN values will
  no longer read in the whole tile as NaNs.

  Fixed segmentation faults that occurred when compressing/decompressing data
  with the PLIO_1 algorithm. [#14252]

- ``Card`` now uses the default Python representation for floating point
  values. [#14508]

- ``ImageHDU`` now properly rejects Numpy scalars, avoiding data corruption. [#14528]

- Fix issues with double quotes in CONTINUE cards. [#14598]

- Fixes an issue where FITS_rec was incorrectly raising a ValueError exception when the heapsize was greater than 2**31
  when the Column type was 'Q' instead of 'P'. [#14810]

astropy.io.misc
^^^^^^^^^^^^^^^

- Columns with big-endian byte ordering (such as those read in from a FITS table) can now be serialized with Parquet. [#14373]

astropy.modeling
^^^^^^^^^^^^^^^^

- Bugfix for using ``getter/setter`` in properties to adjust the internal (computational)
  value of a property vs its external proxy value when the values involve units. [#14512]

- Fix issue with ``filter_non_finite`` option when fitting with ``weights`` via passing
  the ``weights`` through the non-finite-filter alongside the input data. [#14695]

- Fixed an issue with Parameter where a getter could be input without a
  setter (or vice versa). [#14708]

astropy.time
^^^^^^^^^^^^

- Using quantities with units of time for ``Time`` format 'decimalyear' will now
  raise an error instead of converting the quantity to days and then
  interpreting the value as years. An error is raised instead of attempting to
  interpret the unit as years, since the interpretation is ambiguous: in
  'decimaltime' years are equal to 365 or 366 days, while for regular time units
  the year is defined as 365.25 days. [#14566]

astropy.uncertainty
^^^^^^^^^^^^^^^^^^^

- Ensure that ``Distribution`` can be compared with ``==`` and ``!=``
  with regular arrays or scalars, and that inplace operations like
  ``dist[dist<0] *= -1`` work. [#14421]

astropy.units
^^^^^^^^^^^^^

- Modified ``astropy.units.Quantity.__array_ufunc__()`` to return ``NotImplemented`` instead of raising a ``ValueError`` if the inputs are incompatible. [#13977]

- Modified the behavior of ``numpy.array_equal()`` and ``numpy.array_equiv()`` to
  return ``False`` instead of raising an error if their arguments are
  ``astropy.units.Quantity`` instances with incompatible units. [#14163]

- Spaces have been regularized for the ``unicode`` and ``console`` output
  formats: no extraneous spaces in front of the unit, and always a space
  between a possible scale factor and the unit. [#14413]

- Prefixed degrees and arcmin are now typeset without using the symbol in
  ``latex`` and ``unicode`` formats (i.e., ``mdeg`` instead of ``m°``),
  as was already the case for arcsec. [#14419]

- Ensure the unit is kept in ``np.median`` even if the result is a scalar ``nan``
  (the unit was lost for numpy < 1.22). [#14635]

- Ensure that ``Quantity`` with structured dtype can be set using non-structured
  ``Quantity`` (if units match), and that structured dtype names are inferred
  correctly in the creation of ``StructuredUnit``, thus avoiding mismatches
  when setting units. [#14680]

astropy.utils
^^^^^^^^^^^^^

- When using astropy in environments with sparse file systems (e.g., where the temporary directory and astropy data directory resides in different volumes), ``os.rename`` may fail with ``OSError: [Errno 18] Invalid cross-device link``.
  This may affect some clean-up operations executed by the ``data`` module, causing them to fail.
  This patch is to catch ``OSError`` with ``errno == EXDEV`` (i.e., Errno 18) when performing these operations and try to use ``shutil.move`` instead to relocate the data. [#13730]

- Ensure masks are propagated correctly for ``outer`` methods of ufuncs also if
  one of the inputs is not actually masked. [#14624]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- The location of a ``astropy.visualization.wcsaxes.frame.Spine`` in a plot is now
  correctly calculated when the DPI of a figure changes between a WCSAxes being
  created and the figure being drawn. [#13989]

- ``CoordinateHelper.set_ticks()`` now accepts ``number=0``. Previously it errored. [#14160]

- ``WCSAxes.plot_coord`` and ``plot_scatter`` now work correctly for APE 14 compliant WCSes where the units are not always converted to degrees. [#14251]

- Fixed a bug where coordinate overlays did not automatically determine the
  longitude wrap angle or the appropriate units. [#14326]

astropy.wcs
^^^^^^^^^^^

- Fix bugs with high-level WCS API on ``wcs.WCS`` object when using ``-TAB``
  coordinates. [#13571]

- Fixed a bug in how WCS handles ``PVi_ja`` header coefficients when ``CTYPE``
  has ``-SIP`` suffix and in how code detects TPV distortions. [#14295]


Other Changes and Additions
---------------------------

- The minimum supported version of Python is now 3.9, changing from 3.8. [#14286]

- The minimum supported version of Numpy is now 1.21. [#14349]

- The minimum supported version of matplotlib is now 3.3. [#14286, #14321]

- ``astropy`` no longer publishes wheels for i686 architecture. [#14517]

- Added a pre-commit configuration for codespell. [#13985]

- Removed a large fraction of the bundled CFITSIO code and internally refactored
  FITS compression-related code, which has resulted in a speedup when compiling
  astropy from source (40% faster in some cases). [#14252]

- The CFITSIO library is no longer bundled in full with astropy and
  the option to build against an external installation of CFITSIO
  has now been removed, so the ASTROPY_USE_SYSTEM_CFITSIO environment
  variable will be ignored during building. [#14311]

- Updated CDS URL for Sesame look-up as the old URL is deprecated. [#14681]


Version 5.2.2 (2023-03-28)
==========================

Bug Fixes
---------

astropy.io.ascii
^^^^^^^^^^^^^^^^

- CDS and MRT tables with units that contain with multiple divisions, such as
  ``km/s/Mpc`` now parse correctly as being equal to ``km/(s.Mpc)``. [#14369]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fix ``FITSDiff`` when table contains a VLA column with the Q type. [#14539]

astropy.table
^^^^^^^^^^^^^

- Fix a bug when creating a ``QTable`` when a ``Quantity`` input column is present and the
  ``units`` argument modifies the unit of that column. This now works as expected where
  previously this caused an exception. [#14357]

astropy.units
^^^^^^^^^^^^^

- CDS units with multiple divisions, such as ``km/s/Mpc`` now parse
  correctly as being equal to ``km/(s.Mpc)``. [#14369]

astropy.wcs
^^^^^^^^^^^

- Fixed a bug that caused subclasses of BaseHighLevelWCS and HighLevelWCSMixin to
  not work correctly under certain conditions if they did not have ``world_n_dim``
  and ``pixel_n_dim`` defined on them. [#14495]


Version 5.2.1 (2023-01-06)
==========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Fix to ITRS frame ``earth_location`` attribute to give the correct result for
  a topocentric frame. [#14180]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Bounds are no longer passed to the scipy minimizer for methods Brent and
  Golden. The scipy minimizer never used the bounds but silently accepted them.
  In scipy v1.11.0.dev0+ an error is raised, so we now pass None as the bounds
  to the minimizer. Users should not be affected by this change. [#14232]

astropy.io.fits
^^^^^^^^^^^^^^^

- Tables with multidimensional variable length array can now be properly read
  and written. [#13417]

astropy.units
^^^^^^^^^^^^^

- Modified the behavior of ``numpy.histogram()``,
  ``numpy.histogram_bin_edges()``, ``numpy.histogram2d()``, and
  ``numpy.histogramdd()`` so that the ``range`` argument must a compatible
  instance of ``astropy.units.Quantity`` if the other arguments are instances of
  ``astropy.units.Quantity``. [#14213]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Improved the performance of drawing WCSAxes grids by skipping some unnecessary
  computations. [#14164]

- Fixed WCSAxes sometimes triggering a NumPy RuntimeWarning when determining the
  coordinate range of the axes. [#14211]

Other Changes and Additions
---------------------------

- Fix compatibility with Numpy 1.24. [#14193]

Version 5.2 (2022-12-12)
========================

New Features
------------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Adds new topocentric ITRS frame and direct transforms to and from the observed
  frames ``AltAz`` and ``HADec`` with the ability to add or remove refraction
  corrections as required. Since these frames are all within the ITRS, there are
  no corrections applied other than refraction in the transforms. This makes the
  topocentric ITRS frame and these transforms convenient for observers of near
  Earth objects where stellar aberration should be omitted. [#13398]

- Allow comparing ``SkyCoord`` to frames with data. [#13477]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Cosmology instance can be parsed from or converted to a HTML table using
  the new HTML methods in Cosmology's ``to/from_format`` I/O. [#13075]

- A new comparison function has been added -- ``cosmology_equal()`` -- that
  mirrors its ``numpy`` counterpart but allows for the arguments to be converted
  to a ``Cosmology`` and to compare flat cosmologies with their non-flat
  equivalents. [#13104]

- Cosmology equivalence for flat FLRW cosmologies has been generalized to apply
  to all cosmologies using the FlatCosmology mixin. [#13261]

- The cosmological redshift unit now has a physical type of ``"redshift"``. [#13561]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Add ability to read and write a fixed width ASCII table that includes additional
  header rows specifying any or all of the column dtype, unit, format, and
  description. This is available in the ``fixed_width`` and
  ``fixed_width_two_line`` formats via the new ``header_rows`` keyword argument. [#13734]

astropy.io.fits
^^^^^^^^^^^^^^^

- Added support to the ``io.fits`` API for reading and writing file paths of the
  form ``~/file.fits`` or ``~<username>/file.fits``, referring to the home
  directory of the current user or the specified user, respectively. [#13131]

- Added support for opening remote and cloud-hosted FITS files using the
  ``fsspec`` package, which has been added as an optional dependency. [#13238]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Added support in ``io.votable`` for reading and writing file paths of the form
  ``~/file.xml`` or ``~<username>/file.xml``, referring to the home directory of
  the current user or the specified user, respectively. [#13149]

astropy.modeling
^^^^^^^^^^^^^^^^

- Add option to non-linear fitters which enables automatic
  exclusion of non-finite values from the fit data. [#13259]

astropy.nddata
^^^^^^^^^^^^^^

- Modified ``Cutout2D`` to allow objects of type ``astropy.io.fits.Section``
  to be passed to the ``data`` parameter. [#13238]

- Add a PSF image representation to ``astropy.nddata.NDData`` and ``astropy.nddata.CCDData``. [#13743]

astropy.table
^^^^^^^^^^^^^

- An Astropy table can now be converted to a scalar NumPy object array. For NumPy
  >= 1.20, a list of Astropy tables can be converted to an NumPy object array of
  tables. [#13469]

astropy.time
^^^^^^^^^^^^

- Added the ``astropy.time.Time.mean()`` method which also enables the ``numpy.mean()`` function to be used on instances of ``astropy.time.Time``. [#13508]

- Improve the performance of getting the string representation of a large ``Time``
  or ``TimeDelta`` object. This is done via a new ``to_string()`` method that does
  the time string format conversion only for the outputted values. Previously the
  entire array was formatted in advance. [#13555]

astropy.units
^^^^^^^^^^^^^

- It is now possible to use unit format names as string format specifiers for a
  ``Quantity``, e.g. ``f'{1e12*u.m/u.s:latex_inline}'`` now produces the string
  ``'$1 \\times 10^{12} \\; \\mathrm{m\\,s^{-1}}$'``. [#13050]

- Ensure that the ``argmin`` and ``argmax`` methods of ``Quantity`` support the
  ``keepdims`` argument when numpy does (numpy version 1.22 and later). [#13329]

- ``numpy.lib.recfunctions.merge_arrays()`` is registered with numpy overload for
  ``Quantity``. [#13669]

- Added SI prefixes for quecto ("q", :math:`10^{-30}`), ronto ("r",
  :math:`10^{-27}`), ronna ("R", :math:`10^{27}`), and quetta ("Q",
  :math:`10^{30}`). [#14046]

astropy.utils
^^^^^^^^^^^^^

- Added the ``use_fsspec``, ``fsspec_kwargs``, and ``close_files`` arguments
  to ``utils.data.get_readable_fileobj``. [#13238]

- Ensure that the ``argmin`` and ``argmax`` methods of ``Masked`` instances
  support the ``keepdims`` argument when numpy does (numpy version 1.22 and
  later). [#13329]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Add helper functions for WCSAxes instances to draw the instrument beam and a physical scale. [#12102]

- Add a ``scatter_coord`` method to the ``wcsaxes`` functionality based on the
  existing ``plot_coord`` method but that calls ``matplotlib.pyplot.scatter``. [#13562]

- Added a ``sinh`` stretch option to ``simple_norm``. [#13746]

- It is now possible to define "tickable" gridlines for the purpose of placing ticks or tick labels in the interior of WCSAxes plots. [#13829]


API Changes
-----------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Removed deprecated ``MexicanHat1DKernel`` and ``MexicanHat2DKernel``
  classes. Please use ``RickerWavelet1DKernel`` and
  ``RickerWavelet2DKernel`` instead. [#13300]

astropy.units
^^^^^^^^^^^^^

- Multiplying a ``LogQuantity`` like ``Magnitude`` with dimensionless physical
  units by an array will no longer downcast to ``Quantity``. [#12579]

- Quantity normally upcasts integer dtypes to floats, unless the dtype is
  specifically provided.
  Before this happened when ``dtype=None``; now the default has been changed to
  ``dtype=numpy.inexact`` and ``dtype=None`` has the same meaning as in `numpy`. [#12941]

- In "in-place unit changes" of the form ``quantity <<= new_unit``, the result
  will now share memory with the original only if the conversion could be done
  through a simple multiplication with a scale factor. Hence, memory will not be
  shared if the quantity has integer ```dtype``` or is structured, or when the
  conversion is through an equivalency. [#13638]

- When ``Quantity`` is constructed from a structured array and ``unit`` is
  ``None``, the default unit is now structured like the input data. [#13676]

astropy.utils
^^^^^^^^^^^^^

- ``astropy.utils.misc.suppress`` has been removed, use ``contextlib.suppress``
  instead. ``astropy.utils.namedtuple_asdict`` has been removed, instead use
  method ``._asdict`` on a ``namedtuple``. ``override__dir__`` has been deprecated
  and will be removed in a future version, see the docstring for the better
  alternative. [#13636]

- ``astropy.utils.misc.possible_filename`` has been removed. [#13661]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Rename number-of-samples keyword ``nsamples`` in ``ZScaleInterval`` to align
  with the ``n_samples`` keyword used in all other ``Interval`` classes in
  this module. [#13810]


Bug Fixes
---------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Fixed convolution Kernels to ensure the that returned kernels
  are normalized to sum to one (e.g., ``Gaussian1DKernel``,
  ``Gaussian2DKernel``). Also fixed the Kernel ``truncation`` calculation. [#13299]

- Fix import error with setuptools v65.6.0 by replacing
  ``numpy.ctypeslib.load_library`` with Cython to load the C convolution
  extension. [#14035]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- ``BaseCoordinateFrame.get_frame_attr_names()`` had a misleading name,
  because it actually provided a ``dict`` of attribute names and
  their default values. It is now deprecated and replaced by ``BaseCoordinateFrame.get_frame_attr_defaults()``.
  The fastest way to obtain the attribute names is ``BaseFrame.frame_attributes.keys()``. [#13484]

- Fixed bug that caused ``earth_orientation.nutation_matrix()`` to error instead of returning output. [#13572]

- Ensure that ``angle.to_string()`` continues to work after pickling,
  and that units passed on to ``to_string()`` or the ``Angle``
  initializer can be composite units (like ``u.hour**1``), which might
  result from preceding calculations. [#13933]

astropy.io.fits
^^^^^^^^^^^^^^^

- ``report_diff_values()`` have now two new parameters ``rtol`` and ``atol`` to make the
  report consistent with ``numpy.allclose`` results.
  This fixes ``FITSDiff`` with multi-dimensional columns. [#13465]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Fixed two bugs in validator.validator.make_validation_report:
  - ProgressBar iterator was not called correctly.
  - make_validation_report now handles input string urls correctly. [#14102]

astropy.timeseries
^^^^^^^^^^^^^^^^^^

- Fixed a performance regression in ``timeseries.aggregate_downsample``
  introduced in Astropy 5.0 / #11266. [#13069]

astropy.units
^^^^^^^^^^^^^

- Unit changes of the form ``quantity <<= new_unit`` will now work also if the
  quantity is integer. The result will always be float. This means that the result
  will not share memory with the original. [#13638]

- Ensure dimensionless quantities can be added inplace to regular ndarray. [#13913]

astropy.utils
^^^^^^^^^^^^^

- Fixed an incompatibility with latest Python 3.1x versions that kept
  ``astropy.utils.data.download_file`` from switching to TLS+FTP mode. [#14092]

- ``np.quantile`` and ``np.percentile`` can now be used on ``Masked``
   arrays and quantities also with ``keepdims=True``. [#14113]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Significantly improve performance of ``ManualInterval`` when both limits
  are specified manually. [#13898]


Other Changes and Additions
---------------------------

- The deprecated private ``astropy._erfa`` module has been removed. Use
  ``pyerfa``, which is a dependency of ``astropy`` and can be imported directly
  using ``import erfa``. [#13317]

- The minimum version required for numpy is now 1.20 and that for scipy 1.5. [#13885]

- Updated the bundled CFITSIO library to 4.2.0. [#14020]


Version 5.1.1 (2022-10-23)
==========================

API Changes
-----------

astropy.wcs
^^^^^^^^^^^

- The ``pixel`` argument to ``astropy.visualization.wcsaxes.ticklabels.TickLabels.add``
  no longer does anything, is deprecated, and will be removed in a future
  astropy version. It has been replaced by a new required ``data`` argument, which
  should be used to specify the data coordinates of the tick label being added.

  This changes has been made because it is (in general) not possible to correctly
  calculate pixel coordinates before Matplotlib is drawing a figure. [#12630]


Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Fixed a bug that prevented ``SkyOffsetFrame`` instances to be pickled by adding
  a custom ``__reduce__`` method to the class (see issue #9249). [#13305]

- Fixed the check for invalid ``Latitude`` values for float32 values.
  ``Latitude`` now accepts the float32 value of pi/2, which was rejected
  before because a comparison was made using the slightly smaller float64 representation.
  See issue #13708. [#13745]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fixed confusing chained exception messages of ``read()`` function when it fails. [#13170]

- When writing out a :class:`~astropy.table.Table` to HTML format, the
  ``formats`` keyword argument to the :meth:`~astropy.table.Table.write` method
  will now be applied. [#13453]

astropy.io.fits
^^^^^^^^^^^^^^^

- ``heapsize`` is now checked for VLA tables. An error is thrown whether P format is used
  but the heap size is bigger than what can be indexed with a 32 bit signed int. [#13429]

- Fix parsing of ascii TFORM when precision is missing. [#13520]

- A compressed image HDU created from the header of a PRIMARY HDU, now correctly updates
  'XTENSION' and 'SIMPLE' keywords. [#13557]

- Empty variable-length arrays are now properly handled when pathological combinations of
  heapoffset and heapsize are encountered. [#13621]

- ``PCOUNT`` and ``GCOUNT`` keywords are now removed from an uncompressed Primary header,
  for compliance with ``fitsverify`` behavior. [#13753]

astropy.modeling
^^^^^^^^^^^^^^^^

- Bugfix for using ``MagUnit`` units on model parameters. [#13158]

- Fix bug in using non-linear fitters to fit 0-degree polynomials using weights. [#13628]

astropy.table
^^^^^^^^^^^^^

- Fix a problem where accessing one field of a structured column returned a Column
  with the same info as the original column. This resulted in unintuitive behavior
  in general and an exception if the format for the column was set. [#13269]

- Tables with columns with structured data can now be properly stacked and joined. [#13306]

- Update jQuery to 3.6.0, to pick up security fixes. [#13438]

- Fix a Python 3.11 compatibility issue. Ensure that when removing a table column
  that the ``pprint_include_names`` or ``pprint_exclude_names`` attributes get
  updated correctly. [#13639]

- When using ``add_columns`` with same indexes in ``indexes`` option or without
  specifying the option, the order of the new columns will now be kept. [#13783]

- Fix a bug when printing or getting the representation of a multidimensional
  table column that has a zero dimension. [#13838]

- Ensure that mixin columns and their ``info`` are not shared between tables
  even when their underlying data is shared with ``copy=False``. [#13842]

astropy.time
^^^^^^^^^^^^

- Fix ``Time.insert()`` on times which have their ``out_subfmt`` set. [#12732]

- Prevent ``Time()`` from being initialized with an invalid precision
  leading to incorrect results when representing the time as a string. [#13068]

- Fix a bug in Time where a date string like "2022-08-01.123" was being parsed
  as an ISO-format time "2022-08-01 00:00:00.123". The fractional part at the
  end of the string was being taken as seconds. Now this raises an exception
  because the string is not in ISO format. [#13731]

astropy.units
^^^^^^^^^^^^^

- Significantly improved the performance of parsing composite units with the FITS
  format, by ensuring the ``detailed_exception`` argument is properly passed on
  and thus used. [#12699]

- Ensure that ``np.concatenate`` on quantities can take a ``dtype`` argument (added in numpy 1.20). [#13323]

- Ensure that the units of any ``initial`` argument to reductions such as
  ``np.add.reduce`` (which underlies ``np.sum``) are properly taken into account. [#13340]

astropy.utils
^^^^^^^^^^^^^

- Ensure that ``np.concatenate`` on masked data can take a ``dtype`` argument (added in numpy 1.20). [#13323]

- Fix error when suppressing download progress bar while using non-default
  ``sys.stdout`` stream. [#13352]

- Ensure ``str`` and ``repr`` work properly for ``Masked`` versions of
  structured subarrays. [#13404]

- If an attribute is created using ``deprecated_attribute()`` with the
  ``alternative`` argument then getting or setting the value of the deprecated
  attribute now accesses its replacement. [#13824]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Fixed calling ``.tight_layout()`` on a WCSAxes. [#12418]

astropy.wcs
^^^^^^^^^^^

- ``WCS.pixel_to_world`` now creates an ``EarthLocation`` object using ``MJD-AVG``
  if present before falling back to the old behaviour of using ``MJD-OBS``. [#12598]

- The locations of ``WCSAxes`` ticks and tick-labels are now correctly calculated
  when the DPI of a figure changes between a WCSAxes being created and the figure
  being drawn, or when a rasterized artist is added to the WCSAxes. [#12630]

- Fix a bug where ``SlicedLowLevelWCS.world_to_pixel_values`` would break when
  the result of the transform is dependent on the coordinate of a sliced out
  pixel. [#13579]

- Updated bundled WCSLIB version to 7.12. This update includes bug fixes to
  ``wcssub()`` in how it handles temporal axes with -TAB and fixes handling
  of status returns from ``linp2x()`` and ``linx2p()`` relating to distortion
  functions, in particular affecting TPV distortions - see #13509. For a full
  list of changes - see http://www.atnf.csiro.au/people/mcalabre/WCS/CHANGES or
  `astropy/cextern/wcslib/CHANGES <https://github.com/astropy/astropy/blob/24e8730c63902d035cb9110eae2a9ebec12d8905/cextern/wcslib/CHANGES>`_. [#13635]

- Fixed WCS validation not working properly if HDUList is needed
  for multi-extension FITS file. [#13668]


Other Changes and Additions
---------------------------

- Development wheels of astropy should now be installed from
  https://pypi.anaconda.org/astropy/simple instead of from
  https://pkgs.dev.azure.com/astropy-project/astropy/_packaging/nightly/pypi/simple. [#13431]

- Compatibility with Python 3.11, 3.10.7, 3.9.14, 3.8.14 [#13614]


Version 5.1 (2022-05-24)
========================

New Features
------------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- The ephemeris used in ``astropy.coordinates`` can now be set to any version of
  the JPL ephemeris available from https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/. [#12541]

- ``Angle.to_string()`` now accepts the ``'latex_inline'`` unit format. [#13056]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Cosmology instance can be parsed from or converted to a YAML string using
  the new "yaml" format in Cosmology's ``to/from_format`` I/O. [#12279]

- Register "astropy.row" into Cosmology's to/from format I/O, allowing a
  Cosmology instance to be parse from or converted to an Astropy Table Row. [#12313]

- A method ``clone`` has been added to ``Parameter`` to quickly deep copy the
  object and change any constructor argument.
  A supporting equality method is added, and ``repr`` is enhanced to be able to
  roundtrip -- ``eval(repr(Parameter()))`` -- if the Parameter's arguments can
  similarly roundtrip.
  Parameter's arguments are made keyword-only. [#12479]

- Add methods ``Otot`` and ``Otot0`` to FLRW cosmologies to calculate the total
  energy density of the Universe. [#12590]

- Add property ``is_flat`` to cosmologies to calculate the curvature of the Universe.
  ``Cosmology`` is now an abstract class and subclasses must override the
  abstract property ``is_flat``. [#12606]

- For converting a cosmology to a mapping, two new boolean keyword arguments are
  added: ``cosmology_as_str`` for turning the class reference to a string,
  instead of the class object itself, and ``move_from_meta`` to merge the
  metadata with the rest of the returned mapping instead of adding it as a
  nested dictionary. [#12710]

- Register format "astropy.cosmology" with Cosmology I/O. [#12736]

- Cosmological equivalency (``Cosmology.is_equivalent``) can now be extended
  to any Python object that can be converted to a Cosmology, using the new
  keyword argument ``format``.
  This allows e.g. a properly formatted Table to be equivalent to a Cosmology. [#12740]

- The new module ``cosmology/tests/helper.py`` has been added to provide tools
  for testing the cosmology module and related extensions. [#12966]

- A new property ``nonflat`` has been added to flat cosmologies
  (``FlatCosmologyMixin`` subclasses) to get an equivalent cosmology, but of the
  corresponding non-flat class. [#13076]

- ``clone`` has been enhanced to allow for flat cosmologies to clone on the
  equivalent non-flat cosmology. [#13099]

- ``cosmology`` file I/O uses the Unified Table I/O interface, which has added
  support for reading and writing file paths of the form ``~/file.ecsv`` or
  ``~<username>/file.ecsv``, referring to the home directory of the current user
  or the specified user, respectively. [#13129]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Simplify the way that the ``converters`` argument of ``io.ascii.read`` is
  provided. Previously this required wrapping each data type as the tuple returned
  by the ``io.ascii.convert_numpy()`` function and ensuring that the value is a
  ``list``. With this update you can write ``converters={'col1': bool}`` to force
  conversion as a ``bool`` instead of the previous syntax ``converters={'col1':
  [io.ascii.convert_numpy(bool)]}``. Note that this update is back-compatible with
  the old behavior. [#13073]

- Added support in ``io.ascii`` for reading and writing file paths of the form
  ``~/file.csv`` or ``~<username>/file.csv``, referring to the home directory of
  the current user or the specified user, respectively. [#13130]

astropy.io.fits
^^^^^^^^^^^^^^^

- Add option ``unit_parse_strict`` to ``astropy.io.fits.connect.read_table_fits``
  to enable warnings or errors about invalid FITS units when using ``astropy.table.Table.read``.
  The default for this new option is ``"warn"``, which means warnings are now raised for
  columns with invalid units. [#11843]

- Changes default FITS behavior to use buffered I/O
  rather than unbuffered I/O for performance reasons. [#12081]

- ``astropy.io.fits.Header`` now has a method to calculate the size
  (in bytes) of the data portion (with or without padding) following
  that header. [#12110]

astropy.io.misc
^^^^^^^^^^^^^^^

- Allow serialization of model unit equivalencies. [#10198]

- Built-in Cosmology subclasses can now be converted to/from YAML with the
  functions ``dump`` and ``load`` in ``astropy.io.misc.yaml``. [#12279]

- Add asdf support for ``Cosine1D``, ``Tangent1D``, ``ArcSine1D``,
  ``ArcCosine1D``, and ``ArcTangent1D`` models. [#12895]

- Add asdf support for ``Spline1D`` models. [#12897]

astropy.io.registry
^^^^^^^^^^^^^^^^^^^

- Added support to the Unified Table I/O interface for reading and writing file
  paths of the form ``~/file.csv`` or ``~<username>/file.csv``, referring to the
  home directory of the current user or the specified user, respectively. [#13129]

astropy.modeling
^^^^^^^^^^^^^^^^

- Add new fitters based on ``scipy.optimize.least_squares`` method of non-linear
  least-squares optimization: [#12051]

      - ``TRFLSQFitter`` using the Trust Region Reflective algorithm.
      - ``LMLSQFitter`` using the Levenberg-Marquardt algorithm (implemented by ``scipy.optimize.least_squares``).
      - ``DogBoxLSQFitter`` using the dogleg algorithm.

- Enable direct use of the ``ignored`` feature of ``ModelBoundingBox`` by users in
  addition to its use as part of enabling ``CompoundBoundingBox``. [#12384]

- Switch ``modeling.projections`` to use ``astropy.wcs.Prjprm`` wrapper internally and provide access to the ``astropy.wcs.Prjprm`` structure. [#12558]

- Add error to non-finite inputs to the ``LevMarLSQFitter``, to protect against soft scipy failure. [#12811]

- Allow the ``Ellipse2D`` and ``Sersic2D`` theta parameter to be input as
  an angular quantity. [#13030]

- Added ``Schechter1D`` model. [#13116]

astropy.nddata
^^^^^^^^^^^^^^

- Add support for converting between uncertainty types. This uncertainty
  conversion system uses a similar flow to the coordinate subsystem, where
  Cartesian is used as the common system. In this case, variance is used as the
  common system. [#12057]

- The ``as_image_hdu`` option is now available for ``CCDData.to_hdu`` and
  ``CCDData.write``. This option allows the user to get an ``ImageHDU`` as the
  first item of the returned ``HDUList``, instead of the default ``PrimaryHDU``. [#12962]

- File I/O through ``nddata.CCDData`` uses the Unified I/O interface, which has
  added support for reading and writing file paths of the form ``~/file.csv`` or
  ``~<username>/file.csv``, referring to the home directory of the current user
  or the specified user, respectively. [#13129]

astropy.table
^^^^^^^^^^^^^

- A new keyword-only argument ``kind`` was added to the ``Table.sort`` method to
  specify the sort algorithm. [#12637]

- Columns which are ``numpy`` structured arrays are now fully supported,
  effectively allowing tables within tables. This applies to ``Column``,
  ``MaskedColumn``, and ``Quantity`` columns. These structured data columns
  can be stored in ECSV, FITS, and HDF5 formats. [#12644]

- Improve the performance of ``np.searchsorted`` by a factor of 1000 for a
  bytes-type ``Column`` when the search value is ``str`` or an array of ``str``.
  This happens commonly for string data stored in FITS or HDF5 format files. [#12680]

- Add support for using mixin columns in group aggregation operations when the
  mixin supports the specified operation (e.g. ``np.sum`` works for ``Quantity``
  but not ``Time``). In cases where the operation is not supported the code now
  issues a warning and drops the column instead of raising an exception. [#12825]

- Added support to the Unified Table I/O interface for reading and writing file
  paths of the form ``~/file.csv`` or ``~<username>/file.csv``, referring to the
  home directory of the current user or the specified user, respectively. [#13129]

astropy.time
^^^^^^^^^^^^

- Add support for calling ``numpy.linspace()`` with two ``Time`` instances to
  generate a or multiple linearly spaced set(s) of times. [#13132]

astropy.units
^^^^^^^^^^^^^

- ``structured_to_unstructured`` and ``unstructured_to_structured`` in
  ``numpy.lib.recfunctions`` now work with Quantity. [#12486]

- Implement multiplication and division of LogQuantities and numbers [#12566]

- New ``doppler_redshift`` equivalency to convert between
  Doppler redshift and radial velocity. [#12709]

- Added the ``where`` keyword argument to the ``mean()``,``var()``, ``std()`` and ``nansum()`` methods of
  ``astropy.units.Quantity``. Also added the ``initial`` keyword argument to ``astropy.units.Quantity.nansum()``. [#12891]

- Added "Maxwell" as a unit for magnetic flux to the CGS module. [#12975]

- ``Quantity.to_string()`` and ``FunctionUnitBase.to_string()`` now accept the
  ``'latex_inline'`` unit format. The output of ``StructuredUnit.to_string()``
  when called with ``format='latex_inline'`` is now more consistent with the
  output when called with ``format='latex'``. [#13056]

astropy.utils
^^^^^^^^^^^^^

- Added the ``where`` keyword argument to the ``mean()``, ``var()``, ``std()``, ``any()``, and ``all()`` methods of
  ``astropy.utils.masked.MaskedNDArray``. [#12891]

- Improve handling of unavailable IERS-A (predictive future Earth rotation) data
  in two ways. First, allow conversions with degraded accuracy if the IERS-A data
  are missing or do not cover the required time span. This is done with a new
  config item ``conf.iers_degraded_accuracy`` which specifies the behavior when
  times are outside the range of IERS table. The options are 'error' (raise an
  ``IERSRangeError``, default), 'warn' (issue a ``IERSDegradedAccuracyWarning``)
  or 'ignore' (ignore the problem). Second, the logic for auto-downloads was
  changed to guarantee that no matter what happens with the IERS download
  operations, only warnings will be issued. [#13052]

astropy.wcs
^^^^^^^^^^^

- ``astropy.wcs.Celprm`` and ``astropy.wcs.Prjprm`` have been added
  to allow access to lower level WCSLIB functionality and to allow direct
  access to the ``cel`` and ``prj`` members of ``Wcsprm``. [#12514]

- Add ``temporal`` properties for convenient access of/selection of/testing for
  the ``TIME`` axis introduced in WCSLIB version 7.8. [#13094]

API Changes
-----------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- The ``dms_to_degrees`` and ``hms_to_hours`` functions (and implicitly
  tuple-based initialization of ``Angle``) is now deprecated, as it was
  difficult to be sure about the intent of the user for signed values of
  the degrees/hours, minutes, and seconds. [#13162]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- The already deprecated ``Planck18_arXiv_v2`` has been removed.
  Use ``Planck18`` instead [#12354]

- ``default_cosmology.get_cosmology_from_string`` is deprecated and will be
  removed in two minor versions.
  Use ``getattr(astropy.cosmology, <str>)`` instead. [#12375]

- In I/O, conversions of Parameters move more relevant information from the
  Parameter to the Column.
  The default Parameter ``format_spec`` is changed from ``".3g"`` to ``""``. [#12612]

- Units of redshift are added to ``z_reion`` in built-in realizations' metadata. [#12624]

- Cosmology realizations (e.g. ``Planck18``) and parameter dictionaries are now
  lazily loaded from source files. [#12746]

- The Cosmology Parameter argument "fmt" for specifying a format spec
  has been deprecated in favor of using the built-in string representation from
  the Parameter's value's dtype. [#13072]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- When reading an ECSV file, changed the type checking
  to issue an ``InvalidEcsvDatatypeWarning`` instead of raising a ``ValueError``
  exception if the ``datatype`` of a column is not recognized in the ECSV standard.
  This also applies to older versions of ECSV files which used to silently
  proceed but now warn first. [#12841]

astropy.io.fits
^^^^^^^^^^^^^^^

- Removed deprecated ``clobber`` argument from functions in ``astropy.io.fits``. [#12258]

- Add ``-s/--sort`` argument to ``fitsheader`` to sort the fitsort-mode output. [#13106]

astropy.io.misc
^^^^^^^^^^^^^^^

- Deprecate asdf in astropy core in favor of the asdf-astropy package. [#12903, #12930]

astropy.modeling
^^^^^^^^^^^^^^^^

- Made ``astropy.modeling.fitting._fitter_to_model_params`` and ``astropy.modeling.fitting._model_to_fit_params``
  public methods. [#12585]

astropy.table
^^^^^^^^^^^^^

- Change the repr of the Table object to replace embedded newlines and tabs with
  ``r'\n'`` and ``r'\t'`` respectively. This improves the display of such tables. [#12631]

- A new keyword-only argument ``kind`` was added to the ``Table.sort`` method to
  specify the sort algorithm. The signature of ``Table.sort`` was modified so that
  the ``reverse`` argument is now keyword-only. Previously ``reverse`` could be
  specified as the second positional argument. [#12637]

- Changed behavior when a structured ``numpy.ndarray`` is added as a column to a
  ``Table``. Previously this was converted to a ``NdarrayMixin`` subclass of
  ``ndarray`` and added as a mixin column. This was because saving as a file (e.g.
  HDF5, FITS, ECSV) was not supported for structured array columns. Now a
  structured ``numpy.ndarray`` is added to the table as a native ``Column`` and
  saving to file is supported. [#13236]

astropy.tests
^^^^^^^^^^^^^

- Backward-compatible import of ``astropy.tests.disable_internet``
  has been removed; use ``pytest_remotedata.disable_internet``
  from ``pytest-remotedata`` instead. [#12633]

- Backward-compatible import of ``astropy.tests.helper.remote_data``
  has been removed; use ``pytest.mark.remote_data`` from ``pytest-remotedata``
  instead. [#12633]

- The following are deprecated and will be removed in a future release.
  Use ``pytest`` warning and exception handling instead: [#12633]

  * ``astropy.io.ascii.tests.common.raises``
  * ``astropy.tests.helper.catch_warnings``
  * ``astropy.tests.helper.ignore_warnings``
  * ``astropy.tests.helper.raises``
  * ``astropy.tests.helper.enable_deprecations_as_exceptions``
  * ``astropy.tests.helper.treat_deprecations_as_exceptions``

- Backward-compatible plugin ``astropy.tests.plugins.display``
  has been removed; use ``pytest-astropy-header`` instead. [#12633]

astropy.time
^^^^^^^^^^^^

- Creating an `~astropy.time.TimeDelta` object with numerical inputs
  that do not have a unit and without specifying an explicit format,
  for example ``TimeDelta(5)``,
  now results in a `~astropy.time.TimeDeltaMissingUnitWarning`.
  This also affects statements like ``Time("2020-01-01") + 5`` or
  ``Time("2020-01-05") - Time("2020-01-03") < 5``, which implicitly
  transform the right-hand side into an `~astropy.time.TimeDelta` instance. [#12888]

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- The machinery that makes observatory locations available as ``EarthLocation``
  objects is now smarter about processing observatory names from its data files.
  More names are available for use and the empty string is no longer considered
  to be a valid name. [#12721]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fixed ``io.ascii`` read and write functions for most formats to correctly handle
  data fields with embedded newlines for both the fast and pure-Python readers and
  writers. [#12631]

- Fix an issue when writing ``Time`` table columns to a file when the time
  ``format`` is one of ``datetime``, ``datetime64``, or ``ymdhms``. Previously,
  writing a ``Time`` column with one of these formats could result in an exception
  or else an incorrect output file that cannot be read back in. [#12842]

astropy.io.fits
^^^^^^^^^^^^^^^

- Add a ``mask_invalid`` option to ``Table.read`` to allow deactivating the
  masking of NaNs in float columns and empty strings in string columns. This
  option is necessary to allow effective use of memory-mapped reading with
  ``memmap=True``. [#12544]

- Fix ``CompImageHeader.clear()``. [#13102]

astropy.modeling
^^^^^^^^^^^^^^^^

- Bugfix for ``ignore`` functionality failing in ``ModelBoundingBox`` when using
  ``ignore`` option alongside passing bounding box data as tuples. [#13032]

astropy.table
^^^^^^^^^^^^^

- Fixed a bug in ``Table.show_in_browser`` using the ``jsviewer=True`` option
  to display the table with sortable columns. Previously the sort direction arrows
  were not being shown due to missing image files for the arrows. [#12716]

- Fix an issue when writing ``Time`` table columns to a file when the time
  ``format`` is one of ``datetime``, ``datetime64``, or ``ymdhms``. Previously,
  writing a ``Time`` column with one of these formats could result in an exception
  or else an incorrect output file that cannot be read back in. [#12842]

- Fixed a bug where it is not possible to set the ``.info.format`` property of a
  table structured column and get formatted output. [#13233]

- Fixed a bug when adding a masked structured array to a table. Previously this
  was auto-converted to a ``NdarrayMixin`` which loses the mask. With this fix
  the data are added to the table as a ``MaskedColumn`` and the mask is preserved. [#13236]

astropy.time
^^^^^^^^^^^^

- Fix an issue when writing ``Time`` table columns to a file when the time
  ``format`` is one of ``datetime``, ``datetime64``, or ``ymdhms``. Previously,
  writing a ``Time`` column with one of these formats could result in an exception
  or else an incorrect output file that cannot be read back in. [#12842]

astropy.utils
^^^^^^^^^^^^^

- Fixed a bug which caused ``numpy.interp`` to produce incorrect
  results when ``Masked`` arrays were passed. [#12978]

- Fixed HAS_YAML not working as intended. [#13066]

astropy.wcs
^^^^^^^^^^^

- Convert ``NoConvergence`` errors to warnings in ``world_to_pixel_values`` so that callers can work at least with the non converged solution. [#11693]

- Expose the ability to select TIME axis introduced in WCSLIB version 7.8. [#13062]

- Do not call ``wcstab`` on ``wcscopy`` and copy ``wtb`` members from the original WCS. [#13063]

- Updated bundled WCSLIB version to 7.11. This update together with 7.10
  includes bug fixes to ``tabini()`` and ``tabcpy()`` as well as several
  print formatting enhancements. For a full list of
  changes - see http://www.atnf.csiro.au/people/mcalabre/WCS/CHANGES [#13171]

- Fixed error that occurred in ``WCS.world_to_pixel`` for ``WCS`` objects with a
  spectral axis and observer location information when passing a ``SpectralCoord``
  that had missing observer or target information. [#13228]

Version 5.0.4 (2022-03-31)
==========================

Bug Fixes
---------

astropy.modeling
^^^^^^^^^^^^^^^^

- Fixed the ``Gaussian2D`` ``bounding_box`` when ``theta`` is an angular
  ``Quantity``. [#13021]

astropy.utils
^^^^^^^^^^^^^

- Reverted ``astropy.utils.iers.iers.IERS_A_URL`` to ``maia.usno.navy.mil`` domain instead
  of NASA FTP to work around server issues. [#13004]

Other Changes and Additions
---------------------------

- Updated bundled WCSLIB to version 7.9 with several bugfixes and added
  support for time coordinate axes in ``wcsset()`` and ``wcssub()``. The
  four-digit type code for the time axis will have the first digit set to 4,
  i.e., four digit code will be 4xxx where x is a digit 0-9. For a full list of
  bug fixes see https://www.atnf.csiro.au/people/mcalabre/WCS/CHANGES [#12994]


Version 5.0.3 (2022-03-25)
==========================

Bug Fixes
---------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Bugfix in ``astropy.convolution.utils.discretize_model`` which allows the function to handle a ``CompoundModel``.
  Before this fix, ``discretize_model`` was confusing ``CompoundModel`` with a callable function. [#12959]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fix write and read FITS tables with multidimensional items, using ``from_columns``
  without previously defined ``ColDefs`` structure. [#12863]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Fix VOTable linting to avoid use of shell option. [#12985]

astropy.utils
^^^^^^^^^^^^^

- Fix XML linting to avoid use of shell option. [#12985]

Other Changes and Additions
---------------------------

- Updated the bundled CFITSIO library to 4.1.0. [#12967]

Version 5.0.2 (2022-03-10)
==========================

Bug Fixes
---------

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Bugfix to add backwards compatibility for reading ECSV version 0.9 files with
  non-standard column datatypes (such as ``object``, ``str``, ``datetime64``,
  etc.), which would raise a ValueError in ECSV version 1.0. [#12880]

astropy.io.misc
^^^^^^^^^^^^^^^

- Bugfix for ``units_mapping`` schema's property name conflicts. Changes:
      * ``inputs`` to ``unit_inputs``
      * ``outputs`` to ``unit_outputs`` [#12800]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Fixed a bug where ``astropy.io.votable.validate`` was printing output to
  ``sys.stdout`` when the ``output`` parameter was set to ``None``. ``validate``
  now returns a string when ``output`` is set to ``None``, as documented.
  [#12604]

astropy.modeling
^^^^^^^^^^^^^^^^

- Fix handling of units on ``scale`` parameter in BlackBody model. [#12318]

- Indexing on models can now be used with all types of integers
  (like ``numpy.int64``) instead of just ``int``. [#12561]

- Fix computation of the separability of a ``CompoundModel`` where another
  ``CompoundModel`` is on the right hand side of the ``&`` operator. [#12907]

- Provide a hook (``Model._calculate_separability_matrix``) to allow subclasses
  of ``Model`` to define how to compute their separability matrix. [#12900]

astropy.stats
^^^^^^^^^^^^^

- Fixed a bug in which running ``kuiper_false_positive_probability(D,N)`` on
  distributions with many data points could produce NaN values for the false
  positive probability of the Kuiper statistic. [#12896]

astropy.wcs
^^^^^^^^^^^

- Fixed a bug due to which ``naxis``, ``pixel_shape``, and
  ``pixel_bounds`` attributes of ``astropy.wcs.WCS`` were not restored when
  an ``astropy.wcs.WCS`` object was unpickled. This fix also eliminates
  ``FITSFixedWarning`` warning issued during unpiclikng of the WCS objects
  related to the number of axes. This fix also eliminates errors when
  unpickling WCS objects originally created using non-default values for
  ``key``, ``colsel``, and ``keysel`` parameters. [#12844]

Version 5.0.1 (2022-01-26)
==========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Trying to create an instance of ``astropy.coordinates.Distance`` by providing
  both ``z`` and ``parallax`` now raises the expected ``ValueError``. [#12531]

- Fixed a bug where changing the wrap angle of the longitude component of a
  representation could raise a warning or error in certain situations. [#12556]

- ``astropy.coordinates.Distance`` constructor no longer ignores the ``unit``
  keyword when ``parallax`` is provided. [#12569]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- ``astropy.cosmology.utils.aszarr`` can now convert ``Column`` objects. [#12525]

- Reading a cosmology from an ECSV will load redshift and Hubble parameter units
  from the cosmology units module. [#12636]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fix formatting issue in ``_dump_coldefs`` and add tests for ``tabledump`` and
  ``tableload`` convenience functions. [#12526]

astropy.io.misc
^^^^^^^^^^^^^^^

- YAML can now also represent quantities and arrays with structured dtype,
  as well as structured scalars based on ``np.void``. [#12509]

astropy.modeling
^^^^^^^^^^^^^^^^

- Fixes error when fitting multiplication or division based compound models
  where the sub-models have different output units. [#12475]

- Bugfix for incorrectly initialized and filled ``parameters`` data for ``Spline1D`` model. [#12523]

- Bugfix for ``keyerror`` thrown by ``Model.input_units_equivalencies`` when
  used on ``fix_inputs`` models which have no set unit equivalencies. [#12597]

astropy.table
^^^^^^^^^^^^^

- ``astropy.table.Table.keep_columns()`` and
  ``astropy.table.Table.remove_columns()`` now work with generators of column
  names. [#12529]

- Avoid duplicate storage of info in serialized columns if the column
  used to serialize already can hold that information. [#12607]

astropy.timeseries
^^^^^^^^^^^^^^^^^^

- Fixed edge case bugs which emerged when using ``aggregate_downsample`` with custom bins. [#12527]

astropy.units
^^^^^^^^^^^^^

- Structured units can be serialized to/from yaml. [#12492]

- Fix bad typing problems by removing interaction with ``NDArray.__class_getitem__``. [#12511]

- Ensure that ``Quantity.to_string(format='latex')`` properly typesets exponents
  also when ``u.quantity.conf.latex_array_threshold = -1`` (i.e., when the threshold
  is taken from numpy). [#12573]

- Structured units can now be copied with ``copy.copy`` and ``copy.deepcopy``
  and also pickled and unpicked also for ``protocol`` >= 2.
  This does not work for big-endian architecture with older ``numpy<1.21.1``. [#12583]

astropy.utils
^^^^^^^^^^^^^

- Ensure that a ``Masked`` instance can be used to initialize (or viewed
  as) a ``numpy.ma.Maskedarray``. [#12482]

- Ensure ``Masked`` also works with numpy >=1.22, which has a keyword argument
  name change for ``np.quantile``. [#12511]

- ``astropy.utils.iers.LeapSeconds.auto_open()`` no longer emits unnecessary
  warnings when ``astropy.utils.iers.conf.auto_max_age`` is set to ``None``. [#12713]


Version 5.0 (2021-11-15)
========================


New Features
------------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Added dealiasing support to ``convolve_fft``. [#11495]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Added missing coordinate transformations where the starting and ending frames
  are the same (i.e., loopback transformations). [#10909]

- Allow negation, multiplication and division also of representations that
  include a differential (e.g., ``SphericalRepresentation`` with a
  ``SphericalCosLatDifferential``).  For all operations, the outcome is
  equivalent to transforming the representation and differential to cartesian,
  then operating on those, and transforming back to the original representation
  (except for ``UnitSphericalRepresentation``, which will return a
  ``SphericalRepresentation`` if there is a scale change). [#11470]

- ``RadialRepresentation.transform`` can work with a multiplication matrix only.
  All other matrices still raise an exception. [#11576]

- ``transform`` methods are added to ``BaseDifferential`` and ``CartesianDifferential``.
  All transform methods on Representations now delegate transforming differentials
  to the differential objects. [#11654]

- Adds new ``HADec`` built-in frame with transformations to/from ``ICRS`` and ``CIRS``.
  This frame complements ``AltAz`` to give observed coordinates (hour angle and declination)
  in the ``ITRS`` for an equatorially mounted telescope. [#11676]

- ``SkyCoord`` objects now have a ``to_table()`` method, which allows them to be
  converted to a ``QTable``. [#11743]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Cosmologies now store metadata in a mutable parameter ``meta``.
  The initialization arguments ``name`` and ``meta`` are keyword-only. [#11542]

- A new unit, ``redshift``, is defined. It is a dimensionless unit to distinguish
  redshift quantities from other non-redshift values. For compatibility with
  dimensionless quantities the equivalency ``dimensionless_redshift`` is added.
  This equivalency is enabled by default. [#11786]

- Add equality operator for comparing Cosmology instances. Comparison is done on
  all immutable fields (this excludes 'meta').

  Now the following will work:

  .. code-block:: python

      >>> from astropy.cosmology import Planck13, Planck18
      >>> Planck13 == Planck18
      False

      >>> Planck18 == Planck18
      True [#11813]

- Added ``read/write`` methods to Cosmology using the Unified I/O registry.
  Now custom file format readers, writers, and format-identifier functions
  can be registered to read, write, and identify, respectively, Cosmology
  objects. Details are discussed in an addition to the docs. [#11948]

- Added ``to_format/from_format`` methods to Cosmology using the Unified I/O
  registry. Now custom format converters and format-identifier functions
  can be registered to transform Cosmology objects.
  The transformation between Cosmology and dictionaries is pre-registered.
  Details are discussed in an addition to the docs. [#11998]

- Added units module for defining and collecting cosmological units and
  equivalencies. [#12092]

- Flat cosmologies are now set by a mixin class, ``FlatCosmologyMixin`` and its
  FLRW-specific subclass ``FlatFLRWMixin``. All ``FlatCosmologyMixin`` are flat,
  but not all flat cosmologies are instances of ``FlatCosmologyMixin``. As
  example, ``LambdaCDM`` **may** be flat (for the a specific set of parameter
  values),  but ``FlatLambdaCDM`` **will** be flat.

  Cosmology parameters are now descriptors. When accessed from a class they
  transparently stores information, like the units and accepted equivalencies.
  On a cosmology instance, the descriptor will return the parameter value.
  Parameters can have custom ``getter`` methods.

  Cosmological equality is refactored to check Parameters (and the name)
  A new method, ``is_equivalent``, is added to check Cosmology equivalence, so
  a ``FlatLambdaCDM`` and flat ``LambdaCDM`` are equivalent. [#12136]

- Replaced ``z = np.asarray(z)`` with ``z = u.Quantity(z, u.dimensionless_unscaled).value``
  in Cosmology methods. Input of values with incorrect units raises a UnitConversionError
  or TypeError. [#12145]

- Cosmology Parameters allow for custom value setters.
  Values can be set once, but will error if set a second time.
  If not specified, the default setter is used, which will assign units
  using the Parameters ``units`` and ``equivalencies`` (if present).
  Alternate setters may be registered with Parameter to be specified by a str,
  not a decorator on the Cosmology. [#12190]

- Cosmology instance conversion to dict now accepts keyword argument ``cls`` to
  determine dict type, e.g. ``OrderedDict``. [#12209]

- A new equivalency is added between redshift and the Hubble parameter and values
  with units of little-h.
  This equivalency is also available in the catch-all equivalency ``with_redshift``. [#12211]

- A new equivalency is added between redshift and distance -- comoving, lookback,
  and luminosity. This equivalency is also available in the catch-all equivalency
  ``with_redshift``. [#12212]

- Register Astropy Table into Cosmology's ``to/from_format`` I/O, allowing
  a Cosmology instance to be parsed from or converted to a Table instance.
  Also adds the ``__astropy_table__`` method allowing ``Table(cosmology)``. [#12213]

- The WMAP1 and WMAP3 are accessible as builtin cosmologies. [#12248]

- Register Astropy Model into Cosmology's ``to/from_format`` I/O, allowing
  a Cosmology instance to be parsed from or converted to a Model instance. [#12269]

- Register an ECSV reader and writer into Cosmology's I/O, allowing a Cosmology
  instance to be read from from or written to an ECSV file. [#12321]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Added new way to specify the dtype for tables that are read: ``converters``
  can specify column names with wildcards. [#11892]

- Added a new ``astropy.io.ascii.Mrt`` class to write tables in the American
  Astronomical Society Machine-Readable Table format,
  including documentation and tests for the same. [#11897, #12301, #12302]

- When writing, the input data are no longer copied, improving performance.
  Metadata that might be changed, such as format and serialization
  information, is copied, hence users can continue to count on no
  changes being made to the input data. [#11919]

astropy.io.misc
^^^^^^^^^^^^^^^

- Add Parquet serialization of Tables with pyarrow, including metadata support and
  columnar access. [#12215]

astropy.modeling
^^^^^^^^^^^^^^^^

- Added fittable spline models to ``modeling``. [#11634]

- Extensive refactor of ``BoundingBox`` for better usability and maintainability. [#11930]

- Added ``CompoundBoundingBox`` feature to ``~astropy.modeling``, which allows more flexibility in
  defining bounding boxes for models that are applied to images with many slices. [#11942]

- Improved parameter support for ``astropy.modeling.core.custom_model`` created models. [#11984]

- Added the following trigonometric models and linked them to their appropriate inverse models:
    * ``Cosine1D`` [#12158]
    * ``Tangent1D``
    * ``ArcSine1D``
    * ``ArcCosine1D``
    * ``ArcTangent1D`` [#12185]

astropy.table
^^^^^^^^^^^^^

- Added a new method ``Table.update()`` which does a dictionary-style update of a
  ``Table`` by adding or replacing columns. [#11904]

- Masked quantities are now fully supported in tables.  This includes ``QTable``
  automatically converting ``MaskedColumn`` instances to ``MaskedQuantity``,
  and ``Table`` doing the reverse. [#11914]

- Added new keyword arguments ``keys_left`` and ``keys_right`` to the table ``join``
  function to support joining tables on key columns with different names. In
  addition the new keywords can accept a list of column-like objects which are
  used as the match keys. This allows joining on arbitrary data which are not part
  of the tables being joined. [#11954]

- Formatting of any numerical values in the output of ``Table.info()`` and
  ``Column.info()`` has been improved. [#12022]

- It is now possible to add dask arrays as columns in tables
  and have them remain as dask arrays rather than be converted
  to Numpy arrays. [#12219]

- Added a new registry for mixin handlers, which can be used
  to automatically convert array-like Python objects into
  mixin columns when assigned to a table column. [#12219]

astropy.time
^^^^^^^^^^^^

- Adds a new method ``earth_rotation_angle`` to calculate the Local Earth Rotation Angle.
  Also adjusts Local Sidereal Time for the Terrestrial Intermediate Origin (``TIO``)
  and adds a rigorous correction for polar motion. The ``TIO`` adjustment is approximately
  3 microseconds per century from ``J2000`` and the polar motion correction is at most
  about +/-50 nanoseconds. For models ``IAU1982`` and ``IAU1994``, no such adjustments are
  made as they pre-date the TIO concept. [#11680]

astropy.timeseries
^^^^^^^^^^^^^^^^^^

- A custom binning scheme is now available in ``aggregate_downsample``.
  It allows ``time_bin_start`` and ``time_bin_size`` to be arrays, and adds
  an optional ``time_bin_end``.
  This scheme mirrors the API for ``BinnedTimeSeries``. [#11266]

astropy.units
^^^^^^^^^^^^^

- ``Quantity`` gains a ``__class_getitem__`` to create unit-aware annotations
   with the syntax ``Quantity[unit or physical_type, shape, numpy.dtype]``.
   If the python version is 3.9+ or ``typing_extensions`` is installed,
   these are valid static type annotations. [#10662]

- Each physical type is added to ``astropy.units.physical``
  (e.g., ``physical.length`` or ``physical.electrical_charge_ESU``).
  The attribute-accessible names (underscored, without parenthesis) also
  work with ``astropy.units.physical.get_physical_type``. [#11691]

- It is now possible to have quantities based on structured arrays in
  which the unit has matching structure, giving each field its own unit,
  using units constructed like ``Unit('AU,AU/day')``. [#11775]

- The milli- prefix has been added to ``astropy.units.Angstrom``. [#11788]

- Added attributes ``base``, ``coords``, and ``index`` and method ``copy()`` to
  ``QuantityIterator`` to match ``numpy.ndarray.flatiter``. [#11796]

- Added "angular frequency" and "angular velocity" as aliases for the "angular
  speed" physical type. [#11865]

- Add light-second to units of length [#12128]

astropy.utils
^^^^^^^^^^^^^

- The ``astropy.utils.deprecated_renamed_argument()`` decorator now supports
  custom warning messages. [#12305]

- The NaN-aware numpy functions such as ``np.nansum`` now work on Masked
  arrays, with masked values being treated as NaN, but without raising
  warnings or exceptions. [#12454]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Added a feature so that SphericalCircle will accept center parameter as a SkyCoord object. [#11790]

astropy.wcs
^^^^^^^^^^^

- ``astropy.wcs.utils.obsgeo_to_frame`` has been added to convert the obsgeo coordinate
  array on ``astropy.wcs.WCS`` objects to an ``ITRS`` coordinate frame instance. [#11716]

- Updated bundled ``WCSLIB`` to version 7.7 with several bugfixes. [#12034]


API Changes
-----------

astropy.config
^^^^^^^^^^^^^^

- ``update_default_config`` and ``ConfigurationMissingWarning`` are deprecated. [#11502]

astropy.constants
^^^^^^^^^^^^^^^^^

- Removed deprecated ``astropy.constants.set_enabled_constants`` context manager. [#12105]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Positions for the Moon using the 'builtin' ephemeris now use the new
  ``erfa.moon98`` function instead of our own implementation of the Meeus
  algorithm. As this also corrects a misunderstanding of the frame returned by
  the Meeus, this improves the agreement with the JPL ephemeris from about 30 to
  about 6 km rms. [#11753]

- Removed deprecated ``representation`` attribute from
  ``astropy.coordinates.BaseCoordinateFrame`` class. [#12257]

- ``SpectralQuantity`` and ``SpectralCoord`` ``.to_value`` method can now be called without
  ``unit`` argument in order to maintain a consistent interface with ``Quantity.to_value`` [#12440]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- ``z_at_value`` now works with arrays for all arguments (except ``func``,
  ``verbose``, and ``method``). Consequently, ``coordinates.Distance.z`` can
  be used when Distance is an array. [#11778]

- Remove deprecation warning and error remapping in ``Cosmology.clone``.
  Now unknown arguments will raise a ``TypeError``, not an ``AttributeError``. [#11785]

- The ``read/write`` and ``to/from_format`` Unified I/O registries are separated
  and apply only to ``Cosmology``. [#12015]

- Cosmology parameters in ``cosmology.parameters.py`` now have units,
  where applicable. [#12116]

- The function ``astropy.cosmology.utils.inf_like()`` is deprecated. [#12175]

- The function ``astropy.cosmology.utils.vectorize_if_needed()`` is deprecated.
  A new function ``astropy.cosmology.utils.vectorize_redshift_method()`` is added
  as replacement. [#12176]

- Cosmology base class constructor now only accepts arguments ``name`` and ``meta``.
  Subclasses should add relevant arguments and not pass them to the base class. [#12191]

astropy.io
^^^^^^^^^^

- When ``astropy`` raises an ``OSError`` because a file it was told to write
  already exists, the error message now always suggests the use of the
  ``overwrite=True`` argument. The wording is now consistent for all I/O formats. [#12179]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Removed deprecated ``overwrite=None`` option for
  ``astropy.io.ascii.ui.write()``. Overwriting existing files now only happens if
  ``overwrite=True``. [#12171]

astropy.io.fits
^^^^^^^^^^^^^^^

- The internal class _CardAccessor is no longer registered as a subclass of
  the Sequence or Mapping ABCs. [#11923]

- The deprecated ``clobber`` argument will be removed from the
  ``astropy.io.fits`` functions in version 5.1, and the deprecation warnings now
  announce that too. [#12311]

astropy.io.registry
^^^^^^^^^^^^^^^^^^^

- The ``write`` function now is allowed to return possible content results, which
  means that custom writers could, for example, create and return an instance of
  some container class rather than a file on disk. [#11916]

- The registry functions are refactored into a class-based system.
  New Read-only, write-only, and read/write registries can be created.
  All functions accept a new argument ``registry``, which if not specified,
  defaults to the global default registry. [#12015]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Deprecated the ``pedantic`` keyword argument in the
  ``astropy.io.votable.table.parse`` function and the corresponding configuration
  setting. It has been replaced by the ``verify`` option. [#12129]

astropy.modeling
^^^^^^^^^^^^^^^^

- Refactored how ``astropy.modeling.Model`` handles model evaluation in order to better
  organize the code. [#11931]

- Removed the following deprecated modeling features:
      ``astropy.modeling.utils.ExpressionTree`` class,
      ``astropy.modeling.functional_models.MexicanHat1D`` model,
      ``astropy.modeling.functional_models.MexicanHat2D`` model,
      ``astropy.modeling.core.Model.inputs`` setting in model initialize,
      ``astropy.modeling.core.CompoundModel.inverse`` setting in model initialize, and
      ``astropy.modeling.core.CompoundModel.both_inverses_exist()`` method. [#11978]

- Deprecated the ``AliasDict`` class in ``modeling.utils``. [#12411]

astropy.nddata
^^^^^^^^^^^^^^

- Removed ``block_reduce`` and ``block_replicate`` functions from
  ``nddata.utils``. These deprecated functions in ``nddata.utils`` were
  moved to ``nddata.blocks``. [#12288]

astropy.stats
^^^^^^^^^^^^^

- Removed the following deprecated features from ``astropy.stats``:

  * ``conf`` argument for ``funcs.binom_conf_interval()`` and
    ``funcs.binned_binom_proportion()``,
  * ``conflevel`` argument for ``funcs.poisson_conf_interval()``, and
  * ``conf_lvl`` argument for ``jackknife.jackknife_stats()``. [#12200]

astropy.table
^^^^^^^^^^^^^

- Printing a ``Table`` now shows the qualified class name of mixin columns in the
  dtype header row  instead of "object". This applies to all the ``Table`` formatted output
  methods whenever ``show_dtype=True`` is selected. [#11660]

- The 'overwrite' argument has been added to the jsviewer table writer.
  Overwriting an existing file requires 'overwrite' to be True. [#11853]

- The 'overwrite' argument has been added to the pandas table writers.
  Overwriting an existing file requires 'overwrite' to be True. [#11854]

- The table ``join`` function now accepts only the first four arguments ``left``,
  ``right``, ``keys``, and ``join_type`` as positional arguments. All other
  arguments must be supplied as keyword arguments. [#11954]

- Adding a dask array to a Table will no longer convert
  that dask to a Numpy array, so accessing t['dask_column']
  will now return a dask array instead of a Numpy array. [#12219]

astropy.time
^^^^^^^^^^^^

- Along with the new method ``earth_rotation_angle``, ``sidereal_time`` now accepts
  an ``EarthLocation`` as the ``longitude`` argument. [#11680]

astropy.units
^^^^^^^^^^^^^

- Unit ``littleh`` and equivalency ``with_H0`` have been moved to the
  ``cosmology`` module and are deprecated from ``astropy.units``. [#12092]

astropy.utils
^^^^^^^^^^^^^

- ``astropy.utils.introspection.minversion()`` now uses
  ``importlib.metadata.version()``. Therefore, its ``version_path`` keyword is no
  longer used and deprecated. This keyword will be removed in a future release. [#11714]

- Updated ``utils.console.Spinner`` to better resemble the API of
  ``utils.console.ProgressBar``, including an ``update()`` method and
  iterator support. [#11772]

- Removed deprecated ``check_hashes`` in ``check_download_cache()``. The function also
  no longer returns anything. [#12293]

- Removed unused ``download_cache_lock_attempts`` configuration item in
  ``astropy.utils.data``. Deprecation was not possible. [#12293]

- Removed deprecated ``hexdigest`` keyword from ``import_file_to_cache()``. [#12293]

- Setting ``remote_timeout`` configuration item in ``astropy.utils.data`` to 0 will
  no longer disable download from the Internet; Set ``allow_internet`` configuration
  item to ``False`` instead. [#12293]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Removed deprecated ``imshow_only_kwargs`` keyword from ``imshow_norm``. [#12290]

astropy.wcs
^^^^^^^^^^^

- Move complex logic from ``HighLevelWCSMixin.pixel_to_world`` and
  ``HighLevelWCSMixin.world_to_pixel`` into the helper functions
  ``astropy.wcs.wcsapi.high_level_api.high_level_objects_to_values`` and
  ``astropy.wcs.wcsapi.high_level_api.values_to_high_level_objects`` to allow
  reuse in other places. [#11950]


Bug Fixes
---------

astropy.config
^^^^^^^^^^^^^^

- ``generate_config`` no longer outputs wrong syntax for list type. [#12037]

astropy.constants
^^^^^^^^^^^^^^^^^

- Fixed a bug where an older constants version cannot be set directly after
  astropy import. [#12084]

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Passing an ``array`` argument for any Kernel1D or Kernel2D subclasses (with the
  exception of CustomKernel) will now raise a ``TypeError``. [#11969]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- If a ``Table`` containing a ``SkyCoord`` object as a column is written to a
  FITS, ECSV or HDF5 file then any velocity information that might be present
  will be retained. [#11750]

- The output of ``SkyCoord.apply_space_motion()`` now always has the same
  differential type as the ``SkyCoord`` itself. [#11932]

- Fixed bug where Angle, Latitude and Longitude with NaN values could not be printed. [#11943]

- Fixed a bug with the transformation from ``PrecessedGeocentric`` to ``GCRS``
  where changes in ``obstime``, ``obsgeoloc``, or ``obsgeovel`` were ignored.
  This bug would also affect loopback transformations from one ``PrecessedGeocentric``
  frame to another ``PrecessedGeocentric`` frame. [#12152]

- Fixed a bug with the transformations between ``TEME`` and ``ITRS`` or between ``TEME``
  and itself where a change in ``obstime`` was ignored. [#12152]

- Avoid unnecessary transforms through CIRS for AltAz and HADec and
  use ICRS as intermediate frame for these transformations instead. [#12203]

- Fixed a bug where instantiating a representation with a longitude component
  could mutate input provided for that component even when copying is specified. [#12307]

- Wrapping an ``Angle`` array will now ignore NaN values instead of attempting to wrap
  them, which would produce unexpected warnings/errors when working with coordinates
  and representations due to internal broadcasting. [#12317]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Dictionaries for in-built cosmology realizations are not altered by creating
  the realization and are also made immutable. [#12278]

astropy.io.fits
^^^^^^^^^^^^^^^

- Prevent zero-byte writes for FITS binary tables to
  speed up writes on the Lustre filesystem. [#11955]

- Enable ``json.dump`` for FITS_rec with variable length (VLF) arrays. [#11957]

- Add support for reading and writing int8 images [#11996]

- Ensure header passed to ``astropy.io.fits.CompImageHDU`` does not need to contain
  standard cards that can be automatically generated, such as ``BITPIX`` and ``NAXIS``. [#12061]

- Fixed a bug where ``astropy.io.fits.HDUDiff`` would ignore the ``ignore_blank_cards``
  keyword argument. [#12122]

- Open uncompressed file even if extension says it's compressed [#12135]

- Fix the computation of the DATASUM in a ``CompImageHDU`` when the data is >1D. [#12138]

- Reading files where the SIMPLE card is present but with an invalid format now
  issues a warning instead of raising an exception [#12234]

- Convert UNDEFINED to None when iterating over card values. [#12310]

astropy.io.misc
^^^^^^^^^^^^^^^

- Update ASDF tag versions in ExtensionType subclasses to match ASDF Standard 1.5.0. [#11986]

- Fix ASDF serialization of model inputs and outputs and add relevant assertion to
  test helper. [#12381]

- Fix bug preventing ASDF serialization of bounding box for models with only one input. [#12385]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Now accepting UCDs containing phot.color. [#11982]

astropy.modeling
^^^^^^^^^^^^^^^^

- Added ``Parameter`` descriptions to the implemented models which were
  missing. [#11232]

- The ``separable`` property is now correctly set on models constructed with
  ``astropy.modeling.custom_model``. [#11744]

- Minor bugfixes and improvements to modeling including the following:
      * Fixed typos and clarified several errors and their messages throughout
        modeling.
      * Removed incorrect try/except blocks around scipy code in
        ``convolution.py`` and ``functional_models.py``.
      * Fixed ``Ring2D`` model's init to properly accept all combinations
        of ``r_in``, ``r_out``, and ``width``.
      * Fixed bug in ``tau`` validator for the ``Logarithmic1D`` and
        ``Exponential1D`` models when using them as model sets.
      * Fixed ``copy`` method for ``Parameter`` in order to prevent an
        automatic ``KeyError``, and fixed ``bool`` for ``Parameter`` so
        that it functions with vector values.
      * Removed unreachable code from ``Parameter``, the ``_Tabular`` model,
        and the ``Drude1D`` model.
      * Fixed validators in ``Drude1D`` model so that it functions in a
        model set.
      * Removed duplicated code from ``polynomial.py`` for handing of
        ``domain`` and ``window``.
      * Fixed the ``Pix2Sky_HEALPixPolar`` and ``Sky2Pix_HEALPixPolar`` modes
        so that their ``evaluate`` and ``inverse`` methods actually work
        without raising an error. [#12232]

astropy.nddata
^^^^^^^^^^^^^^

- Ensure that the ``wcs=`` argument to ``NDData`` is always parsed into a high
  level WCS object. [#11985]

astropy.stats
^^^^^^^^^^^^^

- Fixed a bug in sigma clipping where the bounds would not be returned for
  completely empty or masked data. [#11994]

- Fixed a bug in ``biweight_midvariance`` and ``biweight_scale`` where
  output data units would be dropped for constant data and where the
  result was a scalar NaN. [#12146]

astropy.table
^^^^^^^^^^^^^

- Ensured that ``MaskedColumn.info`` is propagated in all cases, so that when
  tables are sliced, writing will still be as requested on
  ``info.serialize_method``. [#11917]

- ``table.conf.replace_warnings`` and ``table.jsviewer.conf.css_urls`` configuration
  items now have correct ``'string_list'`` type. [#12037]

- Fixed an issue where initializing from a list of dict-like rows (Mappings) did
  not work unless the row values were instances of ``dict``. Now any object that
  is an instance of the more general ``collections.abc.Mapping`` will work. [#12417]

astropy.uncertainty
^^^^^^^^^^^^^^^^^^^

- Ensure that scalar ``QuantityDistribution`` unit conversion in ufuncs
  works properly again. [#12471]

astropy.units
^^^^^^^^^^^^^

- Add quantity support for ``scipy.special`` dimensionless functions
  erfinv, erfcinv, gammaln and loggamma. [#10934]

- ``VOUnit.to_string`` output is now compliant with IVOA VOUnits 1.0 standards. [#11565]

- Units initialization with unicode has been expanded to include strings such as
  'M☉' and 'e⁻'. [#11827]

- Give a more informative ``NotImplementedError`` when trying to parse a unit
  using an output-only format such as 'unicode' or 'latex'. [#11829]

astropy.utils
^^^^^^^^^^^^^

- Fixed a bug in ``get_readable_fileobj`` that prevented the unified file read
  interface from closing ASCII files. [#11809]

- The function ``astropy.utils.decorators.deprecated_attribute()`` no longer
  ignores its ``message``, ``alternative``, and ``pending`` arguments. [#12184]

- Ensure that when taking the minimum or maximum of a ``Masked`` array,
  any masked NaN values are ignored. [#12454]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- The tick labelling for radians has been fixed to remove a redundant ``.0`` in
  the label for integer multiples of pi at 2pi and above. [#12221]

- Fix a bug where non-``astropy.wcs.WCS`` WCS instances were not accepted in
  ``WCSAxes.get_transform``. [#12286]

- Fix compatibility with Matplotlib 3.5 when using the ``grid_type='contours'``
  mode for drawing grid lines. [#12447]

astropy.wcs
^^^^^^^^^^^

- Enabled ``SlicedLowLevelWCS.pixel_to_world_values`` to handle slices including
  non-``int`` integers, e.g. ``numpy.int64``. [#11980]


Other Changes and Additions
---------------------------

- In docstrings, Sphinx cross-reference targets now use intersphinx, even if the
  target is an internal link (``link`` is now ``'astropy:link``).
  When built in Astropy these links are interpreted as internal links. When built
  in affiliate packages, the link target is set by the key 'astropy' in the
  intersphinx mapping. [#11690]

- Made PyYaml >= 3.13 a strict runtime dependency. [#11903]

- Minimum version of required Python is now 3.8. [#11934]

- Minimum version of required Scipy is now 1.3. [#11934]

- Minimum version of required Matplotlib is now 3.1. [#11934]

- Minimum version of required Numpy is now 1.18. [#11935]

- Fix deprecation warnings with Python 3.10 [#11962]

- Speed up ``minversion()`` in cases where a module with a ``__version__``
  attribute is passed. [#12174]

- ``astropy`` now requires ``packaging``. [#12199]

- Updated the bundled CFITSIO library to 4.0.0. When compiling with an external
  library, version 3.35 or later is required. [#12272]


Version 4.3.1 (2021-08-11)
==========================

Bug Fixes
---------

astropy.io.fits
^^^^^^^^^^^^^^^

- In ``fits.io.getdata`` do not fall back to first non-primary extension when
  user explicitly specifies an extension. [#11860]

- Ensure multidimensional masked columns round-trip properly to FITS. [#11911]

- Ensure masked times round-trip to FITS, even if multi-dimensional. [#11913]

- Raise ``ValueError`` if an ``np.float32`` NaN/Inf value is assigned to a
  header keyword. [#11922]

astropy.modeling
^^^^^^^^^^^^^^^^

- Fixed bug in ``fix_inputs`` handling of bounding boxes. [#11908]

astropy.table
^^^^^^^^^^^^^

- Fix an error when converting to pandas any ``Table`` subclass that
  automatically adds a table index when the table is created. An example is a
  binned ``TimeSeries`` table. [#12018]

astropy.units
^^^^^^^^^^^^^

- Ensure that unpickling quantities and units in new sessions does not change
  hashes and thus cause problems with (de)composition such as getting different
  answers from the ``.si`` attribute. [#11879]

- Fixed cannot import name imperial from astropy.units namespace. [#11977]

astropy.utils
^^^^^^^^^^^^^

- Ensure any ``.info`` on ``Masked`` instances is propagated correctly when
  viewing or slicing. As a consequence, ``MaskedQuantity`` can now be correctly
  written to, e.g., ECSV format with ``serialize_method='data_mask'``. [#11910]


Version 4.3 (2021-07-26)
========================

New Features
------------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Change padding sizes for ``fft_pad`` in ``convolve_fft`` from powers of
  2 only to scipy-optimized numbers, applied separately to each dimension;
  yielding some performance gains and avoiding potential large memory
  impact for certain multi-dimensional inputs. [#11533]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Adds the ability to create topocentric ``CIRS`` frames. Using these,
  ``AltAz`` calculations are now accurate down to the milli-arcsecond
  level. [#10994]

- Adds a direct transformation from ``ICRS`` to ``AltAz`` frames. This
  provides a modest speedup of approximately 10 percent. [#11079]

- Adds new ``WGS84GeodeticRepresentation``, ``WGS72GeodeticRepresentation``,
  and ``GRS80GeodeticRepresentation``. These are mostly for use inside
  ``EarthLocation`` but can also be used to convert between geocentric
  (cartesian) and different geodetic representations directly. [#11086]

- ``SkyCoord.guess_from_table`` now also searches for differentials in the table.
  In addition, multiple regex matches can be resolved when they are exact
  component names, e.g. having both columns “dec” and “pm_dec” no longer errors
  and will be included in the SkyCoord. [#11417]

- All representations now have a ``transform`` method, which allows them to be
  transformed by a 3x3 matrix in a Cartesian basis. By default, transformations
  are routed through ``CartesianRepresentation``. ``SphericalRepresentation`` and
  ``PhysicssphericalRepresentation`` override this for speed and to prevent NaN
  leakage from the distance to the angular components.
  Also, the functions ``is_O3`` and ``is_rotation`` have been added to
  ``matrix_utities`` for checking whether a matrix is in the O(3) group or is a
  rotation (proper or improper), respectively. [#11444]

- Moved angle formatting and parsing utilities to
  ``astropy.coordinates.angle_formats``.
  Added new functionality to ``astropy.coordinates.angle_utilities`` for
  generating points on or in spherical surfaces, either randomly or on a grid. [#11628]

- Added a new method to ``SkyCoord``, ``spherical_offsets_by()``, which is the
  conceptual inverse of ``spherical_offsets_to()``: Given angular offsets in
  longitude and latitude, this method returns a new coordinate with the offsets
  applied. [#11635]

- Refactor conversions between ``GCRS`` and ``CIRS,TETE`` for better accuracy
  and substantially improved speed. [#11069]

- Also refactor ``EarthLocation.get_gcrs`` for an increase in performance of
  an order of magnitude, which enters as well in getting observed positions of
  planets using ``get_body``. [#11073]

- Refactored the usage of metaclasses in ``astropy.coordinates`` to instead use
  ``__init_subclass__`` where possible. [#11090]

- Removed duplicate calls to ```transform_to``` from ```match_to_catalog_sky```
  and ```match_to_catalog_3d```, improving their performance. [#11449]

- The new DE440 and DE440s ephemerides are now available via shortcuts 'de440'
  and 'de440s'.  The DE 440s ephemeris will probably become the default
  ephemeris when choosing 'jpl' in 5.0. [#11601]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Cosmology parameter dictionaries now also specify the Cosmology class to which
  the parameters correspond. For example, the dictionary for
  ``astropy.cosmology.parameters.Planck18`` has the added key-value pair
  ("cosmology", "FlatLambdaCDM"). [#11530]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Added support for reading and writing ASCII tables in QDP (Quick and Dandy
  Plotter) format. [#11256]

- Added support for reading and writing multidimensional column data (masked and
  unmasked) to ECSV. Also added formal support for reading and writing object-type
  column data which can contain items consisting of lists, dicts, and basic scalar
  types. This can be used to store columns of variable-length arrays. Both of
  these features use JSON to convert the object to a string that is stored in the
  ECSV output. [#11569, #11662, #11720]

astropy.io.fits
^^^^^^^^^^^^^^^

- Added ``append`` keyword to append table objects to an existing FITS file [#2632, #11149]

- Check that the SIMPLE card is present when opening a file, to ensure that the
  file is a valid FITS file and raise a better error when opening a non FITS
  one. ``ignore_missing_simple`` can be used to skip this verification. [#10895]

- Expose ``Header.strip`` as a public method, to remove the most common
  structural keywords. [#11174]

- Enable the use of ``os.PathLike`` objects when dealing with (mainly FITS) files. [#11580]

astropy.io.registry
^^^^^^^^^^^^^^^^^^^

- Readers and writers can now set a priority, to assist with resolving which
  format to use. [#11214]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Version 1.4 VOTables now use the VOUnit format specification. [#11032]

- When reading VOTables using the Unified File Read/Write Interface (i.e. using
  the ``Table.read()`` or ``QTable.read()`` functions) it is now possible to
  specify all keyword arguments that are valid for
  ``astropy.io.votable.table.parse()``. [#11643]

astropy.modeling
^^^^^^^^^^^^^^^^

- Added a state attribute to models to allow preventing the syncing of
  constraint values from the constituent models. This syncing can
  greatly slow down fitting if there are large numbers of fit parameters.
  model.sync_constraints = True means check constituent model constraints
  for compound models every time the constraint is accessed, False, do not.
  Fitters that support constraints will set this to False on the model copy
  and then set back to True when the fit is complete before returning. [#11365]

- The ``convolve_models_fft`` function implements model convolution so that one
  insures that the convolution remains consistent across multiple different
  inputs. [#11456]

astropy.nddata
^^^^^^^^^^^^^^

- Prevent unnecessary copies of the data during ``NDData`` arithmetic when units
  need to be added. [#11107]

- NDData str representations now show units, if present. [#11553]

astropy.stats
^^^^^^^^^^^^^

- Added the ability to specify stdfunc='mad_std' when doing sigma clipping,
  which will use a built-in function and lead to significant performance
  improvements if cenfunc is 'mean' or 'median'. [#11664]


- Significantly improved the performance of sigma clipping when cenfunc and
  stdfunc are passed as strings and the ``grow`` option is not used. [#11219]

- Improved performance of ``bayesian_blocks()`` by removing one ``np.log()``
  call [#11356]

astropy.table
^^^^^^^^^^^^^

- Add table attributes to include or exclude columns from the output when
  printing a table. This functionality includes a context manager to
  include/exclude columns temporarily. [#11190]

- Improved the string representation of objects related to ``Table.indices`` so
  they now indicate the object type and relevant attributes. [#11333]

astropy.timeseries
^^^^^^^^^^^^^^^^^^

- An exception is raised when ``n_bins`` is passed as an argument while
  any of the parameters ``time_bin_start`` or ``time_bin_size`` is not
  scalar. [#11463]

astropy.units
^^^^^^^^^^^^^

- The ``physical_type`` attributes of each unit are now objects of the (new)
  ``astropy.units.physical.PhysicalType`` class instead of strings and the
  function ``astropy.units.physical.get_physical_type`` can now translate
  strings to these objects. [#11204]

-  The function ``astropy.units.physical.def_physical_type`` was created to
   either define entirely new physical types, or to add more physical type
   names to an existing physical types. [#11204]

- ``PhysicalType``'s can be operated on using operations multiplication,
  division, and exponentiation are to facilitate dimensional analysis. [#11204]

- It is now possible to define aliases for units using
  ``astropy.units.set_enabled_aliases``. This can be used when reading files
  that have misspelled units. [#11258]

- Add a new "DN" unit, ``units.dn`` or ``units.DN``, representing data number
  for a detector. [#11591]

astropy.utils
^^^^^^^^^^^^^

- Added ``ssl_context`` and ``allow_insecure`` options to ``download_file``,
  as well as the ability to optionally use the ``certifi`` package to provide
  root CA certificates when downloading from sites secured with
  TLS/SSL. [#10434]

- ``astropy.utils.data.get_pkg_data_path`` is publicly scoped (previously the
  private function ``_find_pkg_data_path``) for obtaining file paths without
  checking if the file/directory exists, as long as the package and module
  do. [#11006]

- Deprecated ``astropy.utils.OrderedDescriptor`` and
  ``astropy.utils.OrderedDescriptorContainer``, as new features in Python 3
  make their use less compelling. [#11094, #11099]

- ``astropy.utils.masked`` provides a new ``Masked`` class/factory that can be
  used to represent masked ``ndarray`` and all its subclasses, including
  ``Quantity`` and its subclasses.  These classes can be used inside
  coordinates, but the mask is not yet exposed.  Generally, the interface should
  be considered experimental. [#11127, #11792]

- Add new ``utils.parsing`` module to with helper wrappers around
  ``ply``. [#11227]

- Change the Time and IERS leap second handling so that the leap second table is
  updated only when a Time transform involving UTC is performed. Previously this
  update check was done the first time a ``Time`` object was created, which in
  practice occurred when importing common astropy subpackages like
  ``astropy.coordinates``. Now you can prevent querying internet resources (for
  instance on a cluster) by setting ``iers.conf.auto_download = False``. This
  can  be done after importing astropy but prior to performing any ``Time``
  scale transformations related to UTC. [#11638]


- Added a new module at ``astropy.utils.compat.optional_deps`` to consolidate
  the definition of ``HAS_x`` optional dependency flag variables,
  like ``HAS_SCIPY``. [#11490]

astropy.wcs
^^^^^^^^^^^

- Add IVOA UCD mappings for some FITS WCS keywords commonly used in solar
  physics. [#10965]

- Add ``STOKES`` FITS WCS keyword to the IVOA UCD mapping. [#11236]

- Updated bundled version of WCSLIB to version 7.6. See
  https://www.atnf.csiro.au/people/mcalabre/WCS/CHANGES for a list of
  included changes. [#11549]


API Changes
-----------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- For input to representations, subclasses of the class required for a
  given attribute will now be allowed in. [#11113]

- Except for ``UnitSphericalRepresentation``, shortcuts in representations now
  allow for attached differentials. [#11467]

- Allow coordinate name strings as input to
  ``SkyCoord.is_transformable_to``. [#11552]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Change ``z_at_value`` to use ``scipy.optimize.minimize_scalar`` with default
  method ``Brent`` (other options ``Bounded`` and ``Golden``) and accept
  ``bracket`` option to set initial search region. [#11080]

- Clarified definition of inputs to ``angular_diameter_distance_z1z2``.
  The function now emits ``AstropyUserWarning`` when ``z2`` is less than
  ``z1``. [#11197]

- Split cosmology realizations from core classes, moving the former to new file
  ``realizations``. [#11345]

- Since cosmologies are immutable, the initialization signature and values can
  be stored, greatly simplifying cloning logic and extending it to user-defined
  cosmology classes that do not have attributes with the same name as each
  initialization argument.  [#11515]

- Cloning a cosmology with changed parameter(s) now appends "(modified)" to the
  new instance's name, unless a name is explicitly passed to ``clone``. [#11536]

- Allow ``m_nu`` to be input as any quantity-like or array-like -- Quantity,
  array, float, str, etc. Input is passed to the Quantity constructor and
  converted to eV, still with the prior mass-energy equivalence
  enabled. [#11640]

astropy.io.fits
^^^^^^^^^^^^^^^

- For conversion between FITS tables and astropy ``Table``, the standard mask
  values of ``NaN`` for float and null string for string are now properly
  recognized, leading to a ``MaskedColumn`` with appropriately set mask
  instead of a ``Column`` with those values exposed. Conversely, when writing
  an astropy ``Table`` to a FITS tables, masked values are now consistently
  converted to the standard FITS mask values of ``NaN`` for float and null
  string for string (i.e., not just for tables with ``masked=True``, which no
  longer is guaranteed to signal the presence of ``MaskedColumn``). [#11222]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- The use of ``version='1.0'`` is now fully deprecated in constructing
  a ``astropy.io.votable.tree.VOTableFile``. [#11659]

astropy.modeling
^^^^^^^^^^^^^^^^

- Removed deprecated ``astropy.modeling.blackbody`` module. [#10972]

astropy.table
^^^^^^^^^^^^^

- Added ``Column.value`` as an alias for the existing ``Column.data`` attribute.
  This makes accessing a column's underlying data array consistent with the
  ``.value`` attribute available for ``Time`` and ``Quantity`` objects. [#10962]

- In reading from a FITS tables, the standard mask values of ``NaN`` for float
  and null string for string are properly recognized, leading to a
  ``MaskedColumn`` with appropriately set mask. [#11222]

- Changed the implementation of the ``table.index.Index`` class so instantiating
  from this class now returns an ``Index`` object as expected instead of a
  ``SlicedIndex`` object. [#11333]

astropy.units
^^^^^^^^^^^^^

- The ``physical_type`` attribute of units now returns an instance of
  ``astropy.units.physical.PhysicalType`` instead of a string.  Because
  ``PhysicalType`` instances can be compared to strings, no code changes
  should be necessary when making comparisons.  The string representations
  of different physical types will differ from previous releases. [#11204]

- Calling ``Unit()`` with no argument now returns a dimensionless unit,
  as was documented but not implemented. [#11295]

astropy.utils
^^^^^^^^^^^^^

- Removed deprecated ``utils.misc.InheritDocstrings`` and ``utils.timer``. [#10281]

- Removed usage of deprecated ``ipython`` stream in ``utils.console``. [#10942]

astropy.wcs
^^^^^^^^^^^

- Deprecate ``accuracy`` argument in ``all_world2pix`` which was mistakenly
  *documented*, in the case ``accuracy`` was ever used. [#11055]


Bug Fixes
---------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Fixes for ``convolve_fft`` documentation examples. [#11510]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Allow ``Distance`` instances with negative distance values as input for
  ``SphericalRepresentation``.  This was always noted as allowed in an
  exception message when a negative ``Quantity`` with length units was
  passed in, but was not actually possible to do. [#11113]

- Makes the ``Angle.to_string`` method to follow the format described in the
  docstring with up to 8 significant decimals instead of 4. [#11153]

- Ensure that proper motions can be calculated when converting a ``SkyCoord``
  with cartesian representation to unit-spherical, by fixing the conversion of
  ``CartesianDifferential`` to ``UnitSphericalDifferential``. [#11469]

- When re-representing coordinates from spherical to unit-spherical and vice
  versa, the type of differential will now be preserved. For instance, if only a
  radial velocity was present, that will remain the case (previously, a zero
  proper motion component was added). [#11482]

- Ensure that wrapping of ``Angle`` does not raise a warning even if ``nan`` are
  present.  Also try to make sure that the result is within the wrapped range
  even in the presence of rounding errors. [#11568]

- Comparing a non-SkyCoord object to a ``SkyCoord`` using ``==`` no longer
  raises an error. [#11666]

- Different ``SkyOffsetFrame`` classes no longer interfere with each other,
  causing difficult to debug problems with the ``origin`` attribute. The
  ``origin`` attribute now no longer is propagated, so while it remains
  available on a ``SkyCoord`` that is an offset, it no longer is available once
  that coordinate is transformed to another frame. [#11730] [#11730]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Cosmology instance names are now immutable. [#11535]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fixed bug where writing a table that has comments defined (via
  ``tbl.meta['comments']``) with the 'csv' format was failing. Since the
  formally defined CSV format does not support comments, the comments are now
  just ignored unless ``comment=<comment prefix>`` is supplied to the
  ``write()`` call. [#11475]

- Fixed the issue where the CDS reader failed to treat columns
  as nullable if the ReadMe file contains a limits specifier. [#11531]

- Made sure that the CDS reader does not ignore an order specifier that
  may be present after the null specifier '?'. Also made sure that it
  checks null values only when an '=' symbol is present and reads
  description text even if there is no whitespace after '?'. [#11593]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fix ``ColDefs.add_col/del_col`` to allow in-place addition or removal of
  a column. [#11338]

- Fix indexing of ``fits.Header`` with Numpy integers. [#11387]

- Do not delete ``EXTNAME`` for compressed image header if a default and
  non-default ``EXTNAME`` are present. [#11396]

- Prevent warnings about ``HIERARCH`` with ``CompImageHeader`` class. [#11404]

- Fixed regression introduced in Astropy 4.0.5 and 4.2.1 with verification of
  FITS headers with HISTORY or COMMENT cards with long (> 72 characters)
  values. [#11487]

- Fix reading variable-length arrays when there is a gap between the data and the
  heap. [#11688]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- ``NumericArray`` converter now properly broadcasts scalar mask to array. [#11157]

- VOTables are now written with the correct namespace and schema location
  attributes. [#11659]

astropy.modeling
^^^^^^^^^^^^^^^^

- Fixes the improper propagation of ``bounding_box`` from
  ``astropy.modeling.models`` to their inverses. For cases in which the inverses
  ``bounding_box`` can be determined, the proper calculation has been
  implemented. [#11414]

- Bugfix to allow rotation models to accept arbitrarily-shaped inputs. [#11435]

- Bugfixes for ``astropy.modeling`` to allow ``fix_inputs`` to accept empty
  dictionaries and dictionaries with ``numpy`` integer keys. [#11443]

- Bugfix for how ``SPECIAL_OPERATORS`` are handled. [#11512]

- Fixes ``Model`` crashes when some inputs are scalars and during some types of
  output reshaping. [#11548]

- Fixed bug in ``LevMarLSQFitter`` when using weights and vector inputs. [#11603]

astropy.stats
^^^^^^^^^^^^^

- Fixed a bug with the ``copy=False`` option when carrying out sigma
  clipping - previously if ``masked=False`` this still copied the data,
  but this will now change the array in-place. [#11219]

astropy.table
^^^^^^^^^^^^^

- Ensure that adding a ``Quantity`` or other mixin column to a ``Table``
  does not have side effects, such as creating an associated ``info``
  instance (which would lead to slow-down of, e.g., slicing afterwards). [#11077]

- When writing to a FITS tables, masked values are again always converted to
  the standard FITS mask values of ``NaN`` for float and null string
  for string, not just for table with ``masked=True``. [#11222]

- Using ``Table.to_pandas()`` on an indexed ``Table`` with masked integer values
  now correctly construct the ``pandas.DataFrame``. [#11432]

- Fixed ``Table`` HTML representation in Jupyter notebooks so that it is
  horizontally scrollable within Visual Studio Code. This was done by wrapping
  the ``<table>`` in a ``<div>`` element. [#11476]

- Fix a bug where a string-valued ``Column`` that happened to have a ``unit``
  attribute could not be added to a ``QTable``.  Such columns are now simply
  kept as ``Column`` instances (with a warning). [#11585]

- Fix an issue in ``Table.to_pandas(index=<colname>)`` where the index column name
  was not being set properly for the ``DataFrame`` index. This was introduced by
  an API change in pandas version 1.3.0. Previously when creating a ``DataFrame``
  with the index set to an astropy ``Column``, the ``DataFrame`` index name was
  automatically set to the column name. [#11921]

astropy.time
^^^^^^^^^^^^

- Fix a thread-safety issue with initialization of the leap-second table
  (which is only an issue when ERFA's built-in table is out of date). [#11234]

- Fixed converting a zero-length time object from UTC to
  UT1 when an empty array is passed. [#11516]

astropy.uncertainty
^^^^^^^^^^^^^^^^^^^

- ``Distribution`` instances can now be used as input to ``Quantity`` to
  initialize ``QuantityDistribution``.  Hence, ``distribution * unit``
  and ``distribution << unit`` will work too. [#11210]

astropy.units
^^^^^^^^^^^^^

- Move non-astronomy units from astrophys.py to a new misc.py file. [#11142]

- The physical type of ``astropy.units.mol / astropy.units.m ** 3`` is now
  defined as molar concentration.  It was previously incorrectly defined
  as molar volume. [#11204]

- Make ufunc helper lookup thread-safe. [#11226]

- Make ``Unit`` string parsing (as well as ``Angle`` parsing) thread-safe. [#11227]

- Decorator ``astropy.units.decorators.quantity_input`` now only evaluates
  return type annotations based on ``UnitBase`` or ``FunctionUnitBase`` types.
  Other annotations are skipped over and are not attempted to convert to the
  correct type. [#11506]

astropy.utils
^^^^^^^^^^^^^

- Make ``lazyproperty`` and ``classdecorator`` thread-safe. This should fix a
  number of thread safety issues. [#11224]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Fixed a bug that resulted in some parts of grid lines being visible when they
  should have been hidden. [#11380]

- Fixed a bug that resulted in ``time_support()`` failing for intervals of
  a few months if one of the ticks was the month December. [#11615]

astropy.wcs
^^^^^^^^^^^

- ``fit_wcs_from_points`` now produces a WCS with integer ``NAXIXn``
  values. [#10865]

- Updated bundled version of ``WCSLIB`` to v7.4, fixing a bug that caused
  the coefficients of the TPD distortion function to not be written to the
  header. [#11260]

- Fixed a bug in assigning type when converting ``colsel`` to
  ``numpy.ndarray``. [#11431]

- Added ``WCSCOMPARE_*`` constants to the list of WCSLIB constants
  available/exposed through the ``astropy.wcs`` module. [#11647]

- Fix a bug that caused APE 14 WCS transformations for FITS WCS with ZOPT, BETA,
  VELO, VOPT, or VRAD CTYPE to not work correctly. [#11781]


Other Changes and Additions
---------------------------

- The configuration file is no longer created by default when importing astropy
  and its existence is no longer required. Affiliated packages should update their
  ``__init__.py`` module to remove the block using ``update_default_config`` and
  ``ConfigurationDefaultMissingWarning``. [#10877]

- Replace ``pkg_resources`` (from setuptools) with ``importlib.metadata`` which
  comes from the stdlib, except for Python 3.7 where the backport package is added
  as a new dependency. [#11091]

- Turn on numpydoc's ``numpydoc_xref_param_type``  to create cross-references
  for the parameter types in the Parameters, Other Parameters, Returns and Yields
  sections of the docstrings. [#11118]

- Docstrings across the package are standardized to enable references.
  Also added is an Astropy glossary-of-terms to define standard inputs,
  e.g. ``quantity-like`` indicates an input that can be interpreted by
  ``astropy.units.Quantity``. [#11118]

- Binary wheels are now built to the manylinux2010 specification. These wheels
  should be supported on all versions of pip shipped with Python 3.7+. [#11377]

- The name of the default branch for the astropy git repository has been renamed
  to ``main``, and the documentation and tooling has been updated accordingly.
  If you have made a local clone you may wish to update it following the
  instructions in the repository's README. [#11379]

- Sphinx cross-reference link targets are added for every ``PhysicalType``.
  Now in the parameter types in the Parameters, Other Parameters, Returns and
  Yields sections of the docstring, the physical type of a quantity can be
  annotated in square brackets.
  E.g. `` distance : `~astropy.units.Quantity` ['length'] `` [#11595]

- The minimum supported version of ``ipython`` is now 4.2. [#10942]

- The minimum supported version of ``pyerfa`` is now 1.7.3. [#11637]


Version 4.2.1 (2021-04-01)
==========================

Bug Fixes
---------

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Fixed an issue where specializations of the comoving distance calculation
  for certain cosmologies could not handle redshift arrays. [#10980]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fix bug where manual fixes to invalid header cards were not preserved when
  saving a FITS file. [#11108]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- ``NumericArray`` converter now properly broadcasts scalar mask to array.
  [#11157]

astropy.table
^^^^^^^^^^^^^

- Fix bug when initializing a ``Table`` subclass that uses ``TableAttribute``'s.
  If the data were an instance of the table then attributes provided in the
  table initialization call could be ignored. [#11217]

astropy.time
^^^^^^^^^^^^

- Change epoch of ``TimeUnixTAI`` (``"unix_tai"``) from ``1970-01-01T00:00:00 UTC``
  to ``1970-01-01T00:00:00 TAI`` to match the intended and documented behaviour.
  This essentially changes the resulting times by 8.000082 seconds, the initial
  offset between TAI and UTC. [#11249]

astropy.units
^^^^^^^^^^^^^

- Fixed a bug with the ``quantity_input`` decorator where allowing
  dimensionless inputs for an argument inadvertently disabled any checking of
  compatible units for that argument. [#11283]

astropy.utils
^^^^^^^^^^^^^

- Fix a bug so that ``np.shape``, ``np.ndim`` and ``np.size`` again work on
  classes that use ``ShapedLikeNDArray``, like representations, frames,
  sky coordinates, and times. [#11133]

astropy.wcs
^^^^^^^^^^^

- Fix error when a user defined ``proj_point`` parameter is passed to ``fit_wcs_from_points``. [#11139]


Other Changes and Additions
---------------------------


- Change epoch of ``TimeUnixTAI`` (``"unix_tai"``) from ``1970-01-01T00:00:00 UTC``
  to ``1970-01-01T00:00:00 TAI`` to match the intended and documented behaviour.
  This essentially changes the resulting times by 8.000082 seconds, the initial
  offset between TAI and UTC. [#11249]


Version 4.2 (2020-11-24)
========================

New Features
------------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Methods ``convolve`` and ``convolve_fft`` both now return Quantity arrays
  if user input is given in one. [#10822]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Numpy functions that broadcast, change shape, or index (like
  ``np.broadcast_to``, ``np.rot90``, or ``np.roll``) now work on
  coordinates, frames, and representations. [#10337]

- Add a new science state ``astropy.coordinates.erfa_astrom.erfa_astrom`` and
  two classes ``ErfaAstrom``, ``ErfaAstromInterpolator`` as wrappers to
  the ``pyerfa`` astrometric functions used in the coordinate transforms.
  Using ``ErfaAstromInterpolator``, which interpolates astrometric properties for
  ``SkyCoord`` instances with arrays of obstime, can dramatically speed up
  coordinate transformations while keeping microarcsecond resolution.
  Depending on needed precision and the obstime array in question, speed ups
  reach factors of 10x to >100x. [#10647]

- ``galactocentric_frame_defaults`` can now also be used as a registry, with
  user-defined parameter values and metadata. [#10624]

- Method ``.realize_frame`` from coordinate frames now accepts ``**kwargs``,
  including ``representation_type``. [#10727]

- Avoid an unnecessary call to ``erfa.epv00`` in transformations between
  ``CIRS`` and ``ICRS``, improving performance by 50 %. [#10814]

- A new equatorial coordinate frame, with RA and Dec measured w.r.t to the True
  Equator and Equinox (TETE). This frame is commonly known as "apparent place"
  and is the correct frame for coordinates returned from JPL Horizons. [#10867]

- Added a context manager ``impose_finite_difference_dt`` to the
  ``TransformGraph`` class to override the finite-difference time step
  attribute (``finite_difference_dt``) for all transformations in the graph
  with that attribute. [#10341]

- Improve performance of ``SpectralCoord`` by refactoring internal
  implementation. [#10398]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- The final version of the Planck 2018 cosmological parameters are included
  as the ``Planck18`` object, which is now the default cosmology.  The
  parameters are identical to those of the ``Planck18_arXiv_v2`` object,
  which is now deprecated and will be removed in a future release. [#10915]

astropy.modeling
^^^^^^^^^^^^^^^^

- Added NFW profile and tests to modeling package [#10505]

- Added missing logic for evaluate to compound models [#10002]

- Stop iteration in ``FittingWithOutlierRemoval`` before reaching ``niter`` if
  the masked points are no longer changing. [#10642]

- Keep a (shallow) copy of ``fit_info`` from the last iteration of the wrapped
  fitter in ``FittingWithOutlierRemoval`` and also record the actual number of
  iterations performed in it. [#10642]

- Added attributes for fitting uncertainties (covariance matrix, standard
  deviations) to models. Parameter covariance matrix can be accessed via
  ``model.cov_matrix``, standard deviations by ``model.stds`` or individually
  for each parameter by ``parameter.std``. Currently implemented for
  ``LinearLSQFitter`` and ``LevMarLSQFitter``. [#10552]

- N-dimensional least-squares statistic and specific 1,2,3-D methods [#10670]

astropy.stats
^^^^^^^^^^^^^

- Added ``circstd`` function to obtain a circular standard deviation. [#10690]

astropy.table
^^^^^^^^^^^^^

- Allow initializing a ``Table`` using a list of ``names`` in conjunction with
  a ``dtype`` from a numpy structured array. The list of ``names`` overrides the
  names specified in the ``dtype``. [#10419]

astropy.time
^^^^^^^^^^^^

- Add new ``isclose()`` method to ``Time`` and ``TimeDelta`` classes to allow
  comparison of time objects to within a specified tolerance. [#10646]

- Improve initialization time by a factor of four when creating a scalar ``Time``
  object in a format like ``unix`` or ``cxcsec`` (time delta from a reference
  epoch time). [#10406]

- Improve initialization time by a factor of ~25 or more for large arrays of
  string times in ISO, ISOT or year day-of-year formats. This is done with a new
  C-based time parser that can be adapted for other fixed-format custom time
  formats. [#10360]

- Numpy functions that broadcast, change shape, or index (like
  ``np.broadcast_to``, ``np.rot90``, or ``np.roll``) now work on times.
  [#10337, #10502]

astropy.timeseries
^^^^^^^^^^^^^^^^^^

- Improve memory and speed performance when iterating over the entire time
  column of a ``TimeSeries`` object. Previously this involved O(N^2) operations
  and memory. [#10889]

astropy.units
^^^^^^^^^^^^^

- ``Quantity.to`` has gained a ``copy`` option to allow copies to be avoided
  when the units do not change. [#10517]

- Added the ``spat`` unit of solid angle that represents the full sphere.
  [#10726]

astropy.utils
^^^^^^^^^^^^^

- ``ShapedLikeNDArray`` has gained the capability to use numpy functions
  that broadcast, change shape, or index. [#10337]

- ``get_free_space_in_dir`` now takes a new ``unit`` keyword and
  ``check_free_space_in_dir`` takes ``size`` defined as ``Quantity``. [#10627]

- New ``astropy.utils.data.conf.allow_internet`` configuration item to
  control downloading data from the Internet. Setting ``allow_internet=False``
  is the same as ``remote_timeout=0``. Using ``remote_timeout=0`` to control
  internet access will stop working in a future release. [#10632]

- New ``is_url`` function so downstream packages do not have to secretly use
  the hidden ``_is_url`` anymore. [#10684]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Added the ``Quadrangle`` patch for ``WCSAxes`` for a latitude-longitude
  quadrangle.  Unlike ``matplotlib.patches.Rectangle``, the edges of this
  patch will be rendered as curved lines if appropriate for the WCS
  transformation. [#10862]

- The position of tick labels are now only calculated when needed. If any text
  parameters are changed (color, font weight, size etc.) that don't effect the
  tick label position, the positions are not recomputed, improving performance.
  [#10806]

astropy.wcs
^^^^^^^^^^^

- ``WCS.to_header()`` now appends comments to SIP coefficients. [#10480]

- A new property ``dropped_world_dimensions`` has been added to
  ``SlicedLowLevelWCS`` to record information about any world axes removed by
  slicing a WCS. [#10195]

- New ``WCS.proj_plane_pixel_scales()`` and ``WCS.proj_plane_pixel_area()``
  methods to return pixel scales and area, respectively, as Quantity. [#10872]


API Changes
-----------

astropy.config
^^^^^^^^^^^^^^

- ``set_temp_config`` now preserves the existing cache rather than deleting
  it and relying on reloading it from the previous config file. This ensures
  that any programmatically made changes are preserved as well. [#10474]

- Configuration path detection logic has changed: Now, it looks for ``~`` first
  before falling back to older logic. In addition, ``HOMESHARE`` is no longer
  used in Windows. [#10705]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- The passing of frame classes (as opposed to frame instances) to the
  ``transform_to()`` methods of low-level coordinate-frame classes has been
  deprecated.  Frame classes can still be passed to the ``transform_to()``
  method of the high-level ``SkyCoord`` class, and using ``SkyCoord`` is
  recommended for all typical use cases of transforming coordinates. [#10475]

astropy.stats
^^^^^^^^^^^^^

- Added a ``grow`` parameter to ``SigmaClip``, ``sigma_clip`` and
  ``sigma_clipped_stats``, to allow expanding the masking of each deviant
  value to its neighbours within a specified radius. [#10613]

- Passing float ``n`` to ``poisson_conf_interval`` when using
  ``interval='kraft-burrows-nousek'`` now raises ``TypeError`` as its value
  must be an integer. [#10838]

astropy.table
^^^^^^^^^^^^^

- Change ``Table.columns.keys()`` and ``Table.columns.values()`` to both return
  generators instead of a list. This matches the behavior for Python ``dict``
  objects. [#10543]

- Removed the ``FastBST`` and ``FastRBT`` indexing engines because they depend
  on the ``bintrees`` package, which is no longer maintained and is deprecated.
  Instead, use the ``SCEngine`` indexing engine, which is similar in
  performance and relies on the ``sortedcontainers`` package. [#10622]

- When slicing a mixin column in a table that had indices, the indices are no
  longer copied since they generally are not useful, having the wrong shape.
  With this, the behaviour becomes the same as that for a regular ``Column``.
  (Note that this does not affect slicing of a table; sliced columns in those
  will continue to carry a sliced version of any indices). [#10890]

- Change behavior so that when getting a single item out of a mixin column such
  as ``Time``, ``TimeDelta``, ``SkyCoord`` or ``Quantity``, the ``info``
  attribute is no longer copied. This improves performance, especially when the
  object is an indexed column in a ``Table``. [#10889]

- Raise a TypeError when a scalar column is added to an unsized table. [#10476]

- The order of columns when creating a table from a ``list`` of ``dict`` may be
  changed. Previously, the order was alphabetical because the ``dict`` keys
  were assumed to be in random order. Since Python 3.7, the keys are always in
  order of insertion, so ``Table`` now uses the order of keys in the first row
  to set the column order. To alphabetize the columns to match the previous
  behavior, use ``t = t[sorted(t.colnames)]``. [#10900]

astropy.time
^^^^^^^^^^^^

- Refactor ``Time`` and ``TimeDelta`` classes to inherit from a common
  ``TimeBase`` class. The ``TimeDelta`` class no longer inherits from ``Time``.
  A number of methods that only apply to ``Time`` (e.g. ``light_travel_time``)
  are no longer available in the ``TimeDelta`` class. [#10656]

astropy.units
^^^^^^^^^^^^^

- The ``bar`` unit is no longer wrongly considered an SI unit, meaning that
  SI decompositions like ``(u.kg*u.s**-2* u.sr**-1 * u.nm**-1).si`` will
  no longer include it. [#10586]

astropy.utils
^^^^^^^^^^^^^

- Shape-related items from ``astropy.utils.misc`` -- ``ShapedLikeNDArray``,
  ``check_broadcast``, ``unbroadcast``, and ``IncompatibleShapeError`` --
  have been moved to their own module, ``astropy.utils.shapes``. They remain
  importable from ``astropy.utils``. [#10337]

- ``check_hashes`` keyword in ``check_download_cache`` is deprecated and will
  be removed in a future release. [#10628]

- ``hexdigest`` keyword in ``import_file_to_cache`` is deprecated and will
  be removed in a future release. [#10628]


Bug Fixes
---------

astropy.config
^^^^^^^^^^^^^^

- Fix a few issues with ``generate_config`` when used with other packages.
  [#10893]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Fixed a bug in the coordinate-frame attribute ``CoordinateAttribute`` where
  the internal transformation could behave differently depending on whether
  the input was a low-level coordinate frame or a high-level ``SkyCoord``.
  ``CoordinateAttribute`` now always performs a ``SkyCoord``-style internal
  transformation, including the by-default merging of frame attributes. [#10475]

astropy.modeling
^^^^^^^^^^^^^^^^

- Fixed an issue of ``Model.render`` when the input ``out`` datatype is not
  float64. [#10542]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Fix support for referencing WCSAxes coordinates by their world axes names.
  [#10484]

astropy.wcs
^^^^^^^^^^^

- Objective functions called by ``astropy.wcs.fit_wcs_from_points`` were
  treating longitude and latitude distances equally. Now longitude scaled
  properly. [#10759]


Other Changes and Additions
---------------------------

- Minimum version of required Python is now 3.7. [#10900]

- Minimum version of required Numpy is now 1.17. [#10664]

- Minimum version of required Scipy is now 1.1. [#10900]

- Minimum version of required PyYAML is now 3.13. [#10900]

- Minimum version of required Matplotlib is now 3.0. [#10900]

- The private ``_erfa`` module has been converted to its own package,
  ``pyerfa``, which is a required dependency for astropy, and can be imported
  with ``import erfa``.  Importing ``_erfa`` from ``astropy`` will give a
  deprecation warning.  [#10329]

- Added ``optimize=True`` flag to calls of ``yacc.yacc`` (as already done for
  ``lex.lex``) to allow running in ``python -OO`` session without raising an
  exception in ``astropy.units.format``. [#10379]

- Shortened FITS comment strings for some D2IM and CPDIS FITS keywords to
  reduce the number of FITS ``VerifyWarning`` warnings when working with WCSes
  containing lookup table distortions. [#10513]

- When importing astropy without first building the extension modules first,
  raise an error directly instead of trying to auto-build. [#10883]



Version 4.1 (2020-10-21)
========================

New Features
------------

astropy.config
^^^^^^^^^^^^^^

- Add new function ``generate_config`` to generate the configuration file and
  include it in the documentation. [#10148]

- ``ConfigNamespace.__iter__`` and ``ConfigNamespace.keys`` now yield ``ConfigItem``
  names defined within it. Similarly, ``items`` and ``values`` would yield like a
  Python dictionary would. [#10139]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Added a new ``SpectralCoord`` class that can be used to define spectral
  coordinates and transform them between different velocity frames. [#10185]

- Angle parsing now supports ``cardinal direction`` in the cases
  where angles are initialized as ``string`` instances. eg ``"17°53'27"W"``.[#9859]

- Allow in-place modification of array-valued ``Frame`` and ``SkyCoord`` objects.
  This provides limited support for updating coordinate data values from another
  coordinate object of the same class and equivalent frame attributes. [#9857]

- Added a robust equality operator for comparing ``SkyCoord``, frame, and
  representation objects. A comparison like ``sc1 == sc2`` will now return a
  boolean or boolean array where the objects are strictly equal in all relevant
  frame attributes and coordinate representation values. [#10154]

- Added the True Equator Mean Equinox (TEME) frame. [#10149]

- The ``Galactocentric`` frame will now use the "latest" parameter definitions
  by default. This currently corresponds to the values defined in v4.0, but will
  change with future releases. [#10238]

- The ``SkyCoord.from_name()`` and Sesame name resolving functionality now is
  able to cache results locally and will do so by default. [#9162]

- Allow in-place modification of array-valued ``Representation`` and ``Differential``
  objects, including of representations with attached differentials. [#10210]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Functional Units can now be processed in CDS-tables. [#9971]

- Allow reading in ASCII tables which have duplicate column names. [#9939]

- Fixed failure of ASCII ``fast_reader`` to handle ``names``, ``include_names``,
  ``exclude_names`` arguments for ``RDB`` formatted tables. Homogenised checks
  and exceptions for invalid ``names`` arguments. Improved performance when
  parsing "wide" tables with many columns. [#10306]

- Added type validation of key arguments in calls to ``io.ascii.read()`` and
  ``io.ascii.write()`` functions. [#10005]

astropy.io.misc
^^^^^^^^^^^^^^^
- Added serialization of parameter constraints fixed and bounds.  [#10082]

- Added 'functional_models.py' and 'physical_models.py' to asdf/tags/transform,
  with to allow serialization of all functional and physical models. [#10028, #10293]

- Fix ASDF serialization of circular model inverses, and remove explicit calls
  to ``asdf.yamlutil`` functions that became unnecessary in asdf 2.6.0. [#10189, #10384]

astropy.io.fits
^^^^^^^^^^^^^^^

- Added support for writing Dask arrays to disk efficiently for ``ImageHDU`` and
  ``PrimaryHDU``. [#9742]

- Add HDU name and ver to FITSDiff report where appropriate [#10197]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- New ``exceptions.conf.max_warnings`` configuration item to control the number of times a
  type of warning appears before being suppressed. [#10152]

- No longer ignore attributes whose values were specified as empty
  strings. [#10583]

astropy.modeling
^^^^^^^^^^^^^^^^
- Added Plummer1D model to ``functional_models``. [#9896]

- Added ``UnitsMapping`` model and ``Model.coerce_units`` to support units on otherwise
  unitless models. [#9936]

- Added ``domain`` and ``window`` attributes to ``repr`` and ``str``. Fixed bug with
  ``_format_repr`` in core.py. [#9941]

- Polynomial attributes ``domain`` and ``window`` are now tuples of size 2 and are
  validated. `repr` and `print` show only their non-default values. [#10145]

- Added ``replace_submodel()`` method to ``CompoundModel`` to modify an
  existing instance. [#10176]

- Delay construction of ``CompoundModel`` inverse until property is accessed,
  to support ASDF deserialization of circular inverses in component models. [#10384]

astropy.nddata
^^^^^^^^^^^^^^

- Added support in the ``bitmask`` module for using mnemonic bit flag names
  when specifying the bit flags to be used or ignored when converting a bit
  field to a boolean. [#10095, #10208]

- Added ``reshape_as_blocks`` function to reshape a data array into
  blocks, which is useful to efficiently apply functions on block
  subsets of the data instead of using loops.  The reshaped array is a
  view of the input data array. [#10214]

- Added a ``cache`` keyword option to allow caching for ``CCDData.read`` if
  filename is a URL. [#10265]

astropy.table
^^^^^^^^^^^^^

- Added ability to specify a custom matching function for table joins.  In
  particular this makes it possible to do cross-match table joins on ``SkyCoord``,
  ``Quantity``, or standard columns, where column entries within a specified
  distance are considered to be matched. [#10169]

- Added ``units`` and ``descriptions`` keyword arguments to the Table object
  initialization and ``Table.read()`` methods.  This allows directly setting
  the ``unit`` and ``description`` for the table columns at the time of
  creating or reading the table. [#9671]

- Make table ``Row`` work as mappings, by adding ``.keys()`` and ``.values()``
  methods. With this ``**row`` becomes possible, as does, more simply, turning
  a ``Row`` into a dictionary with ``dict(row)``. [#9712]

- Added two new ``Table`` methods ``.items()`` and ``.values()``, which return
  respectively ``tbl.columns.items()`` (iterator over name, column tuples)  and
  ``tbl.columns.values()`` (list of columns) for a ``Table`` object ``tbl``. [#9780]

- Added new ``Table`` method ``.round()``, which rounds numeric columns to the
  specified number of decimals. [#9862]

- Updated ``to_pandas()`` and ``from_pandas()`` to use and support Pandas
  nullable integer data type for masked integer data. [#9541]

- The HDF5 writer, ``write_table_hdf5()``, now allows passing through
  additional keyword arguments to the ``h5py.Group.create_dataset()``. [#9602]

- Added capability to add custom table attributes to a ``Table`` subclass.
  These attributes are persistent and can be set during table creation. [#10097]

- Added support for ``SkyCoord`` mixin columns in ``dstack``, ``vstack`` and
  ``insert_row`` functions. [#9857]

- Added support for coordinate ``Representation`` and ``Differential`` mixin
  columns. [#10210]

astropy.time
^^^^^^^^^^^^

- Added a new time format ``unix_tai`` which is essentially Unix time but with
  leap seconds included.  More precisely, this is the number of seconds since
  ``1970-01-01 00:00:08 TAI`` and corresponds to the ``CLOCK_TAI`` clock
  available on some linux platforms. [#10081]

astropy.units
^^^^^^^^^^^^^

- Added ``torr`` pressure unit. [#9787]

- Added the ``equal_nan`` keyword argument to ``isclose`` and ``allclose``, and
  updated the docstrings. [#9849]

- Added ``Rankine`` temperature unit. [#9916]

- Added integrated flux unit conversion to ``spectral_density`` equivalency.
  [#10015]

- Changed ``pixel_scale`` equivalency to allow scales defined in any unit.
  [#10123]

- The ``quantity_input`` decorator now optionally allows passing through
  numeric values or numpy arrays with numeric dtypes to arguments where
  ``dimensionless_unscaled`` is an allowed unit. [#10232]

astropy.utils
^^^^^^^^^^^^^

- Added a new ``MetaAttribute`` class to support easily adding custom attributes
  to a subclass of classes like ``Table`` or ``NDData`` that have a ``meta``
  attribute. [#10097]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Added ``invalid`` keyword to ``SqrtStretch``, ``LogStretch``,
  ``PowerStretch``, and ``ImageNormalize`` classes and the
  ``simple_norm`` function.  This keyword is used to replace generated
  NaN values. [#10182]

- Fixed an issue where ticks were sometimes not drawn at the edges of a spherical
  projection on a WCSAxes. [#10442]

astropy.wcs
^^^^^^^^^^^

- WCS objects with a spectral axis will now return ``SpectralCoord``
  objects when calling ``pixel_to_world`` instead of ``Quantity``,
  and can now take either ``Quantity`` or ``SpectralCoord`` as input
  to ``pixel_to_world``. [#10185]

- Implemented support for the ``-TAB`` algorithm (WCS Paper III). [#9641]

- Added an ``_as_mpl_axes`` method to the ``HightLevelWCSWrapper`` class. [#10138]

- Add .upper() to ctype or ctype names to wcsapi/fitwcs.py to mitigate bugs from
  unintended lower/upper case issues [#10557]

API Changes
-----------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- The equality operator for comparing ``SkyCoord``, frame, and representation
  objects was changed. A comparison like ``sc1 == sc2`` was previously
  equivalent to ``sc1 is sc2``. It will now return a boolean or boolean array
  where the objects are strictly equal in all relevant frame attributes and
  coordinate representation values. If the objects have different frame
  attributes or representation types then an exception will be raised. [#10154]

- ```SkyCoord.radial_velocity_correction``` now allows you to pass an ```obstime``` directly
  when the ```SkyCoord``` also has an ```obstime``` set. In this situation, the position of the
  ```SkyCoord``` has space motion applied to correct to the passed ```obstime```. This allows
  mm/s radial velocity precision for objects with large space motion. [#10094]

- For consistency with other astropy classes, coordinate ``Representations``
  and ``Differentials`` can now be initialized with an instance of their own class
  if that instance is passed in as the first argument. [#10210]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Changed the behavior when reading a table where both the ``names`` argument
  is provided (to specify the output column names) and the ``converters``
  argument is provided (to specify column conversion functions). Previously the
  ``converters`` dict names referred to the *input* table column names, but now
  they refer to the *output* table column names. [#9739]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- For FIELDs with datatype="char", store the values as strings instead
  of bytes. [#9505]

astropy.table
^^^^^^^^^^^^^

- ``Table.from_pandas`` now supports a ``units`` dictionary as argument to pass units
  for columns in the ``DataFrame``. [#9472]

astropy.time
^^^^^^^^^^^^

- Require that ``in_subfmt`` and ``out_subfmt`` properties of a ``Time`` object
  have allowed values at the time of being set, either when creating the object
  or when setting those properties on an existing ``Time`` instance.  Previously
  the validation of those properties was not strictly enforced. [#9868]

astropy.utils
^^^^^^^^^^^^^

- Changed the exception raised by ``get_readable_fileobj`` on missing
  compression modules (for ``bz2`` or ``lzma``/``xz`` support) to
  ``ModuleNotFoundError``, consistent with ``io.fits`` file handlers. [#9761]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Deprecated the ``imshow_only_kwargs`` keyword in ``imshow_norm``.
  [#9915]

- Non-finite input values are now automatically excluded in
  ``HistEqStretch`` and ``InvertedHistEqStretch``. [#10177]

- The ``PowerDistStretch`` and ``InvertedPowerDistStretch`` ``a``
  value is restricted to be ``a >= 0`` in addition to ``a != 1``.
  [#10177]

- The ``PowerStretch``, ``LogStretch``, and ``InvertedLogStretch``
  ``a`` value is restricted to be ``a > 0``. [#10177]

- The ``AsinhStretch`` and ``SinhStretch`` ``a`` value is restricted
  to be ``0 < a <= 1``. [#10177]

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Fix a bug where for light deflection by the Sun it was always assumed that the
  source was at infinite distance, which in the (rare and) absolute worst-case
  scenario could lead to errors up to 3 arcsec. [#10666]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- For FIELDs with datatype="char", store the values as strings instead
  of bytes. [#9505]

astropy.table
^^^^^^^^^^^^^

- Fix a bug that prevented ``Time`` columns from being used to sort a table.
  [#10824]

astropy.wcs
^^^^^^^^^^^

- WCS objects with a spectral axis will now return ``SpectralCoord``
  objects when calling ``pixel_to_world`` instead of ``Quantity``
  (note that ``SpectralCoord`` is a sub-class of ``Quantity``). [#10185]

- Add .upper() to ctype or ctype names to wcsapi/fitwcs.py to mitigate bugs from
  unintended lower/upper case issues [#10557]

- Added bounds to ``fit_wcs_from_points`` to ensure CRPIX is on
  input image. [#10346]


Other Changes and Additions
---------------------------

- The way in which users can specify whether to build astropy against
  existing installations of C libraries rather than the bundled one
  has changed, and should now be done via environment variables rather
  than setup.py flags (e.g. --use-system-erfa). The available variables
  are ``ASTROPY_USE_SYSTEM_CFITSIO``, ``ASTROPY_USE_SYSTEM_ERFA``,
  ``ASTROPY_USE_SYSTEM_EXPAT``, ``ASTROPY_USE_SYSTEM_WCSLIB``, and
  ``ASTROPY_USE_SYSTEM_ALL``. These should be set to ``1`` to build
  against the system libraries. [#9730]

- The infrastructure of the package has been updated in line with the
  APE 17 roadmap (https://github.com/astropy/astropy-APEs/blob/master/APE17.rst).
  The main changes are that the ``python setup.py test`` and
  ``python setup.py build_docs`` commands will no longer work. The easiest
  way to replicate these commands is to install the tox
  (https://tox.readthedocs.io) package and run ``tox -e test`` and
  ``tox -e build_docs``. It is also possible to run pytest and sphinx
  directly. Other significant changes include switching to setuptools_scm to
  manage the version number, and adding a ``pyproject.toml`` to opt in to
  isolated builds as described in PEP 517/518. [#9726]

- Bundled ``expat`` is updated to version 2.2.9. [#10038]

- Increase minimum asdf version to 2.6.0. [#10189]

- The bundled version of PLY was updated to 3.11. [#10258]

- Removed dependency on scikit-image. [#10214]

Version 4.0.5 (2021-03-26)
==========================

Bug Fixes
---------

astropy.io.fits
^^^^^^^^^^^^^^^

- Fix bug where manual fixes to invalid header cards were not preserved when
  saving a FITS file. [#11108]

- Fix parsing of RVKC header card patterns that were not recognised
  where multiple spaces were separating field-specifier and value like
  "DP1.AXIS.1:   1". [#11301]

- Fix misleading missing END card error when extra data are found at the end
  of the file. [#11285]

- Fix incorrect wrapping of long card values as CONTINUE cards when some
  words in the value are longer than a single card. [#11304]

astropy.io.misc
^^^^^^^^^^^^^^^

- Fixed problem when writing serialized metadata to HDF5 using h5py >= 3.0.
  With the newer h5py this was writing the metadata table as a variable-length
  string array instead of the previous fixed-length bytes array. Fixed astropy
  to force using a fixed-length bytes array. [#11359]

astropy.modeling
^^^^^^^^^^^^^^^^

- Change ``Voigt1D`` function to use Humlicek's approximation to avoid serious
  inaccuracies + option to use (compiled) ``scipy.special.wofz`` error function
  for yet more accurate results. [#11177]

astropy.table
^^^^^^^^^^^^^

- Fixed bug when initializing a ``Table`` with a column as list of ``Quantity``,
  for example ``Table({'x': [1*u.m, 2*u.m]})``. Previously this resulted in an
  ``object`` dtype with no column ``unit`` set, but now gives a float array with
  the correct unit. [#11329]

- Fixed byteorder conversion in ``to_pandas()``, which had incorrectly
  triggered swapping when native endianness was stored with explicit
  ``dtype`` code ``'<'`` (or ``'>'``) instead of ``'='``. [#11288, #11294]

- Fixed a compatibility issue with numpy 1.21. Initializing a Table with a
  column like ``['str', np.ma.masked]`` was failing in tests due to a change in
  numpy. [#11364]

- Fixed bug when validating the inputs to ``table.hstack``, ``table.vstack``,
  and ``table.dstack``. Previously, mistakenly calling ``table.hstack(t1, t2)``
  (instead of ``table.hstack([t1, t2]))`` would return ``t1`` instead of raising
  an exception. [#11336]

- Fixed byteorder conversion in ``to_pandas()``, which had incorrectly
  triggered swapping when native endianness was stored with explicit
  ``dtype`` code ``'<'`` (or ``'>'``) instead of ``'='``. [#11288]

astropy.time
^^^^^^^^^^^^

- Fix leap second update when using a non english locale. [#11062]

- Fix default assumed location to be the geocenter when transforming times
  to and from solar-system barycenter scales. [#11134]

- Fix inability to write masked times with ``formatted_value``. [#11195]

astropy.units
^^^^^^^^^^^^^

- Ensure ``keepdims`` works for taking ``mean``, ``std``, and ``var`` of
  ``Quantity``. [#11198]

- For ``Quantity.to_string()``, ensure that the precision argument is also
  used when the format is not latex. [#11145]

astropy.wcs
^^^^^^^^^^^

- Allow "un-setting" of auxiliary WCS parameters in the ``aux`` attribute of
  ``Wcsprm``. [#11166]





Version 4.0.4 (2020-11-24)
==========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- The ``norm()`` method for ``RadialDifferential`` no longer requires ``base``
  to be specified.  The ``norm()`` method for other non-Cartesian differential
  classes now gives a clearer error message if ``base`` is not specified. [#10969]

- The transformations between ``ICRS`` and any of the heliocentric ecliptic
  frames (``HeliocentricMeanEcliptic``, ``HeliocentricTrueEcliptic``, and
  ``HeliocentricEclipticIAU76``) now correctly account for the small motion of
  the Sun when transforming a coordinate with velocity information. [#10970]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Partially fixed a performance issue when reading in parallel mode. Parallel
  reading currently has substantially worse performance than the default serial
  reading, so we now ignore the parallel option and fall back to serial reading.
  [#10880]

- Fixed a bug where "" (blank string) as input data for a boolean type column
  was causing an exception instead of indicating a masked value. As a
  consequence of the fix, the values "0" and "1" are now also allowed as valid
  inputs for boolean type columns. These new allowed values apply for both ECSV
  and for basic character-delimited data files ('basic' format with appropriate
  ``converters`` specified). [#10995]

astropy.modeling
^^^^^^^^^^^^^^^^

- Fixed use of weights with ``LinearLSQFitter``. [#10687]

astropy.stats
^^^^^^^^^^^^^

- Fixed an issue in biweight stats when MAD=0 to give the same output
  with and without an input ``axis``. [#10912]

astropy.time
^^^^^^^^^^^^

- Fix a problem with the ``plot_date`` format for matplotlib >= 3.3 caused by
  a change in the matplotlib plot date default reference epoch in that release.
  [#10876]

- Improve initialization time by a factor of four when creating a scalar ``Time``
  object in a format like ``unix`` or ``cxcsec`` (time delta from a reference
  epoch time). [#10406]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Fixed the calculation of the tight bounding box of a ``WCSAxes``. This should
  also significantly improve the application of ``tight_layout()`` to figures
  containing ``WCSAxes``. [#10797]


Version 4.0.3 (2020-10-14)
==========================

Bug Fixes
---------

astropy.table
^^^^^^^^^^^^^

- Fixed a small bug where initializing an empty ``Column`` with a structured dtype
  with a length and a shape failed to give the requested dtype. [#10819]

Other Changes and Additions
---------------------------

- Fixed installation of the source distribution with pip<19. [#10837, #10852]


Version 4.0.2 (2020-10-10)
==========================

New Features
------------

astropy.utils
^^^^^^^^^^^^^

- ``astropy.utils.data.download_file`` now supports FTPS/FTP over TLS. [#9964]

- ``astropy.utils.data`` now uses a lock-free mechanism for caching. This new
  mechanism uses a new cache layout and so ignores caches created using earlier
  mechanisms (which were causing lockups on clusters). The two cache formats can
  coexist but do not share any files. [#10437, #10683]

- ``astropy.utils.data`` now ignores the config item
  ``astropy.utils.data.conf.download_cache_lock_attempts`` since no locking is
  done. [#10437, #10683]

- ``astropy.utils.data.download_file`` and related functions now interpret the
  parameter or config file setting ``timeout=0`` to mean they should make no
  attempt to download files. [#10437, #10683]

- ``astropy.utils.import_file_to_cache`` now accepts a keyword-only argument
  ``replace``, defaulting to True, to determine whether it should replace existing
  files in the cache, in a way as close to atomic as possible. [#10437, #10683]

- ``astropy.utils.data.download_file`` and related functions now treat
  ``http://example.com`` and ``http://example.com/`` as equivalent. [#10631]

astropy.wcs
^^^^^^^^^^^

- The new auxiliary WCS parameters added in WCSLIB 7.1 are now exposed as
  the ``aux`` attribute of ``Wcsprm``. [#10333]

- Updated bundled version of ``WCSLIB`` to v7.3. [#10433]


Bug fixes
---------

astropy.config
^^^^^^^^^^^^^^

- Added an extra fallback to ``os.expanduser('~')`` when trying to find the
  user home directory. [#10570]

astropy.constants
^^^^^^^^^^^^^^^^^

- Corrected definition of parsec to 648 000 / pi AU following IAU 2015 B2 [#10569]

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Fixed a bug where a float-typed integers in the argument ``x_range`` of
  ``astropy.convolution.utils.discretize_oversample_1D`` (and the 2D version as
  well) fails because it uses ``numpy.linspace``, which requires an ``int``.
  [#10696]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Ensure that for size-1 array ``SkyCoord`` and coordinate frames
  the attributes also properly become scalars when indexed with 0.
  [#10113]

- Fixed a bug where ``SkyCoord.separation()`` and ``SkyCoord.separation_3d``
  were not accepting a frame object. [#10332]

- Ensure that the ``lon`` values in ``SkyOffsetFrame`` are wrapped correctly at
  180 degree regardless of how the underlying data is represented. [#10163]

- Fixed an error in the obliquity of the ecliptic when transforming to/from the
  ``*TrueEcliptic`` coordinate frames. The error would primarily result in an
  inaccuracy in the ecliptic latitude on the order of arcseconds. [#10129]

- Fixed an error in the computation of the location of solar system bodies where the
  Earth location of the observer was ignored during the correction for light travel
  time. [#10292]

- Ensure that coordinates with proper motion that are transformed to other
  coordinate frames still can be represented properly. [#10276]

- Improve the error message given when trying to get a cartesian representation
  for coordinates that have both proper motion and radial velocity, but no
  distance. [#10276]

- Fixed an error where ``SkyCoord.apply_space_motion`` would return incorrect
  results when no distance is set and proper motion is high. [#10296]

- Make the parsing of angles thread-safe so that ``Angle`` can be used in
  Python multithreading. [#10556]

- Fixed reporting of ``EarthLocation.info`` which previously raised an exception.
  [#10592]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fixed a bug with the C ``fast_reader`` not correctly parsing newlines when
  ``delimiter`` was also set to ``\n`` or ``\r``; ensured consistent handling
  of input strings without newline characters. [#9929]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fix integer formats of ``TFORMn=Iw`` columns in ASCII tables to correctly read
  values exceeding int32 - setting int16, int32 or int64 according to ``w``. [#9901]

- Fix unclosed memory-mapped FITS files in ``FITSDiff`` when difference found.
  [#10159]

- Fix crash when reading an invalid table file. [#10171]

- Fix duplication issue when setting a keyword ending with space. [#10482]

- Fix ResourceWarning with ``fits.writeto`` and ``pathlib.Path`` object.
  [#10599]

- Fix repr for commentary cards and strip spaces for commentary keywords.
  [#10640]

- Fix compilation of cfitsio with Xcode 12. [#10772]

- Fix handling of 1-dimensional arrays with a single element in ``BinTableHDU`` [#10768]

astropy.io.misc
^^^^^^^^^^^^^^^

- Fix id URL in ``baseframe-1.0.0`` ASDF schema. [#10223]

- Write keys to ASDF only if the value is present, to account
  for a change in behavior in asdf 2.8. [#10674]

astropy.io.registry
^^^^^^^^^^^^^^^^^^^

- Fix ``Table.(read|write).help`` when reader or writer has no docstring. [#10460]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Fixed parsing failure of VOTable with no fields. When detecting a non-empty
  table with no fields, the following warning/exception is issued:
  E25 "No FIELDs are defined; DATA section will be ignored." [#10192]

astropy.modeling
^^^^^^^^^^^^^^^^

- Fixed a problem with mapping ``input_units`` and ``return_units``
  of a ``CompoundModel`` to the units of the constituent models. [#10158]

- Removed hard-coded names of inputs and outputs. [#10174]

- Fixed a problem where slicing a ``CompoundModel`` by name will crash if
  there ``fix_inputs`` operators are present. [#10224]

- Removed a limitation of fitting of data with units with compound models
  without units when the expression involves operators other than addition
  and subtraction. [#10415]

- Fixed a problem with fitting ``Linear1D`` and ``Planar2D`` in model sets. [#10623]

- Fixed reported module name of ``math_functions`` model classes. [#10694]

- Fixed reported module name of ``tabular`` model classes. [#10709]

- Do not create new ``math_functions`` models for ufuncs that are
  only aliases (divide and mod). [#10697]

- Fix calculation of the ``Moffat2D`` derivative with respect to gamma. [#10784]

astropy.stats
^^^^^^^^^^^^^

- Fixed an API regression where ``SigmaClip.__call__`` would convert masked
  elements to ``nan`` and upcast the dtype to ``float64`` in its output
  ``MaskedArray`` when using the ``axis`` parameter along with the defaults
  ``masked=True`` and ``copy=True``. [#10610]

- Fixed an issue where fully masked ``MaskedArray`` input to
  ``sigma_clipped_stats`` gave incorrect results. [#10099]

- Fixed an issue where ``sigma_clip`` and ``SigmaClip.__call__``
  would return a masked array instead of a ``ndarray`` when
  ``masked=False`` and the input was a full-masked ``MaskedArray``.
  [#10099]

- Fixed bug with ``funcs.poisson_conf_interval`` where an integer for N
  with ``interval='kraft-burrows-nousek'`` would throw an error with
  mpmath backend. [#10427]

- Fixed bug in ``funcs.poisson_conf_interval`` with
  ``interval='kraft-burrows-nousek'`` where certain combinations of source
  and background count numbers led to ``ValueError`` due to the choice of
  starting value for numerical optimization. [#10618]

astropy.table
^^^^^^^^^^^^^

- Fixed a bug when writing a table with mixin columns to FITS, ECSV or HDF5.
  If one of the data attributes of the mixin (e.g. ``skycoord.ra``) had the
  same name as one of the table column names (``ra``), the column (``ra``)
  would be dropped when reading the table back. [#10222]

- Fixed a bug when sorting an indexed table on the indexed column after first
  sorting on another column. [#10103]

- Fixed a bug in table argsort when called with ``reverse=True`` for an
  indexed table. [#10103]

- Fixed a performance regression introduced in #9048 when initializing a table
  from Python lists. Also fixed incorrect behavior (for data types other than
  float) when those lists contain ``np.ma.masked`` elements to indicate masked
  data. [#10636]

- Avoid modifying ``.meta`` when serializing columns to FITS. [#10485]

- Avoid crash when reading a FITS table that contains mixin info and PyYAML
  is missing. [#10485]

astropy.time
^^^^^^^^^^^^

- Ensure that for size-1 array ``Time``, the location also properly becomes
  a scalar when indexed with 0. [#10113]

astropy.units
^^^^^^^^^^^^^

- Refined test_parallax to resolve difference between 2012 and 2015 definitions. [#10569]

astropy.utils
^^^^^^^^^^^^^

- The default IERS server has been updated to use the FTPS server hosted by
  CDDIS. [#9964]

- Fixed memory allocation on 64-bit systems within ``xml.iterparse`` [#10076]

- Fix case where ``None`` could be used in a numerical computation. [#10126]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Fixed a bug where the ``ImageNormalize`` ``clip`` keyword was
  ignored when used with calling the object on data. [#10098]

- Fixed a bug where ``axes.xlabel``/``axes.ylabel`` where not correctly set
  nor returned on an ``EllipticalFrame`` class ``WCSAxes`` plot. [#10446]

astropy.wcs
^^^^^^^^^^^

- Handled WCS 360 -> 0 deg crossover in ``fit_wcs_from_points`` [#10155]

- Do not issue ``DATREF`` warning when ``MJDREF`` has default value. [#10440]

- Fixed a bug due to which ``naxis`` argument was ignored if ``header``
  was supplied during the initialization of a WCS object. [#10532]

Other Changes and Additions
---------------------------

- Improved the speed of sorting a large ``Table`` on a single column by a factor
  of around 5. [#10103]

- Ensure that astropy can be used inside Application bundles built with
  pyinstaller. [#8795]

- Updated the bundled CFITSIO library to 3.49. See
  ``cextern/cfitsio/docs/changes.txt`` for additional information.
  [#10256, #10665]

- ``extract_array`` raises a ``ValueError`` if the data type of the
  input array is inconsistent with the ``fill_value``. [#10602]


Version 4.0.1 (2020-03-27)
==========================

Bug fixes
---------

astropy.config
^^^^^^^^^^^^^^

- Fixed a bug where importing a development version of a package that uses
  ``astropy`` configuration system can result in a
  ``~/.astropy/config/package..cfg`` file. [#9975]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Fixed a bug where a vestigal trace of a frame class could persist in the
  transformation graph even after the removal of all transformations involving
  that frame class. [#9815]

- Fixed a bug with ``TransformGraph.remove_transform()`` when the "from" and
  "to" frame classes are not explicitly specified. [#9815]

- Read-only longitudes can now be passed in to ``EarthLocation`` even if
  they include angles outside of the range of -180 to 180 degrees. [#9900]

- ```SkyCoord.radial_velocity_correction``` no longer raises an Exception
  when space motion information is present on the SkyCoord. [#9980]

astropy.io
^^^^^^^^^^

- Fixed a bug that prevented the unified I/O infrastructure from working with
  datasets that are represented by directories rather than files. [#9866]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fixed a bug in the ``fast_reader`` C parsers incorrectly returning entries
  of isolated positive/negative signs as ``float`` instead of ``str``. [#9918]

- Fixed a segmentation fault in the ``fast_reader`` C parsers when parsing an
  invalid file with ``guess=True`` and the file contains inconsistent column
  numbers in combination with a quoted field; e.g., ``"1  2\n 3  4 '5'"``.
  [#9923]

- Magnitude, decibel, and dex can now be stored in ``ecsv`` files. [#9933]

astropy.io.misc
^^^^^^^^^^^^^^^

- Magnitude, decibel, and dex can now be stored in ``hdf5`` files. [#9933]

- Fixed serialization of polynomial models to include non default values of
  domain and window values. [#9956, #9961]

- Fixed a bug which affected overwriting tables within ``hdf5`` files.
  Overwriting an existing path with associated column meta data now also
  overwrites the meta data associated with the table. [#9950]

- Fixed serialization of Time objects with location under time-1.0.0
  ASDF schema. [#9983]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fix regression with ``GroupsHDU`` which needs to modify the header to handle
  invalid headers, and fix accessing ``.data`` for empty HDU. [#9711, #9934]

- Fix ``fitsdiff`` when its arguments are directories that contain other
  directories. [#9711]

- Fix writing noncontiguous data to a compressed HDU. [#9958]

- Added verification of ``disp`` (``TDISP``) keyword to ``fits.Column`` and
  extended tests for ``TFORM`` and ``TDISP`` validation. [#9978]

- Fix checksum verification to process all HDUs instead of only the first one
  because of the lazy loading feature. [#10012]

- Allow passing ``output_verify`` to ``.close`` when using the context manager.
  [#10030]

- Prevent instantiation of ``PrimaryHDU`` and ``ImageHDU`` with a scalar.
  [#10041]

- Fix column access by attribute with FITS_rec: columns with scaling or columns
  from ASCII tables where not properly converted when accessed by attribute
  name. [#10069]

astropy.io.misc
^^^^^^^^^^^^^^^

- Magnitude, decibel, and dex can now be stored in ``hdf5`` files. [#9933]

- Fixed serialization of polynomial models to include non default values of
  domain and window values. [#9956, #9961]

- Fixed a bug which affected overwriting tables within ``hdf5`` files.
  Overwriting an existing path with associated column meta data now also
  overwrites the meta data associated with the table. [#9950]

- Fixed serialization of Time objects with location under time-1.0.0
  ASDF schema. [#9983]

astropy.modeling
^^^^^^^^^^^^^^^^

- Fixed a bug in setting default values of parameters of orthonormal
  polynomials when constructing a model set. [#9987]

astropy.table
^^^^^^^^^^^^^

- Fixed bug in ``Table.reverse`` for tables that contain non-mutable mixin columns
  (like ``SkyCoord``) for which in-place item update is not allowed. [#9839]

- Tables containing Magnitude, decibel, and dex columns can now be saved to
  ``ecsv`` files. [#9933]

- Fixed bug where adding or inserting a row fails on a table with an index
  defined on a column that is not the first one. [#10027]

- Ensured that ``table.show_in_browser`` also worked for mixin columns like
  ``Time`` and ``SkyCoord``. [#10068]

astropy.time
^^^^^^^^^^^^

- Fix inaccuracy when converting between TimeDelta and datetime.timedelta. [#9679]

- Fixed exception when changing ``format`` in the case when ``out_subfmt`` is
  defined and is incompatible with the new format. [#9812]

- Fixed exceptions in ``Time.to_value()``: when supplying any ``subfmt`` argument
  for string-based formats like 'iso', and for ``subfmt='long'`` for the formats
  'byear', 'jyear', and 'decimalyear'. [#9812]

- Fixed bug where the location attribute was lost when creating a new ``Time``
  object from an existing ``Time`` or list of ``Time`` objects. [#9969]

- Fixed a bug where an exception occurred when creating a ``Time`` object
  if the ``val1`` argument was a regular double and the ``val2`` argument
  was a ``longdouble``. [#10034]

astropy.timeseries
^^^^^^^^^^^^^^^^^^

- Fixed issue with reference time for the ``transit_time`` parameter returned by
  the ``BoxLeastSquares`` periodogram. Now, the ``transit_time`` will be within
  the range of the input data and arbitrary time offsets/zero points no longer
  affect results. [#10013]

astropy.units
^^^^^^^^^^^^^

- Fix for ``quantity_input`` annotation raising an exception on iterable
  types that don't define a general ``__contains__`` for checking if ``None``
  is contained (e.g. Enum as of python3.8), by instead checking for instance of
  Sequence. [#9948]

- Fix for ``u.Quantity`` not taking into account ``ndmin`` if constructed from
  another ``u.Quantity`` instance with different but convertible unit [#10066]

astropy.utils
^^^^^^^^^^^^^

- Fixed ``deprecated_renamed_argument`` not passing in user value to
  deprecated keyword when the keyword has no new name. [#9981]

- Fixed ``deprecated_renamed_argument`` not issuing a deprecation warning when
  deprecated keyword without new name is passed in as positional argument.
  [#9985]

- Fixed detection of read-only filesystems in the caching code. [#10007]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Fixed bug from matplotlib >=3.1 where an empty Quantity array is
  sent for unit conversion as an empty list. [#9848]

- Fix bug in ``ZScaleInterval`` to return the array minimum and
  maximum when there are less then ``min_npixels`` in the input array. [#9913]

- Fix a bug in simplifying axis labels that affected non-rectangular frames.
  [#8004, #9991]


Other Changes and Additions
---------------------------

- Increase minimum asdf version to 2.5.2. [#9996, #9819]

- Updated bundled version of ``WCSLIB`` to v7.2. [#10021]



Version 4.0 (2019-12-16)
========================

New Features
------------

astropy.config
^^^^^^^^^^^^^^

- The config and cache directories and the name of the config file are now
  customizable. This allows affiliated packages to put their configuration
  files in locations other than ``CONFIG_DIR/.astropy/``. [#8237]

astropy.constants
^^^^^^^^^^^^^^^^^

- The version of constants can be specified via ScienceState in a way
  that ``constants`` and ``units`` will be consistent. [#8517]

- Default constants now use CODATA 2018 and IAU 2015 definitions. [#8761]

- Constants can be pickled and unpickled. [#9377]

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Fixed a bug [#9168] where having a kernel defined using unitless astropy
  quantity objects would result in a crash [#9300]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Changed ``coordinates.solar_system_ephemeris`` to also accept local files
  as input. The ephemeris can now be selected by either keyword (e.g. 'jpl',
  'de430'), URL or file path. [#8767]

- Added a ``cylindrical`` property to ``SkyCoord`` for shorthand access to a
  ``CylindricalRepresentation`` of the coordinate, as is already available
  for other common representations. [#8857]

- The default parameters for the ``Galactocentric`` frame are now controlled by
  a ``ScienceState`` subclass, ``galactocentric_frame_defaults``. New
  parameter sets will be added to this object periodically to keep up with
  ever-improved measurements of the solar position and motion. [#9346]

- Coordinate frame classes can now have multiple aliases by assigning a list
  of aliases to the class variable ``name``.  Any of the aliases can be used
  for attribute-style access or as the target of ``tranform_to()`` calls.
  [#8834]

- Passing a NaN to ``Distance`` no longer raises a warning. [#9598]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- The pre-publication Planck 2018 cosmological parameters are included as the
  ``Planck2018_arXiv_v2`` object.  Please note that the values are preliminary,
  and when the paper is accepted a final version will be included as
  ``Planck18``. [#8111]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Removed incorrect warnings on ``Overflow`` when reading in
  ``FloatType`` 0.0 with ``use_fast_converter``; synchronised
  ``IntType`` ``Overflow`` warning messages. [#9082]

astropy.io.misc
^^^^^^^^^^^^^^^

- Eliminate deprecated compatibility mode when writing ``Table`` metadata to
  HDF5 format. [#8899]

- Add support for orthogonal polynomial models to ASDF. [#9107]

astropy.io.fits
^^^^^^^^^^^^^^^

- Changed the ``fitscheck`` and ``fitsdiff`` script to use the ``argparse``
  module instead of ``optparse``. [#9148]

- Allow writing of ``Table`` objects with ``Time`` columns that are also table
  indices to FITS files. [#8077]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Support VOTable version 1.4.  The main addition is the new element, TIMESYS,
  which allows defining of metadata for temporal coordinates much like COOSYS
  defines metadata for celestial coordinates. [#9475]

astropy.logger
^^^^^^^^^^^^^^

- Added a configuration option to specify the text encoding of the log file,
  with the default behavior being the platform-preferred encoding. [#9203]

astropy.modeling
^^^^^^^^^^^^^^^^

- Major rework of modeling internals. `See modeling documentation for details.
  <https://docs.astropy.org/en/v4.0.x/modeling/changes_for_4.html>`_ . [#8769]

- Add ``Tabular1D.inverse``. [#9083]

- ``Model.rename`` was changed to add the ability to rename ``Model.inputs``
  and ``Model.outputs``. [#9220]

- New function ``fix_inputs`` to generate new models from others by fixing
  specific inputs variable values to constants. [#9135]

- ``inputs`` and ``outputs`` are now model instance attributes, and ``n_inputs``
  and ``n_outputs`` are class attributes. Backwards compatible default
  values of ``inputs`` and ``outputs`` are generated. ``Model.inputs`` and
  ``Model.outputs`` are now settable which allows renaming them on per user
  case. [#9298]

- Add a new model representing a sequence of rotations in 3D around an
  arbitrary number of axes. [#9369]

- Add many of the numpy ufunc functions as models. [#9401]

- Add ``BlackBody`` model. [#9282]

- Add ``Drude1D`` model. [#9452]

- Added analytical King model (KingProjectedAnalytic1D). [#9084]

- Added Exponential1D and Logarithmic1D models. [#9351]

astropy.nddata
^^^^^^^^^^^^^^

- Add a way for technically invalid but unambiguous units in a fits header
  to be parsed by ``CCDData``. [#9397]

- ``NDData`` now only accepts WCS objects which implement either the high, or
  low level APE 14 WCS API. All WCS objects are converted to a high level WCS
  object, so ``NDData.wcs`` now always returns a high level APE 14 object. Not
  all array slices are valid for wcs objects, so some slicing operations which
  used to work may now fail. [#9067]

astropy.stats
^^^^^^^^^^^^^

- The ``biweight_location``, ``biweight_scale``, and
  ``biweight_midvariance`` functions now allow for the ``axis``
  keyword to be a tuple of integers. [#9309]

- Added an ``ignore_nan`` option to the ``biweight_location``,
  ``biweight_scale``, and ``biweight_midvariance`` functions. [#9457]

- A numpy ``MaskedArray`` can now be input to the ``biweight_location``,
  ``biweight_scale``, and ``biweight_midvariance`` functions. [#9466]

- Removed the warning related to p0 in the Bayesian blocks algorithm. The
  caveat related to p0 is described in the docstring for ``Events``. [#9567]

astropy.table
^^^^^^^^^^^^^

- Improved the implementation of ``Table.replace_column()`` to provide
  a speed-up of 5 to 10 times for wide tables.  The method can now accept
  any input which convertible to a column of the correct length, not just
  ``Column`` subclasses. [#8902]

- Improved the implementation of ``Table.add_column()`` to provide a speed-up
  of 2 to 10 (or more) when adding a column to tables, with increasing benefit
  as the number of columns increases.  The method can now accept any input
  which is convertible to a column of the correct length, not just ``Column``
  subclasses. [#8933]

- Changed the implementation of ``Table.add_columns()`` to use the new
  ``Table.add_column()`` method.  In most cases the performance is similar
  or slightly faster to the previous implementation. [#8933]

- ``MaskedColumn.data`` will now return a plain ``MaskedArray`` rather than
  the previous (unintended) ``masked_BaseColumn``. [#8855]

- Added depth-wise stacking ``dstack()`` in higher level table operation.
  It help will in stacking table column depth-wise. [#8939]

- Added a new table equality method ``values_equal()`` which allows comparison
  table values to another table, list, or value, and returns an
  element-by-element equality table. [#9068]

- Added new ``join_type='cartesian'`` option to the ``join`` operation. [#9288]

- Allow adding a table column as a list of mixin-type objects, for instance
  ``t['q'] = [1 * u.m, 2 * u.m]``. [#9165]

- Allow table ``join()`` using any sortable key column (e.g. Time), not
  just ndarray subclasses. A column is considered sortable if there is a
  ``<column>.info.get_sortable_arrays()`` method that is implemented. [#9340]

- Added ``Table.iterrows()`` for making row-wise iteration faster. [#8969]

- Allow table to be initialized with a list of dict where the dict keys
  are not the same in every row. The table column names are the set of all keys
  found in the input data, and any missing key/value pairs are turned into
  missing data in the table. [#9425]

- Prevent unnecessary ERFA warnings when indexing by ``Time`` columns. [#9545]

- Added support for sorting tables which contain non-mutable mixin columns
  (like ``SkyCoord``) for which in-place item update is not allowed. [#9549]

- Ensured that inserting ``np.ma.masked`` (or any other value with a mask) into
  a ``MaskedColumn`` causes a masked entry to be inserted. [#9623]

- Fixed a bug that caused an exception when initializing a ``MaskedColumn`` from
  another ``MaskedColumn`` that has a structured dtype. [#9651]

astropy.tests
^^^^^^^^^^^^^

- The plugin that handles the custom header in the test output has been
  moved to the ``pytest-astropy-header plugin`` package. `See the README at
  <https://github.com/astropy/pytest-astropy-header>`__ for information about
  using this new plugin. [#9214]

astropy.time
^^^^^^^^^^^^

- Added a new time format ``ymdhms`` for representing times via year, month,
  day, hour, minute, and second attributes. [#7644]

- ``TimeDelta`` gained a ``to_value`` method, so that it becomes easier to
  use it wherever a ``Quantity`` with units of time could be used. [#8762]

- Made scalar ``Time`` and ``TimeDelta`` objects hashable based on JD, time
  scale, and location attributes. [#8912]

- Improved error message when bad input is used to initialize a ``Time`` or
  ``TimeDelta`` object and the format is specified. [#9296]

- Allow numeric time formats to be initialized with numpy ``longdouble``,
  ``Decimal`` instances, and strings.  One can select just one of these
  using ``in_subfmt``.  The output can be similarly set using ``out_subfmt``.
  [#9361]

- Introduce a new ``.to_value()`` method for ``Time`` (and adjusted the
  existing method for ``TimeDelta``) so that one can get values in a given
  ``format`` and possible ``subfmt`` (e.g., ``to_value('mjd', 'str')``. [#9361]

- Prevent unnecessary ERFA warnings when sorting ``Time`` objects. [#9545]

astropy.timeseries
^^^^^^^^^^^^^^^^^^

- Adding ``epoch_phase``, ``wrap_phase`` and ``normalize_phase`` keywords to
  ``TimeSeries.fold()`` to control the phase of the epoch and to return
  normalized phase rather than time for the folded TimeSeries. [#9455]

astropy.uncertainty
^^^^^^^^^^^^^^^^^^^

- ``Distribution`` was rewritten such that it deals better with subclasses.
  As a result, Quantity distributions now behave correctly with ``to`` methods
  yielding new distributions of the kind expected for the starting
  distribution, and ``to_value`` yielding ``NdarrayDistribution`` instances.
  [#9429, #9442]

- The ``pdf_*`` properties that were used to calculate statistical properties
  of ``Distribution`` instances were changed into methods. This allows one
  to pass parameters such as ``ddof`` to ``pdf_std`` and ``pdf_var`` (which
  generally should equal 1 instead of the default 0), and reflects that these
  are fairly involved calculations, not just "properties". [#9613]

astropy.units
^^^^^^^^^^^^^

- Support for unicode parsing. Currently supported are superscripts, Ohm,
  Ångström, and the micro-sign. [#9348]

- Accept non-unit type annotations in @quantity_input. [#8984]

- For numpy 1.17 and later, the new ``__array_function__`` protocol is used to
  ensure that all top-level numpy functions interact properly with
  ``Quantity``, preserving units also in operations like ``np.concatenate``.
  [#8808]

- Add equivalencies for surface brightness units to spectral_density. [#9282]

astropy.utils
^^^^^^^^^^^^^

- ``astropy.utils.data.download_file`` and
  ``astropy.utils.data.get_readable_fileobj`` now provides an ``http_headers``
  keyword to pass in specific request headers for the download. It also now
  defaults to providing ``User-Agent: Astropy`` and ``Accept: */*``
  headers. The default ``User-Agent`` value can be set with a new
  ``astropy.data.conf.default_http_user_agent`` configuration item.
  [#9508, #9564]

- Added a new ``astropy.utils.misc.unbroadcast`` function which can be used
  to return the smallest array that can be broadcasted back to the initial
  array. [#9209]

- The specific IERS Earth rotation parameter table used for time and
  coordinate transformations can now be set, either in a context or per
  session, using ``astropy.utils.iers.earth_rotation_table``. [#9244]

- Added ``export_cache`` and ``import_cache`` to permit transporting
  downloaded data to machines with no Internet connection. Several new
  functions are available to investigate the cache contents; e.g.,
  ``check_download_cache`` can be used to confirm that the persistent
  cache has not become damaged. [#9182]

- A new ``astropy.utils.iers.LeapSeconds`` class has been added to track
  leap seconds. [#9365]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Added a new ``time_support`` context manager/function for making it easy to
  plot and format ``Time`` objects in Matplotlib. [#8782]

- Added support for plotting any WCS compliant with the generalized (APE 14)
  WCS API with WCSAxes. [#8885, #9098]

- Improved display of information when inspecting ``WCSAxes.coords``. [#9098]

- Improved error checking for the ``slices=`` argument to ``WCSAxes``. [#9098]

- Added support for more solar frames in WCSAxes. [#9275]

- Add support for one dimensional plots to ``WCSAxes``. [#9266]

- Add a ``get_format_unit`` to ``wcsaxes.CoordinateHelper``. [#9392]

- ``WCSAxes`` now, by default, sets a default label for plot axes which is the
  WCS physical type (and unit) for that axis. This can be disabled using the
  ``coords[i].set_auto_axislabel(False)`` or by explicitly setting an axis
  label. [#9392]

- Fixed the display of tick labels when plotting all sky images that have a
  coord_wrap less than 360. [#9542]

astropy.wcs
^^^^^^^^^^^

- Added a ``astropy.wcs.wcsapi.pixel_to_pixel`` function that can be used to
  transform pixel coordinates in one dataset with a WCS to pixel coordinates
  in another dataset with a different WCS. This function is designed to be
  efficient when the input arrays are broadcasted views of smaller
  arrays. [#9209]

- Added a ``local_partial_pixel_derivatives`` function that can be used to
  determine a matrix of partial derivatives of each world coordinate with
  respect to each pixel coordinate. [#9392]

- Updated wcslib to v6.4. [#9125]

- Improved the  ``SlicedLowLevelWCS`` class in ``astropy.wcs.wcsapi`` to avoid
  storing chains of nested ``SlicedLowLevelWCS`` objects when applying multiple
  slicing operations in turn. [#9210]

- Added a ``wcs_info_str`` function to ``astropy.wcs.wcsapi`` to show a summary
  of an APE-14-compliant WCS as a string. [#8546, #9207]

- Added two new optional attributes to the APE 14 low-level WCS:
  ``pixel_axis_names`` and ``world_axis_names``. [#9156]

- Updated the WCS class to now correctly take and return ``Time`` objects in
  the high-level APE 14 API (e.g. ``pixel_to_world``. [#9376]

- ``SlicedLowLevelWCS`` now raises ``IndexError`` rather than ``ValueError`` on
  an invalid slice. [#9067]

- Added ``fit_wcs_from_points`` function to ``astropy.wcs.utils``. Fits a WCS
  object to set of matched detector/sky coordinates. [#9469]

- Fix various bugs in ``SlicedLowLevelWCS`` when the WCS being sliced was one
  dimensional. [#9693]


API Changes
-----------

astropy.constants
^^^^^^^^^^^^^^^^^

- Deprecated ``set_enabled_constants`` context manager. Use
  ``astropy.physical_constants`` and ``astropy.astronomical_constants``.
  [#9025]

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Removed the deprecated keyword argument ``interpolate_nan`` from
  ``convolve_fft``. [#9356]

- Removed the deprecated keyword argument ``stddev`` from
  ``Gaussian2DKernel``. [#9356]

- Deprecated and renamed ``MexicanHat1DKernel`` and ``MexicanHat2DKernel``
  to ``RickerWavelet1DKernel`` and ``RickerWavelet2DKernel``. [#9445]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Removed the ``recommended_units`` attribute from Representations; it was
  deprecated since 3.0. [#8892]

- Removed the deprecated frame attribute classes, ``FrameAttribute``,
  ``TimeFrameAttribute``, ``QuantityFrameAttribute``,
  ``CartesianRepresentationFrameAttribute``; deprecated since 3.0. [#9326]

- Removed ``longitude`` and ``latitude`` attributes from ``EarthLocation``;
  deprecated since 2.0. [#9326]

- The ``DifferentialAttribute`` for frame classes now passes through any input
  to the ``allowed_classes`` if only one allowed class is specified, i.e. this
  now allows passing a quantity in for frame attributes that use
  ``DifferentialAttribute``. [#9325]

- Removed the deprecated ``galcen_ra`` and ``galcen_dec`` attributes from the
  ``Galactocentric`` frame. [#9346]

astropy.extern
^^^^^^^^^^^^^^

- Remove the bundled ``six`` module. [#8315]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Masked column handling has changed, see ``astropy.table`` entry below.
  [#8789]

astropy.io.misc
^^^^^^^^^^^^^^^

- Masked column handling has changed, see ``astropy.table`` entry below.
  [#8789]

- Removed deprecated ``usecPickle`` kwarg from ``fnunpickle`` and
  ``fnpickle``. [#8890]

astropy.io.fits
^^^^^^^^^^^^^^^

- Masked column handling has changed, see ``astropy.table`` entry below.
  [#8789]

- ``io.fits.Header`` has been made safe for subclasses for copying and slicing.
  As a result of this change, the private subclass ``CompImageHeader``
  now always should be passed an explicit ``image_header``. [#9229]

- Removed the deprecated ``tolerance`` option in ``fitsdiff`` and
  ``io.fits.diff`` classes. [#9520]

- Removed deprecated keyword arguments for ``CompImageHDU``:
  ``compressionType``, ``tileSize``, ``hcompScale``, ``hcompSmooth``,
  ``quantizeLevel``. [#9520]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Changed ``pedantic`` argument to ``verify`` and change it to have three
  string-based options (``ignore``, ``warn``, and ``exception``) instead of
  just being a boolean. In addition, changed default to ``ignore``, which means
  that warnings will not be shown by default when loading VO tables. [#8715]

astropy.modeling
^^^^^^^^^^^^^^^^

- Eliminates support for compound classes (but not compound instances!) [#8769]

- Slicing compound models more restrictive. [#8769]

- Shape of parameters now includes n_models as dimension. [#8769]

- Parameter instances now hold values instead of models. [#8769]

- Compound model parameters now share instance and value with
  constituent models. [#8769]

- No longer possible to assign slices of parameter values to model parameters
  attribute (it is possible to replace it with a complete array). [#8769]

- Many private attributes and methods have changed (see documentation). [#8769]

- Deprecated ``BlackBody1D`` model and ``blackbody_nu`` and
  ``blackbody_lambda`` functions. [#9282]

- The deprecated ``rotations.rotation_matrix_from_angle`` was removed. [#9363]

- Deprecated and renamed ``MexicanHat1D`` and ``MexicanHat2D``
  to ``RickerWavelet1D`` and ``RickerWavelet2D``. [#9445]

- Deprecated ``modeling.utils.ExpressionTree``. [#9576]

astropy.stats
^^^^^^^^^^^^^

- Removed the ``iters`` keyword from sigma clipping stats functions. [#8948]

- Renamed the ``a`` parameter to ``data`` in biweight stat functions. [#8948]

- Renamed the ``a`` parameter to ``data`` in ``median_absolute_deviation``.
  [#9011]

- Renamed the ``conflevel`` keyword to ``confidence_level`` in
  ``poisson_conf_interval``. Usage of ``conflevel`` now issues
  ``AstropyDeprecationWarning``. [#9408]

- Renamed the ``conf`` keyword to ``confidence_level`` in
  ``binom_conf_interval`` and ``binned_binom_proportion``. Usage of ``conf``
  now issues ``AstropyDeprecationWarning``. [#9408]

- Renamed the ``conf_lvl`` keyword to ``confidence_level`` in
  ``jackknife_stats``. Usage of ``conf_lvl`` now issues
  ``AstropyDeprecationWarning``. [#9408]

astropy.table
^^^^^^^^^^^^^

- The handling of masked columns in the ``Table`` class has changed in a way
  that may impact program behavior. Now a ``Table`` with ``masked=False``
  may contain both ``Column`` and ``MaskedColumn`` objects, and adding a
  masked column or row to a table no longer "upgrades" the table and columns
  to masked.  This means that tables with masked data which are read via
  ``Table.read()`` will now always have ``masked=False``, though specific
  columns will be masked as needed. Two new table properties
  ``has_masked_columns`` and ``has_masked_values`` were added. See the
  `Masking change in astropy 4.0 section within
  <https://docs.astropy.org/en/v4.0.x/table/masking.html>`_ for
  details. [#8789]

- Table operation functions such as ``join``, ``vstack``, ``hstack``, etc now
  always return a table with ``masked=False``, though the individual columns
  may be masked as necessary. [#8957]

- Changed implementation of ``Table.add_column()`` and ``Table.add_columns()``
  methods.  Now it is possible add any object(s) which can be converted or
  broadcasted to a valid column for the table.  ``Table.__setitem__`` now
  just calls ``add_column``. [#8933]

- Changed default table configuration setting ``replace_warnings`` from
  ``['slice']`` to ``[]``.  This removes the default warning when replacing
  a table column that is a slice of another column. [#9144]

- Removed the non-public method
  ``astropy.table.np_utils.recarray_fromrecords``. [#9165]

astropy.tests
^^^^^^^^^^^^^

- In addition to ``DeprecationWarning``, now ``FutureWarning`` and
  ``ImportWarning`` would also be turned into exceptions. [#8506]

- ``warnings_to_ignore_by_pyver`` option in
  ``enable_deprecations_as_exceptions()`` has changed. Please refer to API
  documentation. [#8506]

- Default settings for ``warnings_to_ignore_by_pyver`` are updated to remove
  very old warnings that are no longer relevant and to add a new warning
  caused by ``pytest-doctestplus``. [#8506]

astropy.time
^^^^^^^^^^^^

- ``Time.get_ut1_utc`` now uses the auto-updated ``IERS_Auto`` by default,
  instead of the bundled ``IERS_B`` file. [#9226]

- Time formats that do not use ``val2`` now raise ValueError instead of
  silently ignoring a provided value. [#9373]

- Custom time formats can now accept floating-point types with extended
  precision. Existing time formats raise exceptions rather than discarding
  extended precision through conversion to ordinary floating-point. [#9368]

- Time formats (implemented in subclasses of ``TimeFormat``) now have
  their input and output routines more thoroughly validated, making it more
  difficult to create damaged ``Time`` objects. [#9375]

- The ``TimeDelta.to_value()`` method now can also take the ``format`` name
  as its argument, in which case the value will be calculated using the
  ``TimeFormat`` machinery. For this case, one can also pass a ``subfmt``
  argument to retrieve the value in another form than ``float``. [#9361]

astropy.timeseries
^^^^^^^^^^^^^^^^^^

- Keyword ``midpoint_epoch`` is renamed to ``epoch_time``. [#9455]

astropy.uncertainty
^^^^^^^^^^^^^^^^^^^

- ``Distribution`` was rewritten such that it deals better with subclasses.
  As a result, Quantity distributions now behave correctly with ``to`` methods
  yielding new distributions of the kind expected for the starting distribution,
  and ``to_value`` yielding ``NdarrayDistribution`` instances. [#9442]

astropy.units
^^^^^^^^^^^^^

- For consistency with ``ndarray``, scalar ``Quantity.value`` will now return
  a numpy scalar rather than a python one.  This should help keep track of
  precision better, but may lead to unexpected results for the rare cases
  where numpy scalars behave differently than python ones (e.g., taking the
  square root of a negative number). [#8876]

- Removed the ``magnitude_zero_points`` module, which was deprecated in
  favour of ``astropy.units.photometric`` since 3.1. [#9353]

- ``EquivalentUnitsList`` now has a ``_repr_html_`` method to output a HTML
  table on a call to ``find_equivalent_units`` in Jupyter notebooks. [#9495]

astropy.utils
^^^^^^^^^^^^^

- ``download_file`` and related functions now accept a list of fallback
  sources, and they are able to update the cache at the user's request. [#9182]

- Allow ``astropy.utils.console.ProgressBarOrSpinner.map`` and
  ``.map_unordered`` to take an argument ``multiprocessing_start_method`` to
  control how subprocesses are started; the different methods (``fork``,
  ``spawn``, and ``forkserver``) have different implications in terms of
  security, efficiency, and behavioural anomalies. The option is useful in
  particular for cross-platform testing because Windows supports only ``spawn``
  while Linux defaults to ``fork``. [#9182]

- All operations that act on the astropy download cache now take an argument
  ``pkgname`` that allows one to specify which package's cache to use.
  [#8237, #9182]

- Removed deprecated ``funcsigs`` and ``futures`` from
  ``astropy.utils.compat``. [#8909]

- Removed the deprecated ``astropy.utils.compat.numpy`` module. [#8910]

- Deprecated ``InheritDocstrings`` as it is natively supported by
  Sphinx 1.7 or higher. [#8881]

- Deprecated ``astropy.utils.timer`` module, which has been moved to
  ``astroquery.utils.timer`` and will be part of ``astroquery`` 0.4.0. [#9038]

- Deprecated ``astropy.utils.misc.set_locale`` function, as it is meant for
  internal use only. [#9471]

- The implementation of ``data_info.DataInfo`` has changed (for a considerable
  performance boost). Generally, this should not affect simple subclasses, but
  because the class now uses ``__slots__`` any attributes on the class have to
  be explicitly given a slot. [#8998]

- ``IERS`` tables now use ``nan`` to mark missing values
  (rather than ``1e20``). [#9226]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- The default ``clip`` value is now ``False`` in ``ImageNormalize``. [#9478]

- The default ``clip`` value is now ``False`` in ``simple_norm``.
  [#9698]

- Infinite values are now excluded when calculating limits in
  ``ManualInterval`` and ``MinMaxInterval``.  They were already excluded in
  all other interval classes. [#9480]


Bug Fixes
---------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Fixed ``nan_treatment='interpolate'`` option to ``convolve_fft`` to properly
  take into account ``fill_value``. [#8122]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- The ``QuantityAttribute`` class now supports a None default value if a unit
  is specified. [#9345]

- When ``Representation`` classes with the same name are defined, this no
  longer leads to a ``ValueError``, but instead to a warning and the removal
  of both from the name registry (i.e., one either has to use the class itself
  to set, e.g., ``representation_type``, or refer to the class by its fully
  qualified name). [#8561]

astropy.io.fits
^^^^^^^^^^^^^^^

- Implemented skip (after warning) of header cards with reserved
  keywords in ``table_to_hdu``. [#9390]

- Add ``AstropyDeprecationWarning`` to ``read_table_fits`` when ``hdu=`` is
  selected, but does not match single present table HDU. [#9512]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Address issue #8995 by ignoring BINARY2 null mask bits for string values
  on parsing a VOTable.  In this way, the reader should never create masked
  values for string types. [#9057]

- Corrected a spurious warning issued for the ``value`` attribute of the
  ``<OPTION>`` element in VOTable, as well as a test that erroneously
  treated the warning as acceptable.  [#9470]

astropy.nddata
^^^^^^^^^^^^^^

- ``Cutout2D`` will now get the WCS from its first argument if that argument
  has with WCS property. [#9492]

- ``overlap_slices`` will now raise a ``ValueError`` if the input
  position contains any non-finite values (e.g. NaN or inf). [#9648]

astropy.stats
^^^^^^^^^^^^^

- Fixed a bug where ``bayesian_blocks`` returned a single edge. [#8560]

- Fixed input data type validation for ``bayesian_blocks`` to work int
  arrays. [#9513]

astropy.table
^^^^^^^^^^^^^

- Fix bug where adding a column consisting of a list of masked arrays was
  dropping the masks. [#9048]

- ``Quantity`` columns with custom units can now round-trip via FITS tables,
  as long as the custom unit is enabled during reading (otherwise, the unit
  will become an ``UnrecognizedUnit``). [#9015]

- Fix bug where string values could be truncated when inserting into a
  ``Column`` or ``MaskedColumn``, or when adding or inserting a row containing
  string values. [#9559]

astropy.time
^^^^^^^^^^^^

- Fix bug when ``Time`` object is created with only masked elements. [#9624]

- Fix inaccuracy when converting between TimeDelta and datetime.timedelta.
  [#9679]

astropy.units
^^^^^^^^^^^^^

- Ensure that output from test functions of and comparisons between quantities
  can be stored into pre-allocated output arrays (using ``out=array``) [#9273]

astropy.utils
^^^^^^^^^^^^^

- For the default ``IERS_Auto`` table, which combines IERS A and B values, the
  IERS nutation parameters "dX_2000A" and "dY_2000A" are now also taken from
  the actual IERS B file rather than from the B values stored in the IERS A
  file.  Any differences should be negligible for any practical application,
  but this may help exactly reproducing results. [#9237]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Calling ``WCSAxes.set_axis_off()`` now correctly turns off drawing the Axes.
  [#9411]

- Fix incorrect transformation behavior in ``WCSAxes.plot_coord`` and correctly
  handle when input coordinates are not already in spherical representations.
  [#8927]

- Fixed ``ImageNormalize`` so that when it is initialized without
  ``data`` it will still use the input ``interval`` class. [#9698]

- Fixed ``ImageNormalize`` to handle input data with non-finite
  values. [#9698]

astropy.wcs
^^^^^^^^^^^

- Fix incorrect value returned by
  ``wcsapi.HighLevelWCSWrapper.axis_correlation_matrix``. [#9554]

- Fix NaN-masking of world coordinates when some but not all of the coordinates
  were flagged as invalid by WCSLIB. This occurred for example with WCS with >2
  dimensions where two of the dimensions were celestial coordinates and pixel
  coordinates outside of the 'sky' were converted to world coordinates -
  previously all world coordinates were masked even if uncorrelated with the
  celestial axes, but this is no longer the case. [#9688]

- The default WCS to celestial frame mapping for world coordinate systems that
  specify ``TLON`` and ``TLAT`` coordinates will now return an ITRS frame with
  the representation class set to ``SphericalRepresentation``. This fixes a bug
  that caused ``WCS.pixel_to_world`` to raise an error for such world
  coordinate systems. [#9609]

- ``FITSWCSAPIMixin`` now returns tuples not lists from ``pixel_to_world`` and
  ``world_to_pixel``. [#9678]


Other Changes and Additions
---------------------------

- Versions of Python <3.6 are no longer supported. [#8955]

- Matplotlib 2.1 and later is now required. [#8787]

- Versions of Numpy <1.16 are no longer supported. [#9292]

- Updated the bundled CFITSIO library to 3.470. See
  ``cextern/cfitsio/docs/changes.txt`` for additional information. [#9233]

- The bundled ERFA was updated to version 1.7.0. This is based on
  SOFA 20190722. This includes a fix to avoid precision loss for negative
  JDs, and also includes additional routines to allow updates to the
  leap-second table. [#9323, #9734]

- The default server for the IERS data files has been updated to reflect
  long-term downtime of the canonical USNO server. [#9487, #9508]



Version 3.2.3 (2019-10-27)
==========================

Other Changes and Additions
---------------------------

- Updated IERS A URLs due to USNO prolonged maintenance. [#9443]



Version 3.2.2 (2019-10-07)
==========================

Bug fixes
---------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Fixed a bug in ``discretize_oversample_1D/2D()`` from
  ``astropy.convolution.utils``, which might occasionally introduce unexpected
  oversampling grid dimensions due to a numerical precision issue. [#9293]

- Fixed a bug [#9168] where having a kernel defined using unitless astropy
  quantity objects would result in a crash [#9300]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Fix concatenation of representations for cases where the units were different.
  [#8877]

- Check for NaN values in catalog and match coordinates before building and
  querying the ``KDTree`` for coordinate matching. [#9007]

- Fix sky coordinate matching when a dimensionless distance is provided. [#9008]

- Raise a faster and more meaningful error message when differential data units
  are not compatible with a containing representation's units. [#9064]

- Changed the timescale in ICRS to CIRS from 'tdb' to 'tt' conversion and
  vice-versa, as the erfa function that gets called in the process, pnm06a
  accepts time in TT. [#9079]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fixed the fast reader when used in parallel and with the multiprocessing
  'spawn' method (which is the default on MacOS X with Python 3.8 and later),
  and enable parallel fast reader on Windows. [#8853]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fixes bug where an invalid TRPOS<n> keyword was being generated for FITS
  time column when no location was available. [#8784]

- Fixed a wrong exception when converting a Table with a unit that is not FITS
  compliant and not convertible to a string using ``format='fits'``. [#8906]

- Fixed an issue with A3DTABLE extension that could not be read. [#9012]

- Fixed the update of the header when creating GroupsHDU from data. [#9216]

astropy.nddata
^^^^^^^^^^^^^^

- Fix to ``add_array``, which now accepts ``array_small`` having dimensions
  equal to ``array_large``, instead of only allowing smaller sizes of
  arrays. [#9118]

astropy.stats
^^^^^^^^^^^^^

- Fixed ``median_absolute_deviation`` for the case where ``ignore_nan=True``
  and an input masked array contained both NaNs and infs. [#9307]

astropy.table
^^^^^^^^^^^^^

- Comparisons between ``Column`` instances and ``Quantity`` will now
  correctly take into account the unit (as was already the case for
  regular operations such as addition). [#8904]

astropy.time
^^^^^^^^^^^^

- Allow ``Time`` to be initialized with an empty value for all formats. [#8854]

- Fixed a troubling bug in which ``Time`` could loose precision, with deviations
  of 300 ns. [#9328]

astropy.timeseries
^^^^^^^^^^^^^^^^^^

- Fixed handling of ``Quantity`` input data for all methods of
  ``LombScarge.false_alarm_probabilty``. [#9246]

astropy.units
^^^^^^^^^^^^^

- Allow conversion of ``Column`` with logarithmic units to a suitable
  ``Quantity`` subclass if ``subok=True``. [#9188]

- Ensured that we simplify powers to smaller denominators if that is
  consistent within rounding precision. [#9267]

- Ensured that the powers shown in a unit's repr are always correct,
  not oversimplified. [#9267]

astropy.utils
^^^^^^^^^^^^^

- Fixed ``find_api_page`` access by using custom request headers and HTTPS
  when version is specified. [#9032]

- Make ``download_file`` (and by extension ``get_readable_fileobj`` and others)
  check the size of downloaded files against the size claimed by the server.
  [#9302]

- Fix ``find_current_module`` so that it works properly if astropy is being used
  inside a bundle such as that produced by PyInstaller. [#8845]

- Fix path to renamed classes, which previously included duplicate path/module
  information under certain circumstances. [#8845]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Silence numpy runtime warnings in ``WCSAxes`` when drawing grids. [#8882]

astropy.wcs
^^^^^^^^^^^

- Fixed equality test between ``cunit`` where the first element was equal but
  the following elements differed. [#9154]

- Fixed a crash while loading a WCS from headers containing duplicate SIP
  keywords. [#8893]

- Fixed a possible buffer overflow when using too large negative indices for
  ``cunit`` or ``ctype`` [#9151]

- Fixed reference counting in ``WCSBase.__init__`` [#9166]

- Fix ``SlicedLowLevelWCS`` ``world_to_pixel_values`` and
  ``pixel_to_world_values`` when inputs need broadcasting to the same shape.
  (i.e. when one input is sliced out) [#9250]

- Fixed a bug that caused ``WCS.array_shape``, ``WCS.pixel_shape`` and
  ``WCS.pixel_bounds`` to be incorrect after using ``WCS.sub``. [#9095]


Other Changes and Additions
---------------------------

- Fixed a bug that caused files outside of the astropy module directory to be
  included as package data, resulting in some cases in errors when doing
  repeated builds. [#9039]



Version 3.2.1 (2019-06-14)
==========================

Bug fixes
---------

astropy.io.fits
^^^^^^^^^^^^^^^

- Avoid reporting a warning with ``BinTableHDU.from_columns`` with keywords that
  are not provided by the user.  [#8838]

- Fix ``Header.fromfile`` to work on FITS files. [#8713]

- Fix reading of empty ``BinTableHDU`` when stored in a gzip-compressed file.
  [#8848]

astropy.table
^^^^^^^^^^^^^

- Fix a problem where mask was dropped when creating a ``MaskedColumn``
  from a list of ``MaskedArray`` objects. [#8826]

astropy.wcs
^^^^^^^^^^^

- Added ``None`` to be displayed as a ``world_axis_physical_types`` in
  the ``WCS`` repr, as ``None`` values are now supported in ``APE14``. [#8811]



Version 3.2 (2019-06-10)
========================

New Features
------------

astropy.constants
^^^^^^^^^^^^^^^^^

- Add CODATA 2018 constants but not make them default because the
  redefinition of SI units that will follow has not been implemented
  yet. [#8595]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- New ``BarycentricMeanEcliptic``, ``HeliocentricTrueEcliptic`` and
  ``GeocentricTrueEcliptic`` frames.
  The ecliptic frames are no longer considered experimental. [#8394]

- The default time scale for epochs like 'J2000' or 'B1975' is now "tt",
  which is the correct one for 'J2000' and avoids leap-second warnings
  for epochs in the far future or past. [#8600]

astropy.extern
^^^^^^^^^^^^^^

- Bundled ``six`` now emits ``AstropyDeprecationWarning``. It will be removed
  in 4.0. [#8323]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- IPAC tables now output data types of ``float`` instead of ``double``, or
  ``int`` instead of ``long``, based on the column ``dtype.itemsize``. [#8216]

- Update handling of MaskedColumn columns when using the 'data_mask' serialization
  method.  This can make writing ECSV significantly faster if the data do not
  actually have any masked values. [#8447]

- Fixed a bug that caused newlines to be incorrect when writing out ASCII tables
  on Windows (they were ``\r\r\n`` instead of ``\r\n``). [#8659]

astropy.io.misc
^^^^^^^^^^^^^^^

- Implement serialization of ``TimeDelta`` in ASDF. [#8285]

- Implement serialization of ``EarthLocation`` in ASDF. [#8286]

- Implement serialization of ``SkyCoord`` in ASDF. [#8284]

- Support serialization of Astropy tables with mixin columns in ASDF. [#8337]

- No warnings when reading HDF5 files with only one table and no ``path=``
  argument [#8483]

- The HDF5 writer will now create a default table instead of raising an
  exception when ``path=`` is not specified and when writing to empty/new HDF5
  files. [#8553]

astropy.io.fits
^^^^^^^^^^^^^^^

- Optimize parsing of cards within the ``Header`` class. [#8428]

- Optimize the parsing of headers to get the structural keywords that are
  needed to find extensions. Thanks to this, getting a random HDU from a file
  with many extensions is much faster than before, in particular when the
  extension headers contain many keywords. [#8502]

-  Change behavior of FITS undefined value in ``Header`` such that ``None``
   is used in Python to represent FITS undefined when using dict interface.
   ``Undefined`` can also be assigned and is translated to ``None``.
   Previously setting a header card value to ``None`` resulted in an
   empty string field rather than a FITS undefined value. [#8572]

- Allow ``Header.fromstring`` and ``Card.fromstring`` to accept ``bytes``.
  [#8707]

astropy.io.registry
^^^^^^^^^^^^^^^^^^^

- Implement ``Table`` reader and writer for ``ASDF``. [#8261]

- Implement ``Table`` reader and writer methods to wrap ``pandas`` I/O methods
  for CSV, Fixed width format, HTML, and JSON. [#8381]

- Add ``help()`` and ``list_formats()`` methods to unified I/O ``read`` and
  ``write`` methods. For example ``Table.read.help()`` gives help on available
  ``Table`` read formats and ``Table.read.help('fits')`` gives detailed
  help on the arguments for reading FITS table file. [#8255]

astropy.table
^^^^^^^^^^^^^

- Initializing a table with ``Table(rows=...)``, if the first item is an ``OrderedDict``,
  now uses the column order of the first row. [#8587]

- Added new pprint_all() and pformat_all() methods to Table. These two new
  methods print the entire table by default. [#8577]

- Removed restriction of initializing a Table from a dict with copy=False. [#8541]

- Improved speed of table row access by a factor of about 2-3.  Improved speed
  of Table len() by a factor of around 3-10 (depending on the number of columns).
  [#8494]

- Improved the Table - pandas ``DataFrame`` interface (``to_pandas()`` and
  ``from_pandas()``).  Mixin columns like ``Time`` and ``Quantity`` can now be
  converted to pandas by flattening the columns as necessary to plain
  columns.  ``Time`` and ``TimeDelta`` columns get converted to
  corresponding pandas date or time delta types.  The ``DataFrame``
  index is now handled in the conversion methods. [#8247]

- Added ``rename_columns`` method to rename multiple columns in one call.
  [#5159, #8070]

- Improved Table performance by reducing unnecessary calls to copy and deepcopy,
  especially as related to the table and column ``meta`` attributes.  Changed the
  behavior when slicing a table (either in rows or with a list of column names)
  so now the sliced output gets a light (key-only) copy of ``meta`` instead of a
  deepcopy.  Changed the ``Table.meta`` class-level descriptor so that assigning
  directly to ``meta``, e.g. ``tbl.meta = new_meta`` no longer does a deepcopy
  and instead just directly assigns the ``new_meta`` object reference.  Changed
  Table initialization so that input ``meta`` is copied only if ``copy=True``.
  [#8404]

- Improved Table slicing performance with internal implementation changes
  related to column attribute access and certain input validation. [#8493]

- Added ``reverse`` argument to the ``sort`` and ``argsort`` methods to allow
  sorting in reverse order. [#8528]

- Improved ``Table.sort()`` performance by removing ``self[keys]`` from code
  which is creating deep copies of ``meta`` attribute and adding a new keyword
  ``names`` in ``get_index()`` to get index by using a list or tuple containing
  names of columns. [#8570]

- Expose ``represent_mixins_as_columns`` as a public function in the
  ``astropy.table`` subpackage.  This previously-private function in the
  ``table.serialize`` module is used to represent mixin columns in a Table as
  one or more plain Column objects. [#7729]

astropy.timeseries
^^^^^^^^^^^^^^^^^^

- Added a new astropy.timeseries sub-package to represent and manipulate
  sampled and binned time series. [#8540]

- The ``BoxLeastSquares`` and ``LombScargle`` classes have been moved to
  ``astropy.timeseries.periodograms`` from ``astropy.stats``. [#8591]

- Added the ability to provide absolute ``Time`` objects to the
  ``BoxLeastSquares`` and ``LombScargle`` periodogram classes. [#8599]

- Added model inspection methods (``model_parameters()``, ``design_matrix()``,
  and ``offset()``) to ``astropy.timeseries.LombScargle`` class [#8397].

astropy.units
^^^^^^^^^^^^^

- ``Quantity`` overrides of ``ndarray`` methods such as ``sum``, ``min``,
  ``max``, which are implemented via reductions, have been removed since they
  are dealt with in ``Quantity.__array_ufunc__``. This should not affect
  subclasses, but they may consider doing similarly. [#8316]  Note that this
  does not include methods that use more complicated python code such as
  ``mean``, ``std`` and ``var``. [#8370]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^
- Added ``CompositeStretch``, which inherits from ``CompositeTransform`` and
  also ``BaseStretch`` so that it can be used with ``ImageNormalize``. [#8564]

- Added a ``log_a`` argument to the ``simple_norm`` method. Similar to the
  exposing of the ``asinh_a`` argument for ``AsinhStretch``, the new
  ``log_a`` argument is now exposed for ``LogStretch``. [#8436]

astropy.wcs
^^^^^^^^^^^

- WCSLIB was updated to v 6.2.
  This adds support for time-related WCS keywords (WCS Paper VII).
  FITS headers containing ``Time`` axis are parsed and the axis is included in
  the WCS object. [#8592]

- The ``OBSGEO`` attribute as expanded to 6 members - ``XYZLBH``. [#8592]

- Added a new class ``SlicedLowLevelWCS`` in ``astropy.wcs.wcsapi`` that can be
  used to slice any WCS that conforms to the ``BaseLowLevelWCS`` API. [#8546]

- Updated implementation of ``WCS.__getitem__`` and ``WCS.slice`` to now return
  a ``SlicedLowLevelWCS`` rather than raising an error when reducing the
  dimensionality of the WCS. [#8546]


API Changes
-----------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- ``QuantityAttribute`` no longer has a default value for ``default``.  The
  previous value of None was misleading as it always was an error. [#8450]

- The default J2000 has been changed to use be January 1, 2000 12:00 TT instead
  of UTC.  This is more in line with convention. [#8594]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- IPAC tables now output data types of ``float`` instead of ``double``, or
  ``int`` instead of ``long``, based on the column ``dtype.itemsize``. [#8216]

astropy.io.misc
^^^^^^^^^^^^^^^

- Unit equivalencies can now be serialized to ASDF. [#8252]

astropy.modeling
^^^^^^^^^^^^^^^^

- Composition of model classes is deprecated and will be removed in 4.0.
  Composition of model instances remain unaffected. [#8234, #8408]

astropy.stats
^^^^^^^^^^^^^

- The ``BoxLeastSquares`` and ``LombScargle`` classes have been moved to the
  ``astropy.timeseries.periodograms`` module and will now emit a deprecation
  warning when imported from ``astropy.stats``. [#8591]

astropy.table
^^^^^^^^^^^^^

- Converting an empty table to an array using ``as_array`` method now returns
  an empty array instead of ``None``. [#8647]

- Changed the behavior when slicing a table (either in rows or with a list of column
  names) so now the sliced output gets a light (key-only) copy of ``meta`` instead of
  a deepcopy.  Changed the ``Table.meta`` class-level descriptor so that assigning
  directly to ``meta``, e.g. ``tbl.meta = new_meta`` no longer does a deepcopy
  and instead just directly assigns the ``new_meta`` object reference. Changed
  Table initialization so that input ``meta`` is copied only if ``copy=True``.
  [#8404]

- Added a keyword ``names`` in ``Table.as_array()``.  If provided this specifies
  a list of column names to include for the returned structured array. [#8532]

astropy.tests
^^^^^^^^^^^^^

- Removed ``pytest_plugins`` as they are completely broken for ``pytest>=4``.
  [#7786]

- Removed the ``astropy.tests.plugins.config`` plugin and removed the
  ``--astropy-config-dir`` and ``--astropy-cache-dir`` options from
  testing. Please use caching functionality that is natively in ``pytest``.
  [#7787, #8489]

astropy.time
^^^^^^^^^^^^

- The default time scale for epochs like 'J2000' or 'B1975' is now "tt",
  which is the correct one for 'J2000' and avoids leap-second warnings
  for epochs in the far future or past. [#8600]

astropy.units
^^^^^^^^^^^^^

- Unit equivalencies can now be introspected. [#8252]

astropy.wcs
^^^^^^^^^^^

- The ``world_to_pixel``, ``world_to_array_index*``, ``pixel_to_world*`` and
  ``array_index_to_world*`` methods now all consistently return scalars, arrays,
  or objects not wrapped in a one-element tuple/list when only one scalar,
  array, or object (as was previously already the case for ``WCS.pixel_to_world``
  and ``WCS.array_index_to_world``). [#8663]

astropy.utils
^^^^^^^^^^^^^

- It is now possible to control the number of cores used by ``ProgressBar.map``
  by passing a positive integer as the ``multiprocess`` keyword argument. Use
  ``True`` to use all cores. [#8083]


Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- ``BarycentricTrueEcliptic``, ``HeliocentricTrueEcliptic`` and
  ``GeocentricTrueEcliptic`` now use the correct transformation
  (including nutation), whereas the new ``*MeanEcliptic`` classes
  use the nutation-free transformation. [#8394]

- Representations with ``float32`` coordinates can now be transformed,
  although the output will always be ``float64``. [#8759]

- Fixed bug that prevented using differentials with HCRS<->ICRS
  transformations. [#8794]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fixed a bug where an exception was raised when writing a table which includes
  mixin columns (e.g. a Quantity column) and the output format was specified
  using the ``formats`` keyword. [#8681]

astropy.io.misc
^^^^^^^^^^^^^^^

- Fixed bug in ASDF tag that inadvertently introduced dependency on ``pytest``.
  [#8456]

astropy.modeling
^^^^^^^^^^^^^^^^

- Fixed slowness for certain compound models consisting of large numbers
  of multi-input models [#8338, #8349]

- Fixed bugs in fitting of compound models with units. [#8369]

astropy.nddata
^^^^^^^^^^^^^^

- Fixed bug in reading multi-extension FITS files written by earlier versions
  of ``CCDData``. [#8534]

- Fixed two errors in the way ``CCDData`` handles FITS files with WCS in the
  header. Some of the WCS keywords that should have been removed from the
  header were not, potentially leading to FITS files with inconsistent
  WCS. [#8602]

astropy.table
^^^^^^^^^^^^^

- Fixed a bug when initializing from an empty list: ``Table([])`` no longer
  results in a crash. [#8647]

- Fixed a bug when initializing from an existing ``Table``.  In this case the
  input ``meta`` argument was being ignored.  Now the input ``meta``, if
  supplied, will be used as the ``meta`` for the new ``Table``. [#8404]

- Fix the conversion of bytes values to Python ``str`` with ``Table.tolist``.
  [#8739]

astropy.time
^^^^^^^^^^^^

- Fixed a number of issues to ensure a consistent output type resulting from
  multiplication or division involving a ``TimeDelta`` instance. The output is
  now always a ``TimeDelta`` if the result is a time unit (like u.s or u.d),
  otherwise it will be a ``Quantity``. [#8356]

- Multiplication between two ``TimeDelta`` instances is now possible, resulting
  in a ``Quantity`` with units of time squared (division already correctly
  resulted in a dimensionless ``Quantity``). [#8356]

- Like for comparisons, addition, and subtraction of ``Time`` instances with
  with non-time instances, multiplication and division of ``TimeDelta``
  instances with incompatible other instances no longer immediately raise an
  ``UnitsError`` or ``TypeError`` (depending on the other instance), but
  rather go through the regular Python mechanism of ``TimeDelta`` returning
  ``NotImplemented`` (which will lead to a regular ``TypeError`` unless the
  other instance can handle ``TimeDelta``). [#8356]

- Corrected small rounding errors that could cause the ``jd2`` values in
  ``Time`` to fall outside the range of -0.5 to 0.5. [#8763]

astropy.units
^^^^^^^^^^^^^

- Added a ``Quantity.to_string`` method to add flexibility to the string formatting
  of quantities. It produces unadorned or LaTeX strings, and accepts two different
  sets of delimiters in the latter case: ``inline`` and ``display``. [#8313]

- Ensure classes that mimic quantities by having a ``unit`` attribute and/or
  ``to`` and ``to_value`` methods can be properly used to initialize ``Quantity``
  or set ``Quantity`` instance items. [#8535]

- Add support for ``<<`` to create logarithmic units. [#8290]

- Add support for the ``clip`` ufunc, which in numpy 1.17 is used to implement
  ``np.clip``.  As part of that, remove the ``Quantity.clip`` method under
  numpy 1.17. [#8747]

- Fix parsing of numerical powers in FITS-compatible units. [#8251]

astropy.wcs
^^^^^^^^^^^

- Added a ``PyUnitListProxy_richcmp`` method in ``UnitListProxy`` class to enable
  ``WCS.wcs.cunit`` equality testing. It helps to check whether the two instances of
  ``WCS.wcs.cunit`` are equal or not by comparing the data members of
  ``UnitListProxy`` class [#8480]

- Fixed ``SlicedLowLevelWCS`` when ``array_shape`` is ``None``. [#8649]

- Do not attempt to delete repeated distortion keywords multiple times when
  loading distortions with ``_read_distortion_kw`` and
  ``_read_det2im_kw``. [#8777]


Other Changes and Additions
---------------------------

- Update bundled expat to 2.2.6. [#8343]

- Added instructions for uploading releases to Zenodo. [#8395]

- The bug fixes to the behaviour of ``TimeDelta`` for multiplication and
  division, which ensure that the output is now always a ``TimeDelta`` if the
  result is a time unit (like u.s or u.d) and otherwise a ``Quantity``, imply
  that sometimes the output type will be different than it was before. [#8356]

- For types unrecognized by ``TimeDelta``, multiplication and division now
  will consistently return a ``TypeError`` if the other instance cannot handle
  ``TimeDelta`` (rather than ``UnitsError`` or ``TypeError`` depending on
  presumed abilities of the other instance). [#8356]

- Multiplication between two ``TimeDelta`` instances will no longer result in
  an ``OperandTypeError``, but rather result in a ``Quantity`` with units of
  time squared (division already correctly resulted in a dimensionless
  ``Quantity``). [#8356]

- Made running the tests insensitive to local user configuration when running
  the tests in parallel mode or directly with pytest. [#8727]

- Added a narrative style guide to the documentation for contributor reference.
  [#8588]

- Ensure we call numpy equality functions in a way that reduces the number
  of ``DeprecationWarning``. [#8755]

Installation
^^^^^^^^^^^^

- We now require setuptools 30.3.0 or later to install the core astropy
  package. [#8240]

- We now define groups of dependencies that can be installed with pip, e.g.
  ``pip install astropy[all]`` (to install all optional dependencies). [#8198]



Version 3.1.2 (2019-02-23)
==========================

Bug fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Convert the default of ``QuantityAttribute``, thereby catching the error case
  case of it being set to None at attribute creation, and giving a more useful
  error message in the process. [#8300]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Fix elliptic analytical solution for comoving distance. Only
  relevant for non-flat cosmologies without radiation and ``Om0`` > ``Ode0``.
  [#8391]

astropy.modeling
^^^^^^^^^^^^^^^^

- Fixed slowness for certain compound models consisting of large numbers
  of multi-input models [#8338, #8349]

astropy.visualization.wcsaxes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fix a bug that caused an error when passing an array with all values the same
  to contour or contourf. [#8321]

- Fix a bug that caused contour and contourf to return None instead of the
  contour set. [#8321]


Version 3.1.1 (2018-12-31)
==========================

Bug fixes
---------

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fix error when writing out empty table. [#8279]

astropy.io.fits
^^^^^^^^^^^^^^^

- ``fitsdiff --ignore-hdus`` now prints input filenames in the diff report
  instead of ``<HDUList object at 0x1150f9778>``. [#8295]

astropy.units
^^^^^^^^^^^^^

- Ensure correctness of units when raising to a negative power. [#8263]

- Fix ``with_H0`` equivalency to use the correct direction of
  conversion. [#8292]



Version 3.1 (2018-12-06)
========================

New Features
------------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- ``convolve`` now accepts any array-like input, not just ``numpy.ndarray`` or
  lists. [#7303]

- ``convolve`` Now raises AstropyUserWarning if nan_treatment='interpolate' and
  preserve_nan=False and NaN values are present post convolution. [#8088]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- The ``SkyCoord.from_name`` constructor now has the ability to create
  coordinate objects by parsing object catalogue names that have embedded
  J-coordinates. [#7830]

- The new function ``make_transform_graph_docs`` can be used to create a
  docstring graph from a custom ``TransformGraph`` object. [#7135]

- ``KDTree`` for catalog matching is now built with sliding midpoint rule
  rather than standard.  In code, this means setting ``compact_nodes=False``
  and ``balanced_tree=False`` in ``cKDTree``. The sliding midpoint rule is much
  more suitable for catalog matching, and results in 1000x speedup in some
  cases. [#7324]

- Additional information about a site loaded from the Astropy sites registry is
  now available in ``EarthLocation.info.meta``. [#7857]

- Added a ``concatenate_representations`` function to combine coordinate
  representation data and any associated differentials. [#7922]

- ``BaseCoordinateFrame`` will now check for a method named
  ``_astropy_repr_in_frame`` when constructing the string forms of attributes.
  Allowing any class to control how ``BaseCoordinateFrame`` represents it when
  it is an attribute of a frame. [#7745]

- Some rarely-changed attributes of frame classes are now cached, resulting in
  speedups (up to 50% in some cases) when creating new scalar frame or
  ``SkyCoord`` objects. [#7949, #5952]

- Added a ``directional_offset_by`` method to ``SkyCoord`` that computes a new
  coordinate given a coordinate, position angle, and angular separation [#5727]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- The default cosmology has been changed from ``WMAP9`` to ``Planck15``. [#8123]

- Distance calculations with ``LambaCDM`` with no radiation (T_CMB0=0)
  are now 20x faster by using elliptic integrals for non-flat cases. [#7155]

- Distance calculations with ``FlatLambaCDM`` with no radiation (T_CMB0=0)
  are now 20x faster by using the hypergeometric function solution
  for this special case. [#7087]

- Age calculations with ``FlatLambdaCDM`` with no radiation (Tcmb0=0)
  are now 1000x faster by using analytic solutions instead of integrating.
  [#7117]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Latex reader now ignores ``\toprule``, ``\midrule``, and ``\bottomrule``
  commands. [#7349]

- Added the RST (Restructured-text) table format and the fast version of the
  RDB reader to the set of formats that are guessed by default. [#5578]

- The read trace (used primarily for debugging) now includes guess argument
  sets that were skipped entirely e.g. for not supporting user-supplied kwargs.
  All guesses thus removed from ``filtered_guess_kwargs`` are now listed as
  "Disabled" at the beginning of the trace. [#5578]

- Emit a warning when reading an ECSV file without specifying the ``format``
  and without PyYAML installed.  Previously this silently fell through to
  parsing as a basic format file and the file metadata was lost. [#7580]

- Optionally allow writing masked columns to ECSV with the mask explicitly
  specified as a separate column instead of marking masked elements with ""
  (empty string).  This allows handling the case of a masked string column
  with "" data rows.  [#7481]

astropy.io.misc
^^^^^^^^^^^^^^^

- Added support for saving all representation classes and many coordinate
  frames to the asdf format. [#7079]

- Added support for saving models with units to the asdf format. [#7237]

- Added a new ``character_as_bytes`` keyword to the HDF5 Table reading
  function to control whether byte string columns in the HDF5 file
  are left as bytes or converted to unicode.  The default is to read
  as bytes (``character_as_bytes=True``). [#7024, #8017]

astropy.io.fits
^^^^^^^^^^^^^^^

- ``HDUList.pop()`` now accepts string and tuple extension name
  specifications. [#7236]

- Add an ``ignore_hdus`` keyword to ``FITSDiff`` to allow ignoring HDUs by
  NAME when diffing two FITS files [#7538]

- Optionally allow writing masked columns to FITS with the mask explicitly
  specified as a separate column instead of using the FITS standard of
  certain embedded null values (``NaN`` for float, ``TNULL`` for integers).
  This can be used to work around limitations in the FITS standard. [#7481]

- All time coordinates can now be written to and read from FITS binary tables,
  including those with vectorized locations. [#7430]

- The ``fitsheader`` command line tool now supports a ``dfits+fitsort`` mode,
  and the dotted notation for keywords (e.g. ``ESO.INS.ID``). [#7240]

- Fall back to reading arrays using mode='denywrite' if mode='readonly' fails
  when using memory-mapping. This solves cases on some platforms when the
  available address space was less than the file size (even when using memory
  mapping). [#7926]

astropy.modeling
^^^^^^^^^^^^^^^^

- Add a ``Multiply`` model which preserves unit through evaluate, unlike
  ``Scale`` which is dimensionless. [#7210]

- Add a ``uses_quantity`` property to ``Model`` which allows introspection of if
  the ``Model`` can accept ``Quantity`` objects. [#7417]

- Add a ``separability_matrix`` function which returns the correlation matrix
  of inputs and outputs. [#7803]

- Fixed compatibility of ``JointFitter`` with the latest version of Numpy. [#7984]

- Add ``prior`` and ``posterior`` constraints to modeling parameters. These are
  not used by any current fitters, but are provided to allow user code to
  experiment with Bayesian fitters.  [#7558]

astropy.nddata
^^^^^^^^^^^^^^

- ``NDUncertainty`` objects now have a ``quantity`` attribute for simple
  conversion to quantities. [#7704]

- Add a ``bitmask`` module that provides functions for manipulating bitmasks
  and data quality (DQ) arrays. [#7944]

astropy.stats
^^^^^^^^^^^^^

- Add an ``astropy.stats.bls`` module with an implementation of the "box least
  squares" periodogram that is commonly used for discovering transiting
  exoplanets and eclipsing binaries. [#7391]

astropy.table
^^^^^^^^^^^^^

- Added support for full use of ``Time`` mixin column for join, hstack, and
  vstack table operations. [#6888]

- Added a new table index engine, ``SCEngine``, based on the Sorted Containers
  package. [#7574]

- Add a new keyword argument ``serialize_method`` to ``Table.write`` to
  control how ``Time`` and ``MaskedColumn`` columns are written. [#7481]

- Allow mixin columns to be used in table ``group`` and ``unique``
  functions. This applies to both the key columns and the other data
  columns. [#7712]

- Added support for stacking ``Column``, mixin column (e.g. ``Quantity``,
  ``Time``) or column-like objects. [#7674]

- Added support for inserting a row into a Table that has ``Time`` or
  ``TimeDelta`` column(s). [#7897]

astropy.tests
^^^^^^^^^^^^^

- Added an option ``--readonly`` to the test command to change the
  permissions on the temporary installation location to read-only. [#7598]

astropy.time
^^^^^^^^^^^^

- Allow array-valued ``Time`` object to be modified in place. [#6028]

- Added support for missing values (masking) to the ``Time`` class. [#6028]

- Added supper for a 'local' time scale (for free-running clocks, etc.),
  and round-tripping to the corresponding FITS time scale. [#7122]

- Added `datetime.timedelta` format class for ``TimeDelta``. [#7441]

- Added ``strftime`` and ``strptime`` methods to ``Time`` class.
  These methods are similar to those in the Python standard library
  `time` package and provide flexible input and output formatting. [#7323]

- Added ``datetime64`` format to the ``Time`` class to support working with
  ``numpy.datetime64`` dtype arrays. [#7361]

- Add fractional second support for ``strftime`` and ``strptime`` methods
  of ``Time`` class. [#7705]

- Added an ``insert`` method to allow inserting one or more values into a
  ``Time`` or ``TimeDelta`` object. [#7897]

- Remove timescale from string version of FITS format time string.
  The timescale is not part of the FITS standard and should not be included.
  This change may cause some compatibility issues for code that relies on
  round-tripping a FITS format string with a timescale. Strings generated
  from previous versions of this package are still understood but a
  DeprecationWarning will be issued. [#7870]

astropy.uncertainty
^^^^^^^^^^^^^^^^^^^

- This sub-package was added as a "preview" (i.e. API unstable), containing
  the ``Distribution`` class and associated convenience functions. [#6945]

astropy.units
^^^^^^^^^^^^^

- Add complex numbers support for ``Quantity._repr_latex_``. [#7676]

- Add ``thermodynamic_temperature`` equivalency to convert between
  Jy/sr and "thermodynamic temperature" for cosmology. [#7054]

- Add millibar unit. [#7863]

- Add maggy and nanomaggy unit, as well as associated ``zero_point_flux``
  equivalency. [#7891]

- ``AB`` and ``ST`` are now enabled by default, and have alternate names
  ``ABflux`` and ``STflux``. [#7891]

- Added ``littleh`` unit and associated ``with_H0`` equivalency. [#7970]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Added ``imshow_norm`` function, which combines imshow and creation of a
  ``ImageNormalize`` object. [#7785]

astropy.visualization.wcsaxes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Add support for setting ``set_separator(None)`` in WCSAxes to use default
  separators. [#7570]

- Added two keyword argument options to ``CoordinateHelper.set_format_unit``:
  ``decimal`` can be used to specify whether to use decimal formatting for the
  labels (by default this is False for degrees and hours and True otherwise),
  and ``show_decimal_unit`` can be used to determine whether the units should be
  shown for decimal labels. [#7318]

- Added documentation for ``transform=`` and ``coord_meta=``. [#7698]

- Allow ``coord_meta=`` to optionally include ``format_unit=``. [#7848]

- Add support for more rcParams related to the grid, ticks, and labels, and
  should work with most built-in Matplotlib styles. [#7961]

- Improved rendering of outward-facing ticks. [#7961]

- Add support for ``tick_params`` (which is a standard Matplotlib
  function/method) on both the ``WCSAxes`` class and the individual
  ``CoordinateHelper`` classes. Note that this is provided for compatibility
  with Matplotlib syntax users may be familiar with, but it is not the
  preferred way to change settings. Instead, methods such as ``set_ticks``
  should be preferred. [#7969]

- Moved the argument ``exclude_overlapping`` from ``set_ticks`` to
  ``set_ticklabel``. [#7969]

- Added a ``pad=`` argument to ``set_ticklabel`` to provide a way to control
  the padding between ticks and tick labels. [#7969]

- Added support for setting the tick direction in ``set_ticks`` using the
  ``direction=`` keyword argument. [#7969]

astropy.wcs
^^^^^^^^^^^

- Map ITRS frames to terrestrial WCS coordinates. This will make it possible to
  use WCSAxes to make figures that combine both celestial and terrestrial
  features. An example is plotting the coordinates of an astronomical transient
  over an all- sky satellite image to illustrate the position relative to the
  Earth at the time of the event. The ITRS frame is identified with WCSs that
  use the ``TLON-`` and ``TLAT-`` coordinate types. There are several examples
  of WCSs where this syntax is used to describe terrestrial coordinate systems:
  Section 7.4.1 of `WCS in FITS "Paper II" <https://ui.adsabs.harvard.edu/abs/2002A%26A...395.1077C>`_
  and the `WCSTools documentation <http://tdc-www.harvard.edu/software/wcstools/wcstools.multiwcs.html>`_.
  [#6990]

- Added the abstract base class for the low-level WCS API described in APE 14
  (https://doi.org/10.5281/zenodo.1188875). [#7325]

- Add ``WCS.footprint_contains()`` function to check if the WCS footprint contains a given sky coordinate. [#7273]

- Added the abstract base class for the high-level WCS API described in APE 14
  (https://doi.org/10.5281/zenodo.1188875). [#7325]

- Added the high-level wrapper class for low-level WCS objects as described in
  APE 14 (https://doi.org/10.5281/zenodo.1188875). [#7326]

- Added a new property ``WCS.has_distortion``. [#7326]

- Deprecated ``_naxis1`` and ``_naxis2`` in favor of ``pixel_shape``. [#7973]

- Added compatibility to wcslib version 6. [#8093]


API Changes
-----------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- ``kernel`` can now be a tuple. [#7561]

- Not technically an API changes, however, the doc string indicated that ``boundary=None``
  was the default when actually it is ``boundary='fill'``. The doc string has been corrected,
  however, someone may interpret this as an API change not realising that nothing has actually
  changed. [#7293]

- ``interpolate_replace_nans()`` can no longer accept the keyword argument
  ``preserve_nan``. It is explicitly set to ``False``. [#8088]


astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Fixed ``astropy.coordinates.concatenate`` to include velocity data in the
  concatenation. [#7922]

- Changed the name of the single argument to ``Frame.realize_frame()`` from the
  (incorrect) ``representation_type`` to ``data``. [#7923]

- Negative parallaxes passed to ``Distance()`` now raise an error by default
  (``allow_negative=False``), or are converted to NaN values with a warning
  (``allow_negative=True``). [#7988]

- Negating a ``SphericalRepresentation`` object now changes the angular
  coordinates (by rotating 180º) instead of negating the distance. [#7988]

- Creation of new frames now generally creates copies of frame attributes,
  rather than inconsistently either copying or making references. [#8204]

- The frame class method ``is_equivalent_frame`` now checks for equality of
  components to determine if a frame is the same when it has frame attributes
  that are representations, rather than checking if they are the same
  object. [#8218]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- If a fast reader is explicitly selected (e.g. ``fast_reader='force'``) and
  options which are incompatible with the fast reader are provided
  (e.g. ``quotechar='##'``) then now a ``ParameterError`` exception will be
  raised. [#5578]

- The fast readers will now raise ``InconsistentTableError`` instead of
  ``CParserError`` if the number of data and header columns do not match.
  [#5578]

- Changed a number of ``ValueError`` exceptions to ``InconsistentTableError``
  in places where the exception is related to parsing a table which is
  inconsistent with the specified table format.  Note that
  ``InconsistentTableError`` inherits from ``ValueError`` so no user code
  changes are required. [#7425]

astropy.io.fits
^^^^^^^^^^^^^^^

- The ``fits.table_to_hdu()`` function will translate any column ``format``
  attributes to a TDISPn format string, if possible, and store it as a TDISPn
  keyword in the ``HDU`` header. [#7226]

astropy.modeling
^^^^^^^^^^^^^^^^

- Change the order of the return values from ``FittingWithOutlierRemoval``,
  such that ``fitted_model`` comes first, for consistency with other fitters.
  For the second value, return only a boolean outlier ``mask``, instead of the
  previous ``MaskedArray`` (which included a copy of the input data that was
  both redundant and inadvertently corrupted at masked points). Return a
  consistent type for the second value when ``niter=0``. [#7407]

- Set the minimum value for the ``bolometric_flux`` parameter of the
  ``BlackBody1D`` model to zero. [#7045]

astropy.nddata
^^^^^^^^^^^^^^

- Add two new uncertainty classes, ``astropy.nddata.VarianceUncertainty`` and
  ``astropy.nddata.InverseVariance``. [#6971]

astropy.stats
^^^^^^^^^^^^^

- String values can now be used for the ``cenfunc`` and ``stdfunc``
  keywords in the ``SigmaClip`` class and ``sigma_clip`` and
  ``sigma_clipped_stats`` functions. [#7478]

- The ``SigmaClip`` class and ``sigma_clip`` and
  ``sigma_clipped_stats`` functions now have a ``masked`` keyword,
  which can be used to return either a masked array (default) or an
  ndarray with the min/max values. [#7478]

- The ``iters`` keyword has been renamed (and deprecated) to
  ``maxiters`` in the ``SigmaClip`` class and ``sigma_clip`` and
  ``sigma_clipped_stats`` functions. [#7478]

astropy.table
^^^^^^^^^^^^^

- ``Table.read()`` on a FITS binary table file will convert any TDISPn header
  keywords to a Python formatting string when possible, and store it in the
  column ``format`` attribute. [#7226]

- No values provided to stack will now raise ``ValueError`` rather than
  ``TypeError``. [#7674]

astropy.tests
^^^^^^^^^^^^^

- ``from astropy.tests.helper import *`` no longer includes
  ``quantity_allclose``. However,
  ``from astropy.tests.helper import quantity_allclose`` would still work.
  [#7381]

- ``warnings_to_ignore_by_pyver`` option in
  ``enable_deprecations_as_exceptions()`` now takes ``None`` as key.
  Any deprecation message that is mapped to ``None`` will be ignored
  regardless of the Python version. [#7790]

astropy.time
^^^^^^^^^^^^

- Added the ability to use ``local`` as time scale in ``Time`` and
  ``TimeDelta``. [#6487]

- Comparisons, addition, and subtraction of ``Time`` instances with non-time
  instances will now return ``NotImplemented`` rather than raise the
  ``Time``-specific ``OperandTypeError``.  This will generally lead to a
  regular ``TypeError``.  As a result, ``OperandTypeError`` now only occurs if
  the operation is between ``Time`` instances of incompatible type or scale.
  [#7584]

astropy.units
^^^^^^^^^^^^^

- In ``UnitBase.compose()``, if a sequence (list|tuple) is passed in to
  ``units``, the default for ``include_prefix_units`` is set to
  `True`, so that no units get ignored. [#6957]

- Negative parallaxes are now converted to NaN values when using the
  ``parallax`` equivalency. [#7988]

astropy.utils
^^^^^^^^^^^^^

- ``InheritDocstrings`` now also works on class properties. [#7166]

- ``diff_values()``, ``report_diff_values()``, and ``where_not_allclose()``
  utility functions are moved from ``astropy.io.fits.diff``. [#7444]

- ``invalidate_caches()`` has been removed from the
  ``astropy.utils.compat`` namespace, use it directly from ``importlib``. [#7872]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- In ``ImageNormalize``, the default for ``clip`` is set to ``True``. [#7800]

- Changed ``AsymmetricPercentileInterval`` and ``MinMaxInterval`` to
  ignore NaN values in arrays. [#7360]

- Automatically default to using ``grid_type='contours'`` in WCSAxes when using
  a custom ``Transform`` object if the transform has no inverse. [#7847]


Performance Improvements
------------------------

- Reduced import time by more cautious use of the standard library. [#7647]

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Major performance overhaul to ``convolve()``. [#7293]

- ``convolve()``: Boundaries ``fill``, ``extend``, and ``wrap`` now use a single
  implementation that pads the image with the correct boundary values before convolving.
  The runtimes of these three were significantly skewed. They now have
  equivalent runtimes that are also faster than before due to performant contiguous
  memory access. However, this does increase the memory footprint as an entire
  new image array is required plus that needed for the padded region.[#7293]

- ``convolve()``: Core computation ported from Cython to C. Several optimization
  techniques have been implemented to achieve performance gains, e.g. compiler
  hoisting, and vectorization, etc. Compiler optimization level ``-O2`` required for
  hoisting and ``-O3`` for vectorization. [#7293]

- ``convolve()``: ``nan_treatment=‘interpolate’`` was slow to compute irrespective of
  whether any NaN values exist within the array. The input array is now
  checked for NaN values and interpolation is disabled if non are found. This is a
  significant performance boost for arrays without NaN values. [#7293]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Sped up creating SkyCoord objects by a factor of ~2 in some cases. [#7615]

- Sped up getting xyz vectors from ``CartesianRepresentation`` (which
  is used a lot internally). [#7638]

- Sped up transformations and some representation methods by replacing
  python code with (compiled) ``erfa`` ufuncs. [#7639]

- Sped up adding differential (velocity) data to representations by a factor of
  ~20, which improves the speed of frame and SkyCoord initialization. [#7924]

- Refactored ``SkyCoord`` initializer to improve performance and code clarity.
  [#7958]

- Sped up initialization of ``Longitude`` by ~40%. [#7616]

astropy.stats
^^^^^^^^^^^^^

- The ``SigmaClip`` class and ``sigma_clip`` and
  ``sigma_clipped_stats`` functions are now significantly faster. [#7478]

- A Cython implementation for `astropy.stats.kuiper_two` and a vectorized
  implementation for `astropy.stats.kuiper_false_positive_probability` have
  been added, speeding up both functions.  [#8104]

astropy.units
^^^^^^^^^^^^^

- Sped up creating new composite units, and raising units to some power
  [#7549, #7649]

- Sped up Unit.to when target unit is the same as the original unit. [#7643]

- Lazy-load ``scipy.special`` to shorten ``astropy.units`` import time. [#7636]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Significantly sped up drawing of contours in WCSAxes. [#7568]


Bug Fixes
---------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Fixed bug in ``convolve_fft`` where masked input was copied with
  ``numpy.asarray`` instead of ``numpy.asanyarray``.
  ``numpy.asarray`` removes the mask subclass causing
  ``numpy.ma.ismasked(input)`` to fail, causing ``convolve_fft``
  to ignore all masked input. [#8137]

- Remove function side-effects of input data from ``convolve_fft``.
  It was possible for input data to remain modified if particular exceptions
  were raised. [#8152]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- ``EarthLocation.of_address`` now uses the OpenStreetMap geocoding API by
  default to retrieve coordinates, with the Google API (which now requires an
  API key) as an option. [#7918]

- Fixed a bug that caused frame objects with NaN distances to have NaN sky
  positions, even if valid sky coordinates were specified. [#7988]

- Fixed ``represent_as()`` to not round-trip through cartesian if the same
  representation class as the instance is passed in. [#7988]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fixed a problem when ``guess=True`` that ``fast_reader`` options
  could be dropped after the first fast reader class was tried. [#5578]

- Units in CDS-formatted tables are now parsed correctly by the units
  module. [#7348]

astropy.io.misc
^^^^^^^^^^^^^^^

- Fixed bug when writing a table with masked columns to HDF5. Previously
  the mask was being silently dropped.  If the ``serialize_meta`` option is
  enabled the data mask will now be written as an additional column and the
  masked columns will round-trip correctly. [#7481]

- Fixed a bug where writing to HDF5 failed for for tables with columns of
  unicode strings.  Now those columns are first encoded to UTF-8 and
  written as byte strings. [#7024, #8017]

- Fixed a bug with serializing the bounding_box of models initialized
  with ``Quantities`` . [#8052]

astropy.io.fits
^^^^^^^^^^^^^^^

- Added support for ``copy.copy`` and ``copy.deepcopy`` for ``HDUList``. [#7218]

- Override ``HDUList.copy()`` to return a shallow HDUList instance. [#7218]

astropy.modeling
^^^^^^^^^^^^^^^^

- Fix behaviour of certain models with units, by making certain unit-related
  attributes readonly. [#7210]

- Fixed an issue with validating a ``bounding_box`` whose items are
  ``Quantities``. [#8052]

- Fix ``Moffat1D`` and ``Moffat2D`` derivatives. [#8108]

astropy.nddata
^^^^^^^^^^^^^^

- Fixed rounding behavior in ``overlap_slices`` for even-sized small
  arrays. [#7859]

- Added support for pickling ``NDData`` instances that have an uncertainty.
  [#7383]

astropy.stats
^^^^^^^^^^^^^

- Fix errors in ``kuiper_false_positive_probability``. [#7975]

astropy.tests
^^^^^^^^^^^^^

- Fixing bug that prevented to run the doctests on only a single rst documentation
  file rather than all of them. [#8055]

astropy.time
^^^^^^^^^^^^

- Fix a bug when setting a ``TimeDelta`` array item with plain float value(s).
  This was always interpreted as a JD (day) value regardless of the
  ``TimeDelta`` format. [#7990]

astropy.units
^^^^^^^^^^^^^

- To simplify fast creation of ``Quantity`` instances from arrays, one can now
  write ``array << unit`` (equivalent to ``Quantity(array, unit, copy=False)``).
  If ``array`` is already a ``Quantity``, this will convert the quantity to the
  requested units; in-place conversion can be done with ``quantity <<= unit``.
  [#7734]

astropy.utils
^^^^^^^^^^^^^

- Fixed a bug due to which ``report_diff_values()`` was reporting incorrect
  number of differences when comparing two ``numpy.ndarray``. [#7470]

- The download progress bar is now only displayed in terminals, to avoid
  polluting piped output. [#7577]

- Ignore URL mirror caching when there is no internet. [#8163]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Right ascension coordinates are now shown in hours by default, and the
  ``set_format_unit`` method on ``CoordinateHelper`` now works correctly
  with angle coordinates. [#7215]


Other Changes and Additions
---------------------------

- The documentation build now uses the Sphinx configuration from sphinx-astropy
  rather than from astropy-helpers. [#7139]

- Versions of Numpy <1.13 are no longer supported. [#7058]

- Running tests now suppresses the output of the installation stage by default,
  to allow easier viewing of the test results. To re-enable the output as
  before, use ``python setup.py test --verbose-install``. [#7512]

- The ERFA functions are now wrapped in ufuncs instead of custom C code,
  leading to some speed improvements, and setting the stage for allowing
  overrides with ``__array_ufunc__``. [#7502]

- Updated the bundled CFITSIO library to 3.450. See
  ``cextern/cfitsio/docs/changes.txt`` for additional information. [#8014]

- The ``representation`` keywords in coordinate frames are now deprecated in
  favor of the ``representation_type`` keywords (which are less
  ambiguously named). [#8119]



Version 3.0.5 (2018-10-14)
==========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Fixed bug in which consecutive ``StaticMatrixTransform``'s in a frame
  transform path would be combined in the incorrect order. [#7707]

astropy.tests
^^^^^^^^^^^^^

- Fixing bug that doctests were not picked up from the narrative
  documentation when tests were run for all modules. [#7767]



Version 3.0.4 (2018-08-02)
==========================

API Changes
-----------

astropy.table
^^^^^^^^^^^^^

- The private ``_parent`` attribute in the ``info`` attribute of table
  columns was changed from a direct reference to the parent column to a weak
  reference.  This was in response to a memory leak caused by having a
  circular reference cycle.  This change means that expressions like
  ``col[3:5].info`` will now fail because at the point of the ``info``
  property being evaluated the ``col[3:5]`` weak reference is dead.  Instead
  force a reference with ``c = col[3:5]`` followed by
  ``c.info.indices``. [#6277, #7448]


Bug Fixes
---------

astropy.nddata
^^^^^^^^^^^^^^

- Fixed an bug when creating the ``WCS`` of a cutout (see ``nddata.Cutout2D``)
  when input image's ``WCS`` contains ``SIP`` distortion corrections by
  adjusting the ``crpix`` of the ``astropy.wcs.Sip`` (in addition to
  adjusting the ``crpix`` of the ``astropy.wcs.WCS`` object). This bug
  had the potential to produce large errors in ``WCS`` coordinate
  transformations depending on the position of the cutout relative
  to the input image's ``crpix``. [#7556, #7550]

astropy.table
^^^^^^^^^^^^^

- Fix memory leak where updating a table column or deleting a table
  object was not releasing the memory due to a reference cycle
  in the column ``info`` attributes. [#6277, #7448]

astropy.wcs
^^^^^^^^^^^

- Fixed an bug when creating the ``WCS`` slice (see ``WCS.slice()``)
  when ``WCS`` contains ``SIP`` distortion corrections by
  adjusting the ``WCS.sip.crpix`` in addition to adjusting
  ``WCS.wcs.crpix``. This bug had the potential to produce large errors in
  ``WCS`` coordinate transformations depending on the position of the slice
  relative to ``WCS.wcs.crpix``. [#7556, #7550]


Other Changes and Additions
---------------------------

- Updated bundled wcslib to v 5.19.1 [#7688]


Version 3.0.3 (2018-06-01)
==========================

Bug Fixes
---------

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fix stripping correct (header) comment line from ``meta['comments']``
  in the ``CommentedHeader`` reader for all ``header_start`` settings. [#7508]

astropy.io.fits
^^^^^^^^^^^^^^^

- Raise error when attempting to open gzipped FITS file in 'append' mode.
  [#7473]

- Fix a bug when writing to FITS a table that has a column description
  with embedded blank lines. [#7482]

astropy.tests
^^^^^^^^^^^^^

- Enabling running tests for multiple packages when specified comma
  separated. [#7463]


Version 3.0.2 (2018-04-23)
==========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Computing a 3D separation between two ``SkyCoord`` objects (with the
  ``separation_3d`` method) now works with or without velocity data attached to
  the objects. [#7387]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Fix validate with xmllint=True. [#7255, #7283]

astropy.modeling
^^^^^^^^^^^^^^^^

- ``FittingWithOutlierRemoval`` now handles model sets, as long as the
  underlying fitter supports masked values. [#7199]

- Remove assumption that ``model_set_axis == 0`` for 2D models in
  ``LinearLSQFitter``. [#7317, #7199]

- Fix the shape of the outputs when a model set is evaluated with
  ``model_set_axis=False`` . [#7317]

astropy.stats
^^^^^^^^^^^^^

- Accept a tuple for the ``axis`` parameter in ``sigma_clip``, like the
  underlying ``numpy`` functions and some other functions in ``stats``. [#7199]

astropy.tests
^^^^^^^^^^^^^

- The function ``quantity_allclose`` was moved to the ``units`` package with
  the new, shorter name ``allclose``. This eliminates a runtime dependency on
  ``pytest`` which was causing issues for some affiliated packages. The old
  import will continue to work but may be deprecated in the future. [#7252]

astropy.units
^^^^^^^^^^^^^

- Added a units-aware ``allclose`` function (this was previously available in
  the ``tests`` module as ``quantity_allclose``). To complement ``allclose``,
  a new ``isclose`` function is also added and backported. [#7252]


Version 3.0.1 (2018-03-12)
==========================

Bug Fixes
---------

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fix a unicode decode error when reading a table with non-ASCII characters.
  The fast C reader cannot handle unicode so the code now uses the pure-Python
  reader in this case. [#7103]

astropy.io.fits
^^^^^^^^^^^^^^^

- Updated the bundled CFITSIO library to 3.430. This is to remedy a critical
  security vulnerability that was identified by NASA. See
  ``cextern/cfitsio/docs/changes.txt`` for additional information. [#7274]

astropy.io.misc
^^^^^^^^^^^^^^^

- Make sure that a sufficiently recent version of ASDF is installed when
  running test suite against ASDF tags and schemas. [#7205]

astropy.io.registry
^^^^^^^^^^^^^^^^^^^

- Fix reading files with serialized metadata when using a Table subclass. [#7213]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Fix lookup fields by ID. [#7208]

astropy.modeling
^^^^^^^^^^^^^^^^

- Fix model set evaluation over common input when model_set_axis > 0. [#7222]

- Fixed the evaluation of compound models with units. This required adding the
  ability to have ``input_units_strict`` and ``input_units_allow_dimensionless``
  be dictionaries with input names as keys. [#6952]

astropy.units
^^^^^^^^^^^^^

- ``quantity_helper`` no longer requires ``scipy>=0.18``. [#7219]


Version 3.0 (2018-02-12)
========================

New Features
------------

astropy.constants
^^^^^^^^^^^^^^^^^

- New context manager ``set_enabled_constants`` to temporarily use an older
  version. [#7008]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- The ``Distance`` object now accepts ``parallax`` as a keyword in the
  initializer, and supports retrieving a parallax (as an ``Angle``) via
  the ``.parallax`` attributes. [#6855]

- The coordinate frame classes (subclasses of ``BaseCoordinateFrame``) now
  always have ``.velocity``, ``.proper_motion``, and ``.radial_velocity``
  properties that provide shorthands to the full-space Cartesian velocity as
  a ``CartesianDifferential``, the 2D proper motion as a ``Quantity``, and the
  radial or line-of-sight velocity as a ``Quantity``. [#6869]

- ``SkyCoord`` objects now support storing and transforming differentials - i.e.,
  both radial velocities and proper motions. [#6944]

- All frame classes now automatically get sensible representation mappings for
  velocity components. For example, ``d_x``, ``d_y``, ``d_z`` are all
  automatically mapped to frame component namse ``v_x``, ``v_y``, ``v_z``.
  [#6856]

- ``SkyCoord`` objects now support updating the position of a source given its
  space motion and a new time or time difference. [#6872]

- The frame classes now accept a representation class or differential class, or
  string names for either, through the keyword arguments ``representation_type``
  and ``differential_type`` instead of ``representation`` and
  ``differential_cls``. [#6873]

- The frame classes (and ``SkyCoord``) now give more useful error messages when
  incorrect attribute names are given.  Instead of using the representation
  attribute names, they use the frame attribute names. [#7106]

- ``EarthLocation`` now has a method to compute the  gravitational redshift due
  due to solar system bodies.  [#6861, #6935]

- ``EarthLocation`` now has a ``get_gcrs`` convenience method to get the
  location in GCRS coordinates.  [#6861, #6935]

astropy.io.fits
^^^^^^^^^^^^^^^

- Expanded the FITS ``Column`` interface to accept attributes pertaining to the FITS
  World Coordinate System, which includes spatial(celestial) and time coordinates. [#6359]

- Added ``ver`` attribute to set the ``EXTVER`` header keyword to ``ImageHDU``
  and ``TableHDU``. [#6454]

- The performance for reading FITS tables has been significantly improved,
  in particular for cases where the tables contain one or more string columns
  and when done through ``Table.read``. [#6821]

- The performance for writing tables from ``Table.write`` has now been
  significantly improved for tables containing one or more string columns. [#6920]

- The ``Table.read`` now supports a ``memmap=`` keyword argument to control
  whether or not to use  memory mapping when reading the table. [#6821]

- When reading FITS tables with ``fits.open``, a new keyword argument
  ``character_as_bytes`` can be passed - when set to `True`, character columns
  are returned as Numpy byte arrays (Numpy type S) while when set to `False`,
  the same columns are decoded to Unicode strings (Numpy type U) which uses more
  memory. [#6821]

- The ``table_to_hdu`` function and the ``BinTableHDU.from_columns`` and
  ``FITS_rec.from_columns`` methods now include a ``character_as_bytes``
  keyword argument - if set to `True`, then when string columns are accessed,
  byte columns will be returned, which can provide significantly improved
  performance. [#6920]

- Added support for writing and reading back a table which has "mixin columns"
  such as ``SkyCoord`` or ``EarthLocation`` with no loss of information. [#6912]

- Enable tab-completion for ``FITS_rec`` column names and ``Header`` keywords
  with IPython 5 and later. [#7071]

astropy.io.misc
^^^^^^^^^^^^^^^

- When writing to HDF5 files, the serialized metadata are now saved in a new
  dataset, instead of the HDF5 dataset attributes. This allows for metadata of
  any dimensions. [#6304]

- Added support in HDF5 for writing and reading back a table which has "mixin
  columns" such as ``SkyCoord`` or ``EarthLocation`` with no loss of
  information. [#7007]

- Add implementations of astropy-specific ASDF tag types. [#6790]

- Add ASDF tag and schema for ICRSCoord. [#6904]

astropy.modeling
^^^^^^^^^^^^^^^^

- Add unit support for tabular models. [#6529]

- A ``deepcopy()`` method was added to models. [#6515]

- Added units support to ``AffineTransformation``. [#6853]

- Added ``is_separable`` function to modeling to test the
  separability of a model. [#6746]

- Added ``Model.separable`` property. It returns a boolean value or
  ``None`` if not set. [#6746]

- Support masked array values in ``LinearLSQFitter`` (instead of silently
  ignoring the mask). [#6927]

astropy.stats
^^^^^^^^^^^^^

- Added false alarm probability computation to ``astropy.stats.LombScargle``
  [#6488]

- Implemented Kuiper functions in ``astropy.stats`` [#3724, #6565]

astropy.table
^^^^^^^^^^^^^

- Added support for reading and writing ``astropy.time.Time`` Table columns
  to and from FITS tables, to the extent supported by the FITS standard. [#6176]

- Improved exception handling and error messages when column ``format``
  attribute is incorrect for the column type. [#6385]

- Allow to pass ``htmldict`` option to the jsviewer writer. [#6551]

- Added new table operation ``astropy.table.setdiff`` that returns the set
  difference of table rows for two tables. [#6443]

- Added support for reading time columns in FITS compliant binary tables
  as ``astropy.time.Time`` Table columns. [#6442]

- Allowed to remove table rows through the ``__delitem__`` method. [#5839]

- Added a new ``showtable`` command-line script to view binary or ASCII table
  files. [#6859]

- Added new table property ``astropy.table.Table.loc_indices`` that returns the
  location of rows by indexes. [#6831]

- Allow updating of table by indices through the property ``astropy.table.Table.loc``. [#6831]

- Enable tab-completion for column names with IPython 5 and later. [#7071]

- Allow getting and setting a table Row using multiple column names. [#7107]

astropy.tests
^^^^^^^^^^^^^

- Split pytest plugins into separate modules. Move remotedata, openfiles,
  doctestplus plugins to standalone repositories. [#6384, #6606]

- When testing, astropy (or the package being tested) is now installed to
  a temporary directory instead of copying the build. This allows
  entry points to work correctly. [#6890]

- The tests_require setting in setup.py now works properly when running
  'python setup.py test'. [#6892]

astropy.units
^^^^^^^^^^^^^

- Deprecated conversion of quantities to truth values. Currently, the expression
  ``bool(0 * u.dimensionless_unscaled)`` evaluates to ``True``. In the future,
  attempting to convert a ``Quantity`` to a ``bool`` will raise ``ValueError``.
  [#6580, #6590]

- Modify the ``brightness_temperature`` equivalency to provide a surface
  brightness equivalency instead of the awkward assumed-per-beam equivalency
  that previously existed [#5173, #6663]

- Support was added for a number of ``scipy.special`` functions. [#6852]

astropy.utils
^^^^^^^^^^^^^

- The ``astropy.utils.console.ProgressBar.map`` class method now supports the
  ``ipython_widget`` option. You can now pass it both ``multiprocess=True`` and
  ``ipython_widget=True`` to get both multiprocess speedup and a progress bar
  widget in an IPython Notebook. [#6368]

- The ``astropy.utils.compat.funcsigs`` module has now been deprecated. Use the
  Python 'inspect' module directly instead. [#6598]

- The ``astropy.utils.compat.futures`` module has now been deprecated. Use the
  Python 'concurrent.futures' module directly instead. [#6598]

- ``JsonCustomEncoder`` is expanded to handle ``Quantity`` and ``UnitBase``.
  [#5471]

- Added a ``dcip_xy`` method to IERS that interpolates along the dX_2000A and
  dY_2000A columns of the IERS table.  Hence, the data for the CIP offsets is
  now available for use in coordinate frame conversion. [#5837]

- The functions ``matmul``, ``broadcast_arrays``, ``broadcast_to`` of the
  ``astropy.utils.compat.numpy`` module have been deprecated. Use the
  NumPy functions directly. [#6691]

- The ``astropy.utils.console.ProgressBar.map`` class method now returns
  results in sequential order. Previously, if you set ``multiprocess=True``,
  then the results could arrive in any arbitrary order, which could be a nasty
  shock. Although the function will still be evaluated on the items in
  arbitrary order, the return values will arrive in the same order in which the
  input items were provided. The method is now a thin wrapper around
  ``astropy.utils.console.ProgressBar.map_unordered``, which preserves the old
  behavior. [#6439]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Enable Matplotlib's subtraction shorthand syntax for composing and
  inverting transformations for the ``WCSWorld2PixelTransform`` and
  ``WCSPixel2WorldTransform`` classes by setting ``has_inverse`` to ``True``.
  In order to implement a unit test, also implement the equality comparison
  operator for both classes. [#6531]

- Added automatic hiding of axes labels when no tick labels are drawn on that
  axis. This parameter can be configured with
  ``WCSAxes.coords[*].set_axislabel_visibility_rule`` so that labels are automatically
  hidden when no ticks are drawn or always shown. [#6774]

astropy.wcs
^^^^^^^^^^^

- Added a new function ``celestial_frame_to_wcs`` to convert from
  coordinate frames to WCS (the opposite of what ``wcs_to_celestial_frame``
  currently does. [#6481]

- ``wcslib`` was updated to v 5.18. [#7066]


API Changes
-----------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- ``Gaussian2DKernel`` now accepts ``x_stddev`` in place of ``stddev`` with
  an option for ``y_stddev``, if different. It also accepts ``theta`` like
  ``Gaussian2D`` model. [#3605, #6748]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Deprecated ``recommended_units`` for representations. These were used to
  ensure that any angle was presented in degrees in sky coordinates and
  frames. This is more logically done in the frame itself. [#6858]

- As noted above, the frame class attributes ``representation`` and
  ``differential_cls`` are being replaced by ``representation_type`` and
  ``differential_type``. In the next version, using ``representation`` will raise
  a deprecation warning. [#6873]

- Coordinate frame classes now can't be added to the frame transform graph if
  they have frame attribute names that conflict with any component names. This
  is so ``SkyCoord`` can uniquely identify and distinguish frame attributes from
  frame components. [#6871]

- Slicing and reshaping of ``SkyCoord`` and coordinate frames no longer passes
  the new object through ``__init__``, but directly sets attributes on a new
  instance. This speeds up those methods by an order of magnitude, but means
  that any customization done in ``__init__`` is by-passed. [#6941]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Allow ECSV files to be auto-identified by ``Table.read`` or
  ``Table.write`` based on the ``.ecsv`` file name suffix. In this case it
  is not required to provide the ``format`` keyword. [#6552]

astropy.io.fits
^^^^^^^^^^^^^^^

- Automatically detect and handle compression in FITS files that are opened by
  passing a file handle to ``fits.open`` [#6373]

- Remove the ``nonstandard`` checksum option. [#6571]

astropy.io.misc
^^^^^^^^^^^^^^^

- When writing to HDF5 files, the serialized metadata are now saved in a new
  dataset instead of the HDF5 dataset attributes. This allows for metadata of
  any dimensions. [#6304]

- Deprecated the ``usecPickle`` kwarg of ``fnunpickle`` and ``fnpickle`` as
  it was needed only for Python2 usage. [#6655]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Add handling of ``tree.Group`` elements to ``tree.Resource``.  Unified I/O
  or conversion to astropy tables is not affected. [#6262]

astropy.modeling
^^^^^^^^^^^^^^^^

- Removed deprecated ``GaussianAbsorption1D`` model.
  Use ``Const1D - Gaussian1D`` instead. [#6542]

- Removed the registry from modeling. [#6706]

astropy.table
^^^^^^^^^^^^^

- When setting the column ``format`` attribute the value is now immediately
  validated. Previously one could set to any value and it was only checked
  when actually formatting the column. [#6385]

- Deprecated the ``python3_only`` kwarg of the
  ``convert_bytestring_to_unicode`` and ``convert_unicode_to_bytestring``
  methods it was needed only for Python2 usage. [#6655]

- When reading in FITS tables with ``Table.read``, string columns are now
  represented using Numpy byte (dtype ``S``) arrays rather than Numpy
  unicode arrays (dtype ``U``). The ``Column`` class then ensures the
  bytes are automatically converted to string as needed. [#6821]

- When getting a table row using multiple column names, if one of the
  names is not a valid column name then a ``KeyError`` exception is
  now raised (previously ``ValueError``).  When setting a table row,
  if the right hand side is not a sequence with the correct length
  then a ``ValueError`` is now raised (previously in certain cases
  a ``TypeError`` was raised). [#7107]

astropy.utils
^^^^^^^^^^^^^

- ``download_files_in_parallel`` now always uses ``cache=True`` to make the
  function work on Windows. [#6671]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- The Astropy matplotlib plot style has been deprecated. It will continue to
  work in future but is no longer documented. [#6991]


Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Frame objects now use the default differential even if the representation is
  explicitly provided as long as the representation provided is the same type as
  the default representation. [#6944]

- Coordinate frame classes now raise an error when they are added to the frame
  transform graph if they have frame attribute names that conflict with any
  component names. [#6871]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Added support for reading very large tables in chunks to reduce memory
  usage. [#6458]

- Strip leading/trailing white-space from latex lines to avoid issues when
  matching ``\begin{tabular}`` statements.  This is done by introducing a new
  ``LatexInputter`` class to override the ``BaseInputter``. [#6311]

astropy.io.fits
^^^^^^^^^^^^^^^

- Properly handle opening of FITS files from ``http.client.HTTPResponse`` (i.e.
  it now works correctly when passing the results of ``urllib.request.urlopen``
  to ``fits.open``). [#6378]

- Fix the ``fitscheck`` script for updating invalid checksums, or removing
  checksums. [#6571]

- Fixed potential problems with the compression module [#6732]

- Always use the 'D' format for floating point values in ascii tables. [#6938]

astropy.table
^^^^^^^^^^^^^

- Fix getting a table row when using multiple column names (for example
  ``t[3]['a', 'b', 'c']``).  Also fix a problem when setting an entire row:
  if setting one of the right-hand side values failed this could result in
  a partial update of the referenced parent table before the exception is
  raised. [#7107]

astropy.time
^^^^^^^^^^^^

- Initialization of ``Time`` instances with bytes or arrays with dtype ``S``
  will now automatically attempt to decode as ASCII. This ensures ``Column``
  instances with ASCII strings stored with dtype ``S`` can be used.
  [#6823, #6903]

astropy.units
^^^^^^^^^^^^^

- Fixed a bug that caused PLY files to not be generated correctly in Python 3.
  [#7174]

astropy.utils
^^^^^^^^^^^^^

- The ``deprecated`` decorator applied to a class will now modify the class
  itself, rather than to create a class that just looks and behave like the
  original. This is needed so that the Python 3 ``super`` without arguments
  works for decorated classes. [#6615]

- Fixed ``HomogeneousList`` when setting one item or a slice. [#6773]

- Also check the type when creating a new instance of
  ``HomogeneousList``. [#6773]

- Make ``HomogeneousList`` work with iterators and generators when creating the
  instance, extending it, or using when setting a slice. [#6773]


Other Changes and Additions
---------------------------

- Versions of Python <3.5 are no longer supported. [#6556]

- Versions of Pytest <3.1 are no longer supported. [#6419]

- Versions of Numpy <1.10 are no longer supported. [#6593]

- The bundled CFITSIO was updated to version 3.41 [#6477]

- ``analytic_functions`` sub-package is removed.
  Use ``astropy.modeling.blackbody``. [#6541]

- ``astropy.vo`` sub-package is removed. Use ``astropy.samp`` for SAMP and
  ``astroquery`` for VO cone search. [#6540]

- The guide to setting up Emacs for code development was simplified, and
  updated to recommend ``flycheck`` and ``flake8`` for syntax checks. [#6692]

- The bundled version of PLY was updated to 3.10. [#7174]



Version 2.0.16 (2019-10-27)
===========================

Bug Fixes
---------

astropy.time
^^^^^^^^^^^^

- Fixed a troubling bug in which ``Time`` could loose precision, with deviations
  of 300 ns. [#9328]


Other Changes and Additions
---------------------------

- Updated IERS A URLs due to USNO prolonged maintenance. [#9443]



Version 2.0.15 (2019-10-06)
===========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Fixed a bug where the string representation of a ``BaseCoordinateFrame``
  object could become garbled under specific circumstances when the frame
  defines custom component names via ``RepresentationMapping``. [#8869]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fix uint conversion in ``FITS_rec`` when slicing a table. [#8982]

- Fix reading of unsigned 8-bit integer with compressed fits. [#9219]

astropy.nddata
^^^^^^^^^^^^^^

- Fixed a bug in ``overlap_slices`` where the ``"strict"`` mode was
  too strict for a small array along the upper edge of the large array.
  [#8901]

- Fixed a bug in ``overlap_slices`` where a ``NoOverlapError`` would
  be incorrectly raised for a 0-shaped small array at the origin.
  [#8901]

astropy.samp
^^^^^^^^^^^^

- Fixed a bug that caused an incorrectly constructed warning message
  to raise an error. [#8966]

astropy.table
^^^^^^^^^^^^^

- Fix ``FixedWidthNoHeader`` to pay attention to ``data_start`` keyword when
  finding first data line to split columns [#8485, #8511]

- Fix bug when initializing ``Table`` with ``rows`` as a generator. [#9315]

- Fix ``join`` when there are multiple mixin (Quantity) columns as keys. [#9313]

astropy.units
^^^^^^^^^^^^^

- ``Quantity`` now preserves the ``dtype`` for anything that is floating
  point, including ``float16``. [#8872]

- ``Unit()`` now accepts units with fractional exponents such as ``m(3/2)``
  in the default/``fits`` and ``vounit`` formats that would previously
  have been rejected for containing multiple solidi (``/``). [#9000]

- Fixed the LaTeX representation of units containing a superscript. [#9218]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Fixed compatibility issues with latest versions of Matplotlib. [#8961]


Other Changes and Additions
---------------------------

- Updated required version of Cython to v0.29.13 to make sure that
  generated C files are compatible with the upcoming Python 3.8 release
  as well as earlier supported versions of Python. [#9198]



Version 2.0.14 (2019-06-14)
===========================

Bug Fixes
---------

astropy.io.fits
^^^^^^^^^^^^^^^

- Fix ``Header.update`` which was dropping the comments when passed
  a ``Header`` object. [#8840]

astropy.modeling
^^^^^^^^^^^^^^^^

- ``Moffat1D.fwhm`` and ``Moffat2D.fwhm`` will return a positive value when
  ``gamma`` is negative. [#8801, #8815]

astropy.units
^^^^^^^^^^^^^

- Fixed a bug that prevented ``EarthLocation`` from being initialized with
  numpy >=1.17. [#8849]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Fixed ``quantity_support`` to work around the fact that matplotlib
  does not detect subclasses in its ``units`` framework. With this,
  ``Angle`` and other subclasses work correctly. [#8818]

- Fixed ``quantity_support`` to work properly if multiple context managers
  are nested. [#8844]


Version 2.0.13 (2019-06-08)
===========================

Bug Fixes
---------

astropy.io.fits
^^^^^^^^^^^^^^^
- Fixed bug in ``ColDefs._init_from_array()`` that caused unsigned datatypes
  with the opposite endianness as the host architecture to fail the
  TestColumnFunctions.test_coldefs_init_from_array unit test. [#8460]

astropy.io.misc
^^^^^^^^^^^^^^^

- Explicitly set PyYAML default flow style to None to ensure consistent
  astropy YAML output for PyYAML version 5.1 and later. [#8500]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Block floating-point columns from using repr format when converted to Table
  [#8358]

astropy.stats
^^^^^^^^^^^^^

- Fixed issue in ``bayesian_blocks`` when called with the ``ncp_prior``
  keyword. [#8339]

astropy.units
^^^^^^^^^^^^^

- Fix ``take`` when one gets only a single element from a ``Quantity``,
  ensuring it returns a ``Quantity`` rather than a scalar. [#8617]



Version 2.0.12 (2019-02-23)
===========================

New Features
------------

astropy.utils
^^^^^^^^^^^^^

- The ``deprecated_renamed_argument`` decorator now capable deprecating an
  argument without renaming it. It also got a new ``alternative`` keyword
  argument to suggest alternative functionality instead of the removed
  one. [#8324]


Bug Fixes
---------

astropy.io.fits
^^^^^^^^^^^^^^^

- Fixed bug in ``ColDefs._init_from_array()`` that caused non-scalar unsigned
  entries to not have the correct bzero value set. [#8353]

astropy.modeling
^^^^^^^^^^^^^^^^

- Fixed compatibility of ``JointFitter`` with the latest version of Numpy.
  [#7984]

astropy.table
^^^^^^^^^^^^^

- Fix ``.quantity`` property of ``Column`` class for function-units (e.g.,
  ``dex``). Previously setting this was possible, but getting raised
  an error. [#8425]

- Fixes a bug where initializing a new ``Table`` from the final row of an
  existing ``Table`` failed.  This happened when that row was generated using
  the item index ``[-1]``. [#8422]

astropy.wcs
^^^^^^^^^^^

- Fix bug that caused ``WCS.has_celestial``, ``wcs_to_celestial_frame``, and
  other functionality depending on it to fail in the presence of correlated
  celestial and other axes. [#8420]


Other Changes and Additions
---------------------------

- Fixed ``make clean`` for the documentation on Windows to ensure it
  properly removes the ``api`` and ``generated`` directories. [#8346]

- Updating bundled ``pytest-openfiles`` to v0.3.2. [#8434]

- Making ``ErfaWarning`` and ``ErfaError`` available via
  ``astropy.utils.exceptions``. [#8441]



Version 2.0.11 (2018-12-31)
===========================

Bug Fixes
---------

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fix fast reader C tokenizer to handle double quotes in quoted field.
  [#8283]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fix a bug in ``io.fits`` with writing Fortran-ordered arrays to file
  objects. [#8282]

astropy.units
^^^^^^^^^^^^^

- Add support for ``np.matmul`` as a ``ufunc`` (new in numpy 1.16).
  [#8264, #8305]

astropy.utils
^^^^^^^^^^^^^

- Fix failures caused by IERS_A_URL being unavailable by introducing
  IERS_A_URL_MIRROR. [#8308]



Version 2.0.10 (2018-12-04)
===========================

Bug Fixes
---------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Fix Moffat2DKernel's FWHM computation, which has an influence on the default
  size of the kernel when no size is given. [#8105]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Disable ``of_address`` usage due to Google API now requiring API key. [#7993]

astropy.io.fits
^^^^^^^^^^^^^^^

- ``fits.append`` now correctly handles file objects with valid modes other
  than ``ostream``. [#7856]

astropy.table
^^^^^^^^^^^^^

- Fix ``Table.show_in_notebook`` failure when mixin columns are present. [#8069]

astropy.tests
^^^^^^^^^^^^^

- Explicitly disallow incompatible versions of ``pytest`` when using the test
  runner. [#8188]

astropy.units
^^^^^^^^^^^^^

- Fixed the spelling of the 'luminous emittance/illuminance' physical
  property. [#7942]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Fixed a bug that caused origin to be incorrect if not specified. [#7927]

- Fixed a bug that caused an error when plotting grids multiple times
  with grid_type='contours'. [#7927]

- Put an upper limit on the number of bins in ``hist`` and ``histogram`` and
  factor out calculation of bin edges into public function
  ``calculate_bin_edges``. [#7991]


Other Changes and Additions
---------------------------

- Fixing ``astropy.__citation__`` to provide the full bibtex entry of the 2018
  paper. [#8110]

- Pytest 4.0 is not supported by the 2.0.x LTS releases. [#8173]

- Updating bundled ``pytest-remotedata`` to v0.3.1. [#8174]

- Updating bundled ``pytest-doctestplus`` to v0.2.0. [#8175]

- Updating bundled ``pytest-openfiles`` to v0.3.0. [#8176]

- Adding ``warning_type`` keyword argument to the "deprecated" decorators to
  allow issuing custom warning types instead of the default
  ``AstropyDeprecationWarning``. [#8178]


Version 2.0.9 (2018-10-14)
==========================

Bug Fixes
---------

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fix reading of big files with the fast reader. [#7885]

astropy.io.fits
^^^^^^^^^^^^^^^

- ``HDUList.__contains__()`` now works with ``HDU`` arguments. That is,
  ``hdulist[0] in hdulist`` now works as expected. [#7282]

- ``HDUList`` s can now be written to streams in Python 3 [#7850]

astropy.nddata
^^^^^^^^^^^^^^

- Fixed the bug in CCData.read when the HDU is not specified and the first one
  is empty so the function searches for the first HDU with data which may not
  have an image extension. [#7739]

astropy.stats
^^^^^^^^^^^^^

- Fixed bugs in biweight statistics functions where a constant data
  array (or if using the axis keyword, constant along an axis) would
  return NaN. [#7737]

astropy.table
^^^^^^^^^^^^^

- Fixed a bug in ``to_pandas()`` where integer type masked columns were always
  getting converted to float. This could cause loss of precision. Now this only
  occurs if there are actually masked data values, in which case ``pandas``
  does require the values to be float so that ``NaN`` can be used to mark the
  masked values. [#7741, #7747]

astropy.tests
^^^^^^^^^^^^^

- Change the name of the configuration variable controlling the location of the
  Astropy cache in the Pytest plugin from ``cache_dir`` to
  ``astropy_cache_dir``. The command line flag also changed to
  ``--astropy-cache-dir``.  This prevents a conflict with the ``cache_dir``
  variable provided by pytest itself. Also made similar change to
  ``config_dir`` option as a precaution. [#7721]

astropy.units
^^^^^^^^^^^^^

- ``UnrecognizedUnit`` instances can now be compared to any other object
  without raising `TypeError`. [#7606]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Fix compatibility with Matplotlib 3.0. [#7839]

- Fix an issue that caused a crash when using WCSAxes with a custom Transform
  object and when using ``grid_type='contours'`` to plot a grid. [#7846]

astropy.wcs
^^^^^^^^^^^

- Instead of raising an error ``astropy.wcs`` now returns the input when
  the input has zero size.                                       [#7746]

- Fix ``malloc(0)`` bug in ``pipeline_all_pixel2world()`` and
  ``pipeline_pix2foc()``. They now raise an exception for input with
  zero coordinates, i.e. shape = (0, n). [#7806]

- Fixed an issue with scalar input when WCS.naxis is one. [#7858]

Other Changes and Additions
---------------------------

- Added a new ``astropy.__citation__`` attribute which gives a citation
  for Astropy in bibtex format. Made sure that both this and
  ``astropy.__bibtex__`` works outside the source environment, too. [#7718]



Version 2.0.8 (2018-08-02)
==========================

Bug Fixes
---------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Correct data type conversion for non-float masked kernels. [#7542]

- Fix non-float or masked, zero sum kernels when ``normalize_kernel=False``.
  Non-floats would yield a type error and masked kernels were not being filled.
  [#7541]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Ensure that relative humidities can be given as Quantities, rather than take
  any quantity and just strip its unit. [#7668]

astropy.nddata
^^^^^^^^^^^^^^

- Fixed ``Cutout2D`` output WCS NAXIS values to reflect the cutout
  image size. [#7552]

astropy.table
^^^^^^^^^^^^^

- Fixed a bug in ``add_columns`` method where ``rename_duplicate=True`` would
  cause an error if there were no duplicates. [#7540]

astropy.tests
^^^^^^^^^^^^^

- Fixed bug in ``python setup.py test --coverage`` on Windows machines. [#7673]

astropy.time
^^^^^^^^^^^^

- Avoid rounding errors when converting ``Quantity`` to ``TimeDelta``. [#7625]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Fixed a bug that caused the position of the tick values in decimal mode
  to be incorrectly determined. [#7332]

astropy.wcs
^^^^^^^^^^^

- Fixed a bug that caused ``wcs_to_celestial_frame``, ``skycoord_to_pixel``, and
  ``pixel_to_skycoord`` to raise an error if the axes of the celestial WCS were
  swapped. [#7691]


Version 2.0.7 (2018-06-01)
==========================

Bug Fixes
---------

astropy.modeling
^^^^^^^^^^^^^^^^

- Fixed ``Tabular`` models to not change the shape of data. [#7411]

astropy.stats
^^^^^^^^^^^^^

- In ``freedman_bin_width``, if the data has too small IQR,
  raise ``ValueError``. [#7248, #7402]

astropy.table
^^^^^^^^^^^^^

- Fix a performance issue in ``MaskedColumn`` where initialization was
  extremely slow for large arrays with the default ``mask=None``. [#7422]

- Fix printing table row indexed with unsigned integer. [#7469]

- Fix copy of mask when copying a Table, as this is no more done systematically
  by Numpy since version 1.14. Also fixed a problem when MaskedColumn was
  initialized with ``mask=np.ma.nomask``. [#7486]

astropy.time
^^^^^^^^^^^^

- Fixed a bug in Time that raised an error when initializing a subclass of Time
  with a Time object. [#7453]

astropy.utils
^^^^^^^^^^^^^

- Fixed a bug that improperly handled unicode case of URL mirror in Python 2.
  [#7493]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Fixed a bug that prevented legends from being added to plots done with
  units. [#7510]


Other Changes and Additions
---------------------------

- Bundled ``pytest-remotedata`` plugin is upgraded to 0.3. [#7493]


Version 2.0.6 (2018-04-23)
==========================

Bug Fixes
---------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- convolve(boundary=None) requires the kernel to be smaller than the image.
  This was never actually checked, it now is and an exception is raised.
  [#7313]

astropy.units
^^^^^^^^^^^^^

- ``u.quantity_input`` no longer errors if the return annotation for a
  function is ``None``. [#7336, #7380]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Explicitly default to origin='lower' in WCSAxes. [#7331]

- Lists of units are now converted in the Matplotlib unit converter. This means
  that for Matplotlib versions later than 2.2, more plotting functions now work
  with units (e.g. errorbar). [#7037]


Other Changes and Additions
---------------------------

- Updated the bundled CFITSIO library to 3.44. This is to remedy another
  critical security vulnerability that was identified by NASA. See
  ``cextern/cfitsio/docs/changes.txt`` for additional information. [#7370]


Version 2.0.5 (2018-03-12)
==========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Add a workaround for a bug in the einsum function in Numpy 1.14.0. [#7187]

- Fix problems with printing ``Angle`` instances under numpy 1.14.1. [#7234]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fixed the ``fitsdiff`` script for matching fits file with one in a
  directory path. [#7085]

- Make sure that lazily-loaded ``HDUList`` is automatically loaded when calling
  ``hdulist.pop``. [#7186]

astropy.modeling
^^^^^^^^^^^^^^^^

- Propagate weights to underlying fitter in ``FittingWithOutlierRemoval`` [#7249]

astropy.tests
^^^^^^^^^^^^^

- Support dotted package names as namespace packages when gathering test
  coverage. [#7170]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Matplotlib axes have the ``axisbelow`` property to control the z-order of
  ticks, tick labels, and grid lines. WCSAxes will now respect this property.
  This is useful for drawing scale bars or inset boxes, which should have a
  z-order that places them above all ticks and gridlines. [#7098]


Other Changes and Additions
---------------------------

- Updated the bundled CFITSIO library to 3.430. This is to remedy a critical
  security vulnerability that was identified by NASA. See
  ``cextern/cfitsio/docs/changes.txt`` for additional information. [#7274, #7275]


Version 2.0.4 (2018-02-06)
==========================

Bug Fixes
---------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Fixed IndexError when ``preserve_nan=True`` in ``convolve_fft``. Added
  testing with ``preserve_nan=True``. [#7000]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- The ``sites.json`` file is now parsed explicitly with a UTF-8 encoding. This
  means that future revisions to the file with unicode observatory names can
  be done without breaking the site registry parser.  [#7082]

- Working around a bug in Numpy 1.14.0 that broke some coordinate
  transformations. [#7105]

- Fixed a bug where negative angles could be rounded wrongly when converting
  to a string with seconds omitted. [#7148]

astropy.io.fits
^^^^^^^^^^^^^^^

- When datafile is missing, fits.tabledump uses input file name to build
  output file name. Fixed how it gets input file name from HDUList. [#6976]

- Fix in-place updates to scaled columns. [#6956]

astropy.io.registry
^^^^^^^^^^^^^^^^^^^

- Fixed bug in identifying inherited registrations from multiple ancestors [#7156]

astropy.modeling
^^^^^^^^^^^^^^^^

- Fixed a bug in ``LevMarLSQFitter`` when fitting 2D models with constraints. [#6705]

astropy.utils
^^^^^^^^^^^^^

- ``download_file`` function will check for cache downloaded from mirror URL
  first before attempting actual download if primary URL is unavailable. [#6987]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Fixed test failures for ``astropy.visualization.wcsaxes`` which were due to
  local matplotlibrc files being taken into account. [#7132]


Other Changes and Additions
---------------------------

- Fixed broken links in the documentation. [#6745]

- Substantial performance improvement (potentially >1000x for some cases) when
  converting non-scalar ``coordinates.Angle`` objects to strings. [#7004]


Version 2.0.3 (2017-12-13)
==========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Ecliptic frame classes now support attributes ``v_x``, ``v_y``, ``v_z`` when
  using with a Cartesian representation. [#6569]

- Added a nicer error message when accidentally calling ``frame.representation``
  instead of ``frame.data`` in the context of methods that use ``._apply()``.
  [#6561]

- Creating a new ``SkyCoord`` from a list of multiple ``SkyCoord`` objects now
  yield the correct type of frame, and works at all for non-equatorial frames.
  [#6612]

- Improved accuracy of velocity calculation in ``EarthLocation.get_gcrs_posvel``.
  [#6699]

- Improved accuracy of radial velocity corrections in
  ``SkyCoord.radial_velocity_correction```. [#6861]

- The precision of ecliptic frames is now much better, after removing the
  nutation from the rotation and fixing the computation of the position of the
  Sun. [#6508]

astropy.extern
^^^^^^^^^^^^^^

- Version 0.2.1 of ``pytest-astropy`` is included as an external package.
  [#6918]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fix writing the result of ``fitsdiff`` to file with ``--output-file``. [#6621]

- Fix a minor bug where ``FITS_rec`` instances can not be indexed with tuples
  and other sequences that end up with a scalar. [#6955, #6966]

astropy.io.misc
^^^^^^^^^^^^^^^

- Fix ``ImportError`` when ``hdf5`` is imported first in a fresh Python
  interpreter in Python 3. [#6604, #6610]

astropy.nddata
^^^^^^^^^^^^^^

- Suppress errors during WCS creation in CCDData.read(). [#6500]

- Fixed a problem with ``CCDData.read`` when the extension wasn't given and the
  primary HDU contained no ``data`` but another HDU did. In that case the header
  were not correctly combined. [#6489]

astropy.stats
^^^^^^^^^^^^^

- Fixed an issue where the biweight statistics functions would
  sometimes cause runtime underflow/overflow errors for float32 input
  arrays. [#6905]

astropy.table
^^^^^^^^^^^^^

- Fixed a problem when printing a table when a column is deleted and
  garbage-collected, and the format function caching mechanism happens
  to reuse the same cache key. [#6714]

- Fixed a problem when comparing a unicode masked column (on left side) to
  a bytes masked column (on right side). [#6899]

- Fixed a problem in comparing masked columns in bytes and unicode when the
  unicode had masked entries. [#6899]

astropy.tests
^^^^^^^^^^^^^

- Fixed a bug that causes tests for rst files to not be run on certain
  platforms. [#6555, #6608]

- Fixed a bug that caused the doctestplus plugin to not work nicely with the
  hypothesis package. [#6605, #6609]

- Fixed a bug that meant that the data.astropy.org mirror could not be used when
  using --remote-data=astropy. [#6724]

- Support compatibility with new ``pytest-astropy`` plugins. [#6918]

- When testing, astropy (or the package being tested) is now installed to
  a temporary directory instead of copying the build. This allows
  entry points to work correctly. [#6890]

astropy.time
^^^^^^^^^^^^

- Initialization of Time instances now is consistent for all formats to
  ensure that ``-0.5 <= jd2 < 0.5``. [#6653]

astropy.units
^^^^^^^^^^^^^

- Ensure that ``Quantity`` slices can be set with objects that have a ``unit``
  attribute (such as ``Column``). [#6123]

astropy.utils
^^^^^^^^^^^^^

- ``download_files_in_parallel`` now respects the given ``timeout`` value.
  [#6658]

- Fixed bugs in remote data handling and also in IERS unit test related to path
  URL, and URI normalization on Windows. [#6651]

- Fixed a bug that caused ``get_pkg_data_fileobj`` to not work correctly when
  used with non-local data from inside packages. [#6724]

- Make sure ``get_pkg_data_fileobj`` fails if the URL can not be read, and
  correctly falls back on the mirror if necessary. [#6767]

- Fix the ``finddiff`` option in ``find_current_module`` to properly deal
  with submodules. [#6767]

- Fixed ``pyreadline`` import in ``utils.console.isatty`` for older IPython
  versions on Windows. [#6800]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Fixed the vertical orientation of the ``fits2bitmap`` output bitmap
  image to match that of the FITS image. [#6844, #6969]

- Added a workaround for a bug in matplotlib so that the ``fits2bitmap``
  script generates the correct output file type. [#6969]


Other Changes and Additions
---------------------------

- No longer require LaTeX to build the documentation locally and
  use mathjax instead. [#6701]

- Ensured that all tests use the Astropy data mirror if needed. [#6767]


Version 2.0.2 (2017-09-08)
==========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Ensure transformations via ICRS also work for coordinates that use cartesian
  representations. [#6440]

- Fixed a bug that was preventing ``SkyCoord`` objects made from lists of other
  coordinate objects from being written out to ECSV files. [#6448]

astropy.io.fits
^^^^^^^^^^^^^^^

- Support the ``GZIP_2`` FITS image compression algorithm as claimed
  in docs. [#6486]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Fixed a bug that wrote out VO table as version 1.2 instead of 1.3. [#6521]

astropy.table
^^^^^^^^^^^^^

- Fix a bug when combining unicode columns via join or vstack.  The character
  width of the output column was a factor of 4 larger than needed. [#6459]

astropy.tests
^^^^^^^^^^^^^

- Fixed running the test suite using --parallel. [#6415]

- Added error handling for attempting to run tests in parallel without having
  the ``pytest-xdist`` package installed. [#6416]

- Fixed issue running doctests with pytest>=3.2. [#6423, #6430]

- Fixed issue caused by antivirus software in response to malformed compressed
  files used for testing. [#6522]

- Updated top-level config file to properly ignore top-level directories.
  [#6449]

astropy.units
^^^^^^^^^^^^^

- Quantity._repr_latex_ now respects precision option from numpy
  printoptions. [#6412]

astropy.utils
^^^^^^^^^^^^^

- For the ``deprecated_renamed_argument`` decorator, refer to the deprecation‘s
  caller instead of ``astropy.utils.decorators``, to makes it easier to find
  where the deprecation warnings comes from. [#6422]


Version 2.0.1 (2017-07-30)
==========================

Bug Fixes
---------

astropy.constants
^^^^^^^^^^^^^^^^^

- Fixed Earth radius to be the IAU2015 value for the equatorial radius.
  The polar value had erroneously been used in 2.0. [#6400]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Added old frame attribute classes back to top-level namespace of
  ``astropy.coordinates``. [#6357]

astropy.io.fits
^^^^^^^^^^^^^^^

- Scaling an image always uses user-supplied values when given. Added
  defaults for scaling when bscale/bzero are not present (float images).
  Fixed a small bug in when to reset ``_orig_bscale``. [#5955]

astropy.modeling
^^^^^^^^^^^^^^^^

- Fixed a bug in initializing compound models with units. [#6398]

astropy.nddata
^^^^^^^^^^^^^^

- Updating CCDData.read() to be more flexible with inputs, don't try to
  delete keywords that are missing from the header. [#6388]

astropy.tests
^^^^^^^^^^^^^
- Fixed the test command that is run from ``setuptools`` to allow it to
  gracefully handle keyboard interrupts and pass them on to the ``pytest``
  subprocess. This prompts ``pytest`` to teardown and display useful traceback
  and test information [#6369]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Ticks and tick labels are now drawn in front of, rather than behind,
  gridlines in WCS axes. This improves legibility in situations where
  tick labels may be on the interior of the axes frame, such as the right
  ascension axis of an all-sky Aitoff or Mollweide projection. [#6361]

astropy.wcs
^^^^^^^^^^^

- Fix the missing wcskey part in _read_sip_kw, this will cause error when
  reading sip wcs while there is no default CRPIX1 CRPIX2 keywords and only
  CRPIX1n CRPIX2n in header. [#6372]



Version 2.0 (2017-07-07)
========================

New Features
------------

astropy.constants
^^^^^^^^^^^^^^^^^

- Constants are now organized into version modules, with physical CODATA
  constants in the ``codata2010`` and ``codata2014`` sub-modules,
  and astronomical constants defined by the IAU in the ``iau2012`` and
  ``iau2015`` sub-modules. The default constants in ``astropy.constants``
  in Astropy 2.0 have been updated from ``iau2012`` to ``iau2015`` and
  from ``codata2010`` to ``codata2014``. The constants for 1.3 can be
  accessed in the ``astropyconst13`` sub-module and the constants for 2.0
  (the default in ``astropy.constants``) can also be accessed in the
  ``astropyconst20`` sub-module [#6083]

- The GM mass parameters recommended by IAU 2015 Resolution B 3 have been
  added as ``GM_sun``, ``GM_jup``, and ``GM_earth``, for the Sun,
  Jupiter and the Earth. [#6083]

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Major change in convolution behavior and keyword arguments. Additional
  details are in the API section. [#5782]

- Convolution with un-normalized and un-normalizable kernels is now possible.
  [#5782]

- Add a new argument, ``normalization_rtol``, to ``convolve_fft``, allowing
  the user to specify the relative error tolerance in the normalization of
  the convolution kernel. [#5649, #5177]

- Models can now be convoluted using ``convolve`` or ``convolve_fft``,
  which generates a regular compound model. [#6015]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Frame attributes set on ``SkyCoord`` are now always validated, and any
  ndarray-like operation (like slicing) will also be done on those. [#5751]

- Caching of  all possible frame attributes was implemented. This greatly
  speeds up many ``SkyCoord`` operations. [#5703, #5751]

- A class hierarchy was added to allow the representation layer to store
  differentials (i.e., finite derivatives) of coordinates.  This is intended
  to enable support for velocities in coordinate frames. [#5871]

- ``replicate_without_data`` and ``replicate`` methods were added to
  coordinate frames that allow copying an existing frame object with various
  reference or copy behaviors and possibly overriding frame attributes. [#6182]

- The representation class instances can now contain differential objects.
  This is primarily useful for internal operations that will provide support
  for transforming velocity components in coordinate frames. [#6169]

- ``EarthLocation.to_geodetic()`` (and ``EarthLocation.geodetic``) now return
  namedtuples instead of regular tuples. [#6237]

- ``EarthLocation`` now has ``lat`` and ``lon`` properties (equivalent to, but
  preferred over, the previous ``latitude`` and ``longitude``). [#6237]

- Added a ``radial_velocity_correction`` method to ``SkyCoord`` to do compute
  barycentric and heliocentric velocity corrections. [#5752]

- Added a new ``AffineTransform`` class for coordinate frame transformations.
  This class supports matrix operations with vector offsets in position or
  any differential quantities (so far, only velocity is supported). The
  matrix transform classes now subclass from the base affine transform.
  [#6218]

- Frame objects now have experimental support for velocity components. Most
  frames default to accepting proper motion components and radial velocity,
  and the velocities transform correctly for any transformation that uses
  one of the ``AffineTransform``-type transformations.  For other
  transformations a finite-difference velocity transformation is available,
  although it is not as numerically stable as those that use
  ``AffineTransform``-type transformations. [#6219, #6226]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Allow to specify encoding in ``ascii.read``, only for Python 3 and with the
  pure-Python readers. [#5448]

- Writing latex tables with only a ``tabular`` environment is now possible by
  setting ``latexdict['tabletyle']`` to ``None``. [#6205]

- Allow ECSV format to support reading and writing mixin columns like
  ``Time``, ``SkyCoord``, ``Latitude``, and ``EarthLocation``. [#6181]

astropy.io.fits
^^^^^^^^^^^^^^^

- Checking available disk space before writing out file. [#5550, #4065]

- Change behavior to warn about units that are not FITS-compliant when
  writing a FITS file but not when reading. [#5675]

- Added absolute tolerance parameter when comparing FITS files. [#4729]

- New convenience function ``printdiff`` to print out diff reports. [#5759]

- Allow to instantiate a ``BinTableHDU`` directly from a ``Table`` object.
  [#6139]

astropy.io.misc
^^^^^^^^^^^^^^^

- YAML representer now also accepts numpy types. [#6077]

astropy.io.registry
^^^^^^^^^^^^^^^^^^^

- New functions to unregister readers, writers, and identifiers. [#6217]

astropy.modeling
^^^^^^^^^^^^^^^^

- Added ``SmoothlyBrokenPowerLaw1D`` model. [#5656]

- Add ``n_submodels`` shared method to single and compound models, which
  allows users to get the number of components of a given single (compound)
  model. [#5747]

- Added a ``name`` setter for instances of ``_CompoundModel``. [#5741]

- Added FWHM properties to Gaussian and Moffat models. [#6027]

- Added support for evaluating models and setting the results for inputs
  outside the bounding_box to a user specified ``fill_value``. This
  is controlled by a new optional boolean keyword ``with_bounding_box``. [#6081]

- Added infrastructure support for units on parameters and during
  model evaluation and fitting, added support for units on all
  functional, power-law, polynomial, and rotation models where this
  is appropriate. A new BlackBody1D model has been added.
  [#4855, #6183, #6204, #6235]

astropy.nddata
^^^^^^^^^^^^^^

- Added an image class, ``CCDData``. [#6173]

astropy.stats
^^^^^^^^^^^^^

- Added ``biweight_midcovariance`` function. [#5777]

- Added ``biweight_scale`` and ``biweight_midcorrelation``
  functions. [#5991]

- ``median_absolute_deviation`` and ``mad_std`` have ``ignore_nan`` option
  that will use ``np.ma.median`` with nans masked out or ``np.nanmedian``
  instead of ``np.median`` when computing the median. [#5232]

- Implemented statistical estimators for Ripley's K Function. [#5712]

- Added ``SigmaClip`` class. [#6206]

- Added ``std_ddof`` keyword option to ``sigma_clipped_stats``.
  [#6066, #6207]

astropy.table
^^^^^^^^^^^^^

- Issue a warning when assigning a string value to a column and
  the string gets truncated.  This can occur because numpy string
  arrays are fixed-width and silently drop characters which do not
  fit within the fixed width. [#5624, #5819]

- Added functionality to allow ``astropy.units.Quantity`` to be written
  as a normal column to FITS files. [#5910]

- Add support for Quantity columns (within a ``QTable``) in table
  ``join()``, ``hstack()`` and ``vstack()`` operations. [#5841]

- Allow unicode strings to be stored in a Table bytestring column in
  Python 3 using UTF-8 encoding.  Allow comparison and assignment of
  Python 3 ``str`` object in a bytestring column (numpy ``'S'`` dtype).
  If comparison with ``str`` instead of ``bytes`` is a problem
  (and ``bytes`` is really more logical), please open an issue on GitHub.
  [#5700]

- Added functionality to allow ``astropy.units.Quantity`` to be read
  from and written to a VOtable file. [#6132]

- Added support for reading and writing a table with mixin columns like
  ``Time``, ``SkyCoord``, ``Latitude``, and ``EarthLocation`` via the
  ASCII ECSV format. [#6181]

- Bug fix for ``MaskedColumn`` insert method, where ``fill_value`` attribute
  was not being passed along to the copy of the ``MaskedColumn`` that was
  returned. [#7585]

astropy.tests
^^^^^^^^^^^^^

- ``enable_deprecations_as_exceptions`` function now accepts additional
  user-defined module imports and warning messages to ignore. [#6223, #6334]

astropy.units
^^^^^^^^^^^^^

- The ``astropy.units.quantity_input`` decorator will now convert the output to
  the unit specified as a return annotation under Python 3. [#5606]

- Passing a logarithmic unit to the ``Quantity`` constructor now returns the
  appropriate logarithmic quantity class if ``subok=True``. For instance,
  ``Quantity(1, u.dex(u.m), subok=True)`` yields ``<Dex 1.0 dex(m)>``. [#5928]

- The ``quantity_input`` decorator now accepts a string physical type in
  addition to of a unit object to specify the expected input ``Quantity``'s
  physical type. For example, ``@u.quantity_input(x='angle')`` is now
  functionally the same as ``@u.quantity_input(x=u.degree)``. [#3847]

- The ``quantity_input`` decorator now also supports unit checking for
  optional keyword arguments and accepts iterables of units or physical types
  for specifying multiple valid equivalent inputs. For example,
  ``@u.quantity_input(x=['angle', 'angular speed'])`` or
  ``@u.quantity_input(x=[u.radian, u.radian/u.yr])`` would both allow either
  a ``Quantity`` angle or angular speed passed in to the argument ``x``.
  [#5653]

- Added a new equivalence ``molar_mass_amu`` between g/mol to
  atomic mass units. [#6040, #6113]

- ``Quantity`` has gained a new ``to_value`` method which returns the value
  of the quantity in a given unit. [#6127]

- ``Quantity`` now supports the ``@`` operator for matrix multiplication that
  was introduced in Python 3.5, for all supported versions of numpy. [#6144]

- ``Quantity`` supports the new ``__array_ufunc__`` protocol introduced in
  numpy 1.13.  As a result, operations that involve unit conversion will be
  sped up considerably (by up to a factor of two for costly operations such
  as trigonometric ones). [#2583]

astropy.utils
^^^^^^^^^^^^^

- Added a new ``dataurl_mirror`` configuration item in ``astropy.utils.data``
  that is used to indicate a mirror for the astropy data server. [#5547]

- Added a new convenience method ``get_cached_urls`` to ``astropy.utils.data``
  for getting a list of the URLs in your cache. [#6242]

astropy.wcs
^^^^^^^^^^^

- Upgraded the included wcslib to version 5.16. [#6225]

  The minimum required version of wcslib in is 5.14.


API Changes
-----------

astropy.analytic_functions
^^^^^^^^^^^^^^^^^^^^^^^^^^

- This entire sub-package is deprecated because blackbody has been moved to
  ``astropy.modeling.blackbody``. [#6191]

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Major change in convolution behavior and keyword arguments.
  ``astropy.convolution.convolve_fft`` replaced ``interpolate_nan`` with
  ``nan_treatment``, and ``astropy.convolution.convolve`` received a new
  ``nan_treatment`` argument. ``astropy.convolution.convolve`` also no longer
  double-interpolates interpolates over NaNs, although that is now available
  as a separate ``astropy.convolution.interpolate_replace_nans`` function. See
  `the backwards compatibility note
  <https://docs.astropy.org/en/v2.0.16/convolution/index.html#a-note-on-backward-compatibility-pre-v2-0>`_
  for more on how to get the old behavior (and why you probably don't want to.)
  [#5782]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- The ``astropy.coordinates.Galactic`` frame previously was had the cartesian
  ordering 'w', 'u', 'v' (for 'x', 'y', and 'z', respectively).  This was an
  error and against the common convention.  The 'x', 'y', and 'z' axes now
  map to 'u', 'v', and 'w', following the right-handed ('u' points to
  the Galactic center) convention. [#6330]

- Removed deprecated ``angles.rotation_matrix`` and
  ``angles.angle_axis``. Use the routines in
  ``coordinates.matrix_utilities`` instead. [#6170]

- ``EarthLocation.latitude`` and ``EarthLocation.longitude`` are now
  deprecated in favor of ``EarthLocation.lat`` and ``EarthLocation.lon``.
  They former will be removed in a future version. [#6237]

- The ``FrameAttribute`` class and subclasses have been renamed to just contain
  ``Attribute``. For example, ``QuantityFrameAttribute`` is now
  ``QuantityAttribute``. [#6300]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Cosmological models do not include any contribution from neutrinos or photons
  by default -- that is, the default value of Tcmb0 is 0.  This does not affect
  built in models (such as WMAP or Planck). [#6112]

astropy.io.fits
^^^^^^^^^^^^^^^

- Remove deprecated ``NumCode`` and ``ImgCode`` properties on FITS
  ``_ImageBaseHDU``.  Use module-level constants ``BITPIX2DTYPE`` and
  ``DTYPE2BITPIX`` instead. [#4993]

- ``comments`` meta key (which is ``io.ascii``'s table convention) is output
  to ``COMMENT`` instead of ``COMMENTS`` header. Similarly, ``COMMENT``
  headers are read into ``comments`` meta [#6097]

- Remove compatibility code which forced loading all HDUs on close. The old
  behavior can be used with ``lazy_load_hdus=False``. Because of this change,
  trying to access the ``.data`` attribute from an HDU which is not loaded
  now raises a ``IndexError`` instead of a ``ValueError``. [#6082]

- Deprecated ``clobber`` keyword; use ``overwrite``. [#6203]

- Add EXTVER column to the output of ``HDUList.info()``. [#6124]

astropy.modeling
^^^^^^^^^^^^^^^^

- Removed deprecated ``Redshift`` model; Use ``RedshiftScaleFactor``. [#6053]

- Removed deprecated ``Pix2Sky_AZP.check_mu`` and ``Pix2Sky_SZP.check_mu``
  methods. [#6170]

- Deprecated ``GaussianAbsorption1D`` model, as it can be better represented
  by subtracting ``Gaussian1D`` from ``Const1D``. [#6200]

- Added method ``sum_of_implicit_terms`` to ``Model``, needed when performing
  a linear fit to a model that has built-in terms with no corresponding
  parameters (primarily the ``1*x`` term of ``Shift``). [#6174]

astropy.nddata
^^^^^^^^^^^^^^

- Removed deprecated usage of parameter ``propagate_uncertainties`` as a
  positional keyword. [#6170]

- Removed deprecated ``support_correlated`` attribute. [#6170]

- Removed deprecated ``propagate_add``, ``propagate_subtract``,
  ``propagate_multiply`` and ``propagate_divide`` methods. [#6170]

astropy.stats
^^^^^^^^^^^^^

- Removed the deprecated ``sig`` and ``varfunc`` keywords in the
  ``sigma_clip`` function. [#5715]

- Added ``modify_sample_size`` keyword to ``biweight_midvariance``
  function. [#5991]

astropy.table
^^^^^^^^^^^^^

- In Python 3, when getting an item from a bytestring Column it is now
  converted to ``str``.  This means comparing a single item to a ``bytes``
  object will always fail, and instead one must compare with a ``str``
  object. [#5700]

- Removed the deprecated ``data`` property of Row. [#5729]

- Removed the deprecated functions ``join``, ``hstack``, ``vstack`` and
  ``get_groups`` from np_utils. [#5729]

- Added ``name`` parameter to method ``astropy.table.Table.add_column`` and
  ``names`` parameter to method ``astropy.table.Table.add_columns``, to
  provide the flexibility to add unnamed columns, mixin objects and also to
  specify explicit names. Default names will be used if not
  specified. [#5996]

- Added optional ``axis`` parameter to ``insert`` method for ``Column`` and
  ``MaskedColumn`` classes. [#6092]

astropy.units
^^^^^^^^^^^^^

- Moved ``units.cgs.emu`` to ``units.deprecated.emu`` due to ambiguous
  definition of "emu". [#4918, #5906]

- ``jupiterMass``, ``earthMass``, ``jupiterRad``, and ``earthRad`` no longer
  have their prefixed units included in the standard units.  If needed, they
  can still  be found in ``units.deprecated``. [#5661]

- ``solLum``,``solMass``, and ``solRad`` no longer have  their prefixed units
  included in the standard units.  If needed, they can still be found in
  ``units.required_by_vounit``, and are enabled by default. [#5661]

- Removed deprecated ``Unit.get_converter``. [#6170]

- Internally, astropy replaced use of ``.to(unit).value`` with the new
  ``to_value(unit)`` method, since this is somewhat faster. Any subclasses
  that overwrote ``.to``, should also overwrite ``.to_value`` (or
  possibly just the private ``._to_value`` method.  (If you did this,
  please let us know what was lacking that made this necessary!). [#6137]

astropy.utils
^^^^^^^^^^^^^

- Removed the deprecated compatibility modules for Python 2.6 (``argparse``,
  ``fractions``, ``gzip``, ``odict``, ``subprocess``) [#5975,#6157,#6164]

- Removed the deprecated ``zest.releaser`` machinery. [#6282]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Removed the deprecated ``scale_image`` function. [#6170]

astropy.vo
^^^^^^^^^^

- Cone Search now issues deprecation warning because it is moved to
  Astroquery 0.3.5 and will be removed from Astropy in a future version.
  [#5558, #5904]

- The ``astropy.vo.samp`` package has been moved to ``astropy.samp``, and no
  longer supports HTTPS/SSL. [#6201, #6213]

astropy.wcs
^^^^^^^^^^^

- Removed deprecated ``wcs.rotateCD``. [#6170]


Bug Fixes
---------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Major change in convolution behavior and keyword arguments:
  ``astropy.convolution.convolve`` was not performing normalized convolution
  in earlier versions of astropy. [#5782]

- Direct convolution previously implemented the wrong definition of
  convolution.  This error only affects *asymmetric* kernels. [#6267]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- The ``astropy.coordinates.Galactic`` frame had an incorrect ordering for the
  'u', 'v', and 'w' cartesian coordinates. [#6330]

- The ``astropy.coordinates.search_around_sky``,
  ``astropy.coordinates.search_around_3d``, and ``SkyCoord`` equivalent methods
  now correctly yield an ``astropy.coordinates.Angle`` as the third return type
  even if there are no matches (previously it returned a raw Quantity). [#6347]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fix an issue where the fast C-reader was dropping table comments for a
  table with no data lines. [#8274]

astropy.io.fits
^^^^^^^^^^^^^^^

- ``comments`` meta key (which is ``io.ascii``'s table convention) is output
  to ``COMMENT`` instead of ``COMMENTS`` header. Similarly, ``COMMENT``
  headers are read into ``comments`` meta [#6097]

- Use more sensible fix values for invalid NAXISj header values. [#5935]

- Close file on error to avoid creating a ``ResourceWarning`` warning
  about an unclosed file. [#6168, #6177]

astropy.modeling
^^^^^^^^^^^^^^^^

- Creating a compound model where one of the submodels is
  a compound model whose parameters were changed now uses the
  updated parameters and not the parameters of the original model. [#5741]

- Allow ``Mapping`` and ``Identity`` to be fittable. [#6018]

- Gaussian models now impose positive ``stddev`` in fitting. [#6019]

- OrthoPolynomialBase (Chebyshev2D / Legendre2D) models were being evaluated
  incorrectly when part of a compound model (using the parameters from the
  original model), which in turn caused fitting to fail as a no-op. [#6085]

- Allow ``Ring2D`` to be defined using ``r_out``. [#6192]

- Make ``LinearLSQFitter`` produce correct results with fixed model
  parameters and allow ``Shift`` and ``Scale`` to be fitted with
  ``LinearLSQFitter`` and ``LevMarLSQFitter``. [#6174]

astropy.stats
^^^^^^^^^^^^^

- Allow to choose which median function is used in ``mad_std`` and
  ``median_absolute_deviation``. And allow to use these functions with
  a multi-dimensional ``axis``. [#5835]

- Fixed ``biweight_midvariance`` so that by default it returns a
  variance that agrees with the standard definition. [#5991]

astropy.table
^^^^^^^^^^^^^

- Fix a problem with vstack for bytes columns in Python 3. [#5628]

- Fix QTable add/insert row for multidimensional Quantity. [#6092]

astropy.time
^^^^^^^^^^^^

- Fixed the initial condition of ``TimeFITS`` to allow scale, FITS scale
  and FITS realization to be checked and equated properly. [#6202]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Fixed a bug that caused the default WCS to return coordinates offset by
  one. [#6339]

astropy.vo
^^^^^^^^^^

- Fixed a bug in vo.samp when stopping a hub for which a lockfile was
  not created. [#6211]


Other Changes and Additions
---------------------------

- Numpy 1.7 and 1.8 are no longer supported. [#6006]

- Python 3.3 is no longer supported. [#6020]

- The bundled ERFA was updated to version 1.4.0. [#6239]

- The bundled version of pytest has now been removed, but the
  astropy.tests.helper.pytest import will continue to work properly.
  Affiliated packages should nevertheless transition to importing pytest
  directly rather than from astropy.tests.helper. This also means that
  pytest is now a formal requirement for testing for both Astropy and
  for affiliated packages. [#5694]


Version 1.3.3 (2017-05-29)
==========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Fixed a bug where ``StaticMatrixTransform`` erroneously copied frame
  attributes from the input coordinate to the output frame. In practice, this
  didn't actually affect any transforms in Astropy but may change behavior for
  users who explicitly used the ``StaticMatrixTransform`` in their own code.
  [#6045]

- Fixed ``get_icrs_coordinates`` to loop through all the urls in case one
  raises an exception. [#5864]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fix table header not written out properly when ``fits.writeto()``
  convenience function is used. [#6042]

- Fix writing out read-only arrays. [#6036]

- Extension headers are written out properly when the ``fits.update()``
  convenience function is used. [#6058]

- Angstrom, erg, G, and barn are no more reported as deprecated FITS units.
  [#5929]

astropy.table
^^^^^^^^^^^^^

- Fix problem with Table pprint/pformat raising an exception for
  non-UTF-8 compliant bytestring data. [#6117]

astropy.units
^^^^^^^^^^^^^

- Allow strings 'nan' and 'inf' as Quantity inputs. [#5958]

- Add support for ``positive`` and ``divmod`` ufuncs (new in numpy 1.13).
  [#5998, #6020, #6116]

astropy.utils
^^^^^^^^^^^^^

- On systems that do not have ``pkg_resources`` non-numerical additions to
  version numbers like ``dev`` or ``rc1`` are stripped in ``minversion`` to
  avoid a ``TypeError`` in ``distutils.version.LooseVersion`` [#5944]

- Fix ``auto_download`` setting ignored in ``Time.ut1``. [#6033]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Fix bug in ManualInterval which caused the limits to be returned incorrectly
  if set to zero, and fix defaults for ManualInterval in the presence of NaNs.
  [#6088]

- Get rid of warnings that occurred when slicing a cube due to the tick
  locator trying to find ticks for the sliced axis. [#6104]

- Accept normal Matplotlib keyword arguments in set_xlabel and set_ylabel
  functions. [#5686, #5692, #6060]

- Fix a bug that caused labels to be missing from frames with labels that
  could change direction mid-axis, such as EllipticalFrame. Also ensure
  that empty tick labels do not cause any warnings. [#6063]


Version 1.3.2 (2017-03-30)
==========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Ensure that checking equivalence of ``SkyCoord`` objects works with
  non-scalar attributes [#5884, #5887]

- Ensure that transformation to frames with multi-dimensional attributes
  works as expected [#5890, #5897]

- Make sure all ``BaseRepresentation`` objects can be output as strings.
  [#5889, #5897]

astropy.units
^^^^^^^^^^^^^

- Add support for ``heaviside`` ufunc (new in numpy 1.13). [#5920]

astropy.utils
^^^^^^^^^^^^^

- Fix to allow the C-based _fast_iterparse() VOTable XML parser to
  relloc() its buffers instead of overflowing them. [#5824, #5869]


Other Changes and Additions
---------------------------

- File permissions are revised in the released source distribution. [#5912]


Version 1.3.1 (2017-03-18)
==========================

New Features
------------

astropy.utils
^^^^^^^^^^^^^

- The ``deprecated_renamed_argument`` decorator got a new ``pending``
  parameter to suppress the deprecation warnings. [#5761]

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Changed ``SkyCoord`` so that frame attributes which are not valid for the
  current ``frame`` (but are valid for other frames) are stored on the
  ``SkyCoord`` instance instead of the underlying ``frame`` instance (e.g.,
  setting ``relative_humidity`` on an ICRS ``SkyCoord`` instance.) [#5750]

- Ensured that ``position_angle`` and ``separation`` give correct answers for
  frames with different equinox (see #5722). [#5762]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fix problem with padding bytes written for BinTable columns converted
  from unicode [#5280, #5287, #5288, #5296].

- Fix out-of-order TUNITn cards when writing tables to FITS. [#5720]

- Recognize PrimaryHDU when non boolean values are present for the
  'GROUPS' header keyword. [#5808]

- Fix the insertion of new keywords in compressed image headers
  (``CompImageHeader``). [#5866]

astropy.modeling
^^^^^^^^^^^^^^^^

- Fixed a problem with setting ``bounding_box`` on 1D models. [#5718]

- Fixed a broadcasting problem with weighted fitting of 2D models
  with ``LevMarLSQFitter``. [#5788]

- Fixed a problem with passing kwargs to fitters, specifically ``verblevel``. [#5815]

- Changed FittingWithOutlierRemoval to reject on the residual to the fit [#5831]

astropy.stats
^^^^^^^^^^^^^

- Fix the psd normalization for Lomb-Scargle periodograms in the presence
  of noise. [#5713]

- Fix bug in the autofrequency range when ``minimum_frequency`` is specified
  but ``maximum_frequency`` is not. [#5738]

- Ensure that a masked array is returned when sigma clipping fully masked
  data. [#5711]

astropy.table
^^^^^^^^^^^^^

- Fix problem where key for caching column format function was not
  sufficiently unique. [#5803]

- Handle sorting NaNs and masked values in jsviewer. [#4052, #5572]

- Ensure mixin columns can be added to a table using a scalar value for the
  right-hand side if the type supports broadcasting. E.g., for an existing
  ``QTable``, ``t['q'] = 3*u.m`` will now add a column as expected. [#5820]

- Fixes the bug of setting/getting values from rows/columns of a table using
  numpy array scalars. [#5772]

astropy.units
^^^^^^^^^^^^^

- Fixed problem where IrreducibleUnits could fail to unpickle. [#5868]

astropy.utils
^^^^^^^^^^^^^

- Avoid importing ``ipython`` in ``utils.console`` until it is necessary, to
  prevent deprecation warnings when importing, e.g., ``Column``. [#5755]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Avoid importing matplotlib.pyplot when importing
  astropy.visualization.wcsaxes. [#5680, #5684]

- Ignore Numpy warnings that happen in coordinate transforms in WCSAxes.
  [#5792]

- Fix compatibility issues between WCSAxes and Matplotlib 2.x. [#5786]

- Fix a bug that caused WCSAxes frame visual properties to not be copied
  over when resetting the WCS. [#5791]

astropy.extern
^^^^^^^^^^^^^^

- Fixed a bug where PLY was overwriting its generated files. [#5728]

Other Changes and Additions
---------------------------

- Fixed a deprecation warning that occurred when running tests with
  astropy.test(). [#5689]

- The deprecation of the ``clobber`` argument (originally deprecated in 1.3.0)
  in the ``io.fits`` write functions was changed to a "pending" deprecation
  (without displaying warnings) for now. [#5761]

- Updated bundled astropy-helpers to v1.3.1. [#5880]


Version 1.3 (2016-12-22)
========================

New Features
------------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- The ``convolve`` and ``convolve_fft`` arguments now support a ``mask`` keyword,
  which allows them to also support ``NDData`` objects as inputs. [#5554]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Added an ``of_address`` classmethod to ``EarthLocation`` to enable fast creation of
  ``EarthLocation`` objects given an address by querying the Google maps API [#5154].

- A new routine, ``get_body_barycentric_posvel`` has been added that allows
  one to calculate positions as well as velocities for solar system bodies.
  For JPL kernels, this roughly doubles the execution time, so if one requires
  only the positions, one should use ``get_body_barycentric``. [#5231]

- Transformations between coordinate systems can use the more accurate JPL
  ephemerides. [#5273, #5436]

- Arithmetic on representations, such as addition of two representations,
  multiplication with a ``Quantity``, or calculating the norm via ``abs``,
  has now become possible. Furthermore, there are new methods ``mean``,
  ``sum``, ``dot``, and ``cross``. For all these, the representations are
  treated as vectors in cartesian space (temporarily converting to
  ``CartesianRepresentation`` if necessary).  [#5301]
  has now become possible. Furthermore, there are news methods ``mean``,
  ``sum``, ``dot``, and ``cross`` with obvious meaning. [#5301]
  multiplication with a ``Quantity`` has now become possible. Furthermore,
  there are new methods ``norm``, ``mean``, ``sum``, ``dot``, and ``cross``.
  In all operations, the representations are treated as vectors. They are
  temporarily converted to ``CartesianRepresentation`` if necessary.  [#5301]

- ``CartesianRepresentation`` can be initialized with plain arrays by passing
  in a ``unit``. Furthermore, for input with a vector array, the coordinates
  no longer have to be in the first dimension, but can be at any ``xyz_axis``.
  To complement the latter, a new ``get_xyz(xyz_axis)`` method allows one to
  get a vector array out along a given axis. [#5439]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Files with "Fortran-style" columns (i.e. double-precision scientific notation
  with a character other than "e", like ``1.495978707D+13``) can now be parsed by
  the fast reader natively. [#5552]

- Allow round-tripping masked data tables in most formats by using an
  empty string ``''`` as the default representation of masked values
  when writing. [#5347]

- Allow reading HTML tables with unicode column values in Python 2.7. [#5410]

- Check for self-consistency of ECSV header column names. [#5463]

- Produce warnings when writing an IPAC table from an astropy table that
  contains metadata not supported by the IPAC format. [#4700]

astropy.io.fits
^^^^^^^^^^^^^^^

- "Lazy" loading of HDUs now occurs - when an HDU is requested, the file is
  only read up to the point where that HDU is found.  This can mean a
  substantial speedup when accessing files that have many HDUs. [#5065]

astropy.io.misc
^^^^^^^^^^^^^^^

- Added ``io.misc.yaml`` module to support serializing core astropy objects
  using the YAML protocol. [#5486]

astropy.io.registry
^^^^^^^^^^^^^^^^^^^

- Added ``delay_doc_updates`` contextmanager to postpone the formatting of
  the documentation for the ``read`` and ``write`` methods of the class to
  optionally reduce the import time. [#5275]

astropy.modeling
^^^^^^^^^^^^^^^^

- Added a class to combine astropy fitters and functions to remove outliers
  e. g., sigma clip. [#4760]

- Added a ``Tabular`` model. [#5105]

- Added ``Hermite1D`` and ``Hermite2D`` polynomial models [#5242]

- Added the injection of EntryPoints into astropy.modeling.fitting if
  they inherit from Fitters class. [#5241]

- Added bounding box to ``Lorentz1D`` and ``MexicanHat1D`` models. [#5393]

- Added ``Planar2D`` functional model. [#5456]

- Updated ``Gaussian2D`` to accept no arguments (will use default x/y_stddev
  and theta). [#5537]

astropy.nddata
^^^^^^^^^^^^^^

- Added ``keep`` and ``**kwargs`` parameter to ``support_nddata``. [#5477]

astropy.stats
^^^^^^^^^^^^^

- Added ``axis`` keyword to ``biweight_location`` and
  ``biweight_midvariance``. [#5127, #5158]

astropy.table
^^^^^^^^^^^^^

- Allow renaming mixin columns. [#5469]

- Support generalized value formatting for mixin columns in tables. [#5274]

- Support persistence of table indices when pickling and copying table. [#5468]

astropy.tests
^^^^^^^^^^^^^

- Install both runtime and test dependencies when running the
  ./setup.py test command. These dependencies are specified by the
  install_requires and tests_require keywords via setuptools. [#5092]

- Enable easier subclassing of the TestRunner class. [#5505]

astropy.time
^^^^^^^^^^^^

- ``light_travel_time`` can now use more accurate JPL ephemerides. [#5273, #5436]

astropy.units
^^^^^^^^^^^^^

- Added ``pixel_scale`` and ``plate_scale`` equivalencies. [#4987]

- The ``spectral_density`` equivalency now supports transformations of
  luminosity density. [#5151]

- ``Quantity`` now accepts strings consisting of a number and unit such
  as '10 km/s'. [#5245]

astropy.utils
^^^^^^^^^^^^^

- Added a new decorator: ``deprecated_renamed_argument``. This can be used to
  rename a function argument, while it still allows for the use of the older
  argument name. [#5214]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Added a ``make_lupton_rgb`` function to generate color images from three
  greyscale images, following the algorithm of Lupton et al. (2004). [#5535]

- Added ``data`` and ``interval`` inputs to the ``ImageNormalize``
  class. [#5206]

- Added a new ``simple_norm`` convenience function. [#5206]

- Added a default stretch for the ``Normalization`` class. [#5206].

- Added a default ``vmin/vmax`` for the ``ManualInterval`` class.
  [#5206].

- The ``wcsaxes`` subpackage has now been integrated in astropy as
  ``astropy.visualization.wcsaxes``.  This allows plotting of astronomical
  data/coordinate systems in Matplotlib. [#5496]

astropy.wcs
^^^^^^^^^^^

- Improved ``footprint_to_file``: allow to specify the coordinate system, and
  use by default the one from ``RADESYS``. Overwrite the file instead of
  appending to it. [#5494]


API Changes
-----------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- ``discretize_model`` now raises an exception if non-integer ranges are used.
  Previously it had incorrect behavior but did not raise an exception. [#5538]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- ``SkyCoord``, ``ICRS``, and other coordinate objects, as well as the
  underlying representations such as ``SphericalRepresentation`` and
  ``CartesianRepresentation`` can now be reshaped using methods named like the
  numpy ones for ``ndarray`` (``reshape``, ``swapaxes``, etc.)
  [#4123, #5254, #5482]

- The ``obsgeoloc`` and ``obsgeovel`` attributes of ``GCRS`` and
  ``PrecessedGeocentric`` frames are now stored and returned as
  ``CartesianRepresentation`` objects, rather than ``Quantity`` objects.
  Similarly, ``EarthLocation.get_gcrs_posvel`` now returns a tuple of
  ``CartesianRepresentation`` objects. [#5253]

- ``search_around_3d`` and ``search_around_sky`` now return units
  for the distance matching their input argument when no match is
  found, instead of ``dimensionless_unscaled``. [#5528]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- ASCII writers now accept an 'overwrite' argument.
  The default behavior is changed so that a warning will be
  issued when overwriting an existing file unless ``overwrite=True``.
  In a future version this will be changed from a warning to an
  exception to prevent accidentally overwriting a file. [#5007]

- The default representation of masked values when writing tables was
  changed from ``'--'`` to the empty string ``''``.  Previously any
  user-supplied ``fill_values`` parameter would overwrite the class
  default, but now the values are prepended to the class default. [#5347]

astropy.io.fits
^^^^^^^^^^^^^^^

- The old ``Header`` interface, deprecated since Astropy 0.1 (PyFITS 3.1), has
  been removed entirely. See :ref:`header-transition-guide` for explanations
  on this change and help on the transition. [#5310]

- The following functions, classes and methods have been removed:
  ``CardList``, ``Card.key``, ``Card.cardimage``, ``Card.ascardimage``,
  ``create_card``, ``create_card_from_string``, ``upper_key``,
  ``Header.ascard``, ``Header.rename_key``, ``Header.get_history``,
  ``Header.get_comment``, ``Header.toTxtFile``, ``Header.fromTxtFile``,
  ``new_table``, ``tdump``, ``tcreate``, ``BinTableHDU.tdump``,
  ``BinTableHDU.tcreate``.

- Removed ``txtfile`` argument to the ``Header`` constructor.

- Removed usage of ``Header.update`` with ``Header.update(keyword, value,
  comment)`` arguments.

- Removed ``startColumn`` and ``endColumn`` arguments to the ``FITS_record``
  constructor.

- The ``clobber`` argument in FITS writers has been renamed to
  ``overwrite``. This change affects the following functions and
  methods: ``tabledump``, ``writeto``, ``Header.tofile``,
  ``Header.totextfile``, ``_BaseDiff.report``,
  ``_BaseHDU.overwrite``, ``BinTableHDU.dump`` and
  ``HDUList.writeto``. [#5171]

- Added an optional ``copy`` parameter to ``fits.Header`` which controls if
  a copy is made when creating an ``Header`` from another ``Header``.
  [#5005, #5326]

astropy.io.registry
^^^^^^^^^^^^^^^^^^^

- ``.fts`` and ``.fts.gz`` files will be automatically identified as
  ``io.fits`` files if no explicit ``format`` is given. [#5211]

- Added an optional ``readwrite`` parameter for ``get_formats`` to filter
  formats for read or write. [#5275]

astropy.modeling
^^^^^^^^^^^^^^^^

- ``Gaussian2D`` now raises an error if ``theta`` is set at the same time as
  ``cov_matrix`` (previously ``theta`` was silently ignored). [#5537]

astropy.table
^^^^^^^^^^^^^

- Setting an existing table column (e.g. ``t['a'] = [1, 2, 3]``) now defaults
  to *replacing* the column with a column corresponding to the new value
  (using ``t.replace_column()``) instead of doing an in-place update.  Any
  existing meta-data in the column (e.g. the unit) is discarded.  An
  in-place update is still done when the new value is not a valid column,
  e.g. ``t['a'] = 0``.  To force an in-place update use the pattern
  ``t['a'][:] = [1, 2, 3]``. [#5556]

- Allow ``collections.Mapping``-like ``data`` attribute when initializing a
  ``Table`` object (``dict``-like was already possible). [#5213]

astropy.tests
^^^^^^^^^^^^^

- The inputs to the ``TestRunner.run_tests()`` method now must be
  keyword arguments (no positional arguments).  This applies to the
  ``astropy.test()`` function as well. [#5505]

astropy.utils
^^^^^^^^^^^^^

- Renamed ``ignored`` context manager in ``compat.misc`` to ``suppress``
  to be consistent with https://bugs.python.org/issue19266 . [#5003]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Deprecated the ``scale_image`` function. [#5206]

- The ``mpl_normalize`` module (containing the ``ImageNormalize``
  class) is now automatically imported with the ``visualization``
  subpackage. [#5491]

astropy.vo
^^^^^^^^^^

- The ``clobber`` argument in ``VOSDatabase.to_json()`` has been
  renamed to ``overwrite``. [#5171]

astropy.wcs
^^^^^^^^^^^

- ``wcs.rotateCD()`` was deprecated without a replacement. [#5240]

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Transformations between CIRS and AltAz now correctly account for the
  location of the observer. [#5591]

- GCRS frames representing a location on Earth with multiple obstimes are now
  allowed. This means that the solar system routines ``get_body``,
  ``get_moon`` and ``get_sun`` now work with non-scalar times and a
  non-geocentric observer. [#5253]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fix issue with units or other astropy core classes stored in table meta.
  [#5605]

astropy.io.fits
^^^^^^^^^^^^^^^

- Copying a ``fits.Header`` using ``copy`` or ``deepcopy`` from the ``copy``
  module will use ``Header.copy`` to ensure that modifying the copy will
  not alter the other original Header and vice-versa. [#4990, #5323]

- ``HDUList.info()`` no longer raises ``AttributeError`` in presence of
  ``BZERO``. [#5508]

- Avoid exceptions with numpy 1.10 and up when using scaled integer data
  where ``BZERO`` has float type but integer value. [#4639, #5527]

- Converting a header card to a string now calls ``self.verify('fix+warn')``
  instead of ``self.verify('fix')`` so headers with invalid keywords will
  not raise a ``VerifyError`` on printing. [#887,#5054]

- ``FITS_Record._convert_ascii`` now converts blank fields to 0 when a
  non-blank null column value is set. [#5134, #5394]

astropy.io.registry
^^^^^^^^^^^^^^^^^^^

- ``read`` now correctly raises an IOError if a file with an unknown
  extension can't be found, instead of raising IORegistryError:
  "Format could not be identified." [#4779]

astropy.time
^^^^^^^^^^^^

- Ensure ``Time`` instances holding a single ``delta_ut1_utc`` can be copied,
  flattened, etc. [#5225]

astropy.units
^^^^^^^^^^^^^

- Operations involving ``Angle`` or ``Distance``, or any other
  ``SpecificTypeQuantity`` instance, now also keep return an instance of the
  same type if the instance was the second argument (if the resulting unit
  is consistent with the specific type). [#5327]

- Inplace operations on ``Angle`` and ``Distance`` instances now raise an
  exception if the final unit is not equivalent to radian and meter, resp.
  Similarly, views as ``Angle`` and ``Distance`` can now only be taken
  from quantities with appropriate units, and views as ``Quantity`` can only
  be taken from logarithmic quanties such as ``Magnitude`` if the physical
  unit is dimensionless. [#5070]

- Conversion from quantities to logarithmic units now correctly causes a
  logarithmic quantity such as ``Magnitude`` to be returned. [#5183]


astropy.wcs
^^^^^^^^^^^

- SIP distortion for an alternate WCS is correctly initialized now by
  looking at the "CTYPE" values matching the alternate WCS. [#5443]

Other Changes and Additions
---------------------------

- The bundled ERFA was updated to version 1.3.0.  This includes the
  leap second planned for 2016 Dec 31.

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Initialization of ``Angle`` has been sped up for ``Quantity`` and ``Angle``
  input. [#4970]

- The use of ``np.matrix`` instances in the transformations has been
  deprecated, since this class does not allow stacks of matrices.  As a
  result, the semi-public functions ``angles.rotation_matrix`` and
  ``angles.angle_axis`` are also deprecated, in favour of the new routines
  with the same name in ``coordinates.matrix_utilities``. [#5104]

- A new ``BaseCoordinateFrame.cache`` dictionary has been created to expose
  the internal cache. This is useful when modifying representation data
  in-place without using ``realize_frame``. Additionally, documentation for
  in-place operations on coordinates were added. [#5575]

- Coordinates and their representations are printed with a slightly different
  format, following how numpy >= 1.12 prints structured arrays. [#5423]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- The default cosmological model has been changed to Planck 2015,
  and the citation strings have been updated. [#5372]

astropy.extern
^^^^^^^^^^^^^^

- Updated the bundled ``six`` module to version 1.10.0. [#5521]

- Updated the astropy shipped version of ``PLY`` to version 3.9. [#5526]

- Updated the astropy shipped version of jQuery to v3.3.1, and dataTables
  to v1.10.12. [#5564]

astropy.io.fits
^^^^^^^^^^^^^^^

- Performance improvements for tables with many columns. [#4985]

- Removed obsolete code that was previously needed to properly
  implement the append mode. [#4793]

astropy.io.registry
^^^^^^^^^^^^^^^^^^^

- Reduced the time spent in the ``get_formats`` function. This also reduces
  the time it takes to import astropy subpackages, i.e.
  ``astropy.coordinates``. [#5262]

astropy.units
^^^^^^^^^^^^^

- The functions ``add_enabled_units``, ``set_enabled_equivalencies`` and
  ``add_enabled_equivalencies`` have been sped up by copying the current
  ``_UnitRegistry`` instead of building it from scratch. [#5306]

- To build the documentation, the ``build_sphinx`` command has been deprecated
  in favor of ``build_docs``. [#5179]

- The ``--remote-data`` option to ``python setup.py test`` can now take
  different arguments: ``--remote-data=none`` is the same as not specifying
  ``--remote-data`` (skip all tests that require the internet),
  ``--remote-data=astropy`` skips all tests that need remote data except those
  that require only data from data.astropy.org, and ``--remote-data=any`` is
  the same as ``--remote-data`` (run all tests that use remote data). [#5506]

- The pytest ``recwarn`` fixture has been removed from the tests in favor of
  ``utils.catch_warnings``. [#5489]

- Deprecated escape sequences in strings (Python 3.6) have been removed. [#5489]


Version 1.2.2 (2016-12-22)
==========================

Bug Fixes
---------

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fix a bug where the ``fill_values`` parameter was ignored when writing a
  table to HTML format. [#5379]

astropy.io.fits
^^^^^^^^^^^^^^^

- Handle unicode FITS BinTable column names on Python 2 [#5204, #4805]

- Fix reading of float values from ASCII tables, that could be read as
  float32 instead of float64 (with the E and F formats). These values are now
  always read as float64. [#5362]

- Fixed memoryleak when using the compression module. [#5399, #5464]

- Able to insert and remove lower case HIERARCH keywords in a consistent
  manner [#5313, #5321]

astropy.stats
^^^^^^^^^^^^^

- Fixed broadcasting in ``sigma_clip`` when using negative ``axis``. [#4988]

astropy.table
^^^^^^^^^^^^^

- Assigning a logarithmic unit to a ``QTable`` column that did not have a
  unit yet now correctly turns it into the appropriate function quantity
  subclass (such as ``Magnitude`` or ``Dex``). [#5345]

- Fix default value for ``show_row_index`` in ``Table.show_in_browser``.
  [#5562]

astropy.units
^^^^^^^^^^^^^

- For inverse trig functions that operate on quantities, catch any warnings
  that occur from evaluating the function on the unscaled quantity value
  between __array_prepare__ and __array_wrap__. [#5153]

- Ensure ``!=`` also works for function units such as ``MagUnit`` [#5345]

astropy.wcs
^^^^^^^^^^^

- Fix use of the ``relax`` keyword in ``to_header`` when used to change the
  output precision. [#5164]

- ``wcs.to_header(relax=True)`` adds a "-SIP" suffix to ``CTYPE`` when SIP
  distortion is present in the WCS object. [#5239]

- Improved log messages in ``to_header``. [#5239]

Other Changes and Additions
---------------------------

- The bundled ERFA was updated to version 1.3.0.  This includes the
  leap second planned for 2016 Dec 31.

astropy.stats
^^^^^^^^^^^^^

- ``poisson_conf_interval`` with ``'kraft-burrows-nousek'`` interval is now
  faster and usable with SciPy versions < 0.14. [#5064, #5290]



Version 1.2.1 (2016-06-22)
==========================

Bug Fixes
---------

astropy.io.fits
^^^^^^^^^^^^^^^

- Fixed a bug that caused TFIELDS to not be in the correct position in
  compressed image HDU headers under certain circumstances, which created
  invalid FITS files. [#5118, #5125]

astropy.units
^^^^^^^^^^^^^

- Fixed an  ``ImportError`` that occurred whenever ``astropy.constants`` was
  imported before ``astropy.units``. [#5030, #5121]

- Magnitude zero points used to define ``STmag``, ``ABmag``, ``M_bol`` and
  ``m_bol`` are now collected in ``astropy.units.magnitude_zero_points``.
  They are not enabled as regular units by default, but can be included
  using ``astropy.units.magnitude_zero_points.enable()``. This makes it
  possible to round-trip magnitudes as originally intended.  [#5030]

Version 1.2 (2016-06-19)
========================

General
-------

- Astropy now requires Numpy 1.7.0 or later. [#4784]

New Features
------------

astropy.constants
^^^^^^^^^^^^^^^^^

- Add ``L_bol0``, the luminosity corresponding to absolute bolometric
  magnitude zero. [#4262]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- ``CartesianRepresentation`` now includes a transform() method that can take
  a 3x3 matrix to transform coordinates. [#4860]

- Solar system and lunar ephemerides accessible via ``get_body``,
  ``get_body_barycentric`` and ``get_moon`` functions. [#4890]

- Added astrometric frames (i.e., a frame centered on a particular
  point/object specified in another frame). [#4909, #4941]

- Added ``SkyCoord.spherical_offsets_to`` method. [#4338]

- Recent Earth rotation (IERS) data are now auto-downloaded so that AltAz
  transformations for future dates now use the most accurate available
  rotation values. [#4436]

- Add support for heliocentric coordinate frames. [#4314]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- ``angular_diameter_distance_z1z2`` now supports the computation of
  the angular diameter distance between a scalar and an array like
  argument. [#4593] The method now supports models with negative
  Omega_k0 (positive curvature universes) [#4661] and allows z2 < z1.

astropy.io.ascii
^^^^^^^^^^^^^^^^

- File name could be passed as ``Path`` object. [#4606]

- Check that columns in ``formats`` specifier exist in the output table
  when writing. [#4508, #4511]

- Allow trailing whitespace in the IPAC header lines. [#4758]

- Updated to filter out the default parser warning of BeautifulSoup.
  [#4551]

- Added support for reading and writing reStructuredText simple tables.
  [#4812]

astropy.io.fits
^^^^^^^^^^^^^^^

- File name could be passed as ``Path`` object. [#4606]

- Header allows a dictionary-like cards argument during creation. [#4663]

- New function ``convenience.table_to_hdu`` to allow creating a FITS
  HDU object directly from an astropy ``Table``. [#4778]

- New optional arguments ``ignore_missing`` and ``remove_all`` are added
  to ``astropy.io.fits.header.remove()``. [#5020]

astropy.io.registry
^^^^^^^^^^^^^^^^^^^

- Added custom ``IORegistryError``. [#4833]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- File name could be passed as ``Path`` object. [#4606]

astropy.modeling
^^^^^^^^^^^^^^^^

- Added the fittable=True attribute to the Scale and Shift models with tests. [#4718]

- Added example plots to docstrings for some built-in models. [#4008]

astropy.nddata
^^^^^^^^^^^^^^

- ``UnknownUncertainty`` new subclass of ``NDUncertainty`` that can be used to
  save uncertainties that cannot be used for error propagation. [#4272]

- ``NDArithmeticMixin``: ``add``, ``subtract``, ``multiply`` and ``divide``
  can be used as classmethods but require that two operands are given. These
  operands don't need to be NDData instances but they must be convertible to
  NDData. This conversion is done internally. Using it on the instance does
  not require (but also allows) two operands. [#4272, #4851]

- ``NDDataRef`` new subclass that implements ``NDData`` together with all
  currently available mixins. This class does not implement additional
  attributes, methods or a numpy.ndarray-like interface like ``NDDataArray``.
  attributes, methods or a numpy.ndarray-like interface like ``NDDataArray``.
  [#4797]

astropy.stats
^^^^^^^^^^^^^

- Added ``axis`` keyword for ``mad_std`` function. [#4688, #4689]

- Added Bayesian and Akaike Information Criteria. [#4716]

- Added Bayesian upper limits for Poisson count rates. [#4622]

- Added ``circstats``; a module for computing circular statistics. [#3705, #4472]

- Added ``jackknife`` resampling method. [#3708, #4439]

- Updated ``bootstrap`` to allow bootstrapping statistics with multiple
  outputs. [#3601]

- Added ``LombScargle`` class to compute Lomb-Scargle periodograms [#4811]

astropy.table
^^^^^^^^^^^^^

- ``Table.show_in_notebook`` and ``Table.show_in_browser(jsviewer=True)`` now
  yield tables with an "idx" column, allowing easy identification of the index
  of a row even when the table is re-sorted in the browser. [#4404]

- Added ``AttributeError`` when trying to set mask on non-masked table. [#4637]

- Allow to use a tuple of keys in ``Table.sort``.  [#4671]

- Added ``itercols``; a way to iterate through columns of a table. [#3805,
  #4888]

- ``Table.show_in_notebook`` and the default notebook display (i.e.,
  ``Table._repr_html_``) now use consistent table styles which can be set
  using the ``astropy.table.default_notebook_table_class`` configuration
  item. [#4886]

- Added interface to create ``Table`` directly from any table-like object
  that has an ``__astropy_table__`` method.  [#4885]

astropy.tests
^^^^^^^^^^^^^

- Enable test runner to obtain documentation source files from directory
  other than "docs". [#4748]

astropy.time
^^^^^^^^^^^^

- Added caching of scale and format transformations for improved performance.
  [#4422]

- Recent Earth rotation (IERS) data are now auto-downloaded so that UT1
  transformations for future times now work out of the box. [#4436]

- Add support for barycentric/heliocentric time corrections. [#4314]

astropy.units
^^^^^^^^^^^^^

- The option to use tuples to indicate fractional powers of units,
  deprecated in 0.3.1, has been removed. [#4449]

- Added slug to imperial units. [#4670]

- Added Earth radius (``R_earth``) and Jupiter radius (``R_jup``) to units.
  [#4818]

- Added a ``represents`` property to allow access to the definition of a
  named unit (e.g., ``u.kpc.represents`` yields ``1000 pc``). [#4806]

- Add bolometric absolute and apparent magnitudes, ``M_bol`` and ``m_bol``.
  [#4262]

astropy.utils
^^^^^^^^^^^^^

- ``Path`` object could be passed to ``get_readable_fileobj``. [#4606]

- Implemented a generic and extensible way of merging metadata. [#4459]

- Added ``format_doc`` decorator which allows to replace and/or format the
  current docstring of an object. [#4242]

- Added a new context manager ``set_locale`` to temporarily set the
  current locale. [#4363]

- Added new IERS_Auto class to auto-download recent IERS (Earth rotation)
  data when required by coordinate or time transformations. [#4436]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Add zscale interval based on Numdisplay's implementation. [#4776]

API changes
-----------

astropy.config
^^^^^^^^^^^^^^

- The deprecated ``ConfigurationItem`` and ``ConfigAlias`` classes and the
  ``save_config``, ``get_config_items``, and ``generate_all_config_items``
  functions have now been removed. [#2767, #4446]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Removed compatibility layer for pre-v0.4 API. [#4447]

- Added ``copy`` keyword-only argument to allow initialization without
  copying the (possibly large) input coordinate arrays. [#4883]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Improve documentation of z validity range of cosmology objects [#4882, #4949]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Add a way to control HTML escaping when writing a table as an HTML file. [#4423]

astropy.io.fits
^^^^^^^^^^^^^^^

- Two optional boolean arguments ``ignore_missing`` and ``remove_all`` are
  added to ``Header.remove``. [#5020]

astropy.modeling
^^^^^^^^^^^^^^^^

- Renamed ``Redshift`` model to ``RedshiftScaleFactor``. [#3672]

- Inputs (``coords`` and ``out``) to ``render`` function in ``Model`` are
  converted to float. [#4697]

- ``RotateNative2Celestial`` and ``RotateCelestial2Native`` are now
  implemented as subclasses of ``EulerAngleRotation``. [#4881, #4940]

astropy.nddata
^^^^^^^^^^^^^^

- ``NDDataBase`` does not set the private uncertainty property anymore. This
  only affects you if you subclass ``NDDataBase`` directly. [#4270]

- ``NDDataBase``: the ``uncertainty``-setter is removed. A similar one is
  added in ``NDData`` so this also only affects you if you subclassed
  ``NDDataBase`` directly. [#4270]

- ``NDDataBase``: ``uncertainty``-getter returns ``None`` instead of the
  private uncertainty and is now abstract. This getter is moved to
  ``NDData`` so it only affects direct subclasses of ``NDDataBase``. [#4270]

- ``NDData`` accepts a Quantity-like data and an explicitly given unit.
  Before a ValueError was raised in this case. The final instance will use the
  explicitly given unit-attribute but doesn't check if the units are
  convertible and the data will not be scaled. [#4270]

- ``NDData`` : the given mask, explicit or implicit if the data was masked,
  will be saved by the setter. It will not be saved directly as the private
  attribute. [#4879]

- ``NDData`` accepts an additional argument ``copy`` which will copy every
  parameter before it is saved as attribute of the instance. [#4270]

- ``NDData``: added an ``uncertainty.getter`` that returns the private
  attribute. It is equivalent to the old ``NDDataBase.uncertainty``-getter.
  [#4270]

- ``NDData``: added an ``uncertainty.setter``. It is slightly modified with
  respect to the old ``NDDataBase.uncertainty``-setter. The changes include:

- if the uncertainty has no uncertainty_type an info message is printed
  instead of a TypeError and the uncertainty is saved as
  ``UnknownUncertainty`` except the uncertainty is None. [#4270]

- the requirement that the uncertainty_type of the uncertainty needs to be a
  string was removed. [#4270]

- if the uncertainty is a subclass of NDUncertainty the parent_nddata
  attribute will be set so the uncertainty knows to which data it belongs.
  This is also a Bugfix. [#4152, #4270]

- ``NDData``: added a ``meta``-getter, which will set and return an empty
  OrderedDict if no meta was previously set. [#4509, #4469]

- ``NDData``: added an ``meta``-setter. It requires that the meta is
  dictionary-like (it also accepts Headers or ordered dictionaries and others)
  or None. [#4509, #4469, #4921]

- ``NDArithmeticMixin``: The operand in arithmetic methods (``add``, ...)
  doesn't need to be a subclass of ``NDData``. It is sufficient if it can be
  converted to one. This conversion is done internally. [#4272]

- ``NDArithmeticMixin``: The arithmetic methods allow several new arguments to
  control how or if different attributes of the class will be processed during
  the operation. [#4272]

- ``NDArithmeticMixin``: Giving the parameter ``propagate_uncertainties`` as
  positional keyword is deprecated and will be removed in the future. You now
  need to specify it as keyword-parameter. Besides ``True`` and ``False`` also
  ``None`` is now a valid value for this parameter. [#4272, #4851]

- ``NDArithmeticMixin``: The wcs attribute of the operands is not compared and
  thus raises no ValueError if they differ, except if a ``compare_wcs``
  parameter is specified. [#4272]

- ``NDArithmeticMixin``: The arithmetic operation was split from a general
  ``_arithmetic`` method to different specialized private methods to allow
  subclasses more control on how the attributes are processed without
  overriding ``_arithmetic``. The ``_arithmetic`` method is now used to call
  these other methods. [#4272]

- ``NDSlicingMixin``: If the attempt at slicing the mask, wcs or uncertainty
  fails with a ``TypeError`` a Warning is issued instead of the TypeError. [#4271]

- ``NDUncertainty``: ``support_correlated`` attribute is deprecated in favor of
  ``supports_correlated`` which is a property. Also affects
  ``StdDevUncertainty``. [#4272]

- ``NDUncertainty``: added the ``__init__`` that was previously implemented in
  ``StdDevUncertainty`` and takes an additional ``unit`` parameter. [#4272]

- ``NDUncertainty``: added a ``unit`` property without setter that returns the
  set unit or if not set the unit of the parent. [#4272]

- ``NDUncertainty``: included a ``parent_nddata`` property similar to the one
  previously implemented in StdDevUncertainty. [#4272]

- ``NDUncertainty``: added an ``array`` property with setter. The setter will
  convert the value to a plain numpy array if it is a list or a subclass of a
  numpy array. [#4272]

- ``NDUncertainty``: ``propagate_multiply`` and similar were removed. Before
  they were abstract properties and replaced by methods with the same name but
  with a leading underscore. The entry point for propagation is a method
  called ``propagate``. [#4272]

- ``NDUncertainty`` and subclasses: implement a representation (``__repr__``).
  [#4787]

- ``StdDevUncertainty``: error propagation allows an explicitly given
  correlation factor, which may be a scalar or an array which will be taken
  into account during propagation.
  This correlation must be determined manually and is not done by the
  uncertainty! [#4272]

- ``StdDevUncertainty``: the ``array`` is converted to a plain numpy array
  only if it's a list or a subclass of numpy.ndarray. Previously it was always
  cast to a numpy array but also allowed subclasses. [#4272]

- ``StdDevUncertainty``: setting the ``parent_nddata`` does not compare if the
  shape of it's array is identical to the parents data shape. [#4272]

- ``StdDevUncertainty``: the ``array.setter`` doesn't compare if the array has
  the same shape as the parents data. [#4272]

- ``StdDevUncertainty``: deprecated ``support_correlated`` in favor of
  ``supports_correlated``. [#4272, #4828]

- ``StdDevUncertainty``: deprecated ``propagate_add`` and similar methods in
  favor of ``propagate``. [#4272, #4828]

- Allow ``data`` to be a named argument in ``NDDataArray``. [#4626]

astropy.table
^^^^^^^^^^^^^

- ``operations.unique`` now has a ``keep`` parameter, which allows
  one to select whether to keep the first or last row in a set of
  duplicate rows, or to remove all rows that are duplicates. [#4632]

- ``QTable`` now behaves more consistently by making columns act as a
  ``Quantity`` even if they are assigned a unit after the table is
  created. [#4497, #4884]

astropy.units
^^^^^^^^^^^^^

- Remove deprecated ``register`` argument for Unit classes. [#4448]

astropy.utils
^^^^^^^^^^^^^

- The astropy.utils.compat.argparse module has now been deprecated. Use the
  Python 'argparse' module directly instead. [#4462]

- The astropy.utils.compat.odict module has now been deprecated. Use the
  Python 'collections' module directly instead. [#4466]

- The astropy.utils.compat.gzip module has now been deprecated. Use the
  Python 'gzip' module directly instead. [#4464]

- The deprecated ``ScienceStateAlias`` class has been removed. [#2767, #4446]

- The astropy.utils.compat.subprocess module has now been deprecated. Use the
  Python 'subprocess' module instead. [#4483]

- The astropy.utils.xml.unescaper module now also unescapes ``'%2F'`` to
  ``'/'`` and ``'&&'`` to ``'&'`` in a given URL. [#4699]

- The astropy.utils.metadata.MetaData descriptor has now two optional
  parameters: doc and copy. [#4921]

- The default IERS (Earth rotation) data now is now auto-downloaded via a
  new class IERS_Auto.  When extrapolating UT1-UTC or polar motion values
  outside the available time range, the values are now clipped at the last
  available value instead of being linearly extrapolated. [#4436]

astropy.wcs
^^^^^^^^^^^

- WCS objects can now be initialized with an ImageHDU or
  PrimaryHDU object. [#4493, #4505]

- astropy.wcs now issues an INFO message when the header has SIP coefficients but
  "-SIP" is missing from CTYPE. [#4814]

Bug fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Ameliorate a problem with ``get_sun`` not round-tripping due to
  approximations in the light deflection calculation. [#4952]

- Ensure that ``angle_utilities.position_angle`` accepts floats, as stated
  in the docstring. [#3800]

- Ensured that transformations for ``GCRS`` frames are correct for
  non-geocentric observers. [#4986]

- Fixed a problem with the ``Quantity._repr_latex_`` method causing errors
  when showing an ``EarthLocation`` in a Jupyter notebook. [#4542, #5068]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fix a problem where the fast reader (with use_fast_converter=False) can
  fail on non-US locales. [#4363]

- Fix astropy.io.ascii.read handling of units for IPAC formatted files.
  Columns with no unit are treated as unitless not dimensionless.
  [#4867, #4947]

- Fix problems the header parsing in the sextractor reader. [#4603, #4910]

astropy.io.fits
^^^^^^^^^^^^^^^

- ``GroupsHDU.is_image`` property is now set to ``False``. [#4742]

- Ensure scaling keywords are removed from header when unsigned integer data
  is converted to signed type. [#4974, #5053]

- Made TFORMx keyword check more flexible in test of compressed images to
  enable compatibility of the test with cfitsio 3.380. [#4646, #4653]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- The astropy.io.votable.validator.html module is updated to handle division
  by zero when generating validation report. [#4699]

- KeyError when converting Table v1.2 numeric arrays fixed. [#4782]

astropy.modeling
^^^^^^^^^^^^^^^^

- Refactored ``AiryDisk2D``, ``Sersic1D``, and ``Sersic2D`` models
  to be able to combine them as classes as well as instances. [#4720]

- Modified the "LevMarLSQFitter" class to use the weights in the
  calculation of the Jacobian. [#4751]

astropy.nddata
^^^^^^^^^^^^^^

- ``NDData`` giving masked_Quantities as data-argument will use the
  implicitly passed mask, unit and value. [#4270]

- ``NDData`` using a subclass implementing ``NDData`` with
  ``NDArithmeticMixin`` now allows error propagation. [#4270]

- Fixed memory leak that happened when uncertainty of ``NDDataArray`` was
  set. [#4825, #4862]

- ``StdDevUncertainty``: During error propagation the unit of the uncertainty
  is taken into account. [#4272]

- ``NDArithmeticMixin``: ``divide`` and ``multiply`` yield correct
  uncertainties if only one uncertainty is set. [#4152, #4272]

astropy.stats
^^^^^^^^^^^^^

- Fix ``sigma_clipped_stats`` to use the ``axis`` argument. [#4726, #4808]

astropy.table
^^^^^^^^^^^^^

- Fixed bug where Tables created from existing Table objects were not
  inheriting the ``primary_key`` attribute. [#4672, #4930]

- Provide more detail in the error message when reading a table fails due to a
  problem converting column string values. [#4759]

astropy.units
^^^^^^^^^^^^^

- Exponentiation using a ``Quantity`` with a unit equivalent to dimensionless
  as base and an ``array``-like exponent yields the correct result. [#4770]

- Ensured that with ``spectral_density`` equivalency one could also convert
  between ``photlam`` and ``STmag``/``ABmag``. [#5017]

astropy.utils
^^^^^^^^^^^^^

- The astropy.utils.compat.fractions module has now been deprecated. Use the
  Python 'fractions' module directly instead. [#4463]

- Added ``format_doc`` decorator which allows to replace and/or format the
  current docstring of an object. [#4242]

- Attributes using the astropy.utils.metadata.MetaData descriptor are now
  included in the sphinx documentation. [#4921]

astropy.vo
^^^^^^^^^^

- Relaxed expected accuracy of Cone Search prediction test to reduce spurious
  failures. [#4382]

astropy.wcs
^^^^^^^^^^^

- astropy.wcs.to_header removes "-SIP" from CTYPE when SIP coefficients
  are not written out, i.e. ``relax`` is either ``False`` or ``None``.
  astropy.wcs.to_header appends "-SIP" to CTYPE when SIP coefficients
  are written out, i.e. ``relax=True``. [#4814]

- Made ``wcs.bounds_check`` call ``wcsprm_python2c``, which means it
  works even if ``wcs.set`` has not been called yet. [#4957, #4966].

- WCS objects can no longer be reverse-indexed, which was technically
  permitted but incorrectly implemented previously [#4962]

Other Changes and Additions
---------------------------

- Python 2.6 is no longer supported. [#4486]

- The bundled version of py.test has been updated to 2.8.3. [#4349]

- Reduce Astropy's import time (``import astropy``) by almost a factor 2. [#4649]

- Cython prerequisite for building changed to v0.19 in install.rst [#4705,
  #4710, #4719]

- All astropy.modeling functionality that was deprecated in Astropy 1.0 has
  been removed. [#4857]

- Added instructions for installing Astropy into CASA. [#4840]

- Added an example gallery to the docs demonstrating short
  snippets/examples. [#4734]


Version 1.1.2 (2016-03-10)
==========================

New Features
------------

astropy.wcs
^^^^^^^^^^^

- The ``astropy.wcs`` module now exposes ``WCSHDO_P*`` constants that can be
  used to allow more control over output precision when using the ``relax``
  keyword argument. [#4616]

Bug Fixes
---------

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fixed handling of CDS data file when no description is given and also
  included stripping out of markup for missing value from description. [#4437, #4474]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fixed possible segfault during error handling in FITS tile
  compression. [#4489]

- Fixed crash on pickling of binary table columns with the 'X', 'P', or
  'Q' format. [#4514]

- Fixed memory / reference leak that could occur when copying a ``FITS_rec``
  object (the ``.data`` for table HDUs). [#520]

- Fixed a memory / reference leak in ``FITS_rec`` that occurred in a wide
  range of cases, especially after writing FITS tables to a file, but in
  other cases as well. [#4539]

astropy.modeling
^^^^^^^^^^^^^^^^

- Fix a bug to allow instantiation of a modeling class having a parameter
  with a custom setter that takes two parameters ``(value, model)`` [#4656]

astropy.table
^^^^^^^^^^^^^

- Fixed bug when replacing a table column with a mixin column like
  Quantity or Time. [#4601]

- Disable initial ordering in jsviewer (``show_in_browser``,
  ``show_in_notebook``) to respect the order from the Table. [#4628]

astropy.units
^^^^^^^^^^^^^

- Fixed sphinx issues on plotting quantities. [#4527]

astropy.utils
^^^^^^^^^^^^^

- Fixed latex representation of function units. [#4563]

- The ``zest.releaser`` hooks included in Astropy are now injected locally to
  Astropy, rather than being global. [#4650]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Fixed ``fits2bitmap`` script to allow ext flag to contain extension
  names or numbers. [#4468]

- Fixed ``fits2bitmap`` default output filename generation for
  compressed FITS files. [#4468]

- Fixed ``quantity_support`` to ensure its conversion returns ndarray
  instances (needed for numpy >=1.10). [#4654]

astropy.wcs
^^^^^^^^^^^

- Fixed possible exception in handling of SIP headers that was introduced in
  v1.1.1. [#4492]

- Fixed a bug that caused WCS objects with a high dynamic range of values for
  certain parameters to lose precision when converted to a header. This
  occurred for example in cases of spectral cubes, where a spectral axis in
  Hz might have a CRVAL3 value greater than 1e10 but the spatial coordinates
  would have CRVAL1/2 values 8 to 10 orders of magnitude smaller. This bug
  was present in Astropy 1.1 and 1.1.1 but not 1.0.x. This has now been fixed
  by ensuring that all WCS keywords are output with 14 significant figures by
  default. [#4616]

Other Changes and Additions
---------------------------

- Updated bundled astropy-helpers to v1.1.2. [#4678]

- Updated bundled copy of WCSLIB to 5.14. [#4579]


Version 1.1.1 (2016-01-08)
==========================

New Features
------------

astropy.io.registry
^^^^^^^^^^^^^^^^^^^

- Allow ``pathlib.Path`` objects (available in Python 3.4 and later) for
  specifying the file name in registry read / write functions. [#4405]

astropy.utils
^^^^^^^^^^^^^

- ``console.human_file_size`` now accepts quantities with byte-equivalent
  units [#4373]

Bug Fixes
---------

astropy.analytic_functions
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed the blackbody functions' handling of overflows on some platforms
  (Windows with MSVC, older Linux versions) with a buggy ``expm1`` function.
  [#4393]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fixed an bug where updates to string columns in FITS tables were not saved
  on Python 3. [#4452]

Other Changes and Additions
---------------------------

- Updated bundled astropy-helpers to v1.1.1. [#4413]


Version 1.1 (2015-12-11)
========================

New Features
------------

astropy.config
^^^^^^^^^^^^^^

- Added new tools ``set_temp_config`` and ``set_temp_cache`` which can be
  used either as function decorators or context managers to temporarily
  use alternative directories in which to read/write the Astropy config
  files and download caches respectively.  This is especially useful for
  testing, though ``set_temp_cache`` may also be used as a way to provide
  an alternative (application specific) download cache for large data files,
  rather than relying on the default cache location in users' home
  directories. [#3975]

astropy.constants
^^^^^^^^^^^^^^^^^

- Added the Thomson scattering cross-section. [#3839]

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Added Moffat2DKernel. [#3965]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Added ``get_constellation`` function and ``SkyCoord.get_constellation``
  convenience method to determine the constellation that a coordinate
  is in. [#3758]

- Added ``PrecessedGeocentric`` frame, which is based on GCRS, but precessed
  to a specific requested mean equinox. [#3758]

- Added ``Supergalactic`` frame to support de Vaucouleurs supergalactic
  coordinates. [#3892]

- ``SphericalRepresentation`` now has a ``._unit_representation`` class attribute to specify
  an equivalent UnitSphericalRepresentation. This allows subclasses of
  representations to pair up correctly. [#3757]

- Added functionality to support getting the locations of observatories by
  name. See ``astropy.coordinates.EarthLocation.of_site``. [#4042]

- Added ecliptic coordinates, including ``GeocentricTrueEcliptic``,
  ``BarycentricTrueEcliptic``, and ``HeliocentricTrueEcliptic``. [#3749]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Add Planck 2015 cosmology [#3476]

- Distance calculations now > 20-40x faster for the supplied
  cosmologies due to implementing Cython scalar versions of
  ``FLRW.inv_efunc``.[#4127]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Automatically use ``guess=False`` when reading if the file ``format`` is
  provided and the format parameters are uniquely specified.  This update
  also removes duplicate format guesses to improve performance. [#3418]

- Calls to ascii.read() for fixed-width tables may now omit one of the keyword
  arguments ``col_starts`` or ``col_ends``. Columns will be assumed to begin and
  end immediately adjacent to each other. [#3657]

- Add a function ``get_read_trace()`` that returns a traceback of the
  attempted read formats for the last call to ``astropy.io.ascii.read``. [#3688]

- Supports LZMA decompression via ``get_readable_fileobj`` [#3667]

- Allow ``-`` character is Sextractor format column names. [#4168]

- Improve DAOphot reader to read multi-aperture files [#3535, #4207]

astropy.io.fits
^^^^^^^^^^^^^^^

- Support reading and writing from bzip2 compressed files. i.e. ``.fits.bz2``
  files. [#3789]

- Included a new command-line script called ``fitsinfo`` to display
  a summary of the HDUs in one or more FITS files. [#3677]

astropy.io.misc
^^^^^^^^^^^^^^^

- Support saving all meta information, description and units of tables and columns
  in HDF5 files [#4103]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- A new method was added to ``astropy.io.votable.VOTable``,
  ``get_info_by_id`` to conveniently find an ``INFO`` element by its
  ``ID`` attribute. [#3633]

- Instances in the votable tree now have better ``__repr__`` methods. [#3639]

astropy.logger.py
^^^^^^^^^^^^^^^^^

- Added log levels (e.g., DEBUG, INFO, CRITICAL) to ``astropy.log`` [#3947]

astropy.modeling
^^^^^^^^^^^^^^^^

- Added a new ``Parameter.validator`` interface for setting a validation
  method on individual model parameters.  See the ``Parameter``
  documentation for more details. [#3910]

- The projection classes that are named based on the 3-letter FITS
  WCS projections (e.g. ``Pix2Sky_TAN``) now have aliases using
  longer, more descriptive names (e.g. ``Pix2Sky_Gnomonic``).
  [#3583]

- All of the standard FITS WCS projection types have been
  implemented in ``astropy.modeling.projections`` (by wrapping
  WCSLIB). [#3906]

- Added ``Sersic1D`` and ``Sersic2D`` model classes. [#3889]

- Added the Voigt profile to existing models. [#3901]

- Added ``bounding_box`` property and ``render_model`` function [#3909]

astropy.nddata
^^^^^^^^^^^^^^

- Added ``block_reduce`` and ``block_replicate`` functions. [#3453]

- ``extract_array`` now offers different options to deal with array
  boundaries [#3727]

- Added a new ``Cutout2D`` class to create postage stamp image cutouts
  with optional WCS propagation. [#3823]

astropy.stats
^^^^^^^^^^^^^

- Added ``sigma_lower`` and ``sigma_upper`` keywords to
  ``sigma_clip`` to allow for non-symmetric clipping. [#3595]

- Added ``cenfunc``, ``stdfunc``, and ``axis`` keywords to
  ``sigma_clipped_stats``. [#3792]

- ``sigma_clip`` automatically masks invalid input values (NaNs, Infs) before
  performing the clipping [#4051]

- Added the ``histogram`` routine, which is similar to ``np.histogram`` but
  includes several additional options for automatic determination of optimal
  histogram bins. Associated helper routines include ``bayesian_blocks``,
  ``friedman_bin_width``, ``scott_bin_width``, and ``knuth_bin_width``.
  This functionality was ported from the astroML library. [#3756]

- Added the ``bayesian_blocks`` routine, which implements a dynamic algorithm
  for locating change-points in various time series. [#3756]

- A new function ``poisson_conf_interval()`` was added to allow easy calculation
  of several standard formulae for the error bars on the mean of a Poisson variable
  estimated from a single sample.

astropy.table
^^^^^^^^^^^^^

- ``add_column()`` and ``add_columns()`` now have ``rename_duplicate``
  option to rename new column(s) rather than raise exception when its name
  already exists. [#3592]

- Added ``Table.to_pandas`` and ``Table.from_pandas`` for converting to/from
  pandas dataframes. [#3504]

- Initializing a ``Table`` with ``Column`` objects no longer requires
  that the column ``name`` attribute be defined. [#3781]

- Added an ``info`` property to ``Table`` objects which provides configurable
  summary information about the table and its columns. [#3731]

- Added an ``info`` property to column classes (``Column`` or mixins).  This
  serves a dual function of providing configurable summary information about
  the column, and acting as a manager of column attributes such as
  name, format, or description. [#3731]

- Updated table and column representation to use the ``dtype_info_name``
  function for the dtype value.  Removed the default "masked=False"
  from the table representation. [#3868, #3869]

- Updated row representation to be consistent with the corresponding
  table representation for that row.  Added HTML representation so a
  row displays nicely in IPython notebook.

- Added a new table indexing engine allowing for the creation of
  indices on one or more columns of a table using ``add_index``. These
  indices enable new functionality such as searching for rows by value
  using ``loc`` and ``iloc``, as well as increased performance for
  certain operations. [#3915, #4202]

- Added capability to include a structured array or recarray in a table
  as a mixin column.  This allows for an approximation of nested tables.
  [#3925]

- Added ``keep_byteorder`` option to ``Table.as_array()``.  See the
  "API Changes" section below. [#4080]

- Added a new method ``Table.replace_column()`` to replace an existing
  column with a new data column. [#4090]

- Added a ``tableclass`` option to ``Table.pformat()`` to allow specifying
  a list of CSS classes added to the HTML table. [#4131]

- New CSS for jsviewer table [#2917, #2982, #4174]

- Added a new ``Table.show_in_notebook`` method that shows an interactive view
  of a Table (similar to ``Table.show_in_browser(jsviewer=True)``) in an
  Python/Jupyter notebook. [#4197]

- Added column alignment formatting for better pprint viewing
  experience. [#3644]

astropy.tests
^^^^^^^^^^^^^

- Added new test config options, ``config_dir`` and ``cache_dir``  (these
  can be edited in ``setup.cfg`` or as extra command-line options to
  py.test) for setting the locations to use for the Astropy config files
  and download caches (see also the related ``set_temp_config/cache``
  features added in ``astropy.config``). [#3975]

astropy.time
^^^^^^^^^^^^

- Add support for FITS standard time strings. [#3547]

- Allow the ``format`` attribute to be updated in place to change the
  default representation of a ``Time`` object. [#3673]

- Add support for shape manipulation (reshape, ravel, etc.). [#3224]

- Add argmin, argmax, argsort, min, max, ptp, sort methods. [#3681]

- Add ``Time.to_datetime`` method for converting ``Time`` objects to
  timezone-aware datetimes. [#4119, #4124]

astropy.units
^^^^^^^^^^^^^

- Added furlong to imperial units. [#3529]

- Added mil to imperial units. [#3716]

- Added stone to imperial units. [#4192]

- Added Earth Mass (``M_earth``) and Jupiter mass (``M_jup``) to units [#3907]

- Added support for functional units, in particular the logarithmic ones
  ``Magnitude``, ``Decibel``, and ``Dex``. [#1894]

- Quantities now work with the unit support in matplotlib.  See
  :ref:`plotting-quantities`. [#3981]

- Clarified imperial mass measurements and added pound force (lbf),
  kilopound (kip), and pound per square inch (psi). [#3409]

astropy.utils
^^^^^^^^^^^^^

- Added new ``OrderedDescriptor`` and ``OrderedDescriptorContainer`` utility
  classes that make it easier to implement classes with declarative APIs,
  wherein class-level attributes have an inherit "ordering" to them that is
  specified by the order in which those attributes are defined in the class
  declaration (by defining them using special descriptors that have
  ``OrderedDescriptor`` as a base class).  See the API documentation for
  these classes for more details. Coordinate frames and models now use this
  interface. [#3679]

- The ``get_pkg_data_*`` functions now take an optional ``package`` argument
  which allows specifying any package to read package data filenames or
  content out of, as opposed to only being able to use data from the package
  that the function is called from. [#4079]

- Added function ``dtype_info_name`` to the ``data_info`` module to provide
  the name of a ``dtype`` for human-readable informational purposes. [#3868]

- Added ``classproperty`` decorator--this is to ``property`` as
  ``classmethod`` is to normal instance methods. [#3982]

- ``iers.open`` now handles network URLs, as well as local paths. [#3850]

- The ``astropy.utils.wraps`` decorator now takes an optional
  ``exclude_args`` argument not shared by the standard library ``wraps``
  decorator (as it is unique to the Astropy version's ability of copying
  the wrapped function's argument signature).  ``exclude_args`` allows
  certain arguments on the wrapped function to be excluded from the signature
  of the wrapper function.  This is particularly useful when wrapping an
  instance method as a function (to exclude the ``self`` argument). [#4017]

- ``get_readable_fileobj`` can automatically decompress LZMA ('.xz')
  files using the ``lzma`` module of Python 3.3+ or, when available, the
  ``backports.lzma`` package on earlier versions. [#3667]

- The ``resolve_name`` utility now accepts any number of additional
  positional arguments that are automatically dotted together with the
  first ``name`` argument. [#4083]

- Added ``is_url_in_cache`` for resolving paths to cached files via URLS
  and checking if files exist. [#4095]

- Added a ``step`` argument to the ``ProgressBar.map`` method to give
  users control over the update frequency of the progress bar. [#4191]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Added a function / context manager ``quantity_support`` for enabling
  seamless plotting of ``Quantity`` instances in matplotlib. [#3981]

- Added the ``hist`` function, which is similar to ``plt.hist`` but
  includes several additional options for automatic determination of optimal
  histogram bins. This functionality was ported from the astroML library.
  [#3756]

astropy.wcs
^^^^^^^^^^^

- The included version of wcslib has been upgraded to 5.10. [#3992, #4239]

  The minimum required version of wcslib in the 4.x series remains 4.24.

  The minimum required version of wcslib in the 5.x series is
  5.8.  Building astropy against a wcslib 5.x prior to 5.8
  will raise an ``ImportError`` when ``astropy.wcs`` is imported.

  The wcslib changes relevant to astropy are:

- The FITS headers returned by ``astropy.wcs.WCS.to_header`` and
  ``astropy.wcs.WCS.to_header_string`` now include values with
  more precision.  This will result in numerical differences in
  your results if you convert ``astropy.wcs.WCS`` objects to FITS
  headers and use the results.

- ``astropy.wcs.WCS`` now recognises the ``TPV``, ``TPD``,
  ``TPU``, ``DSS``, ``TNX`` and ``ZPX`` polynomial distortions.

- Added relaxation flags to allow ``PC0i_0ja``, ``PV0j_0ma``, and
  ``PS0j_0ma`` (i.e. with leading zeroes on the index).

- Tidied up error reporting, particularly relating to translating
  status returns from lower-level functions.

- Changed output formatting of floating point values in
  ``to_header``.

- Enhanced text representation of ``WCS`` objects. [#3604]

- The ``astropy.tests.helper`` module is now part of the public API (and has a
  documentation page).  This module was in previous releases of astropy,
  but was not considered part of the public API until now. [#3890]

- There is a new function ``astropy.online_help`` to search the
  astropy documentation and display the result in a web
  browser. [#3642]

API changes
-----------

astropy.cosmology
^^^^^^^^^^^^^^^^^

- ``FLRW._tfunc`` and ``FLRW._xfunc`` are marked as deprecated.  Users
  should use the new public interfaces ``FLRW.lookback_time_integrand``
  and ``FLRW.abs_distance_integrand`` instead. [#3767]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- The default header line processing was made to be consistent with data line
  processing in that it now ignores blank lines that may have whitespace
  characters.  Any code that explicitly specifies a ``header_start`` value
  for parsing a file with blank lines in the header containing whitespace will
  need to be updated. [#2654]

astropy.io.fits
^^^^^^^^^^^^^^^

- The ``uint`` argument to ``fits.open`` is now True by default; that is,
  arrays using the FITS unsigned integer convention will be detected, and
  read as unsigned integers by default.  A new config option for
  ``io.fits``, ``enable_uint``, can be changed to False to revert to the
  original behavior of ignoring the ``uint`` convention unless it is
  explicitly requested with ``uint=True``. [#3916]

- The ``ImageHDU.NumCode`` and ``ImageHDU.ImgCode`` attributes (and same
  for other classes derived from ``_ImageBaseHDU``) are deprecated.  Instead,
  the ``astropy.io.fits`` module-level constants ``BITPIX2DTYPE`` and
  ``DTYPE2BITPIX`` can be used. [#3916]

astropy.modeling
^^^^^^^^^^^^^^^^

- Note: Comparisons of model parameters with array-like values now
  yields a Numpy boolean array as one would get with normal Numpy
  array comparison.  Previously this returned a scalar True or False,
  with True only if the comparison was true for all elements compared,
  which could lead to confusing circumstances. [#3912]

- Using ``model.inverse = None`` to reset a model's inverse to its
  default is deprecated.  In the future this syntax will explicitly make
  a model not have an inverse (even if it has a default).  Instead, use
  ``del model.inverse`` to reset a model's inverse to its default (if it
  has a default, otherwise this just deletes any custom inverse that has
  been assigned to the model and is still equivalent to setting
  ``model.inverse = None``). [#4236]

- Adds a ``model.has_user_inverse`` attribute which indicates whether or not
  a user has assigned a custom inverse to ``model.inverse``.  This is just
  for informational purposes, for example, for software that introspects
  model objects. [#4236]

- Renamed the parameters of ``RotateNative2Celestial`` and
  ``RotateCelestial2Native`` from ``phi``, ``theta``, ``psi`` to
  ``lon``, ``lat`` and ``lon_pole``. [#3578]

- Deprecated the ``Pix2Sky_AZP.check_mu`` and ``Sky2Pix_AZP.check_mu``
  methods (these were obscure "accidentally public" methods that were
  probably not used by anyone). [#3910]

- Added a phase parameter to the Sine1D model. [#3807]

astropy.stats
^^^^^^^^^^^^^

- Renamed the ``sigma_clip`` ``sig`` keyword as ``sigma``. [#3595]

- Changed the ``sigma_clip`` ``varfunc`` keyword to ``stdfunc``. [#3595]

- Renamed the ``sigma_clipped_stats`` ``mask_val`` keyword to
  ``mask_value``. [#3595]

- Changed the default ``iters`` keyword value to 5 in both the
  ``sigma_clip`` and ``sigma_clipped_stats`` functions. [#4067]

astropy.table
^^^^^^^^^^^^^

- ``Table.as_array()`` always returns a structured array with each column in
  the system's native byte order.  The optional ``keep_byteorder=True``
  option will keep each column's data in its original byteorder. [#4080]

- ``Table.simple_table()`` now creates tables with int64 and float64 types
  instead of int32 and float64. [#4114]

- An empty table can now be initialized without a ``names`` argument as long
  as a valid ``dtype`` argument (with names embedded) is supplied. [#3977]

astropy.time
^^^^^^^^^^^^

- The ``astropy_time`` attribute and time format has been removed from the
  public interface.  Existing code that instantiates a new time object using
  ``format='astropy_time'`` can simply omit the ``format``
  specification. [#3857]

astropy.units
^^^^^^^^^^^^^

- Single-item ``Quantity`` instances with record ``dtype`` will now have
  their ``isscalar`` property return ``True``, consistent with behaviour for
  numpy arrays, where ``np.void`` records are considered scalar. [#3899]

- Three changes relating to the FITS unit format [#3993]:

- The FITS unit format will no longer parse an arbitrary number as a
  scale value.  It must be a power of 10 of the form ``10^^k``,
  ``10^k``, ``10+k``, ``10-k`` and ``10(k)``. [#3993]

- Scales that are powers of 10 can be written out.  Previously, any
  non-1.0 scale was rejected.

- The ``*`` character is accepted as a separator between the scale
  and the units.

- Unit formatter classes now require the ``parse`` and ``to_string``
  methods are now required to be classmethods (and the formatter
  classes themselves are assumed to be singletons that are not
  instantiated).  As unit formatters are mostly an internal implementation
  detail this is not likely to affect any users. [#4001]

- CGS E&M units are now defined separately from SI E&M units, and have
  distinct physical types. [#4255, #4355]

astropy.utils
^^^^^^^^^^^^^

- All of the ``get_pkg_data_*`` functions take an optional ``package``
  argument as their second positional argument.  So any code that previously
  passed other arguments to these functions as positional arguments might
  break.  Use keyword argument passing instead to mitigate this. [#4079]

- ``astropy.utils.iers`` now uses a ``QTable`` internally, which means that
  the numerical columns are stored as ``Quantity``, with full support for
  units.  Furthermore, the ``ut1_utc`` method now returns a ``Quantity``
  instead of a float or an array (as did ``pm_xy`` already). [#3223]

- ``astropy.utils.iers`` now throws an ``IERSRangeError``, a subclass
  of ``IndexError``, rather than a raw ``IndexError``.  This allows more
  fine-grained catching of situations where a ``Time`` is beyond the range
  of the loaded IERS tables. [#4302]

astropy.wcs
^^^^^^^^^^^

- When compiled with wcslib 5.9 or later, the FITS headers returned
  by ``astropy.wcs.WCS.to_header`` and
  ``astropy.wcs.WCS.to_header_string`` now include values with more
  precision.  This will result in numerical differences in your
  results if you convert ``astropy.wcs.WCS`` objects to FITS headers
  and use the results.

- If NAXIS1 or NAXIS2 is not passed with the header object to
  WCS.calc_footprint, a ValueError is raised. [#3557]

Bug fixes
---------

astropy.constants
^^^^^^^^^^^^^^^^^

- The constants ``Ry`` and ``u`` are now properly used inside the
  corresponding units.  The latter have changed slightly as a result. [#4229]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Internally, ``coordinates`` now consistently uses the appropriate time
  scales for using ERFA functions. [#4302]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fix a segfault in the fast C parser when one of the column headers
  is empty [#3545].

- Fix several bugs that prevented the fast readers from being used
  when guessing the file format.  Also improved the read trace
  information to better understand format guessing. [#4115]

- Fix an underlying problem that resulted in an uncaught TypeError
  exception when reading a CDS-format file with guessing enabled. [#4120]

astropy.modeling
^^^^^^^^^^^^^^^^

- ``Simplex`` fitter now correctly passes additional keywords arguments to
  the scipy solver. [#3966]

- The keyword ``acc`` (for accuracy) is now correctly accepted by
  ``Simplex``. [#3966]

astropy.units
^^^^^^^^^^^^^

- The units ``Ryd`` and ``u`` are no longer hard-coded numbers, but depend
  on the appropriate values in the ``constants`` module.  As a result, these
  units now imply slightly different conversions.  [#4229]

Other Changes and Additions
---------------------------

- The ``./setup.py test`` command is now implemented in the ``astropy.tests``
  module again (previously its implementation had been moved into
  astropy-helpers).  However, that made it difficult to synchronize changes
  to the Astropy test runner with changes to the ``./setup.py test`` UI.
  astropy-helpers v1.1 and above will detect this implementation of the
  ``test`` command, when present, and use it instead of the old version that
  was included in astropy-helpers (most users will not notice any difference
  as a result of this change). [#4020]

- The repr for ``Table`` no longer displays ``masked=False`` since tables
  are not masked by default anyway. [#3869]

- The version of ``PLY`` that ships with astropy has been updated to 3.6.

- WCSAxes is now required for doc builds. [#4074]

- The migration guide from pre-v0.4 coordinates has been removed to avoid
  cluttering the ``astropy.coordinates`` documentation with increasingly
  irrelevant material.  To see the migration guide, we recommend you simply look
  to the archived documentation for previous versions, e.g.
  https://docs.astropy.org/en/v1.0/coordinates/index.html#migrating-from-pre-v0-4-coordinates
  [#4203]

- In ``astropy.coordinates``, the transformations between GCRS, CIRS,
  and ITRS have been adjusted to more logically reflect the order in
  which they actually apply.  This should not affect most coordinate
  transformations, but may affect code that is especially sensitive to
  machine precision effects that change when the order in which
  transformations occur is changed. [#4255]

- Astropy v1.1.0 will be the last release series to officially support
  Python 2.6.  A deprecation warning will now be issued when using Astropy
  in Python 2.6 (this warning can be disabled through the usual Python warning
  filtering mechanisms). [#3779]


Version 1.0.13 (2017-05-29)
===========================

Bug Fixes
---------

astropy.io.fits
^^^^^^^^^^^^^^^

- Fix use of quantize level parameter for ``CompImageHDU``. [#6029]

- Prevent crash when a header contains non-ASCII (e.g. UTF-8) characters, to
  allow fixing the problematic cards. [#6084]


Version 1.0.12 (2017-03-05)
===========================

Bug Fixes
---------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Fixed bug in ``discretize_integrate_2D`` in which x and y coordinates
  where swapped. [#5634]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Fixed a bug where ``get_transform`` could sometimes produce confusing errors
  because of a typo in the input validation. [#5645]

astropy.io.fits
^^^^^^^^^^^^^^^

- Guard against extremely unlikely problems in compressed images, which
  could lead to memory unmapping errors. [#5775]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Fixed a bug where stdlib ``realloc()`` was used instead of
  ``PyMem_Realloc()`` [#5696, #4739, #2100]

astropy.utils
^^^^^^^^^^^^^

- Fixed ImportError with NumPy < 1.7 and Python 3.x in
  ``_register_patched_dtype_reduce``. [#5848]


Version 1.0.11 (2016-12-22)
===========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Initialising a SkyCoord from a list containing a single SkyCoord no longer removes
  the distance from the coordinate. [#5270]

- Fix errors in the implementation of the conversion to and from FK4 frames
  without e-terms, which will have affected coordinates not on the unit
  sphere (i.e., with distances). [#4293]

- Fix bug where with cds units enabled it was no longer possible to initialize
  an ``Angle``. [#5483]

- Ensure that ``search_around_sky`` and ``search_around_3d`` return
  integer type index arrays for empty (non) matches. [#4877, #5083]

- Return an empty set of matches for ``search_around_sky`` and
  ``search_around_3d`` when one or both of the input coordinate
  arrays is empty. [#4875, #5083]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fix a bug with empty value at end of tab-delimited table on Windows. [#5370]

- Fix reading of big ASCII tables (more than 2Gb) with the fast reader.
  [#5319]

- Fix segfault with FastCsv and row with too many columns. [#5534]

- Fix problem reading an AASTex format table that does not have ``\\``
  at the end of the last table row. [#5427]

astropy.io.fits
^^^^^^^^^^^^^^^

- Removed raising of AssertionError that could occur after closing or
  deleting compressed image data. [#4690, #4694, #4948]

- Fixed bug that caused an ignored exception to be displayed under certain
  conditions when terminating a script after using fits.getdata(). [#4977]

- Fixed usage of inplace operations that were raising an exception with
  recent versions of Numpy due to implicit casting. [#5250]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Fixed bug of ``Resource.__repr__()`` having undefined attributes and
  variables. [#5382]

astropy.modeling
^^^^^^^^^^^^^^^^

- CompoundModel now correctly inherits _n_models, allowing the use of model sets [#5358]

astropy.units
^^^^^^^^^^^^^

- Fixed bug in Ci definition. [#5106]

- Non-ascii cds unit strings are now correctly represented using ``str`` also
  on python2. This solves bugs in parsing coordinates involving strings too.
  [#5355]

- Ensure ``Quantity`` supports ``np.float_power``, which is new in numpy 1.12.
  [#5480]

astropy.utils
^^^^^^^^^^^^^

- Fixed AttributeError when calling ``utils.misc.signal_number_to_name`` with
  Python3 [#5430].

astropy.wcs
^^^^^^^^^^^

- Update the ``_naxis{x}`` attributes when calling ``WCS.slice``. [#5411]


Other Changes and Additions
---------------------------

- The bundled ERFA was updated to version 1.3.0.  This includes the
  leap second planned for 2016 Dec 31. [#5418]

Version 1.0.10 (2016-06-09)
===========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- ``SkyCoord`` objects created before a new frame which has frame attributes
  is created no longer raise ``AttributeError`` when the new attributes are
  accessed [#5021]

- Fix some errors in the implementation of aberration  for ``get_sun``. [#4979]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fix problem reading a zero-length ECSV table with a bool type column. [#5010]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fix convenience functions (``getdata``, ``getheader``, ``append``,
  ``update``) to close files. [#4786]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- The astropy.io.votable.validator.html module is updated to handle division
  by zero when generating validation report. [#4699]

astropy.table
^^^^^^^^^^^^^

- Fixed a bug where ``pprint()`` sometimes raises ``UnicodeDecodeError``
  in Python 2. [#4946]

- Fix bug when doing outer join on multi-dimensional columns. [#4060]

- Fixed bug where Tables created from existing Table objects were not
  inheriting the ``primary_key`` attribute. [#4672]

astropy.tests
^^^^^^^^^^^^^

- Fix coverage reporting in Python 3. [#4822]

astropy.units
^^^^^^^^^^^^^

- Duplicates between long and short names are now removed in the ``names``
  and ``aliases`` properties of units. [#5036]

astropy.utils
^^^^^^^^^^^^^

- The astropy.utils.xml.unescaper module now also unescapes ``'%2F'`` to
  ``'/'`` and ``'&&'`` to ``'&'`` in a given URL. [#4699]

- Fix two problems related to the download cache: clear_download_cache() does
  not work in Python 2.7 and downloading in Python 2.7 and then Python 3
  can result in an exception. [#4810]

astropy.vo
^^^^^^^^^^

- Cache option now properly caches both downloaded JSON database and XML VO
  tables. [#4699]

- The astropy.vo.validator.conf.conesearch_urls listing is updated to reflect
  external changes to some VizieR Cone Search services. [#4699]

- VOSDatabase decodes byte-string to UTF-8 instead of ASCII to avoid
  UnicodeDecodeError for some rare cases. Fixed a Cone Search test that is
  failing as a side-effect of #4699. [#4757]

Other Changes and Additions
---------------------------

- Updated ``astropy.tests`` test runner code to work with Coverage v4.0 when
  generating test coverage reports. [#4176]


Version 1.0.9 (2016-03-10)
==========================

New Features
------------

astropy.nddata
^^^^^^^^^^^^^^

- ``NDArithmeticMixin`` check for matching WCS now works with
  ``astropy.wcs.WCS`` objects [#4499, #4503]

Bug Fixes
---------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Correct a bug in which ``psf_pad`` and ``fft_pad`` would be ignored [#4366]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fixed addition of new line characters after last row of data in
  ascii.latex.AASTex. [#4561]

- Fixed reading of Latex tables where the ``\tabular`` tag is in the first
  line. [#4595]

- Fix use of plain format strings with the fast writer. [#4517]

- Fix bug writing space-delimited file when table has empty fields. [#4417]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fixed possible segfault during error handling in FITS tile
  compression. [#4489]

- Fixed crash on pickling of binary table columns with the 'X', 'P', or
  'Q' format. [#4514]

- Fixed memory / reference leak that could occur when copying a ``FITS_rec``
  object (the ``.data`` for table HDUs). [#520]

- Fixed a memory / reference leak in ``FITS_rec`` that occurred in a wide
  range of cases, especially after writing FITS tables to a file, but in
  other cases as well. [#4539]

astropy.modeling
^^^^^^^^^^^^^^^^

- Fixed display of compound model expressions and components when printing
  compound model instances. [#4414, #4482]

astropy.stats
^^^^^^^^^^^^^

- the input for median_absolute_deviation will not be cast to plain numpy
  arrays when given subclasses of numpy arrays
  (like Quantity, numpy.ma.MaskedArray, etc.) [#4658]

- Fixed incorrect results when using median_absolute_deviation with masked
  arrays. [#4658]

astropy.utils
^^^^^^^^^^^^^

- The ``zest.releaser`` hooks included in Astropy are now injected locally to
  Astropy, rather than being global. [#4650]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Fixed ``fits2bitmap`` script to allow ext flag to contain extension
  names or numbers. [#4468]

- Fixed ``fits2bitmap`` default output filename generation for
  compressed FITS files. [#4468]


Version 1.0.8 (2016-01-08)
==========================

Bug Fixes
---------

astropy.io.fits
^^^^^^^^^^^^^^^

- Fixed an bug where updates to string columns in FITS tables were not saved
  on Python 3. [#4452]

astropy.units
^^^^^^^^^^^^^

- In-place peak-to-peak calculations now work on ``Quantity``. [#4442]

astropy.utils
^^^^^^^^^^^^^

- Fixed ``find_api_page`` to work correctly on python 3.x [#4378, #4379]


Version 1.0.7 (2015-12-04)
==========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Pickling of ``EarthLocation`` instances now also works on Python 2. [#4304]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fix fast writer so bytestring column output is not prefixed by 'b' in
  Python 3. [#4350]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fixed a regression that could cause writes of large FITS files to be
  truncated. [#4307]

- Astropy v1.0.6 included a fix (#4228) for an obscure case where the TDIM
  of a table column is smaller than the repeat count of its data format.
  This updates that fix in such a way that it works with Numpy 1.10 as well.
  [#4266]

astropy.table
^^^^^^^^^^^^^

- Fix a bug when pickling a Table with mixin columns (e.g. Time). [#4098]

astropy.time
^^^^^^^^^^^^

- Fix incorrect ``value`` attribute for epoch formats like "unix"
  when ``scale`` is different from the class ``epoch_scale``. [#4312]

astropy.utils
^^^^^^^^^^^^^

- Fixed an issue where if ipython is installed but ipykernel is not
  installed then importing astropy from the ipython console gave an
  IPython.kernel deprecation warning. [#4279]

- Fixed crash that could occur in ``ProgressBar`` when ``astropy`` is
  imported in an IPython startup script. [#4274]

Other Changes and Additions
---------------------------

- Updated bundled astropy-helpers to v1.0.6. [#4372]


Version 1.0.6 (2015-10-22)
==========================

Bug Fixes
---------

astropy.analytic_functions
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed blackbody analytic functions to properly support arrays of
  temperatures. [#4251]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Fixed errors in transformations for objects within a few AU of the
  Earth.  Included substantive changes to transformation machinery
  that may change distances at levels ~machine precision for other
  objects. [#4254]

astropy.io.fits
^^^^^^^^^^^^^^^

- ``fitsdiff`` and related functions now do a better job reporting differences
  between values that are different types but have the same representation
  (ex: the string '0' versus the number 0). [#4122]

- Miscellaneous fixes for supporting Numpy 1.10. [#4228]

- Fixed an issue where writing a column of unicode strings to a FITS table
  resulted in a quadrupling of size of the column (i.e. the format of the
  FITS column was 4 characters for every one in the original strings).
  [#4228]

- Added support for an obscure case (but nonetheless allowed by the FITS
  standard) where a column has some TDIMn keyword, but a repeat count in
  the TFORMn column greater than the number of elements implied by the
  TDIMn.  For example TFORMn = 100I, but TDIMn = '(5,5)'.  In this case
  the TDIMn implies 5x5 arrays in the column, but the TFORMn implies
  a 100 element 1-D array in the column.  In this case the TDIM takes
  precedence, and the remaining bytes in the column are ignored. [#4228]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Fixed crash with Python compiler optimization level = 2. [#4231]

astropy.vo
^^^^^^^^^^

- Fixed ``check_conesearch_sites`` with ``parallel=True`` on Python >= 3.3
  and on Windows (it was broken in both those cases for separate reasons).
  [#2970]

Other Changes and Additions
---------------------------

- All tests now pass against Numpy v1.10.x. This implies nominal support for
  Numpy 1.10.x moving forward (but there may still be unknown issues). For
  example, there is already a known performance issue with tables containing
  large multi-dimensional columns--for example, tables that contain entire
  images in one or more of their columns.  This is a known upstream issue in
  Numpy. [#4259]


Version 1.0.5 (2015-10-05)
==========================

Bug Fixes
---------

astropy.constants
^^^^^^^^^^^^^^^^^

- Rename units -> unit and error -> uncertainty in the ``repr`` and ``str``
  of constants to match attribute names. [#4147]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Fix string representation of ``SkyCoord`` objects transformed into
  the ``AltAz`` frame [#4055, #4057]

- Fix the ``search_around_sky`` function to allow ``storekdtree`` to be
  ``False`` as was intended. [#4082, #4212]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fix bug when extending one header (without comments) with another
  (with comments). [#3967]

- Somewhat improved resource usage for FITS data--previously a new ``mmap``
  was opened for each HDU of a FITS file accessed through an ``HDUList``.
  Each ``mmap`` used up a single file descriptor, causing problems with
  system resource limits for some users.  Now only a single ``mmap`` is
  opened, and shared for the data of all HDUs.  Note: The problem still
  persists with using the "convenience" functions.  For example using
  ``fits.getdata`` will create one ``mmap`` per HDU read this way (as
  opposed to opening the file with ``fits.open`` and accessing the HDUs
  through the ``HDUList`` object). [#4097]

- Fix bug where reading a file without a newline failed with an
  unrelated / unhelpful exception. [#4160]

astropy.modeling
^^^^^^^^^^^^^^^^

- Cleaned up ``repr`` of models that have no parameters. [#4076]

astropy.nddata
^^^^^^^^^^^^^^

- Initializing ``NDDataArray`` from another instance now sets ``flags`` as
  expected and no longer fails when ``uncertainty`` is set [#4129].
  Initializing an ``NDData`` subclass from a parent instance
  (eg. ``NDDataArray`` from ``NDData``) now sets the attributes other than
  ``data`` as it should [#4130, #4137].

astropy.table
^^^^^^^^^^^^^

- Fix an issue with setting fill value when column dtype is changed. [#4088]

- Fix bug when unpickling a bare Column where the _parent_table
  attribute was not set.  This impacted the Column representation. [#4099]

- Fix issue with the web browser opening with an empty page, and ensure that
  the url is correctly formatted for Windows. [#4132]

- Fix NameError in table stack exception message. [#4213]

astropy.utils
^^^^^^^^^^^^^

- ``resolve_name`` no longer causes ``sys.modules`` to be cluttered with
  additional copies of modules under a package imported like
  ``resolve_name('numpy')``. [#4084]

- ``console`` was updated to support IPython 4.x and Jupyter 1.x.
  This should suppress a ShimWarning that was appearing at
  import of astropy with IPython 4.0 or later. [#4078]

- Temporary downloaded files created by ``get_readable_fileobj`` when passed
  a URL are now deleted immediately after the file is closed. [#4198]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- The color for axes labels was set to white. Since white labels on white
  background are hard to read, the label color has been changed to black.
  [#4143]

- ``ImageNormalize`` now automatically determines ``vmin``/``vmax``
  (via the ``autoscale_None`` method) when they have not been set
  explicitly. [#4117]

astropy.vo
^^^^^^^^^^

- Cone Search validation no longer crashes when the provider gives an
  incomplete test query. It also ensures search radius for a test query
  is not too large to avoid timeout. [#4158, #4159]

Other Changes and Additions
---------------------------

- Astropy now supports Python 3.5. [#4027]

- Updated bundled version of astropy-helpers to 1.0.5. [#4215]

- Updated tests to support py.test 2.7, and upgraded the bundled copy of
  py.test to v2.7.3. [#4027]


Version 1.0.4 (2015-08-11)
==========================

New Features
------------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Modified Cython functions to release the GIL. This enables convolution
  to be parallelized effectively and gives large speedups when used with
  multithreaded task schedulers such as Dask. [#3949]

API Changes
-----------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Some transformations for an input coordinate that's a scalar now correctly
  return a scalar.  This was always the intended behavior, but it may break
  code that has been written to work-around this bug, so it may be viewed as
  an unplanned API change [#3920, #4039]

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- The ``astropy_mpl_style`` no longer sets ``interactive`` to ``True``, but
  instead leaves it at the user preference.  This makes using the style
  compatible with building docs with Sphinx, and other non-interactive
  contexts. [#4030]

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Fix bug where coordinate representation setting gets reset to default value
  when coordinate array is indexed or sliced. [#3824]

- Fixed confusing warning message shown when using dates outside current IERS
  data. [#3844]

- ``get_sun`` now yields a scalar when the input time is a scalar (this was a
  regression in v1.0.3 from v1.0.2) [#3998, #4039]

- Fixed bug where some scalar coordinates were incorrectly being changed to
  length-1 array coordinates after transforming through certain frames.
  [#3920, #4039]

- Fixed bug causing the ``separation`` methods of ``SkyCoord`` and frame
  classes to fail due to infinite recursion [#4033, #4039]

- Made it so that passing in a list of ``SkyCoord`` objects that are in
  UnitSphericalRepresentation to the ``SkyCoord`` constructor appropriately
  yields a new object in UnitSphericalRepresentation [#3938, #4039]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Fixed wCDM to not ignore the Ob0 parameter on initialization. [#3934]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fixed crash when updating data in a random groups HDU opened in update
  mode. [#3730]

- Fixed incorrect checksum / datasum being written when re-writing a scaled
  HDU (i.e. non-trivial BSCALE and/or BZERO) with
  ``do_not_scale_image_data=False``. [#3883]

- Fixed stray deprecation warning in ``BinTableHDU.copy()``. [#3798]

- Better handling of the ``BLANK`` keyword when auto-scaling scaled image
  data.  The ``BLANK`` keyword is now removed from the header after
  auto-scaling is applied, and it is restored properly (with floating point
  NaNs replaced by the filler value) when updating a file opened with the
  ``scale_back=True`` argument.  Invalid usage of the ``BLANK`` keyword is
  also better warned about during validation. [#3865]

- Reading memmaped scaled images won't fail when
  ``do_not_scale_image_data=True`` (that is, since we're just reading the raw
  / physical data there is no reason mmap can't be used). [#3766]

- Fixed a reference cycle that could sometimes cause FITS table-related
  objects (``BinTableHDU``, ``ColDefs``, etc.) to hang around in memory
  longer than expected. [#4012]

astropy.modeling
^^^^^^^^^^^^^^^^

- Improved support for pickling of compound models, including both compound
  model instances, and new compound model classes. [#3867]

- Added missing default values for ``Ellipse2D`` parameters. [#3903]

astropy.time
^^^^^^^^^^^^

- Fixed iteration of scalar ``Time`` objects so that ``iter()`` correctly
  raises a ``TypeError`` on them (while still allowing ``Time`` arrays to be
  iterated). [#4048]

astropy.units
^^^^^^^^^^^^^

- Added frequency-equivalency check when declaring doppler equivalencies
  [#3728]

- Define ``floor_divide`` (``//``) for ``Quantity`` to be consistent
  ``divmod``, such that it only works where the quotient is dimensionless.
  This guarantees that ``(q1 // q2) * q2 + (q1 % q2) == q1``. [#3817]

- Fixed the documentation of supported units to correctly report support for
  SI prefixes.  Previously the table of supported units incorrectly showed
  several derived unit as not supporting prefixes, when in fact they do.
  [#3835]

- Fix a crash when calling ``astropy.units.cds.enable()``.  This will now
  "set" rather than "add" units to the active set to avoid the namespace
  clash with the default units. [#3873]

- Ensure in-place operations on ``float32`` quantities work. [#4007]

astropy.utils
^^^^^^^^^^^^^

- The ``deprecated`` decorator did not correctly wrap classes that have a
  custom metaclass--the metaclass could be dropped from the deprecated
  version of the class. [#3997]

- The ``wraps`` decorator would copy the wrapped function's name to the
  wrapper function even when ``'__name__'`` is excluded from the ``assigned``
  argument. [#4016]

Misc
^^^^

- ``fitscheck`` no longer causes scaled image data to be rescaled when
  adding checksums to existing files. [#3884]

- Fixed an issue where running ``import astropy`` from within the source
  tree did not automatically build the extension modules if the source is
  from a source distribution (as opposed to a git repository). [#3932]

- Fixed multiple instances of a bug that prevented Astropy from being used
  when compiled with the ``python -OO`` flag, due to it causing all
  docstrings to be stripped out. [#3923]

- Removed source code template files that were being installed
  accidentally alongside installed Python modules. [#4014]

- Fixed a bug in the exception logging that caused a crash in the exception
  handler itself on Python 3 when exceptions do not include a message.
  [#4056]


Version 1.0.3 (2015-06-05)
==========================

New Features
------------

astropy.table
^^^^^^^^^^^^^

- Greatly improved the speed of printing a large table to the screen when
  only a few rows are being displayed. [#3796]

astropy.time
^^^^^^^^^^^^

- Add support for the 2015-Jun-30 leap second. [#3794]

API Changes
-----------

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Note that HTML formatted tables will not always be found with guess mode
  unless it passes certain heuristics that strongly suggest the presence of
  HTML in the input.  Code that expects to read tables from HTML should
  specify ``format='html'`` explicitly. See bug fixes below for more
  details. [#3693]

Bug Fixes
---------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Fix issue with repeated normalizations of ``Kernels``. [#3747]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Fixed ``get_sun`` to yield frames with the ``obstime`` set to what's passed into the function (previously it incorrectly always had J2000). [#3750]

- Fixed ``get_sun`` to account for aberration of light. [#3750]

- Fixed error in the GCRS->ICRS transformation that gave incorrect distances. [#3750]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Remove HTML from the list of automatically-guessed formats when reading if
  the file does not appear to be HTML.  This was necessary to avoid a
  commonly-encountered segmentation fault occurring in the libxml parser on
  MacOSX. [#3693]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fixes to support the upcoming Numpy 1.10. [#3419]

astropy.modeling
^^^^^^^^^^^^^^^^

- Polynomials are now scaled when used in a compound model. [#3702]

- Fixed the ``Ellipse2D`` model to be consistent with ``Disk2D`` in
  how pixels are included. [#3736]

- Fixed crash when evaluating a model that accepts no inputs. [#3772]

astropy.testing
^^^^^^^^^^^^^^^

- The Astropy py.test plugins that disable unintentional internet access
  in tests were also blocking use of local UNIX sockets in tests, which
  prevented testing some multiprocessing code--fixed. [#3713]

astropy.units
^^^^^^^^^^^^^

- Supported full SI prefixes for the barn unit ("picobarn", "femtobarn",
  etc.)  [#3753]

- Fix loss of precision when multiplying non-whole-numbered powers
  of units together.  For example, before this change, ``(u.m **
  1.5) ** Fraction(4, 5)`` resulted in an inaccurate floating-point
  power of ``1.2000000000000002``.  After this change, the exact
  rational number of ``Fraction(6, 5)`` is maintained. [#3790]

- Fixed printing of object ndarrays containing multiple Quantity
  objects with differing / incompatible units. Note: Unit conversion errors
  now cause a ``UnitConversionError`` exception to be raised.  However, this
  is a subclass of the ``UnitsError`` exception used previously, so existing
  code that catches ``UnitsError`` should still work. [#3778]

Other Changes and Additions
---------------------------

- Added a new ``astropy.__bibtex__`` attribute which gives a citation
  for Astropy in bibtex format. [#3697]

- The bundled version of ERFA was updated to v1.2.0 to address leapsecond
  updates. [#3802]


Version 0.4.6 (2015-05-29)
==========================

Bug Fixes
---------

astropy.time
^^^^^^^^^^^^

- Fixed ERFA code to handle the 2015-Jun-30 leap second. [#3795]


Version 1.0.2 (2015-04-16)
==========================

New Features
------------

astropy.modeling
^^^^^^^^^^^^^^^^

- Added support for polynomials with degree 0 or degree greater than 15.
  [#3574, 3589]

Bug Fixes
---------

astropy.config
^^^^^^^^^^^^^^

- The pre-astropy-0.4 configuration API has been fixed. It was
  inadvertently broken in 1.0.1. [#3627]

astropy.io.fits
^^^^^^^^^^^^^^^

- Fixed a severe memory leak that occurred when reading tile compressed
  images. [#3680]

- Fixed bug where column data could be unintentionally byte-swapped when
  copying data from an existing FITS file to a new FITS table with a
  TDIMn keyword for that column. [#3561]

- The ``ColDefs.change_attrib``, ``ColDefs.change_name``, and
  ``ColDefs.change_unit`` methods now work as advertised.  It is also
  possible (and preferable) to update attributes directly on ``Column``
  objects (for example setting ``column.name``), and the change will be
  accurately reflected in any associated table data and its FITS header.
  [#3283, #1539, #2618]

- Fixes an issue with the ``FITS_rec`` interface to FITS table data, where a
  ``FITS_rec`` created by copying an existing FITS table but adding new rows
  could not be sliced or masked correctly.  [#3641]

- Fixed handling of BINTABLE with TDIMn of size 1. [#3580]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Loading a ``TABLE`` element without any ``DATA`` now correctly
  creates a 0-row array. [#3636]

astropy.modeling
^^^^^^^^^^^^^^^^

- Added workaround to support inverses on compound models when one of the
  sub-models is itself a compound model with a manually-assigned custom
  inverse. [#3542]

- Fixed instantiation of polynomial models with constraints for parameters
  (constraints could still be assigned after instantiation, but not during).
  [#3606]

- Fixed fitting of 2D polynomial models with the ``LeVMarLSQFitter``. [#3606]

astropy.table
^^^^^^^^^^^^^

- Ensure ``QTable`` can be pickled [#3590]

- Some corner cases when instantiating an ``astropy.table.Table``
  with a Numpy array are handled [#3637]. Notably:

- a zero-length array is the same as passing ``None``

- a scalar raises a ``ValueError``

- a one-dimensional array is treated as a single row of a table.

- Ensure a ``Column`` without units is treated as an ``array``, not as an
  dimensionless ``Quantity``. [#3648]

astropy.units
^^^^^^^^^^^^^

- Ensure equivalencies that do more than just scale a ``Quantity`` are
  properly handled also in ``ufunc`` evaluations. [#2496, #3586]

- The LaTeX representation of the Angstrom unit has changed from
  ``\overset{\circ}{A}`` to ``\mathring{A}``, which should have
  better support across regular LaTeX, MathJax and matplotlib (as of
  version 1.5) [#3617]

astropy.vo
^^^^^^^^^^

- Using HTTPS/SSL for communication between SAMP hubs now works
  correctly on all supported versions of Python [#3613]

astropy.wcs
^^^^^^^^^^^

- When no ``relax`` argument is passed to ``WCS.to_header()`` and
  the result omits non-standard WCS keywords, a warning is
  emitted. [#3652]

Other Changes and Additions
---------------------------

astropy.vo
^^^^^^^^^^

- The number of retries for connections in ``astropy.vo.samp`` can now be
  configured by a ``n_retries`` configuration option. [#3612]

- Testing

- Running ``astropy.test()`` from within the IPython prompt has been
  provisionally re-enabled. [#3184]


Version 1.0.1 (2015-03-06)
==========================

Bug Fixes
---------

astropy.constants
^^^^^^^^^^^^^^^^^

- Ensure constants can be turned into ``Quantity`` safely. [#3537, #3538]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fix a segfault in the fast C parser when one of the column headers
  is empty [#3545].

- Fixed support for reading inf and nan values with the fast reader in
  Windows.  Also fixed in the case of using ``use_fast_converter=True``
  with the fast reader. [#3525]

- Fixed use of mmap in the fast reader on Windows. [#3525]

- Fixed issue where commented header would treat comments defining the table
  (i.e. column headers) as purely information comments, leading to problems
  when trying to round-trip the table. [#3562]

astropy.modeling
^^^^^^^^^^^^^^^^

- Fixed propagation of parameter constraints ('fixed', 'bounds', 'tied')
  between compound models and their components.  There is may still be some
  difficulty defining 'tied' constraints properly for use with compound
  models, however. [#3481]

astropy.nddata
^^^^^^^^^^^^^^

- Restore several properties to the compatibility class ``NDDataArray`` that
  were inadvertently omitted [#3466].

astropy.time
^^^^^^^^^^^^

- Time objects now always evaluate to ``True``, except when empty. [#3530]

Miscellaneous
-------------

- The ERFA wrappers are now written directly in the Python/C API
  rather than using Cython, for greater performance. [#3521]

- Improve import time of astropy [#3488].

Other Changes and Additions
---------------------------

- Updated bundled astropy-helpers version to v1.0.1 to address installation
  issues with some packages that depend on Astropy. [#3541]


Version 1.0 (2015-02-18)
========================

General
-------

- Astropy now requires Numpy 1.6.0 or later.

New Features
------------

astropy.analytic_functions
^^^^^^^^^^^^^^^^^^^^^^^^^^

- The ``astropy.analytic_functions`` was added to contain analytic functions
  useful for astronomy [#3077].

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

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

astropy.cosmology
^^^^^^^^^^^^^^^^^

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

astropy.io.ascii
^^^^^^^^^^^^^^^^

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

astropy.io.fits
^^^^^^^^^^^^^^^

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

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- ``astropy.io.votable.parse`` now takes a ``datatype_mapping``
  keyword argument to map invalid datatype names to valid ones in
  order to support non-compliant files. [#2675]

astropy.modeling
^^^^^^^^^^^^^^^^

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

astropy.nddata
^^^^^^^^^^^^^^

- New array-related utility functions in ``astropy.nddata.utils`` for adding
  and removing arrays from other arrays with different sizes/shapes. [#3201]

- New metaclass ``NDDataBase`` for enforcing the nddata interface in
  subclasses without restricting implementation of the data storage. [#2905]

- New mixin classes ``NDSlicingMixin`` for slicing, ``NDArithmeticMixin``
  for arithmetic operations, and ``NDIOMixin`` for input/output in NDData. [#2905]

- Added a decorator ``support_nddata`` that can be used to write functions
  that can either take separate arguments or NDData objects. [#2855]

astropy.stats
^^^^^^^^^^^^^

- Added ``mad_std()`` function. [#3208]

- Added ``gaussian_fwhm_to_sigma`` and ``gaussian_sigma_to_fwhm``
  constants. [#3208]

- New function ``sigma_clipped_stats`` which can be used to quickly get
  common statistics for an array, using sigma clipping at the same time.
  [#3201]

astropy.table
^^^^^^^^^^^^^

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

astropy.tests
^^^^^^^^^^^^^

- Added a new Quantity-aware ``assert_quantity_allclose``. [#3273]

astropy.time
^^^^^^^^^^^^

- ``Time`` can now handle arbitrary array dimensions, with operations
  following standard numpy broadcasting rules. [#3138]

astropy.units
^^^^^^^^^^^^^

- Support for VOUnit has been updated to be compliant with version
  1.0 of the standard. [#2901]

- Added an ``insert`` method to insert values into a ``Quantity`` object.
  This is similar to the ``numpy.insert`` function. [#3049]

- When viewed in IPython, ``Quantity`` objects with array values now render
  using LaTeX and scientific notation. [#2271]

- Added ``units.quantity_input`` decorator to validate quantity inputs to a
  function for unit compatibility. [#3072]

- Added ``units.astronomical_unit`` as a long form for ``units.au``. [#3303]

astropy.utils
^^^^^^^^^^^^^

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

astropy.visualization
^^^^^^^^^^^^^^^^^^^^^

- Created ``astropy.visualization`` module and added functionality relating
  to image normalization (i.e. stretching and scaling) as well as a new
  script ``fits2bitmap`` that can produce a bitmap image from a FITS file.
  [#3201]

- Added dictionary ``astropy.visualization.mpl_style.astropy_mpl_style``
  which can be used to set a uniform plotstyle specifically for tutorials
  that is improved compared to matplotlib defaults. [#2719, #2787, #3200]

astropy.wcs
^^^^^^^^^^^

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

Misc
^^^^

- ``astropy._erfa`` was added as a new subpackage wrapping the functionality
  of the ERFA library in python.  This is primarily of use for other astropy
  subpackages, but the API may be made more public in the future. [#2992]


API Changes
-----------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

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

astropy.cosmology
^^^^^^^^^^^^^^^^^

- The functional interface to the cosmological routines as well as
  ``set_current`` and ``get_current`` (deprecated in v0.4) have now been
  removed. [#2990]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Added a new argument to ``htmldict`` in the HTML reader named
  ``parser``, which allows the user to specify which parser
  BeautifulSoup should use as a backend. [#2815]

- Add ``FixedWidthTwoLine`` reader to guessing. This will allows to read
  tables that a copied from screen output like ``print my_table`` to be read
  automatically. Discussed in #3025 and #3099 [#3109]

astropy.io.fits
^^^^^^^^^^^^^^^

- A new optional argument ``cache`` has been added to
  ``astropy.io.fits.open()``.  When opening a FITS file from a URL,
  ``cache`` is a boolean value specifying whether or not to save the
  file locally in Astropy's download cache (``True`` by default). [#3041]

astropy.modeling
^^^^^^^^^^^^^^^^

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

astropy.nddata
^^^^^^^^^^^^^^

- ``flags``, ``shape``, ``size``, ``dtype`` and ``ndim`` properties removed
  from ``astropy.nddata.NDData``. [#2905]

- Arithmetic operations, uncertainty propagation, slicing and automatic
  conversion to a numpy array removed from ``astropy.nddata.NDData``. The
  class ``astropy.nddata.NDDataArray`` is functionally equivalent to the
  old ``NDData``.  [#2905]

astropy.table
^^^^^^^^^^^^^

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

astropy.time
^^^^^^^^^^^^

- The ``Time.val`` and ``Time.vals`` properties (deprecated in v0.3) and the
  ``Time.lon``, and ``Time.lat`` properties (deprecated in v0.4) have now
  been removed. [#2990]

- Add ``decimalyear`` format that represents time as a decimal year. [#3265]

astropy.units
^^^^^^^^^^^^^

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

astropy.utils
^^^^^^^^^^^^^

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
  ``interactive`` kwarg. [#2658, #2789]

astropy.wcs
^^^^^^^^^^^

- The ``WCS.calcFootprint`` method (deprecated in v0.4) has now been removed.
  [#2990]

- An invalid unit in a ``CUNITn`` keyword now displays a warning and
  returns a ``UnrecognizedUnit`` instance rather than raising an
  exception [#3190]

Bug Fixes
---------

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- ``astropy.convolution.discretize_model`` now handles arbitrary callables
  correctly [#2274].

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- ``Angle.to_string`` now outputs unicode arrays instead of object arrays.
  [#2981]

- ``SkyCoord.to_string`` no longer gives an error when used with an array
  coordinate with more than one dimension. [#3340]

- Fixed support for subclasses of ``UnitSphericalRepresentation`` and
  ``SphericalRepresentation`` [#3354, #3366]

- Fixed latex display of array angles in IPython notebook. [#3480]

astropy.io.ascii
^^^^^^^^^^^^^^^^

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

astropy.io.fits
^^^^^^^^^^^^^^^

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

astropy.logger
^^^^^^^^^^^^^^

- Fix a bug that occurred when displaying warnings that produced an error
  message ``dictionary changed size during iteration``. [#3353]

astropy.modeling
^^^^^^^^^^^^^^^^

- Fixed a bug in ``SLSQPLSQFitter`` where the ``maxiter`` argument was not
  passed correctly to the optimizer. [#3339]

astropy.table
^^^^^^^^^^^^^

- Fix a problem where ``table.hstack`` fails to stack multiple references to
  the same table, e.g. ``table.hstack([t, t])``. [#2995]

- Fixed a problem where ``table.vstack`` and ``table.hstack`` failed to stack
  a single table, e.g. ``table.vstack([t])``. [#3313]

- Fix a problem when doing nested iterators on a single table. [#3358]

- Fix an error when an empty list, tuple, or ndarray is used for item access
  within a table.  This now returns the table with no rows. [#3442]

astropy.time
^^^^^^^^^^^^

- When creating a Time object from a datetime object the time zone
  info is now correctly used. [#3160]

- For Time objects, it is now checked that numerical input is finite. [#3396]

astropy.units
^^^^^^^^^^^^^

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

astropy.utils
^^^^^^^^^^^^^

- ``treat_deprecations_as_exceptions`` has been fixed to recognize Astropy
  deprecation warnings. [#3015]

- Converted representation of progress bar units without suffix
  from float to int in console.human_file_size. [#2201, #2202, #2721, #3299]

astropy.wcs
^^^^^^^^^^^

- ``astropy.wcs.WCS.sub`` now accepts unicode strings as input on
  Python 2.x [#3356]

Misc
^^^^

- Some modules and tests that would crash upon import when using a non-final
  release of Numpy (e.g. 1.9.0rc1). [#3471]

Other Changes and Additions
---------------------------

- The bundled copy of astropy-helpers has been updated to v1.0. [#3515]

- Updated ``astropy.extern.configobj`` to Version 5. Version 5 uses ``six``
  and the same code covers both Python 2 and Python 3. [#3149]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- The ``repr`` of ``SkyCoord`` and coordinate frame classes now separate
  frame attributes and coordinate information.  [#2704, #2882]

astropy.io.fits
^^^^^^^^^^^^^^^

- Overwriting an existing file using the ``clobber=True`` option no longer
  displays a warning message. [#1963]

- ``fits.open`` no longer catches ``OSError`` exceptions on missing or
  unreadable files-- instead it raises the standard Python exceptions in such
  cases. [#2756, #2785]

astropy.table
^^^^^^^^^^^^^

- Sped up setting of ``Column`` slices by an order of magnitude. [#2994, #3020]

- Updated the bundled ``six`` module to version 1.7.3 and made 1.7.3 the
  minimum acceptable version of ``six``. [#2814]

- The version of ERFA included with Astropy is now v1.1.1 [#2971]

- The code base is now fully Python 2 and 3 compatible and no longer requires
  2to3. [#2033]

- `funcsigs <https://pypi.org/project/funcsigs>`_ is included in
  utils.compat, but defaults to the inspect module components where available
  (3.3+) [#3151].

- The list of modules displayed in the pytest header can now be customized.
  [#3157]

- `jinja2 <http://jinja.pocoo.org/docs/dev/>`_>=2.7 is now required to build the
  source code from the git repository, in order to allow the ERFA wrappers to
  be generated. [#3166]


Version 0.4.5 (2015-02-16)
==========================

Bug Fixes
---------

- Fixed unnecessary attempt to run ``git`` when importing astropy.  In
  particular, fixed a crash in Python 3 that could result from this when
  importing Astropy when the current working directory is an empty git
  repository. [#3475]

Other Changes and Additions
---------------------------

- Updated bundled copy of astropy-helpers to v0.4.6. [#3508]


Version 0.4.4 (2015-01-21)
==========================

Bug Fixes
---------

astropy.vo.samp
^^^^^^^^^^^^^^^

- ``astropy.vo.samp`` is now usable on Python builds that do not
  support the SSLv3 protocol (which depends both on the version of
  Python and the version of OpenSSL or LibreSSL that it is built
  against.) [#3308]

API Changes
-----------

astropy.vo.samp
^^^^^^^^^^^^^^^

- The default SSL protocol used is now determined from the default
  used in the Python ``ssl`` standard library.  This default may be
  different depending on the exact version of Python you are using.
  [#3308]

astropy.wcs
^^^^^^^^^^^

- WCS allows slices of the form slice(None, x, y), which previously resulted
  in an unsliced copy being returned (note: this was previously incorrectly
  reported as fixed in v0.4.3) [#2909]


Version 0.4.3 (2015-01-15)
==========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

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

astropy.cosmology
^^^^^^^^^^^^^^^^^

- The ``ztol`` keyword argument to z_at_value now works correctly [#2993].

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fix a bug in Python 3 when guessing file format using a file object as
  input.  Also improve performance in same situation for Python 2. [#3132]

- Fix a problem where URL was being downloaded for each guess. [#2001]

astropy.io.fits
^^^^^^^^^^^^^^^

- The ``in`` operator now works correctly for checking if an extension
  is in an ``HDUList`` (as given via EXTNAME, (EXTNAME, EXTVER) tuples,
  etc.) [#3060]

- Added workaround for bug in MacOS X <= 10.8 that caused np.fromfile to
  fail. [#3078]

- Added support for the ``RICE_ONE`` compression type synonym. [#3115]

astropy.modeling
^^^^^^^^^^^^^^^^

- Fixed a test failure on Debian/PowerPC and Debian/s390x. [#2708]

- Fixed crash in evaluating models that have more outputs than inputs--this
  case may not be handled as desired for all conceivable models of this
  format (some may have to implement custom ``prepare_inputs`` and
  ``prepare_outputs`` methods).  But as long as all outputs can be assumed
  to have a shape determined from the broadcast of all inputs with all
  parameters then this can be used safely. [#3250]

astropy.table
^^^^^^^^^^^^^

- Fix a bug that caused join to fail for multi-dimensional columns. [#2984]

- Fix a bug where MaskedColumn attributes which had been changed since
  the object was created were not being carried through when slicing. [#3023]

- Fix a bug that prevented initializing a table from a structured array
  with multi-dimensional columns with copy=True. [#3034]

- Fixed unnecessarily large unicode columns when instantiating a table from
  row data on Python 3. [#3052]

- Improved the warning message when unable to aggregate non-numeric
  columns. [#2700]

astropy.units
^^^^^^^^^^^^^

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

astropy.utils
^^^^^^^^^^^^^

- Fixed an issue with the ``deprecated`` decorator on classes that invoke
  ``super()`` in their ``__init__`` method. [#3004]

- Fixed a bug which caused the ``metadata_conflicts`` parameter to be
  ignored in the ``astropy.utils.metadata.merge`` function. [#3294]

astropy.vo
^^^^^^^^^^

- Fixed an issue with reconnecting to a SAMP Hub. [#2674]

astropy.wcs
^^^^^^^^^^^

- Invalid or out of range values passed to ``wcs_world2pix`` will
  now be correctly identified and returned as ``nan``
  values. [#2965]

- Fixed an issue which meant that Python thought ``WCS`` objects were
  iterable. [#3066]

Misc
^^^^

- Astropy will now work if your Python interpreter does not have the
  ``bz2`` module installed. [#3104]

- Fixed ``ResourceWarning`` for ``astropy/extern/bundled/six.py`` that could
  occur sometimes after using Astropy in Python 3.4. [#3156]

Other Changes and Additions
---------------------------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Improved the agreement of the FK5 <-> Galactic conversion with other
  codes, and with the FK5 <-> FK4 <-> Galactic route. [#3107]


Version 0.4.2 (2014-09-23)
==========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- ``Angle`` accepts hours:mins or deg:mins initializers (without
  seconds). In these cases float minutes are also accepted.

- The ``repr`` for coordinate frames now displays the frame attributes
  (ex: ra, dec) in a consistent order.  It should be noted that as part of
  this fix, the ``BaseCoordinateFrame.get_frame_attr_names()`` method now
  returns an ``OrderedDict`` instead of just a ``dict``. [#2845]

astropy.io.fits
^^^^^^^^^^^^^^^

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

astropy.io.misc
^^^^^^^^^^^^^^^

- Fixed a bug that prevented h5py ``Dataset`` objects from being
  automatically recognized by ``Table.read``. [#2831]

astropy.modeling
^^^^^^^^^^^^^^^^

- Make ``LevMarLSQFitter`` work with ``weights`` keyword. [#2900]

astropy.table
^^^^^^^^^^^^^

- Fixed reference cycle in tables that could prevent ``Table`` objects
  from being freed from memory. [#2879]

- Fixed an issue where ``Table.pprint()`` did not print the header to
  ``stdout`` when ``stdout`` is redirected (say, to a file). [#2878]

- Fixed printing of masked values when a format is specified. [#1026]

- Ensured that numpy ufuncs that return booleans return plain ``ndarray``
  instances, just like the comparison operators. [#2963]

astropy.time
^^^^^^^^^^^^

- Ensure bigendian input to Time works on a little-endian machine
  (and vice versa).  [#2942]

astropy.units
^^^^^^^^^^^^^

- Ensure unit is kept when adding 0 to quantities. [#2968]

astropy.utils
^^^^^^^^^^^^^

- Fixed color printing on Windows with IPython 2.0. [#2878]

astropy.vo
^^^^^^^^^^

- Improved error message on Cone Search time out. [#2687]

Other Changes and Additions
---------------------------

- Fixed a couple issues with files being inappropriately included and/or
  excluded from the source archive distributions of Astropy. [#2843, #2854]

- As part of fixing the fact that masked elements of table columns could not be
  printed when a format was specified, the column format string options were
  expanded to allow simple specifiers such as ``'5.2f'``. [#2898]

- Ensure numpy 1.9 is supported. [#2917]

- Ensure numpy master is supported, by making ``np.cbrt`` work with quantities.
  [#2937]

Version 0.4.1 (2014-08-08)
==========================

Bug Fixes
---------

astropy.config
^^^^^^^^^^^^^^

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

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Fixed the multiplication of ``Kernel`` with numpy floats. [#2174]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- ``Distance`` can now take a list of quantities. [#2261]

- For in-place operations for ``Angle`` instances in which the result unit
  is not an angle, an exception is raised before the instance is corrupted.
  [#2718]

- ``CartesianPoints`` are now deprecated in favor of
  ``CartesianRepresentation``. [#2727]

astropy.io.misc
^^^^^^^^^^^^^^^

- An existing table within an HDF5 file can be overwritten without affecting
  other datasets in the same HDF5 file by simultaneously using
  ``overwrite=True`` and ``append=True`` arguments to the ``Table.write``
  method. [#2624]

astropy.logger
^^^^^^^^^^^^^^

- Fixed a crash that could occur in rare cases when (such as in bundled
  apps) where submodules of the ``email`` package are not importable. [#2671]

astropy.nddata
^^^^^^^^^^^^^^

- ``astropy.nddata.NDData()`` no longer raises a ``ValueError`` when passed
  a numpy masked array which has no masked entries. [#2784]

astropy.table
^^^^^^^^^^^^^

- When saving a table to a FITS file containing a unit that is not
  supported by the FITS standard, a warning rather than an exception
  is raised. [#2797]

astropy.units
^^^^^^^^^^^^^

- By default, ``Quantity`` and its subclasses will now convert to float also
  numerical types such as ``decimal.Decimal``, which are stored as objects
  by numpy. [#1419]

- The units ``count``, ``pixel``, ``voxel`` and ``dbyte`` now output
  to FITS, OGIP and VOUnit formats correctly. [#2798]

astropy.utils
^^^^^^^^^^^^^

- Restored missing information from deprecation warning messages
  from the ``deprecated`` decorator. [#2811]

- Fixed support for ``staticmethod`` deprecation in the ``deprecated``
  decorator. [#2811]

astropy.wcs
^^^^^^^^^^^

- Fixed a memory leak when ``astropy.wcs.WCS`` objects are copied
  [#2754]

- Fixed a crash when passing ``ra_dec_order=True`` to any of the
  ``*2world`` methods. [#2791]

Other Changes and Additions
---------------------------

- Bundled copy of astropy-helpers upgraded to v0.4.1. [#2825]

- General improvements to documentation and docstrings [#2722, #2728, #2742]

- Made it easier for third-party packagers to have Astropy use their own
  version of the ``six`` module (so long as it meets the minimum version
  requirement) and remove the copy bundled with Astropy.  See the
  astropy/extern/README file in the source tree.  [#2623]


Version 0.4 (2014-07-16)
========================

New Features
------------

astropy.constants
^^^^^^^^^^^^^^^^^

- Added ``b_wien`` to represent Wien wavelength displacement law constant.
  [#2194]

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Changed the input parameter in ``Gaussian1DKernel`` and
  ``Gaussian2DKernel`` from ``width`` to ``stddev`` [#2085].

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

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

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Added ``z_at_value`` function to find the redshift at which a cosmology
  function matches a desired value. [#1909]

- Added ``FLRW.differential_comoving_volume`` method to give the differential
  comoving volume at redshift z. [#2103]

- The functional interface is now deprecated in favor of the more-explicit
  use of methods on cosmology objects. [#2343]

- Updated documentation to reflect the removal of the functional
  interface. [#2507]

astropy.io.ascii
^^^^^^^^^^^^^^^^

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

astropy.io.fits
^^^^^^^^^^^^^^^

- Included a new command-line script called ``fitsheader`` to display the
  header(s) of a FITS file from the command line. [#2092]

- Added new verification options ``fix+ignore``, ``fix+warn``,
  ``fix+exception``, ``silentfix+ignore``, ``silentfix+warn``, and
  ``silentfix+exception`` which give more control over how to report fixable
  errors as opposed to unfixable errors.

astropy.modeling
^^^^^^^^^^^^^^^^

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

astropy.nddata
^^^^^^^^^^^^^^

- Allow initialization ``NDData`` or ``StdDevUncertainty`` with a
  ``Quantity``. [#2380]

astropy.stats
^^^^^^^^^^^^^

- Added flat prior to binom_conf_interval and binned_binom_proportion

- Change default in ``sigma_clip`` from ``np.median`` to ``np.ma.median``.
  [#2582]

astropy.sphinx
^^^^^^^^^^^^^^

- Note, the following new features are included in astropy-helpers as well:

- The ``automodapi`` and ``automodsumm`` extensions now include sphinx
  configuration options to write out what ``automodapi`` and ``automodsumm``
  generate, mainly for debugging purposes. [#1975, #2022]

- Reference documentation now shows functions/class docstrings at the
  intended user-facing API location rather than the actual file where
  the implementation is found. [#1826]

- The ``automodsumm`` extension configuration was changed to generate
  documentation of class ``__call__`` member functions. [#1817, #2135]

- ``automodapi`` and ``automodsumm`` now have an ``:allowed-package-names:``
  option that make it possible to document functions and classes that
  are in a different namespace.  [#2370]

astropy.table
^^^^^^^^^^^^^

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

astropy.time
^^^^^^^^^^^^

- Mean and apparent sidereal time can now be calculated using the
  ``sidereal_time`` method [#1418].

- The time scale now defaults to UTC if no scale is provided. [#2091]

- ``TimeDelta`` objects can have all scales but UTC, as well as, for
  consistency with time-like quantities, undefined scale (where the
  scale is taken from the object one adds to or subtracts from).
  This allows, e.g., to work consistently in TDB.  [#1932]

- ``Time`` now supports ISO format strings that end in "Z". [#2211, #2203]

astropy.units
^^^^^^^^^^^^^

- Support for the unit format `Office of Guest Investigator Programs (OGIP)
  FITS files
  <https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/>`__
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

astropy.utils
^^^^^^^^^^^^^

- ``timer.RunTimePredictor`` now uses ``astropy.modeling`` in its
  ``do_fit()`` method. [#1896]

astropy.vo
^^^^^^^^^^

- A new sub-package, ``astropy.vo.samp``, is now available (this was
  previously the SAMPy package, which has been refactored for use in
  Astropy). [#1907]

- Enhanced functionalities for ``VOSCatalog`` and ``VOSDatabase``. [#1206]

astropy.wcs
^^^^^^^^^^^

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

Misc
^^^^

- Includes the new astropy-helpers package which separates some of Astropy's
  build, installation, and documentation infrastructure out into an
  independent package, making it easier for Affiliated Packages to depend on
  these features.  astropy-helpers replaces/deprecates some of the submodules
  in the ``astropy`` package (see API Changes below).  See also
  `APE 4 <https://github.com/astropy/astropy-APEs/blob/master/APE4.rst>`_
  for more details on the motivation behind and implementation of
  astropy-helpers.  [#1563]


API Changes
-----------

astropy.config
^^^^^^^^^^^^^^

- The configuration system received a major overhaul, as part of APE3.  It is
  no longer possible to save configuration items from Python, but instead
  users must edit the configuration file directly.  The locations of
  configuration items have moved, and some have been changed to science state
  values.  The old locations should continue to work until astropy 0.5, but
  deprecation warnings will be displayed.  See the `Configuration transition
  <https://docs.astropy.org/en/v0.4/config/config_0_4_transition.html>`_
  docs for a detailed description of the changes and how to update existing
  code. [#2094]

astropy.io.fits
^^^^^^^^^^^^^^^

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

astropy.modeling
^^^^^^^^^^^^^^^^

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

astropy.nddata
^^^^^^^^^^^^^^

- Issue warning if unit is changed from a non-trivial value by directly
  setting ``NDData.unit``. [#2411]

- The ``mask`` and ``flag`` attributes of ``astropy.nddata.NDData`` can now
  be set with any array-like object instead of requiring that they be set
  with a ``numpy.ndarray``. [#2419]

astropy.sphinx
^^^^^^^^^^^^^^

- Use of the ``astropy.sphinx`` module is deprecated; all new development of
  this module is in ``astropy_helpers.sphinx`` which should be used instead
  (therefore documentation builds that made use of any of the utilities in
  ``astropy.sphinx`` now have ``astropy_helpers`` as a documentation
  dependency).

astropy.table
^^^^^^^^^^^^^

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

astropy.time
^^^^^^^^^^^^

- Correct use of UT in TDB calculation [#1938, #1939].

- ``TimeDelta`` objects can have scales other than TAI [#1932].

- Location information should now be passed on via an ``EarthLocation``
  instance or anything that initialises it, e.g., a tuple containing
  either geocentric or geodetic coordinates. [#1928]

astropy.units
^^^^^^^^^^^^^

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

astropy.wcs
^^^^^^^^^^^

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
  <https://ui.adsabs.harvard.edu/abs/2005ASPC..347..491S>`__. [#2360]

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

Misc
^^^^

- The ``astropy.setup_helpers`` and ``astropy.version_helpers`` modules are
  deprecated; any non-critical fixes and development to those modules should
  be in ``astropy_helpers`` instead.  Packages that use these modules in
  their ``setup.py`` should depend on ``astropy_helpers`` following the same
  pattern as in the Astropy package template.


Bug Fixes
---------

astropy.constants
^^^^^^^^^^^^^^^^^

- ``astropy.constants.Constant`` objects can now be deep
  copied. [#2601]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- The distance modulus function in ``astropy.cosmology`` can now handle
  negative distances, which can occur in certain closed cosmologies. [#2008]

- Removed accidental imports of some extraneous variables in
  ``astropy.cosmology`` [#2025]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- ``astropy.io.ascii.read`` would fail to read lists of strings where some of
  the strings consisted of just a newline ("\n"). [#2648]

astropy.io.fits
^^^^^^^^^^^^^^^

- Use NaN for missing values in FITS when using Table.write for float
  columns. Earlier the default fill value was close to 1e20.[#2186]

- Fixes for checksums on 32-bit platforms.  Results may be different
  if writing or checking checksums in "nonstandard" mode.  [#2484]

- Additional minor bug fixes ported from PyFITS.  [#2575]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- It is now possible to save an ``astropy.table.Table`` object as a
  VOTable with any of the supported data formats, ``tabledata``,
  ``binary`` and ``binary2``, by using the ``tabledata_format``
  kwarg. [#2138]

- Fixed a crash writing out variable length arrays. [#2577]

astropy.nddata
^^^^^^^^^^^^^^

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

astropy.sphinx
^^^^^^^^^^^^^^

- Fix crash in smart resolver when the resolution doesn't work. [#2591]

astropy.table
^^^^^^^^^^^^^

- The ``astropy.table.Column`` object can now use both functions and callable
  objects as formats. [#2313]

- Fixed a problem on 64 bit windows that caused errors
  "expected 'DTYPE_t' but got 'long long'" [#2490]

- Fix initialisation of ``TableColumns`` with lists or tuples.  [#2647]

- Fix removal of single column using ``remove_columns``. [#2699]

- Fix a problem that setting a row element within a masked table did not
  update the corresponding table element. [#2734]

astropy.time
^^^^^^^^^^^^

- Correct UT1->UTC->UT1 round-trip being off by 1 second if UT1 is
  on a leap second. [#2077]

astropy.units
^^^^^^^^^^^^^

- ``Quantity.copy`` now behaves identically to ``ndarray.copy``, and thus
  supports the ``order`` argument (for numpy >=1.6). [#2284]

- Composing base units into identical composite units now works. [#2382]

- Creating and composing/decomposing units is now substantially faster [#2544]

- ``Quantity`` objects now are able to be assigned NaN [#2695]

astropy.wcs
^^^^^^^^^^^

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
  to achieve this. [#2468]

- If the C extension for ``astropy.wcs`` was not built or fails to import for
  any reason, ``import astropy.wcs`` will result in an ``ImportError``,
  rather than getting obscure errors once the ``astropy.wcs`` is used.
  [#2061]

- When the C extension for ``astropy.wcs`` is built using a version of
  ``wscslib`` already present in the system, the package does not try
  to install ``wcslib`` headers under ``astropy/wcs/include``. [#2536]

- Fixes an unresolved external symbol error in the
  ``astropy.wcs._wcs`` C extension on Microsoft Windows when built
  with a Microsoft compiler. [#2478]

Misc
^^^^

- Running the test suite with ``python setup.py test`` now works if
  the path to the source contains spaces. [#2488]

- The version of ERFA included with Astropy is now v1.1.0 [#2497]

- Removed deprecated option from travis configuration and force use of
  wheels rather than allowing build from source. [#2576]

- The short option ``-n`` to run tests in parallel was broken
  (conflicts with the distutils built-in option of "dry-run").
  Changed to ``-j``. [#2566]

Other Changes and Additions
---------------------------

- python setup.py test --coverage will now give more accurate
  results, because the coverage analysis will include early imports of
  astropy.  There doesn't seem to be a way to get this to work when
  doing ``import astropy; astropy.test()``, so the ``coverage``
  keyword to ``astropy.test`` has been removed.  Coverage testing now
  depends only on `coverage.py
  <http://coverage.readthedocs.io/en/latest/>`__, not
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


Version 0.3.2 (2014-05-13)
==========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- if ``sep`` argument is specified to be a single character in
  ``sexagisimal_to_string``, it now includes separators only between
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

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Fixed ``format()`` compatibility with Python 2.6. [#2129]

- Be more careful about converting to floating point internally [#1815, #1818]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- The CDS reader in ``astropy.io.ascii`` can now handle multiple
  description lines in ReadMe files. [#2225]

- When reading a table with values that generate an overflow error during
  type conversion (e.g. overflowing the native C long type), fall through to
  using string. Previously this generated an exception [#2234].

- Recognize any string with one to four dashes as null value. [#1335]

astropy.io.fits
^^^^^^^^^^^^^^^

- Allow pickling of ``FITS_rec`` objects. [#1597]

- Improved behavior when writing large compressed images on OSX by removing
  an unnecessary check for platform architecture. [#2345]

- Fixed an issue where Astropy ``Table`` objects containing boolean columns
  were not correctly written out to FITS files. [#1953]

- Several other bug fixes ported from PyFITS v3.2.3 [#2368]

- Fixed a crash on Python 2.x when writing a FITS file directly to a
  ``StringIO.StringIO`` object. [#2463]

astropy.io.registry
^^^^^^^^^^^^^^^^^^^

- Allow readers/writers with the same name to be attached to different
  classes. [#2312]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- By default, floating point values are now written out using
  ``repr`` rather than ``str`` to preserve precision [#2137]

astropy.modeling
^^^^^^^^^^^^^^^^

- Fixed the ``SIP`` and ``InverseSIP`` models both so that they work in the
  first place, and so that they return results consistent with the SIP
  functions in ``astropy.wcs``. [#2177]

astropy.stats
^^^^^^^^^^^^^

- Ensure the ``axis`` keyword in ``astropy.stats.funcs`` can now be used for
  all axes. [#2173]

astropy.table
^^^^^^^^^^^^^

- Ensure nameless columns can be printed, using 'None' for the header. [#2213]

astropy.time
^^^^^^^^^^^^

- Fixed pickling of ``Time`` objects. [#2123]

astropy.units
^^^^^^^^^^^^^

- ``Quantity._repr_latex_()`` returns ``NotImplementedError`` for quantity
  arrays instead of an uninformative formatting exception. [#2258]

- Ensure ``Quantity.flat`` always returns ``Quantity``. [#2251]

- Angstrom unit renders better in MathJax [#2286]

astropy.utils
^^^^^^^^^^^^^

- Progress bars will now be displayed inside the IPython
  qtconsole. [#2230]

- ``data.download_file()`` now evaluates ``REMOTE_TIMEOUT()`` at runtime
  rather than import time. Previously, setting ``REMOTE_TIMEOUT`` after
  import had no effect on the function's behavior. [#2302]

- Progressbar will be limited to 100% so that the bar does not exceed the
  terminal width.  The numerical display can still exceed 100%, however.

astropy.vo
^^^^^^^^^^

- Fixed ``format()`` compatibility with Python 2.6. [#2129]

- Cone Search validation no longer raises ``ConeSearchError`` for positive RA.
  [#2240, #2242]

astropy.wcs
^^^^^^^^^^^

- Fixed a bug where calling ``astropy.wcs.Wcsprm.sub`` with
  ``WCSSUB_CELESTIAL`` may cause memory corruption due to
  underallocation of a temporary buffer. [#2350]

- Fixed a memory allocation bug in ``astropy.wcs.Wcsprm.sub`` and
  ``astropy.wcs.Wcsprm.copy``.  [#2439]

Misc
^^^^

- Fixes for compatibility with Python 3.4. [#1945]

- ``import astropy; astropy.test()`` now correctly uses the same test
  configuration as ``python setup.py test`` [#1811]


Version 0.3.1 (2014-03-04)
==========================

Bug Fixes
---------

astropy.config
^^^^^^^^^^^^^^

- Fixed a bug where ``ConfigurationItem.set_temp()`` does not reset to
  default value when exception is raised within ``with`` block. [#2117]

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- Fixed a bug where ``_truncation`` was left undefined for ``CustomKernel``.
  [#2016]

- Fixed a bug with ``_normalization`` when ``CustomKernel`` input array
  sums to zero. [#2016]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

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

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Allow passing unicode delimiters when reading or writing tables.  The
  delimiter must be convertible to pure ASCII.  [#1949]

- Fix a problem when reading a table and renaming the columns to names that
  already exist. [#1991]

astropy.io.fits
^^^^^^^^^^^^^^^

- Ported all bug fixes from PyFITS 3.2.1.  See the PyFITS changelog at
  https://pyfits.readthedocs.io/en/v3.2.1/ [#2056]

astropy.io.misc
^^^^^^^^^^^^^^^

- Fixed issues in the HDF5 Table reader/writer functions that occurred on
  Windows. [#2099]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

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

astropy.modeling
^^^^^^^^^^^^^^^^

- Fixed bug in computation of model derivatives in ``LinearLSQFitter``.
  [#1903]

- Raise a ``NotImplementedError`` when fitting composite models. [#1915]

- Fixed bug in the computation of the ``Gaussian2D`` model. [#2038]

- Fixed bug in the computation of the ``AiryDisk2D`` model. [#2093]

astropy.sphinx
^^^^^^^^^^^^^^

- Added slightly more useful debug info for AstropyAutosummary. [#2024]

astropy.table
^^^^^^^^^^^^^

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

astropy.time
^^^^^^^^^^^^

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

astropy.units
^^^^^^^^^^^^^

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

astropy.utils
^^^^^^^^^^^^^

- Fixed crash in ``timer.RunTimePredictor.do_fit``. [#1905]

- Fixed ``astropy.utils.compat.argparse`` for Python 3.1. [#2017]

astropy.wcs
^^^^^^^^^^^

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

Misc
^^^^

- There are a number of improvements to make Astropy work better on big
  endian platforms, such as MIPS, PPC, s390x and SPARC. [#1849]

- The test suite will now raise exceptions when a deprecated feature of
  Python or Numpy is used.  [#1948]

Other Changes and Additions
---------------------------

- A new function, ``astropy.wcs.get_include``, has been added to get the
  location of the ``astropy.wcs`` C header files. [#1755]

- The doctests in the ``.rst`` files in the ``docs`` folder are now
  tested along with the other unit tests.  This is in addition to the
  testing of doctests in docstrings that was already being performed.
  See ``docs/development/testguide.rst`` for more information. [#1771]

- Fix a problem where import fails on Python 3 if setup.py exists
  in current directory. [#1877]


Version 0.3 (2013-11-20)
========================

New Features
------------

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
  <https://docs.astropy.org/en/v0.3/development/codeguide.html#unicode-guidelines>`_
  for more information. [#1441]

- ``astropy.utils.misc.find_api_page`` is now imported into the top-level.
  This allows usage like ``astropy.find_api_page(astropy.units.Quantity)``.
  [#1779]

astropy.convolution
^^^^^^^^^^^^^^^^^^^

- New class-based system for generating kernels, replacing ``make_kernel``.
  [#1255] The ``astropy.nddata.convolution`` sub-package has now been moved
  to ``astropy.convolution``. [#1451]

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

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

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Added support for including massive Neutrinos in the cosmology classes. The
  Planck (2013) cosmology has been updated to use this. [#1364]

- Calculations now use and return ``Quantity`` objects where appropriate.
  [#1237]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Added support for writing IPAC format tables [#1152].

astropy.io.fits
^^^^^^^^^^^^^^^

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

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Updated to support the VOTable 1.3 draft. [#433]

- Added the ability to look up and group elements by their utype attribute.
  [#622]

- The format of the units of a VOTable file can be specified using the
  ``unit_format`` parameter.  Note that units are still always written out
  using the CDS format, to ensure compatibility with the standard.

astropy.modeling
^^^^^^^^^^^^^^^^

- Added a new framework for representing and evaluating mathematical models
  and for fitting data to models.  See "What's New in Astropy 0.3" in the
  documentation for further details. [#493]

astropy.stats
^^^^^^^^^^^^^

- Added robust statistics functions
  ``astropy.stats.funcs.median_absolute_deviation``,
  ``astropy.stats.funcs.biweight_location``, and
  ``astropy.stats.funcs.biweight_midvariance``. [#621]

- Added ``astropy.stats.funcs.signal_to_noise_oir_ccd`` for computing the
  signal to noise ratio for source being observed in the optical/IR using a
  CCD. [#870]

- Add ``axis=int`` option to ``stropy.stats.funcs.sigma_clip`` to allow
  clipping along a given axis for multidimensional data. [#1083]

astropy.table
^^^^^^^^^^^^^

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

astropy.time
^^^^^^^^^^^^

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

astropy.units
^^^^^^^^^^^^^

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

astropy.vo
^^^^^^^^^^

- New package added to support Virtual Observatory Simple Cone Search query
  and service validation. [#552]

astropy.wcs
^^^^^^^^^^^

- Fixed attribute error in ``astropy.wcs.Wcsprm`` (lattype->lattyp) [#1463]

- Included a new command-line script called ``wcslint`` and accompanying API
  for validating the WCS in a given FITS file or header. [#580]

- Upgraded included version of WCSLIB to 4.19.

astropy.utils
^^^^^^^^^^^^^

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

astropy.extern.six
^^^^^^^^^^^^^^^^^^

- Added `six <https://pypi.org/project/six/>`_ for python2/python3
  compatibility

- Astropy now uses the ERFA library instead of the IAU SOFA library for
  fundamental time transformation routines.  The ERFA library is derived, with
  permission, from the IAU SOFA library but is distributed under a BSD license.
  See ``license/ERFA.rst`` for details. [#1293]

astropy.logger
^^^^^^^^^^^^^^

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
-----------

- General

- The configuration option ``utils.console.use_unicode`` has been
  moved to the top level and renamed to ``unicode_output``.  It now
  not only affects console widgets, such as progress bars, but also
  controls whether calling `unicode` on certain classes will return a
  string containing unicode characters.

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

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

astropy.cosmology
^^^^^^^^^^^^^^^^^

- The Planck (2013) cosmology will likely give slightly different (and more
  accurate) results due to the inclusion of Neutrino masses. [#1364]

- Cosmology class properties now return ``Quantity`` objects instead of
  simple floating-point values. [#1237]

- The names of cosmology instances are now truly optional, and are set to
  ``None`` rather than the name of the class if the user does not provide
  them.  [#1705]

astropy.io.ascii
^^^^^^^^^^^^^^^^

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

astropy.io.fits
^^^^^^^^^^^^^^^

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

astropy.io.registry
^^^^^^^^^^^^^^^^^^^

- Identifier functions for reading/writing Table and NDData objects should
  now accept ``(origin, *args, **kwargs)`` instead of ``(origin, args,
  kwargs)``. [#591]

- Added a new ``astropy.io.registry.get_formats`` function for listing
  registered I/O formats and details about the their readers/writers. [#1669]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Added a new option ``use_names_over_ids`` option to use when converting
  from VOTable objects to Astropy Tables. This can prevent a situation where
  column names are not preserved when converting from a VOTable. [#609]

astropy.nddata
^^^^^^^^^^^^^^

- The ``astropy.nddata.convolution`` sub-package has now been moved to
  ``astropy.convolution``, and the ``make_kernel`` function has been removed.
  (the kernel classes should be used instead) [#1451]

astropy.stats.funcs
^^^^^^^^^^^^^^^^^^^

- For ``sigma_clip``, the ``maout`` optional parameter has been removed, and
  the function now always returns a masked array.  A new boolean parameter
  ``copy`` can be used to indicated whether the input data should be copied
  (``copy=True``, default) or used by reference (``copy=False``) in the
  output masked array. [#1083]

astropy.table
^^^^^^^^^^^^^

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

astropy.time
^^^^^^^^^^^^

- For consistency with ``Quantity``, the attributes ``val`` and
  ``is_scalar`` have been renamed to ``value`` and ``isscalar``,
  respectively, and the attribute ``vals`` has been dropped. [#767]

- The double-float64 internal representation of time is used more
  efficiently to enable better accuracy. [#366]

- Format and scale arguments are now allowed to be case-insensitive. [#1128]

astropy.units
^^^^^^^^^^^^^

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

astropy.wcs
^^^^^^^^^^^

- For those including the ``astropy.wcs`` C headers in their project, they
  should now include it as:

  #include "astropy_wcs/astropy_wcs_api.h"

  instead of:

  #include "astropy_wcs_api.h"

  [#1631]

- The ``--enable-legacy`` option for ``setup.py`` has been removed. [#1493]

Bug Fixes
---------

astropy.io.ascii
^^^^^^^^^^^^^^^^

- The ``write()`` function was ignoring the ``fill_values`` argument. [#910]

- Fixed an issue in ``DefaultSplitter.join`` where the delimiter attribute
  was ignored when writing the CSV. [#1020]

- Fixed writing of IPAC tables containing null values. [#1366]

- When a table with no header row was read without specifying the format and
  using the ``names`` argument, then the first row could be dropped. [#1692]

astropy.io.fits
^^^^^^^^^^^^^^^

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

astropy.nddata
^^^^^^^^^^^^^^

- Fixed crash when trying to multiple or divide ``NDData`` objects with
  uncertainties. [#1547]

astropy.table
^^^^^^^^^^^^^

- Using a list of strings to index a table now correctly returns a new table
  with the columns named in the list. [#1454]

- Inequality operators now work properly with ``Column`` objects. [#1685]

astropy.time
^^^^^^^^^^^^

- ``Time`` scale and format attributes are now shown when calling ``dir()``
  on a ``Time`` object. [#1130]

astropy.wcs
^^^^^^^^^^^

- Fixed assignment to string-like WCS attributes on Python 3. [#956]

astropy.units
^^^^^^^^^^^^^

- Fixed a bug that caused the order of multiplication/division of plain
  Numpy arrays with Quantities to matter (i.e. if the plain array comes
  first the units were not preserved in the output). [#899]

- Directly instantiated ``CompositeUnits`` were made printable without
  crashing. [#1576]

Misc
^^^^

- Fixed various modules that hard-coded ``sys.stdout`` as default arguments
  to functions at import time, rather than using the runtime value of
  ``sys.stdout``. [#1648]

- Minor documentation fixes and enhancements [#922, #1034, #1210, #1217,
  #1491, #1492, #1498, #1582, #1608, #1621, #1646, #1670, #1756]

- Fixed a crash that could sometimes occur when running the test suite on
  systems with platform names containing non-ASCII characters. [#1698]

Other Changes and Additions
---------------------------

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

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Angles containing out of bounds minutes or seconds (e.g. 60) can be
  parsed--the value modulo 60 is used with carry to the hours/minutes, and a
  warning is issued rather than raising an exception. [#990]

astropy.io.fits
^^^^^^^^^^^^^^^

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

astropy.io.misc
^^^^^^^^^^^^^^^

- The HDF5 reader can now refer to groups in the path as well as datasets;
  if given a group, the first dataset in that group is read. [#1159]

astropy.nddata
^^^^^^^^^^^^^^

- ``NDData`` objects have more helpful, though still rudimentary ``__str__`
  and ``__repr__`` displays. [#1313]

astropy.units
^^^^^^^^^^^^^

- Added 'cycle' unit. [#1160]

- Extended units supported by the CDS formatter/parser. [#1468]

- Added unicode an LaTeX symbols for liter. [#1618]

astropy.wcs
^^^^^^^^^^^

- Redundant SCAMP distortion parameters are removed with SIP distortions are
  also present. [#1278]

- Added iterative implementation of ``all_world2pix`` that can be reliably
  inverted. [#1281]


Version 0.2.5 (2013-10-25)
==========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Fixed incorrect string formatting of Angles using ``precision=0``. [#1319]

- Fixed string formatting of Angles using ``decimal=True`` which ignored the
  ``precision`` argument. [#1323]

- Fixed parsing of format strings using appropriate unicode characters
  instead of the ASCII ``-`` for minus signs. [#1429]

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fixed a crash in the IPAC table reader when the ``include/exclude_names``
  option is set. [#1348]

- Fixed writing AASTex tables to honor the ``tabletype`` option. [#1372]

astropy.io.fits
^^^^^^^^^^^^^^^

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

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Added a warning for when a VOTable 1.2 file contains no ``RESOURCES``
  elements (at least one should be present). [#1337]

- Fixed a test failure specific to MIPS architecture caused by an errant
  floating point warning. [#1179]

astropy.nddata.convolution
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Prevented in-place modification of the input arrays to ``convolve()``.
  [#1153]

astropy.table
^^^^^^^^^^^^^

- Added HTML escaping for string values in tables when outputting the table
  as HTML. [#1347]

- Added a workaround in a bug in Numpy that could cause a crash when
  accessing a table row in a masked table containing ``dtype=object``
  columns. [#1229]

- Fixed an issue similar to the one in #1229, but specific to unmasked
  tables. [#1403]

astropy.units
^^^^^^^^^^^^^

- Improved error handling for unparsable units and fixed parsing CDS units
  without mantissas in the exponent. [#1288]

- Added a physical type for spectral flux density. [#1410]

- Normalized conversions that should result in a scale of exactly 1.0 to
  round off slight floating point imprecisions. [#1407]

- Added support in the CDS unit parser/formatter for unusual unit prefixes
  that are nonetheless required to be supported by that convention. [#1426]

- Fixed the parsing of ``sqrt()`` in unit format strings which was returning
  ``unit ** 2`` instead of ``unit ** 0.5``. [#1458]

astropy.wcs
^^^^^^^^^^^

- When passing a single array to the wcs transformation functions,
  (``astropy.wcs.Wcs.all_pix2world``, etc.), its second dimension must now
  exactly match the number of dimensions in the transformation. [#1395]

- Improved error message when incorrect arguments are passed to
  ``WCS.wcs_world2pix``. [#1394]

- Fixed a crash when trying to read WCS from FITS headers on Python 3.3
  in Windows. [#1363]

- Only headers that are required as part of the WCSLIB C API are installed
  by the package, per request of system packagers. [#1666]

Misc
^^^^

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
---------------------------

- Astropy installation now requests setuptools >= 0.7 during build/installation
  if neither distribute or setuptools >= 0.7 is already installed.  In other
  words, if ``import setuptools`` fails, ``ez_setup.py`` is used to bootstrap
  the latest setuptools (rather than using ``distribute_setup.py`` to bootstrap
  the now obsolete distribute package). [#1197]

- When importing Astropy from a source checkout without having built the
  extension modules first an ``ImportError`` is raised rather than a
  ``SystemExit`` exception. [#1269]


Version 0.2.4 (2013-07-24)
==========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Fixed the angle parser to support parsing the string "1 degree". [#1168]

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Fixed a crash in the ``comoving_volume`` method on non-flat cosmologies
  when passing it an array of redshifts.

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fixed a bug that prevented saving changes to the comment symbol when
  writing changes to a table. [#1167]

astropy.io.fits
^^^^^^^^^^^^^^^

- Added a workaround for a bug in 64-bit OSX that could cause truncation when
  writing files greater than 2^32 bytes in size. [#839]

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Fixed incorrect reading of tables containing multiple ``<RESOURCE>``
  elements. [#1223]

astropy.table
^^^^^^^^^^^^^

- Fixed a bug where ``Table.remove_column`` and ``Table.rename_column``
  could cause a masked table to lose its masking. [#1120]

- Fixed bugs where subclasses of ``Table`` did not preserver their class in
  certain operations. [#1142]

- Fixed a bug where slicing a masked table did not preserve the mask. [#1187]

astropy.units
^^^^^^^^^^^^^

- Fixed a bug where the ``.si`` and ``.cgs`` properties of dimensionless
  ``Quantity`` objects raised a ``ZeroDivisionError``. [#1150]

- Fixed a bug where multiple subsequent calls to the ``.decompose()`` method
  on array quantities applied a scale factor each time. [#1163]

Misc
^^^^

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
---------------------------

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Added a new ``Plank13`` object representing the Plank 2013 results. [#895]

astropy.units
^^^^^^^^^^^^^

- Performance improvements in initialization of ``Quantity`` objects with
  a large number of elements. [#1231]


Version 0.2.3 (2013-05-30)
==========================

Bug Fixes
---------

astropy.time
^^^^^^^^^^^^

- Fixed inaccurate handling of leap seconds when converting from UTC to UNIX
  timestamps. [#1118]

- Tightened required accuracy in many of the time conversion tests. [#1121]

Misc
^^^^

- Fixed a regression that was introduced in v0.2.2 by the fix to issue #992
  that was preventing installation of Astropy affiliated packages that use
  Astropy's setup framework. [#1124]


Version 0.2.2 (2013-05-21)
==========================

Bug Fixes
---------

astropy.io
^^^^^^^^^^

- Fixed issues in both the ``fits`` and ``votable`` sub-packages where array
  byte order was not being handled consistently, leading to possible crashes
  especially on big-endian systems. [#1003]

astropy.io.fits
^^^^^^^^^^^^^^^

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
  properly from unicode strings (so long as they are convertible to ASCII).

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

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- Stopped deprecation warnings from the ``astropy.io.votable`` package that
  could occur during setup. [#970]

- Fixed an issue where INFO elements were being incorrectly dropped when
  occurring inside a TABLE element. [#1000]

- Fixed obscure test failures on MIPS platforms. [#1010]

astropy.nddata.convolution
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fixed an issue in ``make_kernel()`` when using an Airy function kernel.
  Also removed the superfluous 'brickwall' option. [#939]

astropy.table
^^^^^^^^^^^^^

- Fixed a crash that could occur when adding a row to an empty (rowless)
  table with masked columns. [#973]

- Made it possible to assign to one table row from the value of another row,
  effectively making it easier to copy rows, for example. [#1019]

astropy.time
^^^^^^^^^^^^

- Added appropriate ``__copy__`` and ``__deepcopy__`` behavior; this
  omission caused a seemingly unrelated error in FK5 coordinate separation.
  [#891]

astropy.units
^^^^^^^^^^^^^

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

Misc
^^^^

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
---------------------------

- Some performance improvements to the ``astropy.units`` package, in particular
  improving the time it takes to import the sub-package. [#1015]


Version 0.2.1 (2013-04-03)
==========================

Bug Fixes
---------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- Fixed encoding errors that could occur when formatting coordinate objects
  in code using ``from __future__ import unicode_literals``. [#817]

- Fixed a bug where the minus sign was dropped when string formatting dms
  coordinates with -0 degrees. [#875]

astropy.io.fits
^^^^^^^^^^^^^^^

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

astropy.io.votable
^^^^^^^^^^^^^^^^^^

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

astropy.nddata.convolution
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Added better handling of ``inf`` values to the ``convolve_fft`` family of
  functions. [#893]

astropy.table
^^^^^^^^^^^^^

- Fixed silent failure to assign values to a row on multiple columns. [#764]

- Fixed various buggy behavior when viewing a table after sorting by one of
  its columns. [#829]

- Fixed using ``numpy.where()`` with table indexing. [#838]

- Fixed a bug where opening a remote table with ``Table.read()`` could cause
  the entire table to be downloaded twice. [#845]

- Fixed a bug where ``MaskedColumn`` no longer worked if the column being
  masked is renamed. [#916]

astropy.units
^^^^^^^^^^^^^

- Added missing capability for array ``Quantity``\s to be initializable by
  a list of ``Quantity``\s. [#835]

- Fixed the definition of year and lightyear to be in terms of Julian year
  per the IAU definition. [#861]

- "degree" was removed from the list of SI base units. [#863]

astropy.wcs
^^^^^^^^^^^

- Fixed ``TypeError`` when calling ``WCS.to_header_string()``. [#822]

- Added new method ``WCS.all_world2pix`` for converting from world
  coordinates to pixel space, including inversion of the astrometric
  distortion correction. [#1066, #1281]

Misc
^^^^

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
---------------------------

- Added logo and branding for Windows binary installers. [#741]

- Upgraded included version libexpat to 2.1.0. [#781]

- ~25% performance improvement in unit composition/decomposition. [#836]

- Added previously missing LaTeX formatting for ``L_sun`` and ``R_sun``. [#841]

- ConfigurationItem\s now have a more useful and informative __repr__
  and improved documentation for how to use them. [#855]

- Added a friendlier error message when trying to import astropy from a source
  checkout without first building the extension modules inplace. [#864]

- py.test now outputs more system information for help in debugging issues
  from users. [#869]

- Added unit definitions "mas" and "uas" for "milliarcsecond" and
  "microarcsecond" respectively. [#892]


Version 0.2 (2013-02-19)
========================

New Features
------------

astropy.coordinates
^^^^^^^^^^^^^^^^^^^

- This new subpackage contains a representation of celestial coordinates,
  and provides a wide range of related functionality.  While
  fully-functional, it is a work in progress and parts of the API may
  change in subsequent releases.

astropy.cosmology
^^^^^^^^^^^^^^^^^

- Update to include cosmologies with variable dark energy equations of state.
  (This introduces some API incompatibilities with the older Cosmology
  objects).

- Added parameters for relativistic species (photons, neutrinos) to the
  astropy.cosmology classes. The current treatment assumes that neutrinos are
  massless. [#365]

- Add a WMAP9 object using the final (9-year) WMAP parameters from
  Hinshaw et al. 2013. It has also been made the default cosmology.
  [#629, #724]

- astropy.table I/O infrastructure for custom readers/writers
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

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Improved integration with the ``astropy.table`` Table class so that
  table and column metadata (e.g. keywords, units, description,
  formatting) are directly available in the output table object.  The
  CDS, DAOphot, and IPAC format readers now provide this type of
  integrated metadata.

- Changed to using ``astropy.table`` masked tables instead of NumPy
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

astropy.wcs
^^^^^^^^^^^

- From updating the underlying wcslib 4.16:

- When ``astropy.wcs.WCS`` constructs a default coordinate representation
  it will give it the special name "DEFAULTS", and will not report "Found
  one coordinate representation".

Other Changes and Additions
---------------------------

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

- astropy.nddata.convolution removed requirement of PyFFTW3; uses Numpy's
  FFT by default instead with the added ability to specify an FFT
  implementation to use. [#660]


Bug Fixes
---------

astropy.io.ascii
^^^^^^^^^^^^^^^^

- Fixed crash when pprinting a row with INDEF values. [#511]

- Fixed failure when reading DAOphot files with empty keyword values. [#666]

astropy.io.fits
^^^^^^^^^^^^^^^

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

- Fixed some bugs with WCS distortion paper record-valued keyword cards:

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

astropy.io.votable
^^^^^^^^^^^^^^^^^^

- The ``Table`` class now maintains a single array object which is a
  Numpy masked array.  For variable-length columns, the object that
  is stored there is also a Numpy masked array.

- Changed the ``pedantic`` configuration option to be ``False`` by default
  due to the vast proliferation of non-compliant VO Tables. [#296]

- Renamed ``astropy.io.vo`` to ``astropy.io.votable``.

astropy.table
^^^^^^^^^^^^^

- Added a workaround for an upstream bug in Numpy 1.6.2 that could cause
  a maximum recursion depth RuntimeError when printing table rows. [#341]

astropy.wcs
^^^^^^^^^^^

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
  installed on the user's system. [#616, #640]

- Changed use of ``log.warn`` in the logging module to ``log.warning`` since
  the former is deprecated. [#624]


Version 0.1 (2012-06-19)
========================

- Initial release.
