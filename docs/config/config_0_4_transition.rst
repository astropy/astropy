.. _config-0-4-transition:

Configuration transition
========================

This document describes the changes in the configuration system in
astropy 0.4 and how to update code to use it.

For users
---------

The config file
^^^^^^^^^^^^^^^

If you never edited the configuration file in
``~/.astropy/config/astropy.cfg``, there is nothing for you to do.
The first time you import astropy 0.4, it will automatically be
replaced with the configuration file template for astropy 0.4.

If you did edit the configuration file, it will be left untouched.
However, the template for astropy 0.4 will be installed as
``~/.astropy/config/astropy.0.4.cfg``.  You can manually compare your
changes to this file to determine what customizations should be
brought over.

Saving
^^^^^^

Saving configuration items from Python has been completely removed.
Instead, the configuration file must be edited directly.

Renames
^^^^^^^

The location of the configuration parameters have been simplified, so
they always appear in a high-level subpackage of astropy, rather than
low-level file names (which were really an implementation detail that
shouldn't have been exposed to the user).  On the Python side,
configuration items always are referenced through a ``conf`` object at
the root of a subpackage.

Some configuration items that affect the results of science
calculations have been removed as configuration parameters altogether
and converted to science state objects that must be changed from
Python code.

The following table lists all of the moves (in alphabetical order by
original configuration file location).  The old names will continue to
work both from Python and the configuration file for the astropy 0.4
release cycle, and will be removed altogether in astropy 0.5.

.. list-table:: Renamed configuration parameters
   :widths: 20 20 20 20
   :header-rows: 1

   * - Old config file location
     - Old Python location
     - New config file location
     - New Python location
   * - ``[] unicode_output``
     - ``UNICODE_OUTPUT``
     - *unchanged*
     - ``conf.unicode_output``
   * - ``[coordinates.name_resolve] name_resolve_timeout``
     - ``coordinates.name_resolve.NAME_RESOLVE_TIMEOUT``
     - ``[astropy.utils.data] remote_timeout``
     - ``astropy.utils.data.conf.remote_timeout``
   * - ``[coordinates.name_resolve] sesame_url``
     - ``coordinates.name_resolve.SESAME_URL``
     - *removed*
     - ``coordinates.name_resolve.sesame_url.get/set``
   * - ``[coordinates.name_resolve] sesame_database``
     - ``coordinates.name_resolve.SESAME_DATABASE``
     - *removed*
     - ``coordinates.name_resolve.sesame_database.get/set``
   * - ``[cosmology.core] default_cosmology``
     - ``cosmology.core.DEFAULT_COSMOLOGY``
     - *removed*
     - ``cosmology.default_cosmology.get/set``
   * - ``[io.fits] enable_record_valued_keyword_cards``
     - ``io.fits.ENABLE_RECORD_VALUED_KEYWORD_CARDS``
     - *unchanged*
     - ``io.fits.conf.enable_record_valued_keyword_cards``
   * - ``[io.fits] extension_name_case_sensitive``
     - ``io.fits.EXTENSION_NAME_CASE_SENSITIVE``
     - *unchanged*
     - ``io.fits.conf.extension_name_case_sensitive``
   * - ``[io.fits] strip_header_whitespace``
     - ``io.fits.STRIP_HEADER_WHITESPACE``
     - *unchanged*
     - ``io.fits.conf.strip_header_whitespace``
   * - ``[io.fits] use_memmap``
     - ``io.fits.USE_MEMMAP``
     - *unchanged*
     - ``io.fits.conf.use_memmap``
   * - ``[io.votable.table] pedantic``
     - ``io.votable.table.PEDANTIC``
     - ``[io.votable] pedantic``
     - ``io.votable.conf.pedantic``
   * - ``[logger] log_exceptions``
     - ``logger.LOG_EXCEPTIONS``
     - *unchanged*
     - ``logger.conf.log_exceptions``
   * - ``[logger] log_file_format``
     - ``logger.LOG_FILE_FORMAT``
     - *unchanged*
     - ``logger.conf.log_file_format``
   * - ``[logger] log_file_level``
     - ``logger.LOG_FILE_LEVEL``
     - *unchanged*
     - ``logger.conf.log_file_level``
   * - ``[logger] log_file_path``
     - ``logger.LOG_FILE_PATH``
     - *unchanged*
     - ``logger.conf.log_file_path``
   * - ``[logger] log_level``
     - ``logger.LOG_LEVEL``
     - *unchanged*
     - ``logger.conf.log_level``
   * - ``[logger] log_to_file``
     - ``logger.LOG_TO_FILE``
     - *unchanged*
     - ``logger.conf.log_to_file``
   * - ``[logger] log_warnings``
     - ``logger.LOG_WARNINGS``
     - *unchanged*
     - ``logger.conf.log_warnings``
   * - ``[logger] use_color``
     - ``logger.USE_COLOR``
     - ``[] use_color``
     - ``conf.use_color``
   * - ``[nddata.nddata] warn_unsupported_correlated``
     - ``nddata.nddata.WARN_UNSUPPORTED_CORRELATED``
     - ``[nddata] warn_unsupported_correlated``
     - ``nddata.conf.warn_unsupported_correlated``
   * - ``[table.column] auto_colname``
     - ``table.column.AUTO_COLNAME``
     - ``[table] auto_colname``
     - ``table.conf.auto_colname``
   * - ``[table.jsviewer] jquery_url``
     - ``table.jsviewer.JQUERY_URL``
     - *unchanged*
     - ``table.jsviewer.conf.jquery_url``
   * - ``[table.jsviewer] datatables_url``
     - ``table.jsviewer.DATATABLES_URL``
     - *unchanged*
     - ``table.jsviewer.conf.datatables_url``
   * - ``[table.pprint] max_lines``
     - ``table.pprint.MAX_LINES``
     - ``[] max_lines``
     - ``conf.max_lines``
   * - ``[table.pprint] max_width``
     - ``table.pprint.MAX_WIDTH``
     - ``[] max_width``
     - ``conf.max_width``
   * - ``[utils.console] use_color``
     - ``utils.console.USE_COLOR``
     - ``[] use_color``
     - ``conf.use_color``
   * - ``[utils.data] compute_hash_block_size``
     - ``astropy.utils.data.COMPUTE_HASH_BLOCK_SIZE``
     - *unchanged*
     - ``astropy.utils.data.conf.compute_hash_block_size``
   * - ``[utils.data] dataurl``
     - ``astropy.utils.data.DATAURL``
     - *unchanged*
     - ``astropy.utils.data.conf.dataurl``
   * - ``[utils.data] delete_temporary_downloads_at_exit``
     - ``astropy.utils.data.DELETE_TEMPORARY_DOWNLOADS_AT_EXIT``
     - *unchanged*
     - ``astropy.utils.data.conf.delete_temporary_downloads_at_exit``
   * - ``[utils.data] download_cache_block_size``
     - ``astropy.utils.data.DOWNLOAD_CACHE_BLOCK_SIZE``
     - *unchanged*
     - ``astropy.utils.data.conf.download_cache_block_size``
   * - ``[utils.data] download_cache_lock_attempts``
     - ``astropy.utils.data.download_cache_lock_attempts``
     - *unchanged*
     - ``astropy.utils.data.conf.download_cache_lock_attempts``
   * - ``[utils.data] remote_timeout``
     - ``astropy.utils.data.REMOTE_TIMEOUT``
     - *unchanged*
     - ``astropy.utils.data.conf.remote_timeout``
   * - ``[vo.client.conesearch] conesearch_dbname``
     - ``vo.client.conesearch.CONESEARCH_DBNAME``
     - ``[vo] conesearch_dbname``
     - ``vo.conf.conesearch_dbname``
   * - ``[vo.client.vos_catalog] vos_baseurl``
     - ``vo.client.vos_catalog.BASEURL``
     - ``[vo] vos_baseurl``
     - ``vo.conf.vos_baseurl``
   * - ``[vo.samp.utils] use_internet``
     - ``vo.samp.utils.ALLOW_INTERNET``
     - ``[vo.samp] use_internet``
     - ``vo.samp.conf.use_internet``
   * - ``[vo.validator.validate] cs_mstr_list``
     - ``vo.validator.validate.CS_MSTR_LIST``
     - ``[vo.validator] conesearch_master_list``
     - ``vo.validator.conf.conesearch_master_list``
   * - ``[vo.validator.validate] cs_urls``
     - ``vo.validator.validate.CS_URLS``
     - ``[vo.validator] conesearch_urls``
     - ``vo.validator.conf.conesearch_urls``
   * - ``[vo.validator.validate] noncrit_warnings``
     - ``vo.validator.validate.noncrit_warnings``
     - ``[vo.validator] noncritical_warnings``
     - ``vo.validator.conf.noncritical_warnings``

For affiliated package authors
------------------------------

For an affiliated package to support both astropy 0.3 and 0.4,
following the astropy 0.3 config instructions should continue to work.
Note that saving of configuration items has been removed entirely from
astropy 0.4 without a deprecation cycle, so if saving configuration
programmatically is important to your package, you may want to
consider another method to save that state.

However, by the release of astropy 0.5, the astropy 0.3 config API
will no longer work.  The following describes how to transition an
affiliated package written for astropy 0.3 to support astropy 0.4 and
later.  It will not be possible to support astropy 0.3, 0.4 and 0.5
simultaneously.  Below ``pkgname`` is the name of your affiliated
package.

The automatic generation of configuration files from the
``ConfigurationItem`` objects that it finds has been removed.
Instead, the project should include a hard-coded "template"
configuration file in ``pkgname/pkgname.cfg``.  By convention, and to
ease upgrades for end users, all of the values should be commented
out.  For example:

.. code-block:: ini

    [nddata]

    ## Whether to issue a warning if NDData arithmetic is performed with
    ## uncertainties and the uncertainties do not support the propagation of
    ## correlated uncertainties.
    # warn_unsupported_correlated = True

Affiliated packages should transition to using
`astropy.config.ConfigItem` objects as members of
`astropy.config.ConfigNamespace` subclasses.

For example, the following is an example of the astropy 0.3 and
earlier method to define configuration items::

    from astropy.config import ConfigurationItem

    ENABLE_RECORD_VALUED_KEYWORD_CARDS = ConfigurationItem(
        'enabled_record_valued_keyword_cards', True,
        'If True, enable support for record-valued keywords as described by '
        'FITS WCS Paper IV. Otherwise they are treated as normal keywords.')

    EXTENSION_NAME_CASE_SENSITIVE = ConfigurationItem(
        'extension_name_case_sensitive', False,
        'If True, extension names (i.e. the EXTNAME keyword) should be '
       'treated as case-sensitive.')

The above, converted to the new method, looks like::

    from astropy import config as _config

    class Conf(_config.ConfigNamespace):
        """
        Configuration parameters for `astropy.io.fits`.
        """

        enable_record_valued_keyword_cards = _config.ConfigItem(
            True,
            'If True, enable support for record-valued keywords as described by '
            'FITS WCS Paper IV. Otherwise they are treated as normal keywords.',
            aliases=['astropy.io.fits.enabled_record_valued_keyword_cards'])

        extension_name_case_sensitive = _config.ConfigItem(
            False,
            'If True, extension names (i.e. the ``EXTNAME`` keyword) should be '
            'treated as case-sensitive.')
    conf = Conf()


Moving/renaming configuration items in Python
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``ConfigAlias`` objects can be used when a configuration item has been
moved from an astropy 0.3-style ``ConfigurationItem`` to an astropy
0.4-style ``ConfigItem`` inside of a ``ConfigNamespace``.

In the above example, the following adds backward-compatible hooks so
the old Python locations of the configuration items will continue to
work from user code::

    ENABLE_RECORD_VALUED_KEYWORD_CARDS = _config.ConfigAlias(
        '0.4', 'ENABLE_RECORD_VALUED_KEYWORD_CARDS',
        'enable_record_valued_keyword_cards')

Moving/renaming configuration items in the config file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If a configuration item is moved or renamed within the configuration
file, the ``aliases`` kwarg to ``ConfigItem`` can be used so that the
old location will continue to be used as a fallback.  For example, if
the old location of an item was:

.. code-block:: ini

    [coordinates.name_resolve]
    sesame_url = http://somewhere.com

One might want to drop the fact that that is implemented in the module
``name_resolve`` and just store the configuration in ``coordinates``:

.. code-block:: ini

    [coordinates]
    sesame_url = http://somewhere.com

When defining the ``ConfigItem`` for this entry, the ``aliases`` kwarg
can list the old location(s) of the configuration item::

    sesame_url = _config.ConfigItem(
        ["http://somewhere.com"],
        """Docstring""",
        aliases=['astropy.coordinates.name_resolve.sesame_url'])
