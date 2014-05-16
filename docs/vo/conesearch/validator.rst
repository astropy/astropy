.. doctest-skip-all

Using `astropy.vo.validator`
============================

VO services validator is used by STScI to support :ref:`vo-sec-client-scs`.
Currently, only Cone Search services are supported.
A typical user should not need the validator. However, this could be used by
VO service providers to validate their services. Currently, any service
to be validated has to be registered in STScI VAO Registry.

.. _vo-sec-validator-validate:

Validation for Simple Cone Search
---------------------------------

`astropy.vo.validator.validate` validates VO services.
Currently, only Cone Search validation is done using
:func:`~astropy.vo.validator.validate.check_conesearch_sites`,
which utilizes underlying `astropy.io.votable.validator` library.

A master list of all available Cone Search services is obtained from
`astropy.vo.validator.Conf.conesearch_master_list`, which is a URL
query to STScI VAO Registry by default.  However, by default, only the
ones in `astropy.vo.validator.Conf.conesearch_urls` are validated
(also see :ref:`vo-sec-default-scs-services`), while the rest are
skipped. There are also options to validate a user-defined list of
services or all of them.

All Cone Search queries are done using RA, DEC, and SR given by
``<testQuery>`` XML tag in the registry, and maximum verbosity.
In an uncommon case where ``<testQuery>`` is not defined for a service,
it uses a default search for ``RA=0&DEC=0&SR=0.1``.

The results are separated into 4 groups below. Each group
is stored as a JSON file of `~astropy.vo.client.vos_catalog.VOSDatabase`:

#. ``conesearch_good.json`` Passed validation without critical
     warnings and exceptions. This database residing in
     `astropy.vo.Conf.vos_baseurl` is the one used by
     :ref:`vo-sec-client-scs` by default.
#. ``conesearch_warn.json`` Has critical warnings but no
     exceptions. Users can manually set
     `astropy.vo.Conf.conesearch_dbname` to use this at their own
     risk.
#. ``conesearch_exception.json``
     Has some exceptions. *Never* use this.
     For informational purpose only.
#. ``conesearch_error.json``
     Has network connection error. *Never* use this.
     For informational purpose only.

HTML pages summarizing the validation results are stored in
``'results'`` sub-directory, which also contains downloaded XML
files from individual Cone Search queries.

Warnings and Exceptions
^^^^^^^^^^^^^^^^^^^^^^^

A subset of `astropy.io.votable.exceptions` that is considered
non-critical is defined by
`astropy.vo.validator.Conf.noncritical_warnings`, which will not be
flagged as bad by the validator. However, this does not change the
behavior of `astropy.io.votable.Conf.pedantic`, which still needs to
be set to ``False`` for them not to be thrown out by
:func:`~astropy.vo.client.conesearch.conesearch`.  Despite being
listed as non-critical, user is responsible to check whether the
results are reliable; They should not be used blindly.

Some `units recognized by VizieR <http://cdsarc.u-strasbg.fr/vizier/Units.htx>`_
are considered invalid by Cone Search standards. As a result,
they will give the warning ``'W50'``, which is non-critical by default.

User can also modify `astropy.vo.validator.Conf.noncritical_warnings`
to include or exclude any warnings or exceptions, as desired.
However, this should be done with caution. Adding exceptions to
non-critical list is not recommended.

.. _vo-sec-validator-build-db:

Building the Database from Registry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each Cone Search service is a `~astropy.vo.client.vos_catalog.VOSCatalog` in
a `~astropy.vo.client.vos_catalog.VOSDatabase` (see
:ref:`vo-sec-client-cat-manip` and :ref:`vo-sec-client-db-manip`).

In the master registry, there are duplicate catalog titles with
different access URLs, duplicate access URLs with different titles,
duplicate catalogs with slightly different descriptions, etc.

A Cone Search service is really defined by its access URL
regardless of title, description, etc. By default,
:func:`~astropy.vo.client.vos_catalog.VOSDatabase.from_registry` ensures
each access URL is unique across the database.
However, for user-friendly catalog listing, its title will be
the catalog key, not the access URL.

In the case of two different access URLs sharing the same title,
each URL will have its own database entry, with a sequence number
appended to their titles (e.g., 'Title 1' and 'Title 2'). For
consistency, even if the title does not repeat, it will still be
renamed to 'Title 1'.

In the case of the same access URL appearing multiple times in
the registry, the validator will store the first catalog with
that access URL and throw out the rest. However, it will keep
count of the number of duplicates thrown out in the
``'duplicatesIgnored'`` dictionary key of the catalog kept in the
database.

All the existing catalog tags will be copied over as dictionary
keys, except ``'accessURL'`` that is renamed to ``'url'`` for simplicity.
In addition, new keys from validation are added:

* ``validate_expected``
    Expected validation result category, e.g., "good".
* ``validate_network_error``
    Indication for connection error.
* ``validate_nexceptions``
    Number of exceptions found.
* ``validate_nwarnings``
    Number of warnings found.
* ``validate_out_db_name``
    Cone Search database name this entry belongs to.
* ``validate_version``
    Version of validation software.
* ``validate_warning_types``
    List of warning codes.
* ``validate_warnings``
    Descriptions of the warnings.
* ``validate_xmllint``
    Indication of whether ``xmllint`` passed.
* ``validate_xmllint_content``
    Output from ``xmllint``.

Configurable Items
^^^^^^^^^^^^^^^^^^

These parameters are set via :ref:`astropy_config`:

* `astropy.vo.validator.Conf.conesearch_master_list`
    VO registry query URL that should return a VO table with all the desired
    VO services.
* `astropy.vo.validator.Conf.conesearch_urls`
    Subset of Cone Search access URLs to validate.
* `astropy.vo.validator.Conf.noncritical_warnings`
    List of VO table parser warning codes that are considered non-critical.

Also depends on properties in
:ref:`Simple Cone Search Configurable Items <vo-sec-scs-config>`.

.. _vo-sec-validate-examples:

Examples
^^^^^^^^

>>> from astropy.vo.validator import validate

Validate default Cone Search sites with multiprocessing and write results
in the current directory. Reading the master registry can be slow, so the
default timeout is internally set to 60 seconds for it. However,
``astropy.utils.data.REMOTE_TIMEOUT`` should still be set to account for
accessing the individual services (at least 30 seconds is recommended).
In addition, all VO table warnings from the registry are suppressed because
we are not trying to validate the registry itself but the services it contains:

>>> from astropy.utils import data
>>> with data.conf.set_temp('remote_timeout', 30):
...     validate.check_conesearch_sites()
Downloading http://vao.stsci.edu/directory/NVORegInt.asmx/...
|===========================================|  25M/ 25M (100.00%)        00s
INFO: Only 30/11938 site(s) are validated [astropy.vo.validator.validate]
# ...
INFO: good: 14 catalog(s) [astropy.vo.validator.validate]
INFO: warn: 12 catalog(s) [astropy.vo.validator.validate]
INFO: excp: 0 catalog(s) [astropy.vo.validator.validate]
INFO: nerr: 4 catalog(s) [astropy.vo.validator.validate]
INFO: total: 30 out of 30 catalog(s) [astropy.vo.validator.validate]
INFO: check_conesearch_sites took 451.05685997 s on AVERAGE...

Validate only Cone Search access URLs hosted by ``'stsci.edu'`` without verbose
outputs (except warnings that are controlled by :py:mod:`warnings`) or
multiprocessing, and write results in ``'subset'`` sub-directory instead of the
current directory. For this example, we use ``registry_db`` from
:ref:`VO database examples <vo-sec-client-db-manip-examples>`:

>>> urls = registry_db.list_catalogs_by_url(pattern='stsci.edu')
>>> urls
['http://archive.stsci.edu/befs/search.php?',
 'http://archive.stsci.edu/copernicus/search.php?', ...,
 'http://galex.stsci.edu/gxWS/ConeSearch/gxConeSearch.aspx?',
 'http://gsss.stsci.edu/webservices/vo/ConeSearch.aspx?CAT=GSC23&']
>>> with data.conf.set_temp('remote_timeout', 30):
...     validate.check_conesearch_sites(
...         destdir='./subset', verbose=False, parallel=False, url_list=urls)
INFO: check_conesearch_sites took 84.7241549492 s on AVERAGE...

Add ``'W24'`` from `astropy.io.votable.exceptions` to the list of
non-critical warnings to be ignored and re-run default validation.
This is *not* recommended unless you know exactly what you are doing:

>>> from astropy.vo.validator.validate import conf
>>> with conf.set_temp('noncritical_warnings', conf.noncritical_warnings + ['W24']):
...     with data.conf.set_temp('remote_timeout', 30):
...         validate.check_conesearch_sites()

Validate *all* Cone Search services in the master registry
(this will take a while) and write results in ``'all'`` sub-directory:

>>> with data.conf.set_temp('remote_timeout', 30):
...     validate.check_conesearch_sites(destdir='./all', url_list=None)

To look at the HTML pages of the validation results in the current
directory using Firefox browser (images shown are from STScI server
but your own results should look similar)::

    firefox results/index.html

.. image:: images/validator_html_1.png
    :width: 600px
    :alt: Main HTML page of validation results

When you click on 'All tests' from the page above, you will see all the
Cone Search services validated with a summary of validation results:

.. image:: images/validator_html_2.png
    :width: 600px
    :alt: All tests HTML page

When you click on any of the listed URLs from above, you will see
detailed validation warnings and exceptions for the selected URL:

.. image:: images/validator_html_3.png
    :width: 600px
    :alt: Detailed validation warnings HTML page

When you click on the URL on top of the page above, you will see
the actual VO Table returned by the Cone Search query:

.. image:: images/validator_html_4.png
    :width: 600px
    :alt: VOTABLE XML page


.. _vo-sec-validator-inspect:

Inspection of Validation Results
--------------------------------

`astropy.vo.validator.inspect` inspects results from
:ref:`vo-sec-validator-validate`. It reads in JSON files of
`~astropy.vo.client.vos_catalog.VOSDatabase`
residing in ``astropy.vo.Conf.vos_baseurl``, which
can be changed to point to a different location.

Configurable Items
^^^^^^^^^^^^^^^^^^

This parameter is set via :ref:`astropy_config`:

* `astropy.vo.Conf.vos_baseurl`

Examples
^^^^^^^^

>>> from astropy.vo.validator import inspect

Load Cone Search validation results from
``astropy.vo.Conf.vos_baseurl`` (by default, the one used by
:ref:`vo-sec-client-scs`):

>>> r = inspect.ConeSearchResults()
Downloading http://.../conesearch_good.json
|===========================================|  48k/ 48k (100.00%)        00s
Downloading http://.../conesearch_warn.json
|===========================================|  85k/ 85k (100.00%)        00s
Downloading http://.../conesearch_exception.json
|===========================================| 3.0k/3.0k (100.00%)        00s
Downloading http://.../conesearch_error.json
|===========================================| 4.0k/4.0k (100.00%)        00s

Print tally. In this example, there are 13 Cone Search services that
passed validation with non-critical warnings, 14 with critical warnings,
1 with exceptions, and 2 with network error:

>>> r.tally()
good: 13 catalog(s)
warn: 14 catalog(s)
exception: 1 catalog(s)
error: 2 catalog(s)
total: 30 catalog(s)

Print a list of good Cone Search catalogs, each with title, access URL,
warning codes collected, and individual warnings:

>>> r.list_cats('good')
Guide Star Catalog 2.3 1
http://gsss.stsci.edu/webservices/vo/ConeSearch.aspx?CAT=GSC23&
W48,W50
.../vo.xml:136:0: W50: Invalid unit string 'pixel'
.../vo.xml:155:0: W48: Unknown attribute 'nrows' on TABLEDATA
# ...
USNO-A2 Catalogue 1
http://www.nofs.navy.mil/cgi-bin/vo_cone.cgi?CAT=USNO-A2&
W17,W42,W21
.../vo.xml:4:0: W21: vo.table is designed for VOTable version 1.1 and 1.2...
.../vo.xml:4:0: W42: No XML namespace specified
.../vo.xml:15:15: W17: VOTABLE element contains more than one DESCRIPTION...

List Cone Search catalogs with warnings, excluding warnings that were
ignored in `astropy.vo.validator.Conf.noncritical_warnings`, and
writes the output to a file named ``'warn_cats.txt'`` in the current
directory. This is useful to see why the services failed validations:

>>> with open('warn_cats.txt', 'w') as fout:
...     r.list_cats('warn', fout=fout, ignore_noncrit=True)

List the titles of all good Cone Search catalogs:

>>> r.catkeys['good']
[u'Guide Star Catalog 2.3 1',
 u'SDSS DR7 - Sloan Digital Sky Survey Data Release 7 1',
 u'SDSS DR7 - Sloan Digital Sky Survey Data Release 7 2',
 u'SDSS DR7 - Sloan Digital Sky Survey Data Release 7 3', ...,
 u'USNO-A2 Catalogue 1']

Print the details of catalog titled ``'USNO-A2 Catalogue 1'``:

>>> r.print_cat('USNO-A2 Catalogue 1')
{
    "capabilityClass": "ConeSearch",
    "capabilityStandardID": "ivo://ivoa.net/std/ConeSearch",
    "capabilityValidationLevel": "",
    "contentLevel": "#University#Research#Amateur#",
    # ...
    "version": "",
    "waveband": "#Optical#"
}
Found in good

Load Cone Search validation results from a local directory named ``'subset'``.
This is useful if you ran your own :ref:`vo-sec-validator-validate`
and wish to inspect the output databases. This example reads in
validation of STScI Cone Search services done in
:ref:`Validation for Simple Cone Search Examples <vo-sec-validate-examples>`:

>>> from astropy.vo import conf
>>> with conf.set_temp('vos_baseurl', './subset/'):
>>>     r = inspect.ConeSearchResults()
>>> r.tally()
good: 19 catalog(s)
warn: 7 catalog(s)
exception: 2 catalog(s)
error: 0 catalog(s)
total: 28 catalog(s)
>>> r.catkeys['good']
[u'Advanced Camera for Surveys 1',
 u'Berkeley Extreme and Far-UV Spectrometer 1',
 u'Copernicus Satellite 1',
 u'Extreme Ultraviolet Explorer 1', ...,
 u'Wisconsin Ultraviolet Photo-Polarimeter Experiment 1']
