.. _astropy_vo:

*******************************************
Virtual Observatory Access (``astropy.vo``)
*******************************************

Introduction
============

The ``astropy.vo`` subpackage handles simple access for Virtual Observatory
(VO) services.

Currently, only Simple Cone Search Version 1.03 as defined in IVOA
Recommendation (February 22, 2008) is supported. Cone Search queries an
area encompassed by a given radius centered on a given RA and DEC and returns
all the objects found within the area in the given catalog.


.. _vo-sec-default-scs-services:

Default Cone Search Services
----------------------------

Currently, the default Cone Search services used are a subset of those found in
the STScI VAO Registry. They were hand-picked to represent commonly used
catalogs below:

    * 2MASS All-Sky
    * HST Guide Star Catalog
    * SDSS Data Release 7
    * SDSS-III Data Release 8
    * USNO A1
    * USNO A2
    * USNO B1

This subset undergoes daily validations hosted by STScI using
:ref:`vo-sec-validator-validate`. Those that pass without
critical warnings or exceptions are used by :ref:`vo-sec-client-scs` by
default.

If you are a Cone Search service provider and would like to include your
service in the list above, please contact ``help[at]stsci.edu``.


Caching
-------

Caching of downloaded contents is controlled by `astropy.utils.data`.
To *not* use cached data, some functions in this package have a ``cache``
keyword that can be set to `False`.


Getting Started
===============

This package has four main components across two subpackages:

    * ``astropy.vo.client``:
          * :ref:`vo-sec-client-vos` (``astropy.vo.client.vos_catalog``)
          * :ref:`vo-sec-client-scs` (``astropy.vo.client.conesearch``)
    * ``astropy.vo.validator``:
          * :ref:`vo-sec-validator-validate` (``astropy.vo.validator.validate``)
          * :ref:`vo-sec-validator-inspect` (``astropy.vo.validator.inspect``)

They are designed to be used in a work flow as illustrated below:

.. image:: images/astropy_vo_flowchart.png
    :width: 500px
    :alt: VO work flow

The one that a typical user needs is the :ref:`vo-sec-client-scs` component
(see :ref:`Cone Search Examples <vo-sec-scs-examples>`).


Using ``astropy.vo.client``
===========================

This subpackage contains modules supporting VO client-side operations.


.. _vo-sec-client-vos:

General VO Services Access
--------------------------

`astropy.vo.client.vos_catalog` contains common utilities for accessing
simple VO services.

.. _vo-sec-vos-config:

Configurable Items
^^^^^^^^^^^^^^^^^^

These parameters are set via :ref:`astropy_config`:

    * ``astropy.io.votable.table.PEDANTIC``
    * ``astropy.utils.data.REMOTE_TIMEOUT``
    * ``astropy.vo.client.vos_catalog.BASEURL``

Examples
^^^^^^^^

>>> from astropy.vo.client import vos_catalog

Get all catalogs from a database named 'conesearch_good'
(also see :ref:`Cone Search Examples <vo-sec-scs-examples>`):

>>> my_db = vos_catalog.get_remote_catalog_db('conesearch_good')
Downloading http://stsdas.stsci.edu/astrolib/vo_databases/conesearch_good.json
|============================================|  56/ 56k (100.00%)        00s
>>> my_db
<astropy.vo.client.vos_catalog.VOSDatabase at 0x2b0f3d0>
>>> print(str(my_db))
Guide Star Catalog 2.3 1
SDSS DR7 - Sloan Digital Sky Survey Data Release 7 1
SDSS DR7 - Sloan Digital Sky Survey Data Release 7 2
# ...
USNO-A2 Catalogue 1
USNO-A2.0 1

If you get timeout error, you need to use a custom timeout as follows:

>>> from astropy.utils.data import REMOTE_TIMEOUT
>>> with REMOTE_TIMEOUT.set_temp(30):
...     my_db = vos_catalog.get_remote_catalog_db('conesearch_good')

Find catalog names containing 'usno*a2':

>>> my_db.list_catalogs(pattern='usno*a2')
[u'The USNO-A2.0 Catalogue (Monet+ 1998) 1', u'USNO-A2 Catalogue 1']

Get information for a catalog titled 'USNO-A2 Catalogue 1':

>>> my_cat = my_db.get_catalog('USNO-A2 Catalogue 1')
>>> my_cat
<astropy.vo.client.vos_catalog.VOSCatalog at 0x1f78150>
>>> print(str(my_cat))
title: USNO-A2 Catalogue
url: http://www.nofs.navy.mil/cgi-bin/vo_cone.cgi?CAT=USNO-A2&
>>> print(my_cat.dumps())
{
    "capabilityClass": "ConeSearch", 
    "capabilityStandardID": "ivo://ivoa.net/std/ConeSearch", 
    "capabilityValidationLevel": "", 
    "contentLevel": "#University#Research#Amateur#", 
    # ...
    "version": "", 
    "waveband": "#Optical#"
}
>>> my_cat.keys()
[u'validate_network_error',
 u'capabilityClass',
 u'updated',
 # ...
 u'identifier',
 u'validate_xmllint']
>>> my_cat['url']
u'http://www.nofs.navy.mil/cgi-bin/vo_cone.cgi?CAT=USNO-A2&'
>>> my_cat['maxRadius']
1.0

One can also get information for a catalog using its URL.
If a URL yields multiple catalogs (this can happen when the service provider
re-register the URL with a different title), only the first match is returned:

>>> my_cat2 = my_db.get_catalog_by_url(
...     'http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=I/255/out&')
>>> print(my_cat2.dumps())
{
    "capabilityClass": "ConeSearch", 
    "capabilityStandardID": "ivo://ivoa.net/std/ConeSearch", 
    "capabilityValidationLevel": "", 
    "contentLevel": "#Research#", 
    # ...
    "version": "15-Sep-1999", 
    "waveband": "#Optical#"
}

To see validation warnings generated by :ref:`vo-sec-validator-validate`
for the catalog above:

>>> for w in my_cat2['validate_warnings']:
...     print(w)
/.../vo.xml:13:0: W22: The DEFINITIONS element is deprecated in VOTable 1.1...

To get all the matching catalogs by URL:

>>> matched_cats = [cat for key, cat in my_db.get_catalogs_by_url(
...     'http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=I/255/out&')]
>>> for c in matched_cats:
...     print(str(c))
title: The HST Guide Star Catalog, Version GSC-ACT (Lasker+ 1996-99)
url: http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=I/255/out&

To get all catalogs in the database:

>>> all_cats = [cat for key, cat in my_db.get_catalogs()]

To call a given VO service; In this case, a Cone Search
(also see :ref:`Cone Search Examples <vo-sec-scs-examples>`):

>>> result = vos_catalog.call_vo_service(
...     'conesearch_good', pedantic=False,
...     kwargs={'RA':6.088, 'DEC':-72.086, 'SR':0.5})
Trying http://wfaudata.roe.ac.uk/sdssdr7-dsa/DirectCone?DSACAT=SDSS_DR7&...
WARNING: W25: ... failed with: timed out [...]
# ...
Trying http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=I/243/out&
Downloading ... [Done]
WARNING: W22: ... The DEFINITIONS element is deprecated in VOTable 1.1...
>>> result
<astropy.io.votable.tree.Table at 0x6311f50>


.. _vo-sec-client-scs:

Simple Cone Search
------------------

`astropy.vo.client.conesearch` supports VO Simple Cone Search capabilities.

Available databases are generated on the server-side hosted by STScI using
:ref:`vo-sec-validator-validate`. The default database
(``astropy.vo.client.conesearch.CONESEARCH_DBNAME``), which can be changed
in :ref:`vo-sec-scs-config` below, is 'conesearch_good.json'.

In the default setting, it searches the good Cone Search services one by one,
stops at the first one that gives non-zero matches, and returns the results.
Since the list of services are extracted from a Python dictionary, the search
order might differ from call to call. :ref:`vo-sec-scs-examples` below also
show how to use non-default search behaviors.

.. note::

    Most services currently fail to parse when ``pedantic=True``.

.. warning::

    When Cone Search returns warnings, user should decide
    whether the results are reliable by inspecting the
    warning codes in `astropy.io.votable.exceptions`.

.. _vo-sec-scs-config:

Configurable Items
^^^^^^^^^^^^^^^^^^

These parameters are set via :ref:`astropy_config`:

    * ``astropy.utils.data.REMOTE_TIMEOUT``
    * ``astropy.vo.client.conesearch.CONESEARCH_DBNAME``
    * Also depends on
      :ref:`General VO Services Access Configurable Items <vo-sec-vos-config>`

.. _vo-sec-scs-examples:

Examples
^^^^^^^^

>>> from astropy.vo.client import conesearch

Shows a sorted list of Cone Search services to be searched
(to inspect them in detail, see :ref:`vo-sec-client-vos`):

>>> conesearch.list_catalogs()
[u'Guide Star Catalog 2.3 1',
 u'SDSS DR7 - Sloan Digital Sky Survey Data Release 7 1',
 u'SDSS DR7 - Sloan Digital Sky Survey Data Release 7 2',
 u'SDSS DR7 - Sloan Digital Sky Survey Data Release 7 3', ...,
 u'USNO-A2 Catalogue 1']

Perform cone search for :math:`0.5^{\circ}` radius around 47 Tucanae
(RA :math:`6.088^{\circ}`, DEC :math:`-72.086^{\circ}`) with minimum verbosity,
if supported. The first catalog in the database to
successfully return a result is used. If running this for
the first time, a copy of the catalogs database will be
downloaded to local cache. To run this again without
using cached data, set ``cache=False``:

>>> result = conesearch.conesearch(6.088, -72.086, 0.5, pedantic=False)
Trying http://wfaudata.roe.ac.uk/sdssdr7-dsa/DirectCone?DSACAT=SDSS_DR7&...
WARNING: W25: ... failed with: timed out [...]
# ...
Trying http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=I/243/out&
Downloading ... [Done]
WARNING: W22: ... The DEFINITIONS element is deprecated in VOTable 1.1...

To run the command above using custom timeout of
30 seconds for each Cone Search service query:

>>> from astropy.utils.data import REMOTE_TIMEOUT
>>> with REMOTE_TIMEOUT.set_temp(30):
...     result = conesearch.conesearch(6.088, -72.086, 0.5, pedantic=False)

Extract Numpy array containing the matched objects. See
`numpy` for available operations:

>>> cone_arr = result.array.data
>>> cone_arr
array([(0.499227, 4.465375, -72.106348, '0150-00089898'),
       (0.499716, 4.468892, -72.050834, '0150-00089994'),
       (0.498522, 4.469573, -72.120406, '0150-00090016'), ...,
       (0.499753, 7.709364, -72.059581, '0150-00229066'),
       (0.499739, 7.710762, -72.067145, '0150-00229140'),
       (0.49931, 7.71077, -72.107684, '0150-00229141')],
      dtype=[('_r', '<f8'), ('_RAJ2000', '<f8'), ('_DEJ2000', '<f8'),
             ('USNO-A1.0', '|S13')])
>>> cone_arr.dtype.names
('_r', '_RAJ2000', '_DEJ2000', 'USNO-A1.0')
>>> cone_arr.size
36386
>>> ra_list = cone_arr['_RAJ2000']
>>> ra_list
array([ 4.465375,  4.468892,  4.469573, ...,  7.709364,  7.710762,  7.71077 ])
>>> cone_arr[0]  # First row
(0.499227, 4.465375, -72.106348, '0150-00089898')
>>> cone_arr[-1]  # Last row
(0.49931, 7.71077, -72.107684, '0150-00229141')
>>> cone_arr[:10]  # First 10 rows
array([(0.499227, 4.465375, -72.106348, '0150-00089898'),
       (0.499716, 4.468892, -72.050834, '0150-00089994'),
       (0.498522, 4.469573, -72.120406, '0150-00090016'),
       (0.497264, 4.471928, -72.077298, '0150-00090076'),
       (0.49791, 4.472681, -72.059362, '0150-00090094'),
       (0.496641, 4.473367, -72.085139, '0150-00090125'),
       (0.496347, 4.474137, -72.092803, '0150-00090146'),
       (0.49673, 4.474284, -72.072045, '0150-00090153'),
       (0.4963, 4.474306, -72.094837, '0150-00090154'),
       (0.496058, 4.475092, -72.094639, '0150-00090173')],
      dtype=[('_r', '<f8'), ('_RAJ2000', '<f8'), ('_DEJ2000', '<f8'),
             ('USNO-A1.0', '|S13')])

Sort the matched objects by angular separation in ascending order:

>>> import numpy as np
>>> sep = cone_arr['_r']
>>> i_sorted = np.argsort(sep)
>>> cone_arr[i_sorted]
array([(0.069463, 6.311239, -72.09662, '0150-00163075'),
       (0.074441, 6.328667, -72.093992, '0150-00163906'),
       (0.074641, 6.330117, -72.09117, '0150-00163971'), ...,
       (0.499953, 6.769789, -72.541034, '0150-00186705'),
       (0.499986, 7.664375, -71.970134, '0150-00227015'),
       (0.499991, 7.211106, -72.4507, '0150-00207032')],
      dtype=[('_r', '<f8'), ('_RAJ2000', '<f8'), ('_DEJ2000', '<f8'),
             ('USNO-A1.0', '|S13')])

Result can also be manipulated as :ref:`astropy-io-votable`
and its unit can be manipulated as :ref:`astropy-units`.
In this example, we convert RA values from degree to arcsec:

>>> from astropy import units as u
>>> ra_field = result.get_field_by_id('_RAJ2000')
>>> ra_field.unit
Unit("deg")
>>> ra_field.unit.to(u.arcsec) * ra_list
array([ 16075.35  ,  16088.0112,  16090.4628, ...,  27753.7104,
        27758.7432,  27758.772 ])

Perform the same Cone Search as above but asynchronously using
`~astropy.vo.client.conesearch.AsyncConeSearch`.
Queries to individual Cone Search services are still governed by
``astropy.utils.data.REMOTE_TIMEOUT``. Although Cone Search is forced
to run in silent mode asynchronously, you will still see some
'Downloading ...' messages generated by
:func:`astropy.utils.data.download_file`:

>>> async_search = conesearch.AsyncConeSearch(
...     6.088, -72.086, 0.5, pedantic=False)

Check asynchronous search status:

>>> async_search.running()
True
>>> async_search.done()
False

Get search results after a 30-second wait (not to be
confused with ``astropy.utils.data.REMOTE_TIMEOUT`` that
governs individual Cone Search queries). If search is still not
done after 30 seconds, ``TimeoutError`` is raised. Otherwise,
Cone Search result is returned and can be manipulated as
above. If no ``timeout`` keyword given, it waits until
completion:

>>> async_result = async_search.get(timeout=30)
>>> cone_arr = async_result.array.data
>>> cone_arr.size
36386

Estimate the execution time and the number of objects for
the Cone Search service URL from above. The prediction naively
assumes a linear model, which might not be accurate for some cases.
It also uses the normal :func:`~astropy.vo.client.conesearch.conesearch`,
not the asynchronous version. This example uses a custom
timeout of 30 seconds:

>>> result.url
u'http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=I/243/out&'
>>> with REMOTE_TIMEOUT.set_temp(30):
...     t_est, n_est = conesearch.predict_search(
...         result.url, 6.088, -72.086, 0.5, pedantic=False, plot=True)
>>> t_est  # Predicted execution time
15.109336715274392
>>> n_est  # Predicted number of objects
33277

.. image:: images/client_predict_search_t.png
    :width: 450px
    :alt: Example plot from conesearch.predict_search (t)

.. image:: images/client_predict_search_n.png
    :width: 450px
    :alt: Example plot from conesearch.predict_search (n)

For debugging purpose, one can obtain the actual execution time
and number of objects, and compare them with the predicted values
above. However, keep in mind that running this for every prediction
would defeat the purpose of the prediction itself:

>>> t_real, tab = conesearch.conesearch_timer(
...     6.088, -72.086, 0.5, catalog_db=result.url, pedantic=False)
Trying http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=I/243/out&
Downloading ...
INFO: conesearch_timer took 9.04488611221 s on AVERAGE for 1 call(s). [...]
>>> t_real  # Actual execution time
9.044886112213135
>>> tab.array.size  # Actual number of objects
36386

For better control of which Cone Search service to use, one can
utilize the ``catalog_db`` keyword in
:func:`~astropy.vo.client.conesearch.conesearch`.
In this example, we look for those containing 'guide*star' in their titles
and only perform Cone Search using those services. As the first catalog in
the list to successfully return non-zero results is used, the order of
catalog names given in ``catalog_db`` is important:

>>> gsc_cats = conesearch.list_catalogs(pattern='guide*star')
>>> gsc_cats
[u'Guide Star Catalog 2.3 1',
 u'The HST Guide Star Catalog, Version 1.1 (Lasker+ 1992) 1',
 u'The HST Guide Star Catalog, Version 1.2 (Lasker+ 1996) 1',
 u'The HST Guide Star Catalog, Version GSC-ACT (Lasker+ 1996-99) 1']
>>> gsc_result = conesearch.conesearch(
...     6.088, -72.086, 0.5, catalog_db=gsc_cats, pedantic=False)
Trying http://gsss.stsci.edu/webservices/vo/ConeSearch.aspx?CAT=GSC23&
Trying http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=I/220/out&
Downloading ...
WARNING: W22: ... The DEFINITIONS element is deprecated in VOTable 1.1...
>>> gsc_result.array.size
2985
>>> gsc_result.url
u'http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=I/220/out&'

To repeat the Cone Search above with the services listed in a
different order:

>>> gsc_cats_reordered = [gsc_cats[i] for i in (3, 1, 2, 0)]
>>> gsc_cats_reordered
[u'The HST Guide Star Catalog, Version GSC-ACT (Lasker+ 1996-99) 1',
 u'The HST Guide Star Catalog, Version 1.1 (Lasker+ 1992) 1',
 u'The HST Guide Star Catalog, Version 1.2 (Lasker+ 1996) 1',
 u'Guide Star Catalog 2.3 1']
>>> gsc_result = conesearch.conesearch(
...     6.088, -72.086, 0.5, catalog_db=gsc_cats_reordered, pedantic=False)
Trying http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=I/255/out&
Downloading ...
WARNING: W22: ... The DEFINITIONS element is deprecated in VOTable 1.1...
>>> gsc_result.array.size
2985
>>> gsc_result.url
u'http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=I/255/out&'

One can also specify exactly which Cone Search service to use instead
of using a list. In this example, we only use the Cone Search service
titled 'The HST Guide Star Catalog, Version 1.2 (Lasker+ 1996) 1':

>>> gsc_result = conesearch.conesearch(
...     6.088, -72.086, 0.5,
...     catalog_db='The HST Guide Star Catalog, Version 1.2 (Lasker+ 1996) 1',
...     pedantic=False)
Trying http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=I/254/out&
Downloading ...
WARNING: W22: ... The DEFINITIONS element is deprecated in VOTable 1.1...
>>> gsc_result.array.size
2986
>>> gsc_result.url
u'http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=I/254/out&'

If one is unable to obtain any results using the default
Cone Search database, 'conesearch_good.json', that only contains
sites that cleanly passed validation, one can use :ref:`astropy_config`
to use another database, 'conesearch_warn.json', containing sites with
validation warnings. One should use these sites with caution:

>>> conesearch.CONESEARCH_DBNAME.set('conesearch_warn')
>>> conesearch.list_catalogs()
Downloading http://stsdas.stsci.edu/astrolib/vo_databases/conesearch_warn.json
|============================================|  94/ 94k (100.00%)        00s
[u'2MASS All-Sky Catalog of Point Sources (Cutri+ 2003) 1',
 u'2MASS All-Sky Point Source Catalog 1',
 u'Data release 7 of Sloan Digital Sky Survey catalogs 1',
 u'Data release 7 of Sloan Digital Sky Survey catalogs 2',
 u'Data release 7 of Sloan Digital Sky Survey catalogs 3', ...,
 u'USNO-B1 Catalogue 1']


Using ``astropy.vo.validator``
==============================

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

A master list of all available Cone Search services is
obtained from ``astropy.vo.validator.validate.CS_MSTR_LIST``, which
is a URL query to STScI VAO Registry by default.
However, by default, only the ones in ``astropy.vo.validator.validate.CS_URLS``
are validated (also see :ref:`vo-sec-default-scs-services`),
while the rest are skipped. There are also options to validate
a user-defined list of services or all of them.

All Cone Search queries are done using RA, DEC, and SR given by
``<testQuery>`` XML tag in the registry, and maximum verbosity.

The results are separated into 4 groups below. Each group
is stored as a JSON database:

    #. *conesearch_good.json*
           Passed validation without critical warnings and
           exceptions. This database residing in
           ``astropy.vo.client.vos_catalog.BASEURL`` is the one used
           by :ref:`vo-sec-client-scs` by default.
    #. *conesearch_warn.json*
           Has critical warnings but no exceptions. Users
           can manually set
           ``astropy.vo.client.conesearch.CONESEARCH_DBNAME``
           to use this at their own risk.
    #. *conesearch_exception.json*
           Has some exceptions. *Never* use this.
           For informational purpose only.
    #. *conesearch_error.json*
           Has network connection error. *Never* use this.
           For informational purpose only.

HTML pages summarizing the validation results are stored in
'results' sub-directory, which also contains downloaded XML
files from individual Cone Search queries.

Warnings and Exceptions
^^^^^^^^^^^^^^^^^^^^^^^

A subset of `astropy.io.votable.exceptions` that is considered
non-critical is defined by ``astropy.vo.validator.validate.NONCRIT_WARNINGS``,
which will not be flagged as bad by the validator. However,
this does not change the behavior of ``astropy.io.votable.table.PEDANTIC``,
which still needs to be set to `False` for them not to be thrown out
by :func:`~astropy.vo.client.conesearch.conesearch`.
Despite being listed as non-critical, user is responsible
to check whether the results are reliable; They should not be
used blindly.

Some
`units recognized by VizieR <http://cdsarc.u-strasbg.fr/vizier/Units.htx>`_
are considered invalid by Cone Search standards. As a result,
they will give the warning 'W50', which is non-critical by default.

User can also modify ``astropy.vo.validator.validate.NONCRIT_WARNINGS`` to
include or exclude any warnings or exceptions, as desired.
However, this should be done with caution. Adding exceptions
to non-critical list is not recommended.

Building the Database
^^^^^^^^^^^^^^^^^^^^^

Each Cone Search service is a catalog in the JSON database,
which is represented by a nested dictionary::

    {
        "__version__": 1,
        "catalogs": {
            "catalog 1": {
                "some key": "some value",
                # ...
            },
            "catalog 2": {
                "some key": "some value",
                # ...
            },
            # ...
        }
    }

In the master registry, there are duplicate catalog titles with
different access URLs, duplicate access URLs with different titles,
duplicate catalogs with slightly different descriptions, etc.

A Cone Search service is really defined by its access URL
regardless of title, description, etc. The validator ensures
each access URL is unique across all the output databases.
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
'duplicatesIgnored' dictionary key of the catalog kept in the
database.

All the existing catalog tags will be copied over as dictionary
keys, except 'accessURL' that is renamed to 'url' for simplicity.
In addition, new keys named 'validate_xxx' are added; 'xxx' will
be the original attribute names of
`astropy.io.votable.validator.result.Result`.

Configurable Items
^^^^^^^^^^^^^^^^^^

These parameters are set via :ref:`astropy_config`:

    * ``astropy.utils.data.REMOTE_TIMEOUT``
    * ``astropy.vo.validator.validate.CS_MSTR_LIST``
    * ``astropy.vo.validator.validate.CS_URLS``
    * ``astropy.vo.validator.validate.NONCRIT_WARNINGS``
    * Also depends on properties in
      :ref:`Simple Cone Search Configurable Items <vo-sec-scs-config>`

.. _vo-sec-validate-examples:

Examples
^^^^^^^^

>>> from astropy.vo.validator import validate

Validate default Cone Search sites with multiprocessing
and write results in the current directory. Reading the
master registry can be slow, so setting timeout to at least
30 seconds is recommended:

>>> from astropy.utils.data import REMOTE_TIMEOUT
>>> with REMOTE_TIMEOUT.set_temp(30):
...     validate.check_conesearch_sites()
Downloading http://vao.stsci.edu/directory/NVORegInt.asmx/...
WARNING: W20: None:2:0: W20: No version number specified in file...
# ...
INFO: Only 31/11144 site(s) are validated [astropy.vo.server.validate]
Downloading http://nvo.stsci.edu/vor10/getRecord.aspx?...
# ...
INFO: warn: 15 catalog(s) [astropy.vo.server.validate]
INFO: good: 15 catalog(s) [astropy.vo.server.validate]
INFO: nerr: 1 catalog(s) [astropy.vo.server.validate]
INFO: excp: 0 catalog(s) [astropy.vo.server.validate]
INFO: total: 31 catalog(s) [astropy.vo.server.validate]
INFO: Validation of 31 site(s) took 129.094 s [astropy.vo.server.validate]

From the master registry, select Cone Search access URLs
hosted by 'stsci.edu':

>>> import numpy as np
>>> from astropy.io.votable import parse_single_table
>>> from astropy.utils.data import get_readable_fileobj
>>> with REMOTE_TIMEOUT.set_temp(30):
...     with get_readable_fileobj(validate.CS_MSTR_LIST(),
...                               encoding='binary') as fd:
...         tab_all = parse_single_table(fd, pedantic=False)
Downloading http://vao.stsci.edu/directory/NVORegInt.asmx/...
WARNING: W20: None:2:0: W20: No version number specified in file...
# ...
>>> arr = tab_all.array.data[
...     np.where(tab_all.array['capabilityClass'] == b'ConeSearch')]
>>> urls = [s for s in arr['accessURL'] if b'stsci.edu' in s]
>>> urls
['http://archive.stsci.edu/hst/search.php?sci_data_set_name=Y*&amp;',
 'http://archive.stsci.edu/tues/search.php?',
 'http://archive.stsci.edu/hst/search.php?sci_data_set_name=J*&amp;',
 'http://archive.stsci.edu/hut/search.php?', ...,
 'http://archive.stsci.edu/kepler/kepler_fov/search.php?',
 'http://archive.stsci.edu/kepler/confirmed_planets/search.php?']

Validate only the URLs found above without verbose
outputs or multiprocessing, and write results in
'subset' sub-directory instead of the current directory:

>>> with REMOTE_TIMEOUT.set_temp(30):
...     validate.check_conesearch_sites(
...         destdir='./subset', verbose=False, parallel=False, url_list=urls)
Downloading http://vao.stsci.edu/directory/NVORegInt.asmx/...
WARNING: W49: ... Empty cell illegal for integer fields...
# ...
Downloading http://nvo.stsci.edu/vor10/getRecord.aspx?...
# ...

Add 'W24' from `astropy.io.votable.exceptions` to the list of
non-critical warnings to be ignored and re-run default validation.
This is *not* recommended unless you know exactly what you are doing:

>>> validate.NONCRIT_WARNINGS.set(validate.NONCRIT_WARNINGS() + ['W24'])
>>> with REMOTE_TIMEOUT.set_temp(30):
...     validate.check_conesearch_sites()

Reset the list of ignored warnings back to default value.
Validate *all* Cone Search services in the master registry
(this will take a while) and write results in 'all' sub-directory:

>>> validate.NONCRIT_WARNINGS.set(validate.NONCRIT_WARNINGS.defaultvalue)
>>> with REMOTE_TIMEOUT.set_temp(30):
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
:ref:`vo-sec-validator-validate`. It reads in JSON databases
residing in ``astropy.vo.client.vos_catalog.BASEURL``, which
can be changed to point to a different location.

Configurable Items
^^^^^^^^^^^^^^^^^^

This parameter is set via :ref:`astropy_config`:

    * ``astropy.vo.client.vos_catalog.BASEURL``

Examples
^^^^^^^^

Load Cone Search validation results from
``astropy.vo.client.vos_catalog.BASEURL``
(by default, the one used by :ref:`vo-sec-client-scs`):

>>> from astropy.vo.validator import inspect
>>> r = inspect.ConeSearchResults()
Downloading .../conesearch_good.json
|============================================|  56/ 56k (100.00%)        00s
Downloading .../conesearch_warn.json
|============================================|  94/ 94k (100.00%)        00s
Downloading .../conesearch_exception.json
|============================================|  45/ 45  (100.00%)        00s
Downloading .../conesearch_error.json
|============================================|   1/  1k (100.00%)        00s

Print tally. In this example, there are 15 Cone Search services that
passed validation with non-critical warnings, 15 with critical warnings,
none with exceptions, and 1 with network error:

>>> r.tally()
good: 15 catalog(s)
warn: 15 catalog(s)
exception: 0 catalog(s)
error: 1 catalog(s)
total: 31 catalog(s)

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

List Cone Search catalogs with warnings, excluding warnings that
were ignored in ``astropy.vo.validator.validate.NONCRIT_WARNINGS``,
and writes the output to a file named 'warn_cats.txt' in the current
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

Print the details of catalog titled 'USNO-A2 Catalogue 1':

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

Load Cone Search validation results from a local directory named 'subset'.
This is useful if you ran your own :ref:`vo-sec-validator-validate`
and wish to inspect the output databases. This example reads in
validation of STScI Cone Search services done in
:ref:`Validation for Simple Cone Search Examples <vo-sec-validate-examples>`:

>>> from astropy.vo.client.vos_catalog import BASEURL
>>> with BASEURL.set_temp('./subset/'):
>>>     r = inspect.ConeSearchResults()
>>> r.tally()
good: 21 catalog(s)
warn: 7 catalog(s)
exception: 0 catalog(s)
error: 0 catalog(s)
total: 28 catalog(s)
>>> r.catkeys['good']
[u'Advanced Camera for Surveys 1',
 u'Berkeley Extreme and Far-UV Spectrometer 1',
 u'Copernicus Satellite 1',
 u'Extreme Ultraviolet Explorer 1', ...,
 u'Wisconsin Ultraviolet Photo-Polarimeter Experiment 1']


See Also
========

- `NVO Directory <http://nvo.stsci.edu/vor10/index.aspx>`_

- `Simple Cone Search Version 1.03, IVOA Recommendation (22 February 2008) <http://www.ivoa.net/Documents/REC/DAL/ConeSearch-20080222.html>`_

- `STScI VAO Registry <http://vao.stsci.edu/directory/NVORegInt.asmx?op=VOTCapabilityPredOpt>`_

- `STScI VO Databases <http://stsdas.stsci.edu/astrolib/vo_databases/>`_


Reference/API
=============

.. automodapi:: astropy.vo.client.vos_catalog
   :no-inheritance-diagram:

.. automodapi:: astropy.vo.client.conesearch
   :no-inheritance-diagram:

.. automodapi:: astropy.vo.validator.validate
   :no-inheritance-diagram:

.. automodapi:: astropy.vo.validator.inspect
   :no-inheritance-diagram:
