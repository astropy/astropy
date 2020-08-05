.. _utils-data:

***************************************************
Downloadable Data Management (`astropy.utils.data`)
***************************************************

Introduction
============

A number of Astropy's tools work with data sets that are either awkwardly
large (e.g., `~astropy.coordinates.solar_system_ephemeris`) or
regularly updated (e.g., `~astropy.utils.iers.IERS_B`) or both
(e.g., `~astropy.utils.iers.IERS_A`). This kind of
data - authoritative data made available on the Web, and possibly updated
from time to time - is reasonably common in astronomy. The Astropy Project therefore
provides some tools for working with such data.

The primary tool for this is the ``astropy`` *cache*. This is a repository of
downloaded data, indexed by the URL where it was obtained. The tool
`~astropy.utils.data.download_file` and various other things built upon it can
use this cache to request the contents of a URL, and (if they choose to use the
cache) the data will only be downloaded if it is not already present in the
cache. The tools can be instructed to obtain a new copy of data
that is in the cache but has been updated online.

The ``astropy`` cache is stored in a centralized place (on Linux machines by
default it is ``$HOME/.astropy/cache``; see :ref:`astropy_config` for
more details).  You can check its location on your machine::

   >>> import astropy.config.paths
   >>> astropy.config.paths.get_cache_dir()  # doctest: +SKIP
   '/home/burnell/.astropy/cache'

This centralization means that the cache is persistent and shared between all
``astropy`` runs in any virtualenv by one user on one machine (possibly more if
your home directory is shared between multiple machines).  This can
dramatically accelerate ``astropy`` operations and reduce the load on servers,
like those of the IERS, that were not designed for heavy Web traffic. If you
find the cache has corrupted or outdated data in it, you can remove an entry or
clear the whole thing with `~astropy.utils.data.clear_download_cache`.

The files in the cache directory are named according to a cryptographic hash of
their URL (currently MD5, so in principle malevolent entities can cause
collisions, though the security risks this poses are marginal at most). The
modification times on these files normally indicate when they were last
downloaded from the Internet.

Usage Within Astropy
====================

For the most part, you can ignore the caching mechanism and rely on
``astropy`` to have the correct data when you need it. For example, precise
time conversions and sky locations need measured tables of the Earth's
rotation from the IERS. The table `~astropy.utils.iers.IERS_Auto` provides
the infrastructure for many of these calculations. It makes available
Earth rotation parameters, and if you request them for a time more recent
than its tables cover, it will download updated tables from the IERS. So
for example asking what time it is in UT1 (a timescale that reflects the
irregularity of the Earth's rotation) probably triggers a download of the
IERS data::

   >>> from astropy.time import Time
   >>> Time.now().ut1  # doctest: +SKIP
   Downloading https://maia.usno.navy.mil/ser7/finals2000A.all
   |============================================| 3.2M/3.2M (100.00%)         1s
   <Time object: scale='ut1' format='datetime' value=2019-09-22 08:39:03.812731>

But running it a second time does not require any new download::

   >>> Time.now().ut1  # doctest: +SKIP
   <Time object: scale='ut1' format='datetime' value=2019-09-22 08:41:21.588836>

Some data is also made available from the `Astropy data server`_ either
for use within ``astropy`` or for your convenience. These are available more
conveniently with the ``get_pkg_data_*`` functions::

   >>> from astropy.utils.data import get_pkg_data_contents
   >>> print(get_pkg_data_contents("coordinates/sites-un-ascii"))  # doctest: +SKIP
   # these are all mappings from the name in sites.json (which is ASCII-only) to the "true" unicode names
   TUBITAK->TÜBİTAK

Usage From Outside Astropy
==========================

Users of ``astropy`` can also make use of ``astropy``'s caching and downloading
mechanism. In its simplest form, this amounts to using
`~astropy.utils.data.download_file` with the ``cache=True`` argument to obtain
their data, from the cache if the data is there::

   >>> from astropy.utils.iers import IERS_B_URL, IERS_B
   >>> from astropy.utils.data import download_file
   >>> IERS_B.open(download_file(IERS_B_URL, cache=True))["year","month","day"][-3:]  # doctest: +SKIP
    <IERS_B length=3>
    year month  day
   int64 int64 int64
   ----- ----- -----
    2019     8     4
    2019     8     5
    2019     8     6

If users want to update the cache to a newer version of the
data (note that here the data was already up to date; users
will have to decide for themselves when to obtain new versions),
they can use the ``cache='update'`` argument::

   >>> IERS_B.open(download_file(IERS_B_URL,
   ...                           cache='update')
   ... )["year","month","day"][-3:]  # doctest: +SKIP
   Downloading http://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now
   |=========================================| 3.2M/3.2M (100.00%)         0s
   <IERS_B length=3>
    year month  day
   int64 int64 int64
   ----- ----- -----
    2019     8    18
    2019     8    19
    2019     8    20

If they are concerned that the primary source of the data may be
overloaded or unavailable, they can use the ``sources`` argument
to provide a list of sources to attempt downloading from, in order.
This need not include the original source. Regardless, the data
will be stored in the cache under the original URL requested::

   >>> f = download_file("ftp://ssd.jpl.nasa.gov/pub/eph/planets/bsp/de405.bsp",
   ...     cache=True,
   ...     sources=['https://data.nanograv.org/static/data/ephem/de405.bsp',
   ...              'ftp://ssd.jpl.nasa.gov/pub/eph/planets/bsp/de405.bsp'])  # doctest: +SKIP
   Downloading ftp://ssd.jpl.nasa.gov/pub/eph/planets/bsp/de405.bsp from https://data.nanograv.org/static/data/ephem/de405.bsp
   |========================================|  65M/ 65M (100.00%)        19s

.. _Astropy data server: https://www.astropy.org/astropy-data/

Cache Management
================

Because the cache is persistent, it is possible for it to become inconveniently
large, or become filled with irrelevant data. While it is simply a directory on
disk, each file is supposed to represent the contents of a URL, and many URLs
do not make acceptable on-disk filenames (for example, containing troublesome
characters like ":" and "~"). There is reason to worry that multiple
``astropy`` processes accessing the cache simultaneously might lead to cache
corruption. The data is therefore stored in a subdirectory named after the hash
of the URL, and write access is handled in a way that is resistant to
concurrency problems. So access to the cache is more convenient with a few
helpers provided by `~astropy.utils.data`.

If your cache starts behaving oddly you can use
`~astropy.utils.data.check_download_cache` to examine your cache contents and
raise an exception if it finds any anomalies.  If a single file is undesired or
damaged, it can be removed by calling
`~astropy.utils.data.clear_download_cache` with an argument that is the URL it
was obtained from, the filename of the downloaded file, or the hash of its
contents. Should the cache ever become badly corrupted,
`~astropy.utils.data.clear_download_cache` with no arguments will simply delete
the whole directory, freeing the space and removing any inconsistent data. Of
course, if you remove data using either of these tools, any processes currently
using that data may be disrupted (or, under Windows, deleting the cache may not
be possible until those processes terminate). So use
`~astropy.utils.data.clear_download_cache` with care.

To check the total space occupied by the cache, use
`~astropy.utils.data.cache_total_size`. The contents of the cache can be
listed with `~astropy.utils.data.get_cached_urls`, and the presence of a
particular URL in the cache can be tested with
`~astropy.utils.data.is_url_in_cache`. More general manipulations can be
carried out using `~astropy.utils.data.cache_contents`, which returns a
`~dict` mapping URLs to on-disk filenames of their contents.

If you want to transfer the cache to another computer, or preserve its contents
for later use, you can use the functions `~astropy.utils.data.export_download_cache` to
produce a ZIP file listing some or all of the cache contents, and
`~astropy.utils.data.import_download_cache` to load the ``astropy`` cache from such a
ZIP file.

The Astropy cache has changed format - once in the Python 2 to Python
3 transition, and again before Astropy version 4.0.2 to resolve some
concurrency problems that arose on some compute clusters. Each version of the
cache is in its own subdirectory, so the old versions do not interfere with the
new versions and vice versa, but their contents are not used by this version
and are not cleared by `~astropy.utils.data.clear_download_cache`. To remove
these old cache directories, you can run::

   >>> from shutil import rmtree
   >>> from os.path import join
   >>> from astropy.config.paths import get_cache_dir
   >>> rmtree(join(get_cache_dir(), 'download', 'py2'), ignore_errors=True)  # doctest: +SKIP
   >>> rmtree(join(get_cache_dir(), 'download', 'py3'), ignore_errors=True)  # doctest: +SKIP

Using Astropy With Limited or No Internet Access
================================================

You might want to use ``astropy`` on a telescope control machine behind a strict
firewall. Or you might be running continuous integration (CI) on your ``astropy``
server and want to avoid hammering astronomy servers on every pull request for
every architecture. Or you might not have access to US government or military
web servers. Whichever is the case, you may need to avoid ``astropy`` needing data
from the Internet. There is no simple and complete solution to this problem at
the moment, but there are tools that can help.

Exactly which external data your project depends on will depend on what parts
of ``astropy`` you use and how. The most general solution is to use a computer that
can access the Internet to run a version of your calculation that pulls in all of
the data files you will require, including sufficiently up-to-date versions of
files like the IERS data that update regularly. Then once the cache on this
connected machine is loaded with everything necessary, transport the cache
contents to your target machine by whatever means you have available, whether
by copying via an intermediate machine, portable disk drive, or some other
tool. The cache directory itself is somewhat portable between machines of the
same UNIX flavour; this may be sufficient if you can persuade your CI system to
cache the directory between runs. For greater portability, though, you can
simply use `~astropy.utils.data.export_download_cache` and
`~astropy.utils.data.import_download_cache`, which are portable and will allow
adding files to an existing cache directory.

If your application needs IERS data specifically, you can download the
appropriate IERS table, covering the appropriate time span, by any means you
find convenient. You can then load this file into your application and use the
resulting table rather than `~astropy.utils.iers.IERS_Auto`. In fact, the IERS
B table is small enough that a version (not necessarily recent) is bundled with
``astropy`` as ``astropy.utils.iers.IERS_B_FILE``. Using a specific non-automatic
table also has the advantage of giving you control over exactly which version
of the IERS data your application is using. See also :ref:`iers-working-offline`.

If your issue is with certain specific servers, even if they are the ones
``astropy`` normally uses, if you can anticipate exactly which files will be needed
(or just pick up after ``astropy`` fails to obtain them) and make those files
available somewhere else, you can request they be downloaded to the cache
using `~astropy.utils.data.download_file` with the ``sources`` argument set
to locations you know do work. You can also set ``sources`` to an empty list
to ensure that `~astropy.utils.data.download_file` does not attempt to use
the Internet at all.

If you have a particular URL that is giving you trouble, you can download it
using some other tool (e.g., ``wget``), possibly on another machine, and
then use `~astropy.utils.data.import_file_to_cache`.

Astropy Data and Clusters
=========================

Astronomical calculations often require the use of a large number of different
processes on different machines with a shared home filesystem. This can pose
certain complexities. In particular, if the many different processes attempt to
download a file simultaneously this can overload a server or trigger security
systems. The parallel access to the home directory can also trigger concurrency
problems in the Astropy data cache, though we have tried to minimize these. We
therefore recommend the following guidelines:

 * Write a simple script that sets ``astropy.utils.iers.conf.auto_download = True``
   and then accesses all cached resources your code will need, including source name
   lookups and IERS tables. Run it on the head node from time to time (frequently
   enough to beat the timeout ``astropy.utils.iers.conf.auto_max_age``, which
   defaults to 30 days) to ensure all data is up to date.

 * Make an Astropy config file (see :ref:`astropy_config`) that sets
   ``astropy.utils.iers.conf.auto_download = False`` so that the worker jobs will
   not suddenly notice an out-of-date table all at once and frantically attempt
   to download it.

 * Optionally, in this file, set ``astropy.utils.data.conf.remote_timeout = 0`` to
   prevent any attempt to download any file from the worker nodes; if you do this,
   you will need to override this setting in your script that does the actual
   downloading.

Now your worker nodes should not need to obtain anything from the Internet and
all should run smoothly.
