About
-----

This branch contains the configuration for ReadTheDocs to build archived
versions of the documentation that are otherwise no longer build on RTD.

Details
-------

`ReadTheDocs <http://readthedocs.org>`_ (RTD) allows documentation to be built
for specific tags as well as branches in a git repository. However, as of
October 2018, there is no way to include a custom ``robots.txt`` on the RTD
server, so old versions of the Astropy docs are still indexed by search engines.
Another solution apart from ``robots.txt`` is to add::

    <meta name="robots" content="noindex, nofollow">

to the header of each page, and this is what we have now done in the latest
version of Astropy (see
`astropy/astropy#7909 <https://github.com/astropy/astropy/pull/7909>`_). This
is now in effect for versions of the docs built from November 2018 onwards, but
there is no way to add this to existing RTD builds.

The solution we have adopted is to instead store the pre-October 2018 tagged
versions of the documentation in a dedicated repository
(https://github.com/astropy/archived-documentation), then to create this
``docs-archive`` branch and copy the archived versions of the docs into it, so
that old tagged versions can be accessed as e.g.
http://docs.astropy.org/en/docs-archive/v1.0.1/.
