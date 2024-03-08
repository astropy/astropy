.. currentmodule:: astropy.io.fits

..  _fits_io_cloud:

Obtaining subsets from cloud-hosted FITS files
**********************************************

Astropy offers support for extracting data from FITS files stored in the cloud.
Specifically, the `astropy.io.fits.open` function accepts the ``use_fsspec``
and ``fsspec_kwargs`` parameters, which allow remote files to be accessed in an
efficient way using the |fsspec| package.

``fsspec`` is an optional dependency of Astropy which supports reading
files from a range of remote and distributed storage backends, such as Amazon
and Google Cloud Storage.  This chapter explains its use.

.. note::

    The examples in this chapter require ``fsspec`` which is an optional
    dependency of Astropy.  See :ref:`installing-astropy` for details on
    installing optional dependencies.


Subsetting FITS files hosted on an HTTP web server
==================================================

A common use case for ``fsspec`` is to read subsets of FITS data from a web
server which supports serving partial files via the
`Range Requests <https://en.wikipedia.org/wiki/Byte_serving>`__ feature of the
HTTP protocol.  Most web servers support serving portions of files in this way.

For example, let's assume you want to retrieve data from a large image obtained
by the Hubble Space Telescope available at the following url::

    >>> # Download link for a large Hubble archive image (213 MB)
    >>> url = "https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:HST/product/j8pu0y010_drc.fits"

This file can be opened by passing the url to `astropy.io.fits.open`.
By default, Astropy will download the entire file to local disc before opening
it.  This works fine for small files but tends to require a lot of time and
memory for large files.

You can improve the performance for large files by passing the parameter
``use_fsspec=True`` to `open`.  This will make Astropy use ``fsspec``
to download only the necessary parts of the FITS file.
For example:

.. doctest-requires:: fsspec

    >>> from astropy.io import fits
    ...
    >>> # `fits.open` will download the primary header
    >>> with fits.open(url, use_fsspec=True) as hdul:  # doctest: +REMOTE_DATA
    ...
    ...     # Download a single header
    ...     header = hdul[1].header
    ...
    ...     # Download a single data array
    ...     image = hdul[1].data
    ...
    ...     # Download a 10-by-20 pixel cutout by using .section
    ...     cutout = hdul[2].section[10:20, 30:50]

The example above requires less time and memory than would be required to
download the entire file. This is because ``fsspec`` is able to leverage
two *lazy data loading* features available in Astropy:

1. The ``lazy_load_hdus`` parameter offered by `open` takes care of loading HDU
   header and data attributes on demand rather than reading all HDUs at once.
   This parameter is set to ``True`` by default.  You do not need to pass it
   explicitly, unless you changed its default value in the
   :ref:`astropy:astropy_config`.
2. The `ImageHDU.section` and `CompImageHDU.section` properties enables a
   subset of a data array to be read into memory without downloading the entire
   image or cube. See the :ref:`astropy:data-sections` part of the documentation
   for more details.

Additional tips for achieving good performance when working with remote files
are provided in the :ref:`astropy:optimizing_fsspec` section further down
this page.

.. note::

    The `ImageHDU.section` and `CompImageHDU.section` feature is only efficient
    for files that are not externally compressed (such as ``.fits.gz`` files).
    Files that are compressed using internal tile compression should work properly.
    Use ``.section`` on an externally compressed image will cause the whole FITS
    file to be downloaded.


Subsetting FITS files hosted in Amazon S3 cloud storage
=======================================================

The FITS file used in the example above also happens to be available via
Amazon cloud storage, where it is stored in a `public S3 bucket
<https://registry.opendata.aws/hst/>`__ at the following location::

    >>> s3_uri = "s3://stpubdata/hst/public/j8pu/j8pu0y010/j8pu0y010_drc.fits"

With ``use_fsspec`` enabled, you can obtain a small cutout from a file stored
in Amazon S3 cloud storage in the same way as above.  When opening paths with
prefix ``s3://`` (Amazon S3 Storage) or ``gs://`` (Google Cloud Storage),
`open` will automatically default to ``use_fsspec=True`` for convenience.
For example:

.. doctest-requires:: fsspec

    >>> # Download a small 10-by-20 pixel cutout from a FITS file stored in Amazon S3
    >>> with fits.open(s3_uri, fsspec_kwargs={"anon": True}) as hdul:  # doctest: +REMOTE_DATA
    ...     cutout = hdul[1].section[10:20, 30:50]

Obtaining cutouts from Amazon S3 in this way may be particularly performant if
your code is running on a server in the same Amazon cloud region as the data.

.. note::

    To open paths with prefix ``s3://``, fsspec requires an optional dependency
    called |s3fs|.  A ``ModuleNotFoundError`` will be raised if this dependency
    is missing. See :ref:`installing-astropy` for details on installing optional
    dependencies.


Working with Amazon S3 access credentials
-----------------------------------------

In the example above, we passed ``fsspec_kwargs={"anon": True}`` to enable the
data to be retrieved in an anonymous way without providing Amazon cloud access
credentials.  This is possible because the data is located in a public S3
bucket which has been configured to allow anonymous access.

In some cases you may want to access data stored in an Amazon S3 data bucket
that is private or uses the "Requester Pays" feature. You will have to provide
a secret access key in this case to avoid encountering a ``NoCredentialsError``.
You can use the ``fsspec_kwargs`` parameter to pass extra arguments, such as
access keys, to the `fsspec.open` function as follows:

.. doctest-skip::

    >>> fsspec_kwargs = {"key": "YOUR-SECRET-KEY-ID",
    ...                  "secret": "YOUR-SECRET-KEY"}
    >>> with fits.open(s3_uri, fsspec_kwargs=fsspec_kwargs) as hdul:
    ...     cutout = hdul[2].section[10:20, 30:50]

.. warning::

    Including secret access keys inside Python code is dangerous because you
    may accidentally end up revealing your keys when you share your code with
    others. A better practice is to store your access keys via a configuration
    file or environment variables. See the |s3fs| documentation for guidance.


Using :class:`~astropy.nddata.Cutout2D` with cloud-hosted FITS files
====================================================================

The examples above used the `ImageHDU.section` feature to download
small cutouts given a set of pixel coordinates. For astronomical images it is
often more convenient to obtain cutouts based on a sky position and angular
size rather than array coordinates. For this reason, Astropy provides the
`astropy.nddata.Cutout2D` tool which makes it easy to obtain cutouts informed
by an image's World Coordinate System (`~astropy.wcs.WCS`).

This cutout tool can be used in combination with ``fsspec`` and ``.section``.
For example, assume you happen to know that the image we opened above contains
a nice edge-on galaxy at the following position::

    >>> # Approximate location of an edge-on galaxy
    >>> from astropy.coordinates import SkyCoord
    >>> position = SkyCoord('10h01m41.13s 02d25m20.58s')

We also know that the radius of the galaxy is approximately 5 arcseconds::

    >>> # Approximate size of the galaxy
    >>> from astropy import units as u
    >>> size = 5*u.arcsec

Given this sky position and radius, we can use `~astropy.nddata.Cutout2D`
in combination with ``use_fsspec=True`` and ``.section`` as follows:

.. doctest-requires:: fsspec

    >>> from astropy.nddata import Cutout2D
    >>> from astropy.wcs import WCS
    ...
    >>> with fits.open(s3_uri, use_fsspec=True, fsspec_kwargs={"anon": True}) as hdul:  # doctest: +REMOTE_DATA
    ...     wcs = WCS(hdul[1].header)
    ...     cutout = Cutout2D(hdul[1].section,  # use `.section` rather than `.data`!
    ...                       position=position,
    ...                       size=size,
    ...                       wcs=wcs)

See :ref:`cutout_images` for more details on this feature.


..  _optimizing_fsspec:

Performance improvement tips for subsetting remote FITS files
=============================================================

In the examples above we explained that it is important to use the
``use_fsspec=True`` feature in combination with the ``lazy_load_hdus=True``
parameter and the ``ImageHDU.section`` feature to obtain good performance.

There are two additional factors which significantly impact the performance
you will encounter, namely: (i) the structure of the FITS file, and (ii) the caching
and block size configuration of ``fsspec``.  The remainder of this section
briefly explains these two factors.

Matching the FITS file structure to the data slicing patterns
-------------------------------------------------------------

The order in which multi-dimensional data is organized inside FITS files plays
a major role in the subsetting performance.

Astropy uses the row-major order for indexing FITS data. This means that the
right-most axis is the one that varies the fastest inside the file.
Put differently, the data for the right-most dimension tends to be located in
contiguous regions of the file and is therefore the easiest to extract.

For example, in the case of a 2D image, the slice ``.section[0, :]`` can be
obtained by downloading one contiguous region of bytes from the file.
In contrast, the slice ``.section[:, 0]`` requires accessing bytes spread
across the entire image array. The same is true for higher dimensions,
for example, obtaining the slice ``.section[0, :, :]`` from a 3D cube
will tend to be much faster than requesting ``.section[:, :, 0]``.

Obtaining slices of data that are well matched to the internal layout of
the FITS file generally yields the best performance.
If subsetting performance is important to you, you may have to consider
modifying your FITS files to ensure that the ordering of the dimensions
is well-matched to your data slicing patterns.

Configuring the ``fsspec`` block size and download strategy
-----------------------------------------------------------

The ``fsspec`` package supports different data reading and caching strategies
which aim to find a balance between the number of network requests on one hand
and the total amount of data transferred on the other hand.  By default, fsspec
will attempt to download data in large contiguous blocks using a buffered
*read ahead* strategy, similar to the strategy that is employed when operating
systems load local files into memory.

You can tune the performance of fsspec's buffering strategy by passing custom
``block_size`` and ``cache_type`` parameters to `fsspec.open`.  You can pass
these parameters via the ``fsspec_kwargs`` argument of `astropy.io.fits.open`.
For example, we can configure fsspec to make buffered reads with a minimum
``block_size`` of 1 MB as follows:

.. doctest-requires:: fsspec

    >>> fsspec_kwargs = {"block_size": 1_000_000, "cache_type": "bytes"}
    >>> with fits.open(url, use_fsspec=True, fsspec_kwargs=fsspec_kwargs) as hdul:  # doctest: +REMOTE_DATA
    ...     cutout = hdul[1].section[10:20, 30:50]

The ideal configuration will depend on the latency and throughput of the
network, as well as the exact shape and volume of the data you seek to obtain.

See the |fsspec| documentation for more information on its options.
