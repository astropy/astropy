.. _cosmology_io:

***********************
Cosmology I/O & Convert
***********************

.. _cosmology_io_introduction:

Introduction
============

The ``astropy.cosmology.io`` package provides a unified interface for reading, writing,
and converting |Cosmology| objects. Most of the features of this package are accessible
through the |Cosmology| class, which provides the methods |Cosmology.read|,
|Cosmology.write|, |Cosmology.to_format|, and |Cosmology.from_format| for reading,
writing, and converting |Cosmology| objects, respectively.


Getting Started: Reading and Writing
====================================

The |Cosmology| class includes two methods, |Cosmology.read| and |Cosmology.write|, that
make it possible to read from and write to files.

The registered ``read`` / ``write`` formats include "ascii.ecsv" and "ascii.html", like
for Table. Also, custom ``read`` / ``write`` formats may be registered into the Astropy
Cosmology I/O framework. For more information on the built-in formats, see
:ref:`cosmology_io_builtin`, or :ref:`cosmology_io_custom` for information on
registering custom formats.

Writing a cosmology instance requires only the file location and optionally, or if the
file format cannot be inferred, a keyword argument "format". Additional positional
arguments and keyword arguments are passed to the relevant writer methods.

.. doctest-skip::

    >>> from astropy.cosmology import Planck18
    >>> Planck18.write("example_cosmology.ecsv", format="ascii.ecsv")


Reading back the cosmology is done from |Cosmology|.

.. doctest-skip::

    >>> from astropy.cosmology import Cosmology
    >>> cosmo = Cosmology.read("example_cosmology.ecsv", format="ascii.ecsv")
    >>> cosmo == Planck18
    True

More specific cosmology classes can be used, providing default values for missing
information (see :ref:`cosmology_io_details`). However using the base class is safest as
it provides no default information and therefore requires the file to have all necessary
information to describe a cosmology.

To see a list of the available read/write file formats:

.. code-block:: python

    >>> from astropy.cosmology import Cosmology
    >>> Cosmology.read.list_formats()
      Format   Read Write Auto-identify
    ---------- ---- ----- -------------
    ascii.ecsv  Yes   Yes           Yes
    ascii.html  Yes   Yes           Yes

    >>> Cosmology.write.list_formats()
       Format   Read Write Auto-identify
    ----------- ---- ----- -------------
     ascii.ecsv  Yes   Yes           Yes
     ascii.html  Yes   Yes           Yes
    ascii.latex   No   Yes           Yes

This list will include both built-in and registered 3rd-party formats.


Getting Started: Converting Formats
===================================

Reading and writing |Cosmology| objects go through intermediate representations, often a
dict or |QTable| instance. These intermediate representations are accessible through the
methods |Cosmology.to_format| / |Cosmology.from_format|.

|Cosmology.to_format| / |Cosmology.from_format| parse a Cosmology to/from
another python object. This can be useful for e.g., iterating through an MCMC
of cosmological parameters or printing out a cosmological model to a journal
format, like latex or HTML. When 3rd party cosmology packages register with
Astropy's Cosmology I/O, ``to/from_format`` can be used to convert cosmology
instances between packages!

.. EXAMPLE START: Planck18 to QTable and back

Another pre-registered format is "table", for converting a |Cosmology| to and
from a |QTable|.

.. code-block::

    >>> from astropy.cosmology import Planck18
    >>> ct = Planck18.to_format("astropy.table")
    >>> ct
    <QTable length=1>
      name        H0        Om0    Tcmb0    Neff      m_nu      Ob0
             km / (Mpc s)            K                 eV
      str8     float64    float64 float64 float64  float64[3] float64
    -------- ------------ ------- ------- ------- ----------- -------
    Planck18        67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.04897

Now this |QTable| can be used to load a new cosmological instance identical to
the |Planck18| cosmology from which it was created.

.. code-block::

    >>> cosmo = Cosmology.from_format(ct, format="astropy.table")
    >>> print(cosmo)
    FlatLambdaCDM(name="Planck18", H0=67.66 km / (Mpc s), Om0=0.30966,
                  Tcmb0=2.7255 K, Neff=3.046, m_nu=[0. 0. 0.06] eV, Ob0=0.04897)

Perhaps most usefully, |QTable| itself has ``read/write`` methods with numerous options,
e.g. FITS, that now work with |Cosmology|.

.. EXAMPLE END

To see the a list of the available conversion formats:

.. code-block:: python

    >>> from astropy.cosmology import Cosmology
    >>> Cosmology.to_format.list_formats()
          Format      Read Write Auto-identify
    ----------------- ---- ----- -------------
    astropy.cosmology  Yes   Yes           Yes
        astropy.model  Yes   Yes           Yes
          astropy.row  Yes   Yes           Yes
        astropy.table  Yes   Yes           Yes
              mapping  Yes   Yes           Yes
                 yaml  Yes   Yes            No

This list will include both built-in and registered 3rd-party formats.


Using Cosmology I/O
===================

More details are provided in the following pages:

.. toctree::
   :maxdepth: 2

   details.rst
   builtin.rst
   custom.rst


Reference/API
=============

.. automodapi:: astropy.cosmology.connect
