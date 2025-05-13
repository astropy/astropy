.. _cosmology_io_details:

*********************
Cosmology I/O Details
*********************


.. _cosmology_io_details_pickle:

Cosmology in a Pickle
=====================

For *temporary* storage an easy means to serialize and deserialize a Cosmology
object is using the :mod:`pickle` module. This is good for e.g. passing a
|Cosmology| between threads.

.. doctest-skip::

   >>> import pickle
   >>> from astropy.cosmology import Planck18
   >>> with open("planck18.pkl", mode="wb") as file:
   ...     # use protocol 5 to ensure byteorder is preserved (not switched to native)
   ...     pickle.dump(Planck18, file, protocol=5)
   >>> # and to read back
   >>> with open("planck18.pkl", mode="rb") as file:
   ...     cosmo = pickle.load(file)
   >>> cosmo
   FlatLambdaCDM(name="Planck18", ...

However this method has all the attendant drawbacks of :mod:`pickle` — security
vulnerabilities and non-human-readable files. Pickle files just generally don't
make for good persistent storage.

Solving both these issues, ``astropy`` provides a unified interface for reading
and writing data in different formats.


.. _cosmology_io_renaming_fields:

Renaming Fields
===============

Many I/O methods in :mod:`~astropy.cosmology` support renaming fields of the
|Cosmology| class when converting to a different format. This is done by
passing a ``rename`` dictionary to the ``to_format`` method.
Similarly, when converting from a different format, a ``rename`` dictionary
can be passed to the ``from_format`` method, mapping the fields of the input
to the fields of the |Cosmology| class.

For example, to rename the ``H0`` field to ``Hubble`` when converting to a table
format::

    >>> from astropy.cosmology import Cosmology, Planck18
    >>> renamed_table = Planck18.to_format("astropy.table", rename={"H0": "Hubble"})
    >>> renamed_table
    <QTable length=1>
        name      Hubble      Om0    Tcmb0    Neff      m_nu      Ob0
                km / (Mpc s)            K                 eV
        str8     float64    float64 float64 float64  float64[3] float64
    -------- ------------ ------- ------- ------- ----------- -------
    Planck18        67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.04897

    >>> cosmo = Cosmology.from_format(renamed_table, format="astropy.table",
    ...                               rename={"Hubble": "H0"})
    >>> cosmo == Planck18
    True


.. _cosmology_io_subclasses:

I/O from Subclasses
===================

When a subclass of |Cosmology| is used to read a file, the subclass will provide
a keyword argument ``cosmology=<class>`` to the registered read method. The
method uses this cosmology class, regardless of the class indicated in the
file, and sets parameters' default values from the class' signature.

.. doctest-skip::

    >>> from astropy.cosmology import FlatLambdaCDM
    >>> cosmo = FlatLambdaCDM.read('<file name>')
    >>> cosmo == Planck18
    True
