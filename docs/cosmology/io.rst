.. _read_write_cosmologies:

Read, Write, and Convert Cosmology Objects
******************************************

An easy way to serialize and deserialize a Cosmology object is using the
:mod:`pickle` module.

.. doctest-skip::

   >>> import pickle
   >>> from astropy.cosmology import Planck18
   >>> with open("planck18.pkl", mode="wb") as file:
   ...     pickle.dump(Planck18, file)
   >>> # and to read back
   >>> with open("planck18.pkl", mode="rb") as file:
   ...     cosmo = pickle.load(file)
   >>> cosmo
   FlatLambdaCDM(name="Planck18", ...

However this method has all the attendant drawbacks of :mod:`pickle` — security
vulnerabilities and non-human-readable files.

Solving both these issues, ``astropy`` provides a unified interface for reading
and writing data in different formats.


Getting Started
===============

The :class:`~astropy.cosmology.Cosmology` class includes two methods,
:meth:`~astropy.cosmology.Cosmology.read` and
:meth:`~astropy.cosmology.Cosmology.write`, that make it possible to read from
and write to files.

There are currently no built-in Cosmology readers nor writers, but custom
``read`` / ``write`` formats may be registered into the Astropy Cosmology I/O
framework.

Writing a cosmology instance requires only the file location and optionally,
if the file format cannot be inferred, a keyword argument "format". Additional
positional arguments and keyword arguments are passed to the reader methods.

.. doctest-skip::

    >>> from astropy.cosmology import Planck18
    >>> Planck18.write('<file name>')

Reading back the cosmology is most safely done from ``Cosmology``, the base
class, as it provides no default information and therefore requires the file
to have all necessary information to describe a cosmology.

.. doctest-skip::

    >>> from astropy.cosmology import Cosmology
    >>> cosmo = Cosmology.read('<file name>')
    >>> cosmo == Planck18
    True

When a subclass of ``Cosmology`` is used to read a file, the subclass will
provide a keyword argument ``cosmology=<class>`` to the registered read
method. The method uses this cosmology class, regardless of the class
indicated in the file, and sets parameters' default values from the class'
signature.

.. doctest-skip::

    >>> from astropy.cosmology import FlatLambdaCDM
    >>> cosmo = FlatLambdaCDM.read('<file name>')
    >>> cosmo == Planck18
    True
    
Reading and writing :class:`~astropy.cosmology.Cosmology` objects go through
intermediate representations, often a dict or `~astropy.table.QTable` instance.
These intermediate representations are accessible through the methods
:meth:`~astropy.cosmology.Cosmology.to_format` /
:meth:`~astropy.cosmology.Cosmology.from_format`.

To see the a list of the available formats:

    >>> from astropy.cosmology import Cosmology
    >>> Cosmology.to_format.list_formats()
    Format Read Write Auto-identify
    ------- ---- ----- -------------
    mapping  Yes   Yes           Yes

This list will include both built-in and registered 3rd-party formats.

:meth:`~astropy.cosmology.Cosmology.to_format` /
:meth:`~astropy.cosmology.Cosmology.from_format` parse a Cosmology to/from
another python object. This can be useful for e.g., iterating through an MCMC
of cosmological parameters or printing out a cosmological model to a journal
format, like latex or HTML. When 3rd party cosmology packages register with
Astropy' Cosmology I/O, ``to/from_format`` can be used to convert cosmology
instances between packages!

.. EXAMPLE START: Planck18 to mapping and back

.. code-block::

    >>> from astropy.cosmology import Planck18
    >>> cm = Planck18.to_format("mapping")
    >>> cm
    {'cosmology': <class 'astropy.cosmology.flrw.FlatLambdaCDM'>,
     'name': 'Planck18',
     'H0': <Quantity 67.66 km / (Mpc s)>,
     'Om0': 0.30966,
     ...

Now this dict can be used to load a new cosmological instance identical
to the ``Planck18`` cosmology from which it was created.

.. code-block::

    >>> from astropy.cosmology import Cosmology
    >>> cosmo = Cosmology.from_format(cm, format="mapping")
    >>> cosmo == Planck18
    True

.. EXAMPLE END



.. _custom_cosmology_converters:

Custom Cosmology To/From Formats
================================

Custom representation formats may also be registered into the Astropy Cosmology
I/O framework for use by these methods. For details of the framework see
:ref:`io_registry`.

.. EXAMPLE START : custom to/from format

As an example, the following is an implementation of an `astropy.table.Row`
converter. Note that we can use other registered parsers -- here "mapping"
-- to make the implementation much simpler.

We start by defining the function to parse a `astropy.table.Row` into a
`~astropy.cosmology.Cosmology`. This function should take 1 positional
argument, the row object, and 2 keyword arguments, for how to handle
extra metadata and which Cosmology class to use. Details about metadata
treatment are in ``Cosmology.from_format.help("mapping")``.

.. code-block:: python
    :emphasize-lines: 12,13

    >>> import copy
    >>> from astropy.cosmology import Cosmology
    
    >>> def from_table_row(row, *, move_to_meta=False, cosmology=None):
    ...     # get name from column
    ...     name = row['name'] if 'name' in row.columns else None
    ...     meta = copy.deepcopy(row.meta)
    ...     # turn row into mapping (dict of the arguments)
    ...     mapping = dict(row)
    ...     mapping["meta"] = meta
    ...     # build cosmology from map
    ...     return Cosmology.from_format(mapping,
    ...                                  move_to_meta=move_to_meta, cosmology=cosmology)

The next step is a function to perform the reverse operation: parse a
`~astropy.cosmology.Cosmology` into a `~astropy.table.Row`. This function
requires only the cosmology object and a ``*args`` to absorb unneeded
information passed by `astropy.io.registry.UnifiedReadWrite` (which
implements `astropy.cosmology.Cosmology.to_format`).

.. code-block:: python
    :emphasize-lines: 4

    >>> from astropy.table import QTable
    
    >>> def to_table_row(cosmology, *args):
    ...     p = cosmology.to_format("mapping")
    ...     p["cosmology"] = p["cosmology"].__qualname__  # as string
    ...     meta = p.pop("meta")
    ...     # package parameters into lists for Table parsing
    ...     params = {k: [v] for k, v in p.items()}
    ...     return QTable(params, meta=meta)[0]  # return row

Last we write a function to help with format auto-identification and then
register everything into `astropy.io.registry`.

.. code-block:: python
    :emphasize-lines: 11, 12, 13

    >>> from astropy.cosmology import Cosmology
    >>> from astropy.io import registry as io_registry
    >>> from astropy.table import Row

    >>> def row_identify(origin, format, *args, **kwargs):
    ...     """Identify if object uses the Table format."""
    ...     if origin == "read":
    ...         return isinstance(args[1], Row) and (format in (None, "row"))
    ...     return False

    >>> io_registry.register_reader("row", Cosmology, from_table_row)
    >>> io_registry.register_writer("row", Cosmology, to_table_row)
    >>> io_registry.register_identifier("row", Cosmology, row_identify)

Now the registered functions can be used in
:meth:`astropy.cosmology.Cosmology.from_format` and
:meth:`astropy.cosmology.Cosmology.to_format`.

.. code-block:: python

    >>> from astropy.cosmology import Planck18
    >>> row = Planck18.to_format("row")
    >>> row
    <Row index=0>
      cosmology     name        H0        Om0    Tcmb0    Neff    m_nu [3]    Ob0  
                           km / (Mpc s)            K                 eV            
        str13       str8     float64    float64 float64 float64   float64   float64
    ------------- -------- ------------ ------- ------- ------- ----------- -------
    FlatLambdaCDM Planck18        67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.04897

    >>> cosmo = Cosmology.from_format(row)
    >>> cosmo == Planck18  # test it round-trips
    True


.. doctest::
    :hide:

    >>> io_registry.unregister_reader("row", Cosmology)
    >>> io_registry.unregister_writer("row", Cosmology)
    >>> io_registry.unregister_identifier("row", Cosmology)
    >>> try:
    ...     io_registry.get_reader("row", Cosmology)
    ... except io_registry.IORegistryError:
    ...     pass

.. EXAMPLE END


.. _custom_cosmology_readers_writers:

Custom Cosmology Readers/Writers
================================

Custom ``read`` / ``write`` formats may be registered into the Astropy Cosmology
I/O framework. For details of the framework see :ref:`io_registry`.

.. EXAMPLE START : custom read/write

As an example, in the following we will fully work out a
`~astropy.cosmology.Cosmology` <-> JSON (de)serializer.
Note that we can use other registered parsers -- here "mapping"
-- to make the implementation much simpler.

We start by defining the function to parse JSON into a
`~astropy.cosmology.Cosmology`. This function should take 1 positional
argument, the file object or file path. We will also pass kwargs through to
:meth:`~astropy.cosmology.Cosmology.from_format`, which handles metadata
and which Cosmology class to use. Details of are in
``Cosmology.from_format.help("mapping")``.

.. code-block:: python
    :emphasize-lines: 9,12,17

    >>> import json, os
    >>> import astropy.units as u
    >>> from astropy.cosmology import Cosmology

    >>> def read_json(filename, **kwargs):
    ...     # read file, from path-like or file-like
    ...     if isinstance(filename, (str, bytes, os.PathLike)):
    ...         with open(filename, "r") as file:
    ...             data = file.read()
    ...     else:  # file-like : this also handles errors in dumping
    ...         data = filename.read()
    ...     mapping = json.loads(data)  # parse json mappable to dict
    ...     # deserialize Quantity
    ...     for k, v in mapping.items():
    ...         if isinstance(v, dict) and "value" in v and "unit" in v:
    ...             mapping[k] = u.Quantity(v["value"], v["unit"])
    ...     for k, v in mapping.get("meta", {}).items():  # also the metadata
    ...         if isinstance(v, dict) and "value" in v and "unit" in v:
    ...             mapping["meta"][k] = u.Quantity(v["value"], v["unit"])
    ...     return Cosmology.from_format(mapping, **kwargs)


The next step is a function to write a `~astropy.cosmology.Cosmology` to
JSON. This function requires the cosmology object and a file object/path.
We also require the boolean flag "overwrite" to set behavior for existing
files. Note that `~astropy.units.Quantity` is not natively compatible with
JSON. In both the ``write`` and ``read`` methods we have to create custom
parsers.

.. code-block:: python
    :emphasize-lines: 2,19

    >>> def write_json(cosmology, file, *, overwrite=False, **kwargs):
    ...    data = cosmology.to_format("mapping")  # start by turning into dict
    ...    data["cosmology"] = data["cosmology"].__name__  # change class field to str
    ...    # serialize Quantity
    ...    for k, v in data.items():
    ...        if isinstance(v, u.Quantity):
    ...            data[k] = {"value": v.value.tolist(), "unit": str(v.unit)}
    ...    for k, v in data.get("meta", {}).items():  # also serialize the metadata
    ...        if isinstance(v, u.Quantity):
    ...            data["meta"][k] = {"value": v.value.tolist(), "unit": str(v.unit)}
    ...
    ...    if isinstance(file, (str, bytes, os.PathLike)):
    ...        # check that file exists and whether to overwrite.
    ...        if os.path.exists(file) and not overwrite:
    ...            raise IOError(f"{file} exists. Set 'overwrite' to write over.")
    ...        with open(file, "w") as write_file:
    ...            json.dump(data, write_file)
    ...    else:
    ...        json.dump(data, file)

Last we write a function to help with format auto-identification and then
register everything into `astropy.io.registry`.

.. code-block:: python
    :emphasize-lines: 7,8,9

    >>> from astropy.io import registry as io_registry

    >>> def json_identify(origin, filepath, fileobj, *args, **kwargs):
    ...     """Identify if object uses the JSON format."""
    ...     return filepath is not None and filepath.endswith(".json")

    >>> io_registry.register_reader("json", Cosmology, read_json)
    >>> io_registry.register_writer("json", Cosmology, write_json)
    >>> io_registry.register_identifier("json", Cosmology, json_identify)
    
Now the registered functions can be used in
:meth:`astropy.cosmology.Cosmology.read` and
:meth:`astropy.cosmology.Cosmology.write`.

.. doctest-skip:: win32

    >>> import tempfile
    >>> from astropy.cosmology import Planck18
    >>>
    >>> file = tempfile.NamedTemporaryFile()
    >>> Planck18.write(file.name, format="json", overwrite=True)
    >>> with open(file.name) as f: f.readlines()
    ['{"cosmology": "FlatLambdaCDM", "name": "Planck18",
       "H0": {"value": 67.66, "unit": "km / (Mpc s)"}, "Om0": 0.30966,
       ...
    >>>
    >>> cosmo = Cosmology.read(file.name, format="json")
    >>> file.close()
    >>> cosmo == Planck18  # test it round-trips
    True


.. doctest::
    :hide:

    >>> io_registry.unregister_reader("json", Cosmology)
    >>> io_registry.unregister_writer("json", Cosmology)
    >>> io_registry.unregister_identifier("json", Cosmology)
    >>> try:
    ...     io_registry.get_reader("json", Cosmology)
    ... except io_registry.IORegistryError:
    ...     pass

.. EXAMPLE END


Reference/API
=============

.. automodapi:: astropy.cosmology.connect

.. automodapi:: astropy.cosmology.io.mapping
