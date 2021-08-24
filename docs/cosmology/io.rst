.. doctest-skip-all

.. _read_write_cosmologies:

Read, Write, and Convert Cosmology Objects
******************************************

An easy way to serialize and deserialize a Cosmology object is using the
:mod:`pickle` module.

.. code-block:: python
   :emphasize-lines: 1,4,8

   >>> import pickle
   >>> from astropy.cosmology import Planck18
   >>> with open("planck18.pkl", mode="wb") as file:
   ...     pickle.dump(file, Planck18)

   And to read the file back
   >>> with open("planck18.pkl", mode="rb") as file:
   >>>     cosmo = pickle.load(file)
   >>> cosmo

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
::

    >>> from astropy.cosmology import Planck18
    >>> Planck18.write('[file name]')

Reading back the cosmology is most safely done from ``Cosmology``, the base
class, as it provides no default information and therefore requires the file
to have all necessary information to describe a cosmology.

::

    >>> from astropy.cosmology import Cosmology
    >>> cosmo = Cosmology.read('[file name]')
    >>> cosmo == Planck18
    True

When a subclass of ``Cosmology`` is used to read a file, the subclass will
provide a keyword argument ``cosmology=[class]`` to the registered read
method. The method uses this cosmology class, regardless of the class
indicated in the file, and sets parameters' default values from the class'
signature.

::

    >>> from astropy.cosmology import FlatLambdaCDM
    >>> cosmo = FlatLambdaCDM.read('[file name]')
    >>> cosmo == Planck18
    True
    
Reading and writing :class:`~astropy.cosmology.Cosmology` objects go through
intermediate representations, often a dict or `~astropy.table.QTable` instance.
These intermediate representations are accessible through the methods
``to_format`` / ``from_format``.

.. EXAMPLE START: Planck18 to mapping and back

    >>> from astropy.cosmology import Planck18
    >>> cm = Planck18.to_format("mapping")
    >>> cm
    {'cosmology': astropy.cosmology.core.FlatLambdaCDM,
     'name': 'Planck18',
     'H0': <Quantity 67.66 km / (Mpc s)>,
     'Om0': 0.30966,
     ...

    Now this dict can be used to load a new cosmological instance identical
    to the ``Planck18`` cosmology from which it was created.

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
    converter. Note that we can use other registered parsers -- here ``mapping``
    -- to make the implementation much simpler.

    We start by defining the function to parse a `astropy.table.Row` into a
    `~astropy.cosmology.Cosmology`. This function should take 1 positional
    argument, the row object, and 2 keyword arguments, for how to handle
    extra metadata and which Cosmology class to use. Details of each are in
    ``from_mapping`` or ``Cosmology.from_format.help("mapping")``.

    .. code-block:: python
    
        from astropy.cosmology.io.mapping import from_mapping, to_mapping
    
        def from_table_row(row, *, move_to_meta=False, cosmology=None):
            # get name from column
            name = row['name'] if 'name' in row.columns else None
            meta = copy.deepcopy(row.meta)
            # turn row into mapping (dict of the arguments)
            mapping = dict(row)
            mapping["cosmology"] = meta.pop("cosmology")
            mapping["meta"] = meta
            # build cosmology from map
            return from_mapping(mapping, move_to_meta=move_to_meta, cosmology=cosmology)

    The next step is a function to perform the reverse operation: parse a
    `~astropy.cosmology.Cosmology` into a `astropy.table.Row`. This function
    only the cosmology object and also ``*args`` to absorb unneeded information
    passed by `astropy.io.registry.UnifiedReadWrite`, which implements
    `astropy.cosmology.Cosmology.to_format`.

    .. code-block:: python

        from astropy.cosmology.io.mapping import to_mapping

        def to_table_row(cosmology, *args):
            # start by getting a map representation
            p = to_mapping(cosmology)
            # create metadata from mapping
            meta = p.pop("meta")
            meta["cosmology"] = p.pop("cosmology").__name__
            # package parameters into lists for Table parsing
            params = {k: [v] for k, v in p.items()}
            return QTable(params, meta=meta)[0]  # return row

    Last we write a function to help with format auto-identification and then
    register everything into `astropy.io.registry`.

    .. code-block:: python

        froma astropy.cosmology import Cosmology
        from astropy.io import registry as io_registry
        from astropy.table import Row


        def row_identify(origin, format, *args, **kwargs):
            """Identify if object uses the Table format."""
            if origin == "read":
                return isinstance(args[1], Row) and (format in (None, "row"))
            return False
    
        # register the methods
        io_registry.register_reader("row", Cosmology, from_table_row)
        io_registry.register_writer("row", Cosmology, to_table_row)
        io_registry.register_identifier("row", Cosmology, row_identify)

.. EXAMPLE END


.. _custom_cosmology_readers_writers:

Custom Cosmology Readers/Writers
================================

Custom ``read`` / ``write`` formats may be registered into the Astropy Cosmology
I/O framework. For details of the framework see :ref:`io_registry`.
As a quick example, outlining how to make a
`~astropy.cosmology.Cosmology` <-> JSON serializer:

.. code-block:: python

   import json
   import os
   from astropy.cosmology import Cosmology
   from astropy.io import registry as io_registry

   def read_json(filename, **kwargs):
       # read file, from path-like or file-like
       if isinstance(filename, (str, bytes, os.PathLike)):  # pathlike
           with open(filename, "r") as file:
               data = file.read()
       else:  # file-like
           data = filename.read()

       mapping = json.loads(data)  # parse json to dict
       return Cosmology.from_format(mapping, move_to_meta=move_to_meta)

   def write_json(cosmology, file, overwrite=False, **kwargs):
       data = cosmology.to_format("mapping")  # start by turning into dict
       data["cosmology"] = data["cosmology"].__name__  # change class to str

       if isinstance(file, _all_pathlike):  # pathlike
           # check that file exists and whether to overwrite.
           if os.path.exists(file) and not overwrite:
               raise IOError(f"{file} exists. Set 'overwrite' to write over.")
           with open(file, "w") as write_file:
               json.dump(data, write_file)
       else:  # file-like or error (this handles errors in dumping)
           json.dump(data, file)

   def json_identify(origin, filepath, fileobj, *args, **kwargs):
       """Identify if object uses the JSON format."""
       return filepath is not None and filepath.endswith(".json")

   # register the read/write methods
   io_registry.register_reader("json", Cosmology, read_json)
   io_registry.register_writer("json", Cosmology, write_json)
   io_registry.register_identifier("json", Cosmology, json_identify)


Reference/API
=============

.. automodapi:: astropy.cosmology.connect
