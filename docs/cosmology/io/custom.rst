.. _cosmology_io_custom:

****************************
Custom Cosmology I/O formats
****************************

.. _custom_cosmology_converters:

Custom Cosmology To/From Formats
================================

Custom representation formats may also be registered into the Astropy Cosmology
I/O framework for use by these methods. For details of the framework see
:ref:`io_registry`. Note |Cosmology| ``to/from_format`` uses a custom registry,
available at ``Cosmology.<to/from>_format.registry``.

.. EXAMPLE START : custom to/from format

As an example, the following is an implementation of an |Row| converter. We can
and should use inbuilt parsers, like |QTable|, but to show a more complete
example we limit ourselves to only the "mapping" parser.

We start by defining the function to parse a |Row| into a |Cosmology|. This
function should take 1 positional argument, the row object, and 2 keyword
arguments, for how to handle extra metadata and which Cosmology class to use.
Details about metadata treatment are in
``Cosmology.from_format.help("mapping")``.

.. code-block:: python

    >>> import copy
    >>> from astropy.cosmology import Cosmology

    >>> def from_table_row(row, *, move_to_meta=False, cosmology=None):
    ...     # get name from column
    ...     name = row['name'] if 'name' in row.columns else None
    ...     meta = copy.deepcopy(row.meta)
    ...     # turn row into mapping (dict of the arguments)
    ...     mapping = dict(row)
    ...     mapping['name'] = name
    ...     mapping.setdefault("cosmology", meta.pop("cosmology", None))
    ...     mapping["meta"] = meta
    ...     # build cosmology from map
    ...     return Cosmology.from_format(mapping, move_to_meta=move_to_meta,
    ...                                  cosmology=cosmology)

.. code-block:: python

    convert_registry.register_reader("astropy.row", Cosmology, from_table_row)


The next step is a function to perform the reverse operation: parse a
|Cosmology| into a |Row|. This function requires only the cosmology object and
a ``*args`` to absorb unneeded information passed by
:class:`astropy.io.registry.UnifiedReadWrite` (which implements
|Cosmology.to_format|).

.. code-block:: python

    >>> from astropy.table import QTable

    >>> def to_table_row(cosmology, *args):
    ...     p = cosmology.to_format("mapping", cosmology_as_str=True)
    ...     meta = p.pop("meta")
    ...     # package parameters into lists for Table parsing
    ...     params = {k: [v] for k, v in p.items()}
    ...     return QTable(params, meta=meta)[0]  # return row

.. code-block:: python

    convert_registry.register_writer("astropy.row", Cosmology, to_table_row)


Last we write a function to help with format auto-identification and then
register everything into `astropy.io.registry`.

.. code-block:: python

    >>> from astropy.cosmology import Cosmology
    >>> from astropy.table import Row

    >>> def row_identify(origin, format, *args, **kwargs):
    ...     """Identify if object uses the Table format."""
    ...     if origin == "read":
    ...         return isinstance(args[1], Row) and (format in (None, "astropy.row"))
    ...     return False

.. code-block:: python

    convert_registry.register_identifier("astropy.row", Cosmology, row_identify)


Now the registered functions can be used in |Cosmology.from_format| and
|Cosmology.to_format|.

.. code-block:: python

    >>> from astropy.cosmology import Planck18
    >>> row = Planck18.to_format("astropy.row")
    >>> row
    <Row index=0>
      cosmology     name        H0        Om0    Tcmb0    Neff      m_nu      Ob0
                           km / (Mpc s)            K                 eV
        str13       str8     float64    float64 float64 float64  float64[3] float64
    ------------- -------- ------------ ------- ------- ------- ----------- -------
    FlatLambdaCDM Planck18        67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.04897

    >>> cosmo = Cosmology.from_format(row)
    >>> cosmo == Planck18  # test it round-trips
    True

.. EXAMPLE END


.. _custom_cosmology_readers_writers:

Custom Cosmology Readers/Writers
================================

Custom ``read`` / ``write`` formats may be registered into the Astropy
Cosmology I/O framework. For details of the framework see :ref:`io_registry`.
Note |Cosmology| ``read/write`` uses a custom registry, available at
``Cosmology.<read/write>.registry``.

.. EXAMPLE START : custom read/write

As an example, in the following we will fully work out a |Cosmology| <-> JSON
(de)serializer. Note that we can use other registered parsers -- here "mapping"
-- to make the implementation much simpler.

We start by defining the function to parse JSON into a |Cosmology|. This
function should take 1 positional argument, the file object or file path. We
will also pass kwargs through to |Cosmology.from_format|, which handles
metadata and which Cosmology class to use. Details of which are in
``Cosmology.from_format.help("mapping")``.

.. code-block:: python

    >>> import json, os
    >>> import astropy.units as u
    >>> from astropy.cosmology import Cosmology
    >>> from astropy.cosmology.io import readwrite_registry

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

    >>> readwrite_registry.register_reader("json", Cosmology, read_json)


The next step is a function to write a |Cosmology| to JSON. This function
requires the cosmology object and a file object/path. We also require the
boolean flag "overwrite" to set behavior for existing files. Note that
|Quantity| is not natively compatible with JSON. In both the ``write`` and
``read`` methods we have to create custom parsers.

.. code-block:: python

    >>> def write_json(cosmology, file, *, overwrite=False, **kwargs):
    ...    data = cosmology.to_format("mapping", cosmology_as_str=True)  # start by turning into dict
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

    >>> readwrite_registry.register_writer("json", Cosmology, write_json)

Last we write a function to help with format auto-identification and then
register everything into :mod:`astropy.io.registry`.

.. code-block:: python

    >>> def json_identify(origin, filepath, fileobj, *args, **kwargs):
    ...     """Identify if object uses the JSON format."""
    ...     return filepath is not None and filepath.endswith(".json")

    >>> readwrite_registry.register_identifier("json", Cosmology, json_identify)


Now the registered functions can be used in |Cosmology.read| and
|Cosmology.write|.

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

    >>> from astropy.io.registry import IORegistryError
    >>> readwrite_registry.unregister_reader("json", Cosmology)
    >>> readwrite_registry.unregister_writer("json", Cosmology)
    >>> readwrite_registry.unregister_identifier("json", Cosmology)
    >>> try:
    ...     readwrite_registry.get_reader("json", Cosmology)
    ... except IORegistryError:
    ...     pass

.. EXAMPLE END
