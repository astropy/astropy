.. doctest-skip-all

.. _read_write_cosmologies:

Reading and Writing Cosmology Objects
*************************************

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


.. _custom_cosmology_readers_writers:

Custom Cosmology Readers/Writers
================================

Custom ``read`` / ``write`` formats may be registered into the Astropy Cosmology
I/O framework. For details of the framework see :ref:`io_registry`.
As a quick example, outlining how to make a
`~astropy.cosmology.Cosmology` <-> JSON serializer:

.. code-block:: python

   import json
   from astropy.cosmology import Cosmology
   from astropy.io import registry as io_registry

   def read_json(filename, **kwargs):
       cosmology = ...  # read and parse file
       return cosmology

   def write_json(cosmology, file, **kwargs):
       data = ...  # parse cosmology to dict
       with open(file, "w") as write_file:
            json.dump(data, write_file)

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
