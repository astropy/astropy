.. _cosmology_io_builtin:

**********************************
Built-in Cosmology To/From Formats
**********************************

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


.. _cosmology_io_builtin-cosmology:

Cosmology
*********

.. automodule:: astropy.cosmology._io.cosmology


.. _cosmology_io_builtin-mapping:

Mapping
*******

.. automodule:: astropy.cosmology._io.mapping


.. _cosmology_io_builtin-table:

Table
*****

.. automodule:: astropy.cosmology._io.table

.. _cosmology_io_builtin-model:


Model
*****

.. automodule:: astropy.cosmology._io.model



.. _cosmology_io_builtin-yaml:

YAML
****

.. automodule:: astropy.cosmology._io.yaml


.. _cosmology_io_builtin-row:

Row
***

.. automodule:: astropy.cosmology._io.row



*************************************
Built-in Cosmology Read/Write Formats
*************************************

To see a list of the available read/write file formats:

    >>> from astropy.cosmology import Cosmology
    >>> Cosmology.write.list_formats()
       Format   Read Write Auto-identify
    ----------- ---- ----- -------------
     ascii.ecsv  Yes   Yes           Yes
     ascii.html  Yes   Yes           Yes
    ascii.latex   No   Yes           Yes


.. _cosmology_io_builtin-ecsv:

ECSV
****

.. automodule:: astropy.cosmology._io.ecsv


.. _cosmology_io_builtin-latex:

LaTeX
*****

.. automodule:: astropy.cosmology._io.latex


.. _cosmology_io_builtin-html:

HTML
****

.. automodule:: astropy.cosmology._io.html
