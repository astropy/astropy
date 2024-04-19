.. _cosmology_io-builtin:

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


.. _cosmology_io-builtin-cosmology:

Cosmology
*********

.. automodule:: astropy.cosmology.io._builtin.cosmology


.. _cosmology_io-builtin-mapping:

Mapping
*******

.. automodule:: astropy.cosmology.io._builtin.mapping


.. _cosmology_io-builtin-table:

Table
*****

.. automodule:: astropy.cosmology.io._builtin.table

.. _cosmology_io-builtin-model:


Model
*****

.. automodule:: astropy.cosmology.io._builtin.model



.. _cosmology_io-builtin-yaml:

YAML
****

.. automodule:: astropy.cosmology.io._builtin.yaml


.. _cosmology_io-builtin-row:

Row
***

.. automodule:: astropy.cosmology.io._builtin.row


.. _cosmology_io_builtin_readwrite:

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


.. _cosmology_io-builtin-ecsv:

ECSV
****

.. automodule:: astropy.cosmology.io._builtin.ecsv


.. _cosmology_io-builtin-latex:

LaTeX
*****

.. automodule:: astropy.cosmology.io._builtin.latex


.. _cosmology_io-builtin-html:

HTML
****

.. automodule:: astropy.cosmology.io._builtin.html
