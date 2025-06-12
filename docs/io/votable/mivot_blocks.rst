.. doctest-skip-all

Reading and writing VO model annotations (MIVOT)
------------------------------------------------

A ``RESOURCE`` element can be a ``MIVOT`` block since VOTable version 1.5.

Introduction
^^^^^^^^^^^^
Model Instances in VOTables (`MIVOT <https://ivoa.net/documents/MIVOT/20230620/REC-mivot-1.0.pdf>`_)
defines a syntax to map VOTable data to any model serialised in VO-DML (Virtual Observatory Data Modeling Language).
This annotation schema operates as a bridge between data and the models. It associates both column/param metadata and data
from the VOTable to the data model elements (class, attributes, types, etc.). It also brings up VOTable data or
metadata that were possibly missing in the table, e.g., coordinate system description, or curation tracing.
The data model elements are grouped in an independent annotation block complying with the MIVOT XML schema which
is added as an extra resource above the table element.
The MIVOT syntax allows to describe a data structure as a hierarchy of classes.
It is also able to represent relations and compositions between them. It can moreover build up data model objects by
aggregating instances from different tables of the VOTable.

Astropy implementation
^^^^^^^^^^^^^^^^^^^^^^
The purpose of Astropy is not to process VO annotations.
It is just to allow related packages to get and set MIVOT blocks from/into VOTables.
For this reason, in this implementation MIVOT annotations are both imported and exported as strings.
The current implementation prevents client code from injecting into VOTables strings
that are not MIVOT serializations.

MivotBlock implementation:

- MIVOT blocks are handled by the :class:`astropy.io.votable.tree.MivotBlock` class.
- A MivotBlock instance can only be carried by a resource with "type=meta".
- This instance holds the XML mapping block as a string.
- MivotBlock objects are instanced by the Resource parser.
- The MivotBlock class has its own logic that operates both parsing and IO functionalities.

Example
"""""""

.. code-block:: xml

       <VOTABLE xmlns="http://www.ivoa.net/xml/VOTable/v1.3"
         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" version="1.3">
         <RESOURCE>
           <RESOURCE type="meta">
             <VODML xmlns="http://www.ivoa.net/xml/mivot">
              ...
             </VODML>
           </RESOURCE>
           <TABLE name="myDataTable">
            ....
           </TABLE>
         </RESOURCE>
    </VOTABLE>

Reading a VOTable containing a MIVOT block
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To read in a VOTable file containing or not a MIVOT Resource, pass a file path to `~astropy.io.votable.parse`:

.. code-block:: python

   >>> from astropy.io.votable import parse
   >>> from astropy.utils.data import get_pkg_data_filename
   >>> votable = parse(get_pkg_data_filename("data/test.order.xml", package="astropy.io.votable.tests"))
   <VODML xmlns="http://www.ivoa.net/xml/mivot">
   </VODML>

   <VODML xmlns="http://www.ivoa.net/xml/mivot">
   </VODML>

The parse function will call the MIVOT parser if it detects a MIVOT block.

Building a Resource containing a MIVOT block
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Construct the MIVOT block by passing the XML block as a parameter:

.. code-block:: python

   >>> from astropy.io.votable import tree
   >>> from astropy.io.votable.tree import MivotBlock, Resource, VOTableFile
   >>> mivot_block = MivotBlock("""
   <VODML xmlns="http://www.ivoa.net/xml/mivot" >
      <REPORT status="OK">Unit test mapping block</REPORT>
      <GLOBALS>  </GLOBALS>
   </VODML>
   """)

Build a new resource:

.. code-block:: python

   >>> mivot_resource = Resource()

Give it the type meta:

.. code-block:: python

   >>> mivot_resource.type = "meta"

Then add it the MIVOT block:

.. code-block:: python

   >>> mivot_resource.mivot_block = mivot_block

Now you have a MIVOT resource that you can add to an object Resource creating a new Resource:

.. code-block:: python

   >>> votable = VOTableFile()
   >>> r1 = Resource()
   >>> r1.type = "results"
   >>> r1.resources.append(mivot_resource)

You can add an `astropy.io.votable.tree.TableElement` to the resource:

.. code-block:: python

   >>> table = tree.TableElement(votable)
   >>> r1.tables.append(t1)
   >>> votable.resources.append(r1)
   >>> for resource in votable.resources:
   ...     print(resource.mivot_block.content)
   <VODML xmlns="http://www.ivoa.net/xml/mivot" >
      <REPORT status="OK">Unit test mapping block</REPORT>
      <GLOBALS>  </GLOBALS>
   </VODML>
