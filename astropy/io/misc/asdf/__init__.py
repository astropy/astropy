# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
The **asdf** subpackage contains code that is used to serialize astropy types
so that they can be represented and stored using the Advanced Scientific Data
Format (ASDF). This subpackage defines classes, referred to as **tags**, that
implement the logic for serialization and deserialization.

ASDF makes use of abstract data type definitons called **schemas**. The tags
provided here are simply specific implementations of particular schemas.
Currently astropy only implements tags for a subset of schemas that are defined
externally by the ASDF Standard. However, it is likely that astropy will
eventually define schemas of its own.

Astropy currently has no ability to read or write ASDF files itself. In order
to process ASDF files it is necessary to make use of the standalone **asdf**
package. Users should never need to refer to tag implementations directly.
Their presence should be entirely transparent when processing ASDF files.

If both **asdf** and **astropy** are installed, no futher configuration is
required in order to process ASDF files. The **asdf** package has been designed
to automatically detect the presence of the tags defined by **astropy**.

Documentation on the ASDF Standard can be found `here
<https://asdf-standard.readthedocs.io>`__. Documentation on the ASDF Python
module can be found `here <https://asdf.readthedocs.io>`__.
"""
