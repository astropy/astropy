# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
The **asdf** subpackage contains code that is used to serialize astropy types
so that they can be represented and stored using the Advanced Scientific Data
Format (ASDF). This subpackage defines classes, referred to as **tags**, that
implement the logic for serialization and deserialization.

ASDF makes use of abstract data type definitions called **schemas**. The tag
classes provided here are specific implementations of particular schemas. Some
of the tags in Astropy (e.g., those related to transforms) implement schemas
that are defined by the ASDF Standard. In other cases, both the tags and
schemas are defined within Astropy (e.g., those related to many of the
coordinate frames).

Astropy currently has no ability to read or write ASDF files itself. In order
to process ASDF files it is necessary to make use of the standalone **asdf**
package. Users should never need to refer to tag implementations directly.
Their presence should be entirely transparent when processing ASDF files.

If both **asdf** and **astropy** are installed, no further configuration is
required in order to process ASDF files. The **asdf** package has been designed
to automatically detect the presence of the tags defined by **astropy**.

Documentation on the ASDF Standard can be found `here
<https://asdf-standard.readthedocs.io>`__. Documentation on the ASDF Python
module can be found `here <https://asdf.readthedocs.io>`__.
"""
