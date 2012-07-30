****************************************************
A description of the package (`astropy.packagename`)
****************************************************

When creating a new subpackage's docuementation, this file should be
copied to a file "index.rst" in a directory corresponding to the name of
the package. E.g., ``docs/packagename/index.rst``. And don't forget to
delete this paragraph.

Introduction
============

Include general content that might be useful for understanding the
package here, as well as general scientific or mathematical background
that might be necessary for the "big-picture" of this package.


Getting Started
===============

Short tutorial-like examples of how to do common-tasks - should be
fairly quick, with any more detailed examples in the next section.


Using `packagename`
===================

.. THIS SECTION SHOULD BE EITHER


This section is for the detailed documentation.  For simpler packages, this
should either by paragraphs or sub-divided into sub-sections like:

Sub-topic 1
-----------

Content if needed

A Complex example
-----------------

Content if needed

Sub-sub topic 1
^^^^^^^^^^^^^^^^

Content if needed (note the use of ^^^^ at this level).

Sub-sub-sub topic 1
"""""""""""""""""""

Content if needed (note the use of """"" at this level).
This is probably the deepest level that is practical.  However,
just in case, the next levels of detail should use the +, :, and ~
characters respectively.


.. OR IF MORE COMPLICATED,

For more complicated packages that require multiple documents, this
should just be a table of contents referencing those documents:

.. toctree::
    subdoc1
    subdoc2
    subdoc3


Either a toctree or sub-sections should be used, *not* both.

For example, if your toctree looks like the above example, this document
should be ``docs/packagename/index.rst``, and the other documents should
be ``docs/packagename/subdoc1.rst``, ``docs/packagename/subdoc2.rst``,
and ``docs/packagename/subdoc3.rst``.

In the "more complicated" case of using ``subdoc.rst`` files, each of those
should likewise use the section character header order of ``* = - ^ " + : ~``.


See Also (optional)
===================

Include here any references to related packages, articles, or texts.


Reference/API
=============

.. automodapi:: packagename


Acknowledgments and Licenses (optional)
=======================================

Any acknowledgements or licenses needed for this package - remove the
section if none are necessary.
