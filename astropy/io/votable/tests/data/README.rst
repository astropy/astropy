VOTable test data
=================

Listing where some of the data files used in the tests come from.

``unicode-example.vot``, ``unicode-example-binary.vot``, ``unicode-example-binary2.vot``
    VOTable 1.6 examples by Mark Taylor, copied from the top level
    https://github.com/ivoa-std/VOTable at commit 2836c14. They hold the
    same table in TABLEDATA, BINARY and BINARY2. Used by
    ``test_unicode_examples.py``.

``utf8_char.xml``
    A small VOTable 1.6 file with UTF-8 ``char`` columns, one of them with a
    ``width``. Written for ``test_schema_versions.py``.
