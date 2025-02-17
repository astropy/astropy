.. _table_io_jsviewer:

JSViewer
--------

Provides an interactive HTML export of a Table, like the
:class:`~astropy.io.ascii.HTML` writer but using the DataTables_ library, which
allow to visualize interactively an HTML table (with columns sorting, search,
and pagination).

Example
^^^^^^^

..
  EXAMPLE START
  JSViewer to Provide an Interactive HTML Export of a Table

To write a table ``t`` to a new file::

    >>> t.write('new_table.html', format='jsviewer')

Several additional parameters can be used:

- *table_id*: the HTML ID of the ``<table>`` tag, defaults to ``'table{id}'``
  where ``id`` is the ID of the Table object.
- *max_lines*: maximum number of lines.
- *table_class*: HTML classes added to the ``<table>`` tag, can be useful to
  customize the style of the table.
- *jskwargs*: additional arguments passed to :class:`~astropy.table.JSViewer`.
- *css*: CSS style, default to ``astropy.table.jsviewer.DEFAULT_CSS``.
- *htmldict*: additional arguments passed to :class:`~astropy.io.ascii.HTML`.

.. _Datatables: https://www.datatables.net/

..
  EXAMPLE END
