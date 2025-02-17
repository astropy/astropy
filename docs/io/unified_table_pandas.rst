

.. _table_io_pandas:

Pandas
------

.. _pandas: https://pandas.pydata.org/pandas-docs/stable/index.html

``astropy`` `~astropy.table.Table` supports the ability to read or write tables
using some of the `I/O methods <https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html>`_
available within pandas_. This interface thus provides convenient wrappers to
the following functions / methods:

.. csv-table::
    :header: "Format name", "Data Description", "Reader", "Writer"
    :widths: 25, 25, 25, 25

    ``pandas.csv``,`CSV <https://en.wikipedia.org/wiki/Comma-separated_values>`__,`read_csv() <https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#io-read-csv-table>`_,`to_csv() <https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#io-store-in-csv>`_
    ``pandas.json``,`JSON <http://www.json.org/>`__,`read_json() <https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#io-json-reader>`_,`to_json() <https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#io-json-writer>`_
    ``pandas.html``,`HTML <https://en.wikipedia.org/wiki/HTML>`__,`read_html() <https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#io-read-html>`_,`to_html() <https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#io-html>`_
    ``pandas.fwf``,Fixed Width,`read_fwf() <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_fwf.html#pandas.read_fwf>`_,

**Notes**:

- There is no fixed-width writer in pandas_.
- Reading HTML requires `BeautifulSoup4 <https://pypi.org/project/beautifulsoup4/>`_ and
  `html5lib <https://pypi.org/project/html5lib/>`_ to be installed.

When reading or writing a table, any keyword arguments apart from the
``format`` and file name are passed through to pandas, for instance:

.. doctest-skip::

  >>> t.write('data.csv', format='pandas.csv', sep=' ', header=False)
  >>> t2 = Table.read('data.csv', format='pandas.csv', sep=' ', names=['a', 'b', 'c'])
