# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

import numpy as np

from ...tests.helper import pytest
from ... import table
from ...table import Table
from ...extern.six import PY3
from ...utils import console

BIG_WIDE_ARR = np.arange(2000, dtype=np.float64).reshape(100, 20)
SMALL_ARR = np.arange(18, dtype=np.int64).reshape(6, 3)


@pytest.mark.usefixtures('table_type')
class TestMultiD():

    def test_multidim(self, table_type):
        """Test printing with multidimensional column"""
        arr = [np.array([[1, 2],
                         [10, 20]], dtype=np.int64),
               np.array([[3, 4],
                         [30, 40]], dtype=np.int64),
               np.array([[5, 6],
                         [50, 60]], dtype=np.int64)]
        t = table_type(arr)
        lines = t.pformat()
        assert lines == ['col0 [2] col1 [2] col2 [2]',
                         '-------- -------- --------',
                         '  1 .. 2   3 .. 4   5 .. 6',
                         '10 .. 20 30 .. 40 50 .. 60']

        lines = t.pformat(html=True)
        assert lines == ['<table id="table{id}">'.format(id=id(t)),
                         '<thead><tr><th>col0 [2]</th><th>col1 [2]</th><th>col2 [2]</th></tr></thead>',
                         '<tr><td>1 .. 2</td><td>3 .. 4</td><td>5 .. 6</td></tr>',
                         '<tr><td>10 .. 20</td><td>30 .. 40</td><td>50 .. 60</td></tr>',
                         '</table>']
        assert t._repr_html_().splitlines() == [
            '&lt;{0} masked={1} length=2&gt;'.format(table_type.__name__, t.masked),
            '<table id="table{tid}">'.format(tid=id(t)),
            '<thead><tr><th>col0 [2]</th><th>col1 [2]</th><th>col2 [2]</th></tr></thead>',
            '<thead><tr><th>int64</th><th>int64</th><th>int64</th></tr></thead>',
            '<tr><td>1 .. 2</td><td>3 .. 4</td><td>5 .. 6</td></tr>',
            '<tr><td>10 .. 20</td><td>30 .. 40</td><td>50 .. 60</td></tr>',
            '</table>']

        t = table_type([arr])
        lines = t.pformat()
        assert lines == ['col0 [2,2]',
                         '----------',
                         '   1 .. 20',
                         '   3 .. 40',
                         '   5 .. 60']

    def test_fake_multidim(self, table_type):
        """Test printing with 'fake' multidimensional column"""
        arr = [np.array([[(1,)],
                         [(10,)]], dtype=np.int64),
               np.array([[(3,)],
                         [(30,)]], dtype=np.int64),
               np.array([[(5,)],
                         [(50,)]], dtype=np.int64)]
        t = table_type(arr)
        lines = t.pformat()
        assert lines == ['col0 [1,1] col1 [1,1] col2 [1,1]',
                         '---------- ---------- ----------',
                         '         1          3          5',
                         '        10         30         50']

        lines = t.pformat(html=True)
        assert lines == ['<table id="table{id}">'.format(id=id(t)),
                         '<thead><tr><th>col0 [1,1]</th><th>col1 [1,1]</th><th>col2 [1,1]</th></tr></thead>',
                         '<tr><td>1</td><td>3</td><td>5</td></tr>',
                         '<tr><td>10</td><td>30</td><td>50</td></tr>',
                         '</table>']
        assert t._repr_html_().splitlines() == [
            '&lt;{0} masked={1} length=2&gt;'.format(table_type.__name__, t.masked),
            '<table id="table{id}">'.format(id=id(t)),
            '<thead><tr><th>col0 [1,1]</th><th>col1 [1,1]</th><th>col2 [1,1]</th></tr></thead>',
            '<thead><tr><th>int64</th><th>int64</th><th>int64</th></tr></thead>',
            '<tr><td>1</td><td>3</td><td>5</td></tr>', u'<tr><td>10</td><td>30</td><td>50</td></tr>',
            '</table>']

        t = table_type([arr])
        lines = t.pformat()
        assert lines == ['col0 [2,1,1]',
                         '------------',
                         '     1 .. 10',
                         '     3 .. 30',
                         '     5 .. 50']


def test_html_escaping():
    t = table.Table([(str('<script>alert("gotcha");</script>'), 2, 3)])
    assert t._repr_html_().splitlines() == [
        '&lt;Table length=3&gt;',
        '<table id="table{id}">'.format(id=id(t)),
        '<thead><tr><th>col0</th></tr></thead>',
        '<thead><tr><th>str33</th></tr></thead>',
        '<tr><td>&lt;script&gt;alert(&quot;gotcha&quot;);&lt;/script&gt;</td></tr>',
        '<tr><td>2</td></tr>',
        '<tr><td>3</td></tr>',
        '</table>']

@pytest.mark.usefixtures('table_type')
class TestPprint():

    def _setup(self, table_type):
        self.tb = table_type(BIG_WIDE_ARR)
        self.tb['col0'].format = 'e'
        self.tb['col1'].format = '.6f'

        self.tb['col0'].unit = 'km**2'
        self.tb['col19'].unit = 'kg s m**-2'
        self.ts = table_type(SMALL_ARR)

    def test_empty_table(self, table_type):
        t = table_type()
        lines = t.pformat()
        assert lines == ['<No columns>']
        c = repr(t)
        assert c.splitlines() == ['<{0} masked={1} length=0>'.format(table_type.__name__, t.masked),
                                  '<No columns>']

    def test_format0(self, table_type):
        """Try getting screen size but fail to defaults because testing doesn't
        have access to screen (fcntl.ioctl fails).
        """
        self._setup(table_type)
        arr = np.arange(4000, dtype=np.float64).reshape(100, 40)
        lines = table_type(arr).pformat()
        nlines, width = console.terminal_size()
        assert len(lines) == nlines
        for line in lines[:-1]:  # skip last "Length = .. rows" line
            assert (len(line) > width - 10 and
                    len(line) <= width)

    def test_format1(self, table_type):
        """Basic test of formatting, unit header row included"""
        self._setup(table_type)
        lines = self.tb.pformat(max_lines=8, max_width=40)
        assert lines == ['    col0         col1    ...   col19  ',
                         '    km2                  ... kg s / m2',
                         '------------ ----------- ... ---------',
                         '0.000000e+00    1.000000 ...      19.0',
                         '         ...         ... ...       ...',
                         '1.960000e+03 1961.000000 ...    1979.0',
                         '1.980000e+03 1981.000000 ...    1999.0',
                         'Length = 100 rows']

    def test_format2(self, table_type):
        """Basic test of formatting, unit header row excluded"""
        self._setup(table_type)
        lines = self.tb.pformat(max_lines=8, max_width=40, show_unit=False)
        assert lines == ['    col0         col1    ... col19 ',
                         '------------ ----------- ... ------',
                         '0.000000e+00    1.000000 ...   19.0',
                         '2.000000e+01   21.000000 ...   39.0',
                         '         ...         ... ...    ...',
                         '1.960000e+03 1961.000000 ... 1979.0',
                         '1.980000e+03 1981.000000 ... 1999.0',
                         'Length = 100 rows']

    def test_format3(self, table_type):
        """Include the unit header row"""
        self._setup(table_type)
        lines = self.tb.pformat(max_lines=8, max_width=40, show_unit=True)

        assert lines == ['    col0         col1    ...   col19  ',
                         '    km2                  ... kg s / m2',
                         '------------ ----------- ... ---------',
                         '0.000000e+00    1.000000 ...      19.0',
                         '         ...         ... ...       ...',
                         '1.960000e+03 1961.000000 ...    1979.0',
                         '1.980000e+03 1981.000000 ...    1999.0',
                         'Length = 100 rows']

    def test_format4(self, table_type):
        """Do not include the name header row"""
        self._setup(table_type)
        lines = self.tb.pformat(max_lines=8, max_width=40, show_name=False)
        assert lines == ['    km2                  ... kg s / m2',
                         '------------ ----------- ... ---------',
                         '0.000000e+00    1.000000 ...      19.0',
                         '2.000000e+01   21.000000 ...      39.0',
                         '         ...         ... ...       ...',
                         '1.960000e+03 1961.000000 ...    1979.0',
                         '1.980000e+03 1981.000000 ...    1999.0',
                         'Length = 100 rows']

    def test_noclip(self, table_type):
        """Basic table print"""
        self._setup(table_type)
        lines = self.ts.pformat(max_lines=-1, max_width=-1)
        assert lines == ['col0 col1 col2',
                         '---- ---- ----',
                         '   0    1    2',
                         '   3    4    5',
                         '   6    7    8',
                         '   9   10   11',
                         '  12   13   14',
                         '  15   16   17']

    def test_clip1(self, table_type):
        """max lines below hard limit of 8
        """
        self._setup(table_type)
        lines = self.ts.pformat(max_lines=3, max_width=-1)
        assert lines == ['col0 col1 col2',
                         '---- ---- ----',
                         '   0    1    2',
                         '   3    4    5',
                         '   6    7    8',
                         '   9   10   11',
                         '  12   13   14',
                         '  15   16   17']

    def test_clip2(self, table_type):
        """max lines below hard limit of 8 and output longer than 8
        """
        self._setup(table_type)
        lines = self.ts.pformat(max_lines=3, max_width=-1, show_unit=True, show_dtype=True)
        assert lines == [' col0  col1  col2',
                         '                 ',
                         'int64 int64 int64',
                         '----- ----- -----',
                         '    0     1     2',
                         '  ...   ...   ...',
                         '   15    16    17',
                         'Length = 6 rows']


    def test_clip3(self, table_type):
        """Max lines below hard limit of 8 and max width below hard limit
        of 10
        """
        self._setup(table_type)
        lines = self.ts.pformat(max_lines=3, max_width=1, show_unit=True)
        assert lines == ['col0 ...',
                         '     ...',
                         '---- ...',
                         '   0 ...',
                         ' ... ...',
                         '  12 ...',
                         '  15 ...',
                         'Length = 6 rows']

    def test_clip4(self, table_type):
        """Test a range of max_lines"""
        self._setup(table_type)
        for max_lines in (0, 1, 4, 5, 6, 7, 8, 100, 101, 102, 103, 104, 130):
            lines = self.tb.pformat(max_lines=max_lines, show_unit=False)
            assert len(lines) == max(8, min(102, max_lines))



@pytest.mark.usefixtures('table_type')
class TestFormat():

    def test_column_format(self, table_type):
        t = table_type([[1, 2], [3, 4]], names=('a', 'b'))
        # default (format=None)
        assert str(t['a']) == ' a \n---\n  1\n  2'

        # just a plain format string
        t['a'].format = '5.2f'
        assert str(t['a']) == '  a  \n-----\n 1.00\n 2.00'

        #  Old-style that is almost new-style
        t['a'].format = '{ %4.2f }'
        assert str(t['a']) == '   a    \n--------\n{ 1.00 }\n{ 2.00 }'

        #  New-style that is almost old-style
        t['a'].format = '%{0:}'
        assert str(t['a']) == ' a \n---\n %1\n %2'

        #  New-style with extra spaces
        t['a'].format = ' {0:05d} '
        assert str(t['a']) == '   a   \n-------\n 00001 \n 00002 '

        #  New-style has precedence
        t['a'].format = '%4.2f {0:}'
        assert str(t['a']) == '   a   \n-------\n%4.2f 1\n%4.2f 2'

        #  Invalid format spec
        t['a'].format = 'fail'
        with pytest.raises(ValueError):
            str(t['a'])

    def test_column_format_with_threshold(self, table_type):
        from ... import conf
        with conf.set_temp('max_lines', 8):
            t = table_type([np.arange(20)], names=['a'])
            t['a'].format = '%{0:}'
            assert str(t['a']).splitlines() == [' a ',
                                                '---',
                                                ' %0',
                                                ' %1',
                                                '...',
                                                '%18',
                                                '%19',
                                                'Length = 20 rows']
            t['a'].format = '{ %4.2f }'
            assert str(t['a']).splitlines() == ['    a    ',
                                                '---------',
                                                ' { 0.00 }',
                                                ' { 1.00 }',
                                                '      ...',
                                                '{ 18.00 }',
                                                '{ 19.00 }',
                                                'Length = 20 rows']

    def test_column_format_func(self, table_type):
        # run most of functions twice
        # 1) astropy.table.pprint._format_funcs gets populated
        # 2) astropy.table.pprint._format_funcs gets used

        t = table_type([[1., 2.], [3, 4]], names=('a', 'b'))

        # mathematical function
        t['a'].format = lambda x: str(x * 3.)
        assert str(t['a']) == ' a \n---\n3.0\n6.0'
        assert str(t['a']) == ' a \n---\n3.0\n6.0'

    def test_column_format_callable(self, table_type):
        # run most of functions twice
        # 1) astropy.table.pprint._format_funcs gets populated
        # 2) astropy.table.pprint._format_funcs gets used

        t = table_type([[1., 2.], [3, 4]], names=('a', 'b'))

        # mathematical function
        class format(object):
            def __call__(self, x):
                return str(x * 3.)
        t['a'].format = format()
        assert str(t['a']) == ' a \n---\n3.0\n6.0'
        assert str(t['a']) == ' a \n---\n3.0\n6.0'

    def test_column_format_func_wrong_number_args(self, table_type):
        t = table_type([[1., 2.], [3, 4]], names=('a', 'b'))

        # function that expects wrong number of arguments
        def func(a, b):
            pass

        t['a'].format = func
        with pytest.raises(ValueError):
            str(t['a'])

    def test_column_format_func_multiD(self, table_type):
        arr = [np.array([[1, 2],
                         [10, 20]])]
        t = table_type(arr, names=['a'])

        # mathematical function
        t['a'].format = lambda x: str(x * 3.)
        outstr = '   a [2]    \n------------\n  3.0 .. 6.0\n30.0 .. 60.0'
        assert str(t['a']) == outstr
        assert str(t['a']) == outstr

    def test_column_format_func_not_str(self, table_type):
        t = table_type([[1., 2.], [3, 4]], names=('a', 'b'))

        # mathematical function
        t['a'].format = lambda x: x * 3
        with pytest.raises(ValueError):
            str(t['a'])

    def test_column_alignment(self, table_type):
        t = table_type([[1], [2], [3], [4]],
                       names=('long title a', 'long title b',
                              'long title c', 'long title d'))
        t['long title a'].format = '<'
        t['long title b'].format = '^'
        t['long title c'].format = '>'
        t['long title d'].format = '0='
        assert str(t['long title a']) == 'long title a\n------------\n1           '
        assert str(t['long title b']) == 'long title b\n------------\n     2      '
        assert str(t['long title c']) == 'long title c\n------------\n           3'
        assert str(t['long title d']) == 'long title d\n------------\n000000000004'


class TestFormatWithMaskedElements():

    def test_column_format(self):
        t = Table([[1, 2, 3], [3, 4, 5]], names=('a', 'b'), masked=True)
        t['a'].mask = [True, False, True]
        # default (format=None)
        assert str(t['a']) == ' a \n---\n --\n  2\n --'

        # just a plain format string
        t['a'].format = '5.2f'
        assert str(t['a']) == '  a  \n-----\n   --\n 2.00\n   --'

        #  Old-style that is almost new-style
        t['a'].format = '{ %4.2f }'
        assert str(t['a']) == '   a    \n--------\n      --\n{ 2.00 }\n      --'

        #  New-style that is almost old-style
        t['a'].format = '%{0:}'
        assert str(t['a']) == ' a \n---\n --\n %2\n --'

        #  New-style with extra spaces
        t['a'].format = ' {0:05d} '
        assert str(t['a']) == '   a   \n-------\n     --\n 00002 \n     --'

        #  New-style has precedence
        t['a'].format = '%4.2f {0:}'
        assert str(t['a']) == '   a   \n-------\n     --\n%4.2f 2\n     --'

    def test_column_format_with_threshold(self, table_type):
        from ... import conf
        with conf.set_temp('max_lines', 8):
            t = table_type([np.arange(20)], names=['a'])
            t['a'].format = '%{0:}'
            t['a'].mask[0] = True
            t['a'].mask[-1] = True
            assert str(t['a']).splitlines() == [' a ',
                                                '---',
                                                ' --',
                                                ' %1',
                                                '...',
                                                '%18',
                                                ' --',
                                                'Length = 20 rows']
            t['a'].format = '{ %4.2f }'
            assert str(t['a']).splitlines() == ['    a    ',
                                                '---------',
                                                '       --',
                                                ' { 1.00 }',
                                                '      ...',
                                                '{ 18.00 }',
                                                '       --',
                                                'Length = 20 rows']

    def test_column_format_func(self):
        # run most of functions twice
        # 1) astropy.table.pprint._format_funcs gets populated
        # 2) astropy.table.pprint._format_funcs gets used

        t = Table([[1., 2., 3.], [3, 4, 5]], names=('a', 'b'), masked=True)
        t['a'].mask = [True, False, True]
        # mathematical function
        t['a'].format = lambda x: str(x * 3.)
        assert str(t['a']) == ' a \n---\n --\n6.0\n --'
        assert str(t['a']) == ' a \n---\n --\n6.0\n --'

    def test_column_format_func_with_special_masked(self):
        # run most of functions twice
        # 1) astropy.table.pprint._format_funcs gets populated
        # 2) astropy.table.pprint._format_funcs gets used

        t = Table([[1., 2., 3.], [3, 4, 5]], names=('a', 'b'), masked=True)
        t['a'].mask = [True, False, True]
        # mathematical function
        def format_func(x):
            if x is np.ma.masked:
                return '!!'
            else:
                return str(x * 3.)
        t['a'].format = format_func
        assert str(t['a']) == ' a \n---\n !!\n6.0\n !!'
        assert str(t['a']) == ' a \n---\n !!\n6.0\n !!'

    def test_column_format_callable(self):
        # run most of functions twice
        # 1) astropy.table.pprint._format_funcs gets populated
        # 2) astropy.table.pprint._format_funcs gets used

        t = Table([[1., 2., 3.], [3, 4, 5]], names=('a', 'b'), masked=True)
        t['a'].mask = [True, False, True]

        # mathematical function
        class format(object):
            def __call__(self, x):
                return str(x * 3.)
        t['a'].format = format()
        assert str(t['a']) == ' a \n---\n --\n6.0\n --'
        assert str(t['a']) == ' a \n---\n --\n6.0\n --'

    def test_column_format_func_wrong_number_args(self):
        t = Table([[1., 2.], [3, 4]], names=('a', 'b'), masked=True)
        t['a'].mask = [True, False]

        # function that expects wrong number of arguments
        def func(a, b):
            pass

        t['a'].format = func
        with pytest.raises(ValueError):
            str(t['a'])

        # but if all are masked, it never gets called
        t['a'].mask = [True, True]
        assert str(t['a']) == ' a \n---\n --\n --'

    def test_column_format_func_multiD(self):
        arr = [np.array([[1, 2],
                         [10, 20]])]
        t = Table(arr, names=['a'], masked=True)
        t['a'].mask[0,1] = True
        t['a'].mask[1,1] = True
        # mathematical function
        t['a'].format = lambda x: str(x * 3.)
        outstr = '  a [2]   \n----------\n 3.0 .. --\n30.0 .. --'
        assert str(t['a']) == outstr
        assert str(t['a']) == outstr


def test_pprint_npfloat32():
    """
    Test for #148, that np.float32 cannot by itself be formatted as float,
    but has to be converted to a python float.
    """
    dat = np.array([1., 2.], dtype=np.float32)
    t = Table([dat], names=['a'])
    t['a'].format = '5.2f'
    assert str(t['a']) == '  a  \n-----\n 1.00\n 2.00'


def test_pprint_py3_bytes():
    """
    Test for #1346.  Make sure a bytestring (dtype=S<N>) in Python 3 is printed
    correctly (without the "b" prefix like b'string').
    """
    val = bytes('val', encoding='utf-8') if PY3 else 'val'
    dat = np.array([(val,)], dtype=[(str('col'), 'S3')])
    t = table.Table(dat)
    assert t['col'].pformat() == ['col', '---', 'val']


def test_pprint_nameless_col():
    """Regression test for #2213, making sure a nameless column can be printed
    using None as the name.
    """
    col = table.Column([1.,2.])
    assert str(col).startswith('None')
