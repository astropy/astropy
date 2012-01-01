import numpy as np
from .. import Column


class TestColumn():

    def test_1(self):
        Column('a')

    def test_subclass(self):
        c = Column('a')
        assert isinstance(c, np.ndarray)
        c2 = c * 2
        assert isinstance(c2, Column)
        assert isinstance(c2, np.ndarray)

    def test_numpy_ops(self):
        """Show that basic numpy operations with Column behave sensibly"""

        arr = np.array([1, 2, 3])
        c = Column('a', arr)
        eq = c == arr
        assert np.all(eq)
        assert len(eq) == 3
        assert type(eq) == Column
        assert eq.dtype.str == '|b1'
        eq = arr == c
        assert np.all(eq)

        lt = c - 1 < arr
        assert np.all(lt)

    def test_view(self):
        c = np.array([1, 2, 3]).view(Column)
        assert repr(c) == 'array([1, 2, 3])'

    def test_format(self):
        """Show that the formatted output from str() works"""
        c1 = Column(name='a', data=np.arange(2000), dtype=float,
                    format='%6.2f')
        assert str(c1).startswith('  0.00,   1.00,   2.00')
        assert str(c1).endswith(', 1999.00')
        np.set_printoptions(threshold=5)
        c1.format = '%6.1f'
        assert str(c1) == '   0.0,    1.0,    2.0, ..., 1998.0, 1999.0'


class TestAttrEqual():
    """Bunch of tests originally from ATpy that test the attrs_equal method."""

    def test_5(self):
        c1 = Column(name='a', dtype=int, units='mJy')
        c2 = Column(name='a', dtype=int, units='mJy')
        assert c1.attrs_equal(c2)

    def test_6(self):
        c1 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        assert c1.attrs_equal(c2)

    def test_7(self):
        c1 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='b', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        assert not c1.attrs_equal(c2)

    def test_8(self):
        c1 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='a', dtype=float, units='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        assert not c1.attrs_equal(c2)

    def test_9(self):
        c1 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='a', dtype=int, units='ergs/cm^2/s/Hz', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        assert not c1.attrs_equal(c2)

    def test_10(self):
        c1 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='a', dtype=int, units='mJy', format='%g',
                    description='test column', meta={'c': 8, 'd': 12})
        assert not c1.attrs_equal(c2)

    def test_11(self):
        c1 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='another test column', meta={'c': 8, 'd': 12})
        assert not c1.attrs_equal(c2)

    def test_12(self):
        c1 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'e': 8, 'd': 12})
        assert not c1.attrs_equal(c2)

    def test_13(self):
        c1 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'c': 9, 'd': 12})
        assert not c1.attrs_equal(c2)

