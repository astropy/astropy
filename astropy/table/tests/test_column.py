from .. import Column
import pytest


class TestColumn():

    def test_1(self):
        Column()

    def test_2(self):
        c = Column()
        with pytest.raises(Exception):
            c.dtype = int  # can't set dtype directly

    def test_3(self):
        Column(dtype=int)

    def test_4(self):
        c1 = Column()
        c2 = Column()
        assert c1 == c2

    def test_5(self):
        c1 = Column(name='a', dtype=int, units='mJy')
        c2 = Column(name='a', dtype=int, units='mJy')
        assert c1 == c2

    def test_6(self):
        c1 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        assert c1 == c2

    def test_7(self):
        c1 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='b', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        assert c1 != c2

    def test_8(self):
        c1 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='a', dtype=float, units='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        assert c1 != c2

    def test_9(self):
        c1 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='a', dtype=int, units='ergs/cm^2/s/Hz', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        assert c1 != c2

    def test_10(self):
        c1 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='a', dtype=int, units='mJy', format='%g',
                    description='test column', meta={'c': 8, 'd': 12})
        assert c1 != c2

    def test_11(self):
        c1 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='another test column', meta={'c': 8, 'd': 12})
        assert c1 != c2

    def test_12(self):
        c1 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'e': 8, 'd': 12})
        assert c1 != c2

    def test_13(self):
        c1 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'c': 8, 'd': 12})
        c2 = Column(name='a', dtype=int, units='mJy', format='%i',
                    description='test column', meta={'c': 9, 'd': 12})
        assert c1 != c2
