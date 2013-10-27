import numpy as np

from ..metadata import MetaData
from ..compat.odict import OrderedDict
from ...tests.helper import pytest, raises
from ...io import fits


class OrderedDictSubclass(OrderedDict):
    pass


class Data(object):
    meta = MetaData()

    def __init__(self, meta=None):
        self.meta = meta


def test_none():
    d = Data()
    assert isinstance(d.meta, OrderedDict)
    assert len(d.meta) == 0


@pytest.mark.parametrize(('meta'), ([dict([('a', 1)]),
                                     OrderedDict([('a', 1)]),
                                     OrderedDictSubclass([('a', 1)])]))
def test_mapping_init(meta):
    d = Data(meta=meta)
    assert type(d.meta) == type(meta)
    assert d.meta['a'] == 1


@pytest.mark.parametrize(('meta'), (["ceci n'est pas un meta", 1.2, [1, 2, 3]]))
def test_non_mapping_init(meta):
    with pytest.raises(TypeError):
        d = Data(meta=meta)


@pytest.mark.parametrize(('meta'), ([dict([('a', 1)]),
                                     OrderedDict([('a', 1)]),
                                     OrderedDictSubclass([('a', 1)])]))
def test_mapping_set(meta):
    d = Data(meta=meta)
    assert type(d.meta) == type(meta)
    assert d.meta['a'] == 1


@pytest.mark.parametrize(('meta'), (["ceci n'est pas un meta", 1.2, [1, 2, 3]]))
def test_non_mapping_set(meta):
    with pytest.raises(TypeError):
        d = Data(meta=meta)


def test_meta_fits_header():

    header = fits.header.Header()
    header.set('observer', 'Edwin Hubble')
    header.set('exptime', '3600')

    d = Data(meta=header)

    assert d.meta['OBSERVER'] == 'Edwin Hubble'
