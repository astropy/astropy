import abc

from ..metadata import MetaData
from ..compat.odict import OrderedDict
from ...tests.helper import pytest
from ...io import fits


class OrderedDictSubclass(OrderedDict):
    pass


class MetaBaseTest(object):

    __metaclass__ = abc.ABCMeta

    def test_none(self):
        d = self.test_class(*self.args)
        assert isinstance(d.meta, OrderedDict)
        assert len(d.meta) == 0


    @pytest.mark.parametrize(('meta'), ([dict([('a', 1)]),
                                         OrderedDict([('a', 1)]),
                                         OrderedDictSubclass([('a', 1)])]))
    def test_mapping_init(self, meta):
        d = self.test_class(*self.args, meta=meta)
        assert type(d.meta) == type(meta)
        assert d.meta['a'] == 1


    @pytest.mark.parametrize(('meta'), (["ceci n'est pas un meta", 1.2, [1, 2, 3]]))
    def test_non_mapping_init(self, meta):
        with pytest.raises(TypeError):
            self.test_class(*self.args, meta=meta)


    @pytest.mark.parametrize(('meta'), ([dict([('a', 1)]),
                                         OrderedDict([('a', 1)]),
                                         OrderedDictSubclass([('a', 1)])]))
    def test_mapping_set(self, meta):
        d = self.test_class(*self.args, meta=meta)
        assert type(d.meta) == type(meta)
        assert d.meta['a'] == 1


    @pytest.mark.parametrize(('meta'), (["ceci n'est pas un meta", 1.2, [1, 2, 3]]))
    def test_non_mapping_set(self, meta):
        with pytest.raises(TypeError):
            d = self.test_class(*self.args, meta=meta)


    def test_meta_fits_header(self):

        header = fits.header.Header()
        header.set('observer', 'Edwin Hubble')
        header.set('exptime', '3600')

        d = self.test_class(*self.args, meta=header)

        assert d.meta['OBSERVER'] == 'Edwin Hubble'


class ExampleData(object):
    meta = MetaData()

    def __init__(self, meta=None):
        self.meta = meta


class TestMetaExampleData(MetaBaseTest):
    test_class = ExampleData
    args = ()
