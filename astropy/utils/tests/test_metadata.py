import abc

from collections import OrderedDict

from ..metadata import MetaData, MergeConflictError, merge, enable_merge_strategies
from ...utils import metadata
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


def test_metadata_merging_conflict_exception():
    """Regression test for issue #3294.

    Ensure that an exception is raised when a metadata conflict exists
    and ``metadata_conflicts='error'`` has been set.
    """
    data1 = ExampleData()
    data2 = ExampleData()
    data1.meta['somekey'] = {'x': 1, 'y': 1}
    data2.meta['somekey'] = {'x': 1, 'y': 999}
    with pytest.raises(MergeConflictError):
        merge(data1.meta, data2.meta, metadata_conflicts='error')


import numpy as np

def test_metadata_merging():
    # Recursive merge
    meta1 = {'k1': {'k1': [1, 2],
                    'k2': 2},
             'k2': 2,
             'k4': (1, 2)}
    meta2 = {'k1': {'k1': [3]},
             'k3': 3,
             'k4': (3,)}
    out = merge(meta1, meta2, metadata_conflicts='error')
    assert out == {'k1': {'k2': 2,
                          'k1': [1, 2, 3]},
                   'k2': 2,
                   'k3': 3,
                   'k4': (1, 2, 3)}

    # Merge two ndarrays
    meta1 = {'k1': np.array([1, 2])}
    meta2 = {'k1': np.array([3])}
    out = merge(meta1, meta2, metadata_conflicts='error')
    assert np.all(out['k1'] == np.array([1, 2, 3]))

    # Merge list and np.ndarray
    meta1 = {'k1': [1, 2]}
    meta2 = {'k1': np.array([3])}
    assert np.all(out['k1'] == np.array([1, 2, 3]))

    # Can't merge two scalar types
    meta1 = {'k1': 1}
    meta2 = {'k1': 2}
    with pytest.raises(MergeConflictError):
        merge(meta1, meta2, metadata_conflicts='error')

    # Conflicting shape
    meta1 = {'k1': np.array([1, 2])}
    meta2 = {'k1': np.array([[3]])}
    with pytest.raises(MergeConflictError):
        merge(meta1, meta2, metadata_conflicts='error')

    # Conflicting array type
    meta1 = {'k1': np.array([1, 2])}
    meta2 = {'k1': np.array(['3'])}
    with pytest.raises(MergeConflictError):
        merge(meta1, meta2, metadata_conflicts='error')

    # Conflicting array type with 'silent' merging
    meta1 = {'k1': np.array([1, 2])}
    meta2 = {'k1': np.array(['3'])}
    out = merge(meta1, meta2, metadata_conflicts='silent')
    assert np.all(out['k1'] == np.array(['3']))


def test_metadata_merging_new_strategy():
    original_merge_strategies = list(metadata.MERGE_STRATEGIES)

    class MergeNumbersAsList(metadata.MergeStrategy):
        """
        Scalar float or int values are joined in a list.
        """
        types = ((int, float), (int, float))

        @classmethod
        def merge(cls, left, right):
            return [left, right]

    class MergeConcatStrings(metadata.MergePlus):
        """
        Scalar string values are concatenated
        """
        types = (str, str)
        enabled = False

    # Normally can't merge two scalar types
    meta1 = {'k1': 1, 'k2': 'a'}
    meta2 = {'k1': 2, 'k2': 'b'}

    # Enable new merge strategy
    with enable_merge_strategies(MergeNumbersAsList, MergeConcatStrings):
        assert MergeNumbersAsList.enabled
        assert MergeConcatStrings.enabled
        out = merge(meta1, meta2, metadata_conflicts='error')
    assert out['k1'] == [1, 2]
    assert out['k2'] == 'ab'
    assert not MergeNumbersAsList.enabled
    assert not MergeConcatStrings.enabled

    # Confirm the default enabled=False behavior
    with pytest.raises(MergeConflictError):
        merge(meta1, meta2, metadata_conflicts='error')

    # Enable all MergeStrategy subclasses
    with enable_merge_strategies(metadata.MergeStrategy):
        assert MergeNumbersAsList.enabled
        assert MergeConcatStrings.enabled
        out = merge(meta1, meta2, metadata_conflicts='error')
    assert out['k1'] == [1, 2]
    assert out['k2'] == 'ab'
    assert not MergeNumbersAsList.enabled
    assert not MergeConcatStrings.enabled


    metadata.MERGE_STRATEGIES = original_merge_strategies
