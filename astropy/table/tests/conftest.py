# Licensed under a 3-clause BSD style license - see LICENSE.rst

from ...tests.helper import pytest
from ... import table


@pytest.fixture(params=[table.Column, table.MaskedColumn])
def Column(request):
    # Fixture to run all the Column tests for both an unmasked (ndarray)
    # and masked (MaskedArray) column.
    return request.param


class MaskedTable(table.Table):
    def __init__(self, *args, **kwargs):
        kwargs['masked'] = True
        table.Table.__init__(self, *args, **kwargs)


# Fixture to run all the Column tests for both an unmasked (ndarray)
# and masked (MaskedArray) column.
@pytest.fixture(params=[False, True])
def table_types(request):
    class TableTypes:
        def __init__(self, request):
            self.Table = MaskedTable if request.param else table.Table
            self.Column = table.MaskedColumn if request.param else table.Column
    return TableTypes(request)


# Fixture to run all the Column tests for both an unmasked (ndarray)
# and masked (MaskedArray) column.
@pytest.fixture(params=[False, True])
def table_type(request):
    return MaskedTable if request.param else table.Table


# Fixture to run all the Column tests for both an unmasked (ndarray)
# and masked (MaskedArray) column.
@pytest.fixture(params=[False, True])
def table_data(request):
    class TableData:
        def __init__(self, request):
            self.Table = MaskedTable if request.param else table.Table
            self.Column = table.MaskedColumn if request.param else table.Column
            self.COLS = [
                self.Column(name='a', data=[1, 2, 3], description='da',
                            format='fa', meta={'ma': 1}, unit='ua'),
                self.Column(name='b', data=[4, 5, 6], description='db',
                            format='fb', meta={'mb': 1}, unit='ub'),
                self.Column(name='c', data=[7, 8, 9], description='dc',
                            format='fc', meta={'mc': 1}, unit='ub')]
            self.DATA = self.Table(self.COLS)
    return TableData(request)


class SubclassTable(table.Table):
    pass


@pytest.fixture(params=[True, False])
def tableclass(request):
    return table.Table if request.param else SubclassTable


@pytest.fixture(params=[0, 1, -1])
def protocol(request):
    """
    Fixture to run all the tests for protocols 0 and 1, and -1 (most advanced).
    """
    return request.param


# Fixture to run all tests for both an unmasked (ndarray) and masked
# (MaskedArray) column.
@pytest.fixture(params=[False, True])
def table_type(request):
    # return MaskedTable if request.param else table.Table
    try:
        request.param
        return MaskedTable
    except AttributeError:
        return table.Table
