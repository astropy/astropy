import dask.array as da

from astropy.utils.data_info import ParentDtypeInfo

__all__ = ['as_dask_column']


class DaskInfo(ParentDtypeInfo):
    @staticmethod
    def default_format(val):
        return f'{val.compute()}'


class DaskColumn(da.Array):

    info = DaskInfo()

    def copy(self):
        # Array hard-codes the resulting copied array as Array, so need to
        # overload this since Table tries to copy the array.
        return as_dask_column(self)

    def __getitem__(self, item):
        return as_dask_column(super().__getitem__(item))


def as_dask_column(array):
    return DaskColumn(array.dask, array.name, array.chunks, meta=array)
