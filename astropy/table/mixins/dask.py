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
        return as_dask_column(self, info=self.info)

    def __getitem__(self, item):
        result = super().__getitem__(item)
        if isinstance(item, int):
            return result
        else:
            return as_dask_column(result, info=self.info)

    def insert(self, obj, values, axis=0):
        return as_dask_column(da.insert(self, obj, values, axis=axis),
                              info=self.info)


def as_dask_column(array, info=None):
    result = DaskColumn(array.dask, array.name, array.chunks, meta=array)
    if info is not None:
        result.info = info
    return result
