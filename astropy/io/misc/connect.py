# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This file connects any readers/writers defined in io.misc to the
# astropy.table.Table class

from astropy.io.misc.pyarrow.csv import register_pyarrow_csv_table

from . import hdf5, parquet

hdf5.register_hdf5()
parquet.register_parquet()
parquet.register_parquet_votable()
register_pyarrow_csv_table()
