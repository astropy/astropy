# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.table import QTable

__all__ = ['BaseTimeSeries']


class BaseTimeSeries(QTable):

    _required_columns = None

    def add_columns(self, cols, indexes=None, names=None, **kwargs):

        if names is None:
            new_names = [col.info.name for col in cols]
        else:
            new_names = names

        existing_and_new = set(new_names) | set(self.colnames)

        if self._required_columns is not None:

            if len(set(self._required_columns) & existing_and_new) == len(self._required_columns):
                # Existing and new columns include all required columns, so we are
                # good to go
                self._required_columns = None
            elif len(existing_and_new - set(self._required_columns)) == 0:
                # All existing and new columns are a strict subset of required
                # columns - it doesn't mean that all required columns are present,
                # but no non-required columns are. So we are good to go.
                pass
            else:
                # There are columns from the required ones missing but some non-
                # required columns are being added
                for colname in self._required_columns:
                    if colname not in existing_and_new:
                        raise ValueError("{0} requires a column called '{1}' to be set "
                                         "before data can be added"
                                         .format(self.__class__.__name__, colname))

        return super().add_columns(cols, indexes=indexes, names=names, **kwargs)
