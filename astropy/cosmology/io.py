# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy
import inspect
import json
import os

import numpy as np

from astropy.io.registry import UnifiedReadWriteMethod
from astropy.io import registry
from astropy.table import QTable
from astropy.utils.shapes import IncompatibleShapeError

# from astropy.table.connect import TableRead, TableWrite

__all__ = ["CosmologyRead", "CosmologyWrite"]
__doctest_skip__ = ["CosmologyRead", "CosmologyWrite"]


_VAR_POSITIONAL = inspect.Parameter.VAR_POSITIONAL
_VAR_KEYWORD = inspect.Parameter.VAR_KEYWORD


class CosmologyRead(registry.UnifiedReadWrite):
    """Read and parse a data table and to a `~astropy.cosmology.Cosmology`.

    This function provides the Cosmology interface to the astropy unified I/O
    layer. This allows easily reading a file in many supported data formats
    using syntax such as::

        >>> from astropy.cosmology import Cosmology
        >>> cosmo1 = Cosmology.read('cosmo1.dat', format='ascii.ecsv')
        >>> cosmo2 = Cosmology.read('cosmo2.json', format='json')

    Get help on the available readers using the``help()`` method::

      >>> Cosmology.read.help()  # Get help reading and list supported formats
      >>> Cosmology.read.help('fits')  # Get detailed help on FITS reader
      >>> Cosmology.read.list_formats()  # Print list of available formats

    See also: https://docs.astropy.org/en/stable/io/unified.html

    Parameters
    ----------
    *args : (optional)
        Positional arguments passed through to data reader. If supplied the
        first argument is typically the input filename.
        If the first argument is a Mapping all other arguments are ignored.
    format : str (optional, keyword-only)
        File format specifier.
    units : list, dict (optional, keyword-only)
        List or dict of units to apply to columns
    descriptions : list, dict (optional, keyword-only)
        List or dict of descriptions to apply to columns
    **kwargs
        Keyword arguments passed through to data reader.

    Returns
    -------
    out : `~astropy.cosmology.Cosmology`
        `~astropy.cosmology.Cosmology` corresponding to file contents.
        If the table in the file held multiple rows, the row index must be
        specified.

    Methods
    -------
    from_mapping
    from_table

    Notes
    -----
    `~astropy.cosmology.Cosmology` is an abstract base class and cannot be
    instantiated. The cosmology class may be specified by name in a column
    titled ``cosmology`` (case insensitive). As a classmethod ``read`` may
    be called from any subclass of `~astropy.cosmology.Cosmology`. In that case
    the cosmology-specifying column is not required. If present, the class will
    be checked for compatibility (being a subclass) with the calling class.
    Like the base class, any Cosmology class's ``read`` method can instantiate
    any of its subclasses.

    The other parameters are passed to the init signature. Parameters not in the
    signature will be moved to the metadata, except if the key is already
    present, in which case a `TypeError` will be raised.

    metadata can be included as a column -- called "meta" -- containing a
    Mapping, like a `dict`. Alternatively, for file formats like ECSV that can
    hold metadata, any key of same name as the cosmology instance (held in the
    column "name") with a Mapping value will also be taken as metadata. If both
    exist, they will be merged, with the column taking higher priority.

    Many formats do not allow unit information, so make sure the values are in
    the correct units for initializing the cosmology. For this reason ECSV
    is the recommended file format for Cosmologies.

    Examples
    --------
    TODO

    Raises
    ------
    TypeError
        If a parameter is not in the init signature and already has a value in
        the metadata.
    """

    def __init__(self, instance, cls):
        super().__init__(instance, cls, "read")

    def __call__(self, *args, index=None, **kwargs):
        # shortcuts for determining format
        format = kwargs.get("format", None)
        fmt = str(format).lower()
        fend = str(args[0]).lower().endswith

        # need to determine the format: is it JSON or ECSV?
        # JSON?
        if (fmt == "json") or (format is None and fend(".json")):
            cosmos = self._read_json(*args, key=index, **kwargs)
        # ECSV?
        elif (fmt == "ascii.ecsv") or (format is None and fend(".ecsv")):
            table = QTable.read(*args, **kwargs)
            cosmos = self.from_table(table, index=index)
        else:
            raise ValueError(
                "`format` must be JSON or 'ascii.ecsv' or None. "
                "If format is None, `file` must have suffix '.json' or '.ecsv'."
            )

        return cosmos

    def _read_json(self, file, key=None, **kwargs):
        # read file, from path-like or file-like
        if isinstance(file, (str, bytes, os.PathLike)):  # pathlike
            with open(file, "r") as read_file:
                data = read_file.read()
        else:  # file-like : this also handles errors in dumping
            data = file.read()

        # parse file
        mapping = json.loads(data)

        return self.from_mapping(mapping, key=key)

    def from_mapping(self, mapping, key=None):
        """Load `~astropy.cosmology.Cosmology` from mapping object.

        Parameters
        ----------
        mapping : mapping
        key : hashable or None, optional
            Any valid key into `mapping`.
            If not None, this is used to select a cosmology if mapping is a
            set of cosmologies.

        Returns
        -------
        `~astropy.cosmology.Cosmology`

        """
        # get params subset from mpping
        params = mapping[key] if key is not None else mapping

        # determine class
        subclasses = {c.__qualname__: c for c in
                      self._cls.__subclasses__(deep=True)}
        cls_name = params.pop("cosmology", self._cls.__qualname__)
        try:
            cosmo_cls = subclasses[cls_name]
        except KeyError:
            raise KeyError(f"{cls_name} is not a subclass of {cls.__qualname__}.")

        # match keys to arguments
        inspect_params = cosmo_cls._init_signature.parameters.values()
        expected_params = set(
            [
                p.name
                for p in inspect_params
                if not ((p.kind == _VAR_POSITIONAL) or (p.kind == _VAR_KEYWORD))
            ]
        )
        HAVE_KWARGS = any([p.kind == _VAR_KEYWORD for p in inspect_params])

        # move non-parameters to metadata, if no **kwargs in signature
        meta = params.pop("meta", None) or {}
        for k, v in tuple(params.items()):
            if k not in expected_params and not HAVE_KWARGS:
                v = params.pop(k)  # make sure not in params
                # TODO! raise an error / warning if repeat k?
                meta.setdefault(k, v)

            # TODO! parse list str when Quantity supports it
            # s = repr([2, 3, 4] * u.kpc)[10:-1]
            # num, unit = s.split("]")
            # eval(num + "]") * u.Unit(unit)
        meta = meta if meta else None  # empty meta -> None

        # make cosmo instance
        ba = cosmo_cls._init_signature.bind_partial(**params, meta=meta)
        cosmo = cosmo_cls(*ba.args, **ba.kwargs)
        return cosmo

    def from_table(self, table, index=None):
        """Instantiate a `~astropy.cosmology.Cosmology` from a |QTable|.

        Parameters
        ----------
        table : `~astropy.QTable`
        index : int or None, optional
            The row from table.
            TODO! allow for index to be a string -- the name fromo the row.

        Returns
        -------
        `~astropy.cosmology.Cosmology`

        """
        # get 1-row of table
        row = table.iloc[index] if index is not None else table
        if len(row) != 1:  # check that it's one row
            raise IncompatibleShapeError("Table must be 1-D.")

        # special values
        # get name from 1) column 2) metadata 3) it's None
        name = row.columns.get("name", [row.meta.pop("name", None)])[0]
        # metadata is a mix of the general metadata and row-specific data
        meta = {k: v for k, v in row.meta.items() if not isinstance(k, int)}
        meta.update(row.meta.get(index, {}))
        # the cosmology class is in the table's metadata
        cosmo_cls = meta.pop("cosmology", None)

        # turn row into mapping (dict of the arguments)
        mapping = {k: v[0] for k, v in zip(row.colnames, row.values())}
        mapping["name"] = name
        mapping["meta"] = {**meta, **mapping.get("meta", {})}
        mapping.setdefault("cosmology", cosmo_cls)

        # build cosmology from map
        cosmo = self.from_mapping(mapping)
        return cosmo


class CosmologyWrite(registry.UnifiedReadWrite):
    """Write this Cosmology object out in the specified format.

    This function provides the Table interface to the astropy unified I/O
    layer.  This allows easily writing a file in many supported data formats
    using syntax such as::

      >>> from astropy.table import Table
      >>> dat = Table([[1, 2], [3, 4]], names=('a', 'b'))
      >>> dat.write('table.dat', format='ascii')

    Get help on the available writers for ``Table`` using the``help()`` method::

      >>> Table.write.help()  # Get help writing Table and list supported formats
      >>> Table.write.help('fits')  # Get detailed help on Table FITS writer
      >>> Table.write.list_formats()  # Print list of available formats

    The ``serialize_method`` argument is explained in the section on
    `Table serialization methods
    <https://docs.astropy.org/en/latest/io/unified.html#table-serialization-methods>`_.
    See also: https://docs.astropy.org/en/stable/io/unified.html

    Parameters
    ----------
    *args : tuple, optional
        Positional arguments passed through to data writer. If supplied the
        first argument is the output filename.
    format : str
        File format specifier.
    serialize_method : str, dict, optional
        Serialization method specifier for columns.
    **kwargs : dict, optional
        Keyword arguments passed through to data writer.

    """

    def __init__(self, instance, cls):
        super().__init__(instance, cls, "write")

    def __call__(self, *args, serialize_method=None, **kwargs):
        # shortcuts for determining format
        format = kwargs.get("format", None)
        fmt = str(format).lower()
        fend = str(args[0]).lower().endswith
        # need to determine the format: is it JSON or ECSV?
        # JSON?
        if (fmt == "json") or (format is None and fend(".json")):
            self._write_json(args[0], **kwargs)
        # ECSV?
        elif (fmt == "ascii.ecsv") or (format is None and fend(".ecsv")):
            table = self.to_table()
            table.write(*args, serialize_method=serialize_method, **kwargs)
        else:
            raise ValueError(
                "`format` must be JSON or 'ascii.ecsv' or None. "
                "If format is None, `file` must have suffix '.json' or '.ecsv'."
            )

    def _write_json(self, file, mode="w", overwrite=False, **kwargs):
        data = self.to_mapping()
        data["cosmology"] = data["cosmology"].__qualname__

        if isinstance(file, (str, bytes, os.PathLike)):  # pathlike
            if os.path.exists(file) and not overwrite:
                raise IOError("overwrite")
            with open(file, "w") as write_file:
                json.dump(data, write_file)
        else:  # this also handles errors in dumping
            json.dump(data, file)

    def to_mapping(self):
        """Return the Cosmology class, inputs, and metadata as a dict.

        Has key-values:
        - 'cosmology' : the cosmology's class
        - 'meta' : the contents of the cosmology's metadata attribute
        - keys : values from initialization

        """
        d = {}
        # start with the cosmology class
        d["cosmology"] = self._cls
        # get all the immutable inputs
        d.update({k: v for k, v in self._instance._init_arguments.items()
                  if k != "meta"})
        # add the mutable metadata
        d["meta"] = copy.deepcopy(self._instance.meta)

        return d

    def to_table(self):
        """Return the Cosmology a `~astropy.table.QTable`.

        Has metadata
        Has columns

        """
        # start with name
        params = {"name": [self._instance.name]}
        # get all the immutable inputs, wrapped into lists for QTable
        params.update({k: [v] for k, v in self._instance._init_arguments.items()
                       if k not in ("name", "meta")})
        # get the mutable metadata
        meta = copy.deepcopy(self._instance.meta)
        # and document the cosmology class
        meta["cosmology"] = self._cls.__qualname__

        table = QTable(params, meta=meta)
        return table
