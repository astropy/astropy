# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module provides functions to help with testing against iraf tasks.
"""

import numpy as np

from astropy.logger import log

iraf_models_map = {1.0: "Chebyshev", 2.0: "Legendre", 3.0: "Spline3", 4.0: "Spline1"}


def get_records(fname):
    """
    Read the records of an IRAF database file into a python list.

    Parameters
    ----------
    fname : str
           name of an IRAF database file

    Returns
    -------
        A list of records
    """
    f = open(fname)
    dtb = f.read()
    f.close()
    recs = dtb.split("begin")[1:]
    records = [Record(r) for r in recs]
    return records


def get_database_string(fname):
    """
    Read an IRAF database file.

    Parameters
    ----------
    fname : str
          name of an IRAF database file

    Returns
    -------
        the database file as a string
    """
    f = open(fname)
    dtb = f.read()
    f.close()
    return dtb


class Record:
    """
    A base class for all records - represents an IRAF database record.

    Attributes
    ----------
    recstr: string
            the record as a string
    fields: dict
            the fields in the record
    taskname: string
            the name of the task which created the database file
    """

    def __init__(self, recstr):
        self.recstr = recstr
        self.fields = self.get_fields()
        self.taskname = self.get_task_name()

    def aslist(self):
        reclist = self.recstr.split("\n")
        reclist = [entry.strip() for entry in reclist]
        [reclist.remove(entry) for entry in reclist if len(entry) == 0]
        return reclist

    def get_fields(self):
        # read record fields as an array
        fields = {}
        flist = self.aslist()
        numfields = len(flist)
        for i in range(numfields):
            line = flist[i]
            if line and line[0].isalpha():
                field = line.split()
                if i + 1 < numfields:
                    if not flist[i + 1][0].isalpha():
                        fields[field[0]] = self.read_array_field(
                            flist[i : i + int(field[1]) + 1]
                        )
                    else:
                        fields[field[0]] = " ".join(s for s in field[1:])
                else:
                    fields[field[0]] = " ".join(s for s in field[1:])
            else:
                continue
        return fields

    def get_task_name(self):
        try:
            return self.fields["task"]
        except KeyError:
            return None

    def read_array_field(self, fieldlist):
        # Turn an iraf record array field into a numpy array
        fieldline = [entry.split() for entry in fieldlist[1:]]
        # take only the first 3 columns
        # identify writes also strings at the end of some field lines
        xyz = [entry[:3] for entry in fieldline]
        try:
            farr = np.array(xyz)
        except Exception:
            log.debug(f"Could not read array field {fieldlist[0].split()[0]}")
        return farr.astype(np.float64)


class IdentifyRecord(Record):
    """
    Represents a database record for the onedspec.identify task.

    Attributes
    ----------
    x: array
       the X values of the identified features
       this represents values on axis1 (image rows)
    y: int
       the Y values of the identified features
       (image columns)
    z: array
       the values which X maps into
    modelname: string
        the function used to fit the data
    nterms: int
        degree of the polynomial which was fit to the data
        in IRAF this is the number of coefficients, not the order
    mrange: list
        the range of the data
    coeff: array
        function (modelname) coefficients
    """

    def __init__(self, recstr):
        super().__init__(recstr)
        self._flatcoeff = self.fields["coefficients"].ravel()
        self.x = self.fields["features"][:, 0]
        self.y = self.get_ydata()
        self.z = self.fields["features"][:, 1]
        self.modelname = self.get_model_name()
        self.nterms = self.get_nterms()
        self.mrange = self.get_range()
        self.coeff = self.get_coeff()

    def get_model_name(self):
        return iraf_models_map[self._flatcoeff[0]]

    def get_nterms(self):
        return self._flatcoeff[1]

    def get_range(self):
        low = self._flatcoeff[2]
        high = self._flatcoeff[3]
        return [low, high]

    def get_coeff(self):
        return self._flatcoeff[4:]

    def get_ydata(self):
        image = self.fields["image"]
        left = image.find("[") + 1
        right = image.find("]")
        section = image[left:right]
        if "," in section:
            yind = image.find(",") + 1
            return int(image[yind:-1])
        else:
            return int(section)


class FitcoordsRecord(Record):
    """
    Represents a database record for the longslit.fitccords task.

    Attributes
    ----------
    modelname: string
        the function used to fit the data
    xorder: int
        number of terms in x
    yorder: int
        number of terms in y
    xbounds: list
        data range in x
    ybounds: list
        data range in y
    coeff: array
        function coefficients

    """

    def __init__(self, recstr):
        super().__init__(recstr)
        self._surface = self.fields["surface"].ravel()
        self.modelname = iraf_models_map[self._surface[0]]
        self.xorder = self._surface[1]
        self.yorder = self._surface[2]
        self.xbounds = [self._surface[4], self._surface[5]]
        self.ybounds = [self._surface[6], self._surface[7]]
        self.coeff = self.get_coeff()

    def get_coeff(self):
        return self._surface[8:]


class IDB:
    """
    Base class for an IRAF identify database.

    Attributes
    ----------
    records: list
             a list of all `IdentifyRecord` in the database
    numrecords: int
             number of records
    """

    def __init__(self, dtbstr):
        self.records = [IdentifyRecord(rstr) for rstr in self.aslist(dtbstr)]
        self.numrecords = len(self.records)

    def aslist(self, dtb):
        # return a list of records
        # if the first one is a comment remove it from the list
        rl = dtb.split("begin")
        try:
            rl0 = rl[0].split("\n")
        except Exception:
            return rl
        if len(rl0) == 2 and rl0[0].startswith("#") and not rl0[1].strip():
            return rl[1:]
        else:
            return rl


class ReidentifyRecord(IDB):
    """
    Represents a database record for the onedspec.reidentify task.
    """

    def __init__(self, databasestr):
        super().__init__(databasestr)
        self.x = np.array([r.x for r in self.records])
        self.y = self.get_ydata()
        self.z = np.array([r.z for r in self.records])

    def get_ydata(self):
        y = np.ones(self.x.shape)
        y = y * np.array([r.y for r in self.records])[:, np.newaxis]
        return y
