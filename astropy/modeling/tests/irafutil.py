# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module provides functions to help with testing against iraf tasks.
"""

import numpy as np

from astropy.logger import log


class IdentifyRecord:
    """
    Represents a database record for the onedspec.identify task.

    Attributes
    ----------
    recstr: string
            the record as a string
    fields: dict
            the fields in the record
    x: array
       the X values of the identified features
       this represents values on axis1 (image rows)
    y: int
       the Y values of the identified features
       (image columns)
    z: array
       the values which X maps into
    coeff: array
        function (modelname) coefficients
    """

    def __init__(self, recstr):
        self.recstr = recstr
        self.fields = self.get_fields()
        self._flatcoeff = self.fields["coefficients"].ravel()
        self.x = self.fields["features"][:, 0]
        self.y = self.get_ydata()
        self.z = self.fields["features"][:, 1]
        self.coeff = self._flatcoeff[4:]

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

    def get_range(self):
        low = self._flatcoeff[2]
        high = self._flatcoeff[3]
        return [low, high]

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
