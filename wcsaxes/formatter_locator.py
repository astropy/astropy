import re

import numpy as np

from astropy import units as u
from astropy.coordinates import Angle


DMS_RE = re.compile('^dd(:mm(:ss(.(s)+)?)?)?$')
HMS_RE = re.compile('^hh(:mm(:ss(.(s)+)?)?)?$')
DDEC_RE = re.compile('^d(.(d)+)?$')


class AngleFormatterLocator(object):
    """
    A joint formatter/locator
    """

    def __init__(self, values=None, number=None, spacing=None, format='dd:mm:ss'):
        self._values = values
        self._number = number
        self._spacing = spacing
        self.format = format

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, values):
        self._number = None
        self._spacing = None
        self._values = values

    @property
    def number(self):
        return self._number

    @number.setter
    def number(self, number):
        self._number = number
        self._spacing = None
        self._values = None

    @property
    def spacing(self):
        return self._spacing

    @spacing.setter
    def spacing(self, spacing):
        if not isinstance(spacing, u.Quantity):
            raise TypeError("spacing should be a quantity")
        self._number = None
        self._spacing = spacing
        self._values = None

    @property
    def format(self):
        return self._format

    @format.setter
    def format(self, value):
        self._format = value
        if DMS_RE.match(value) is not None:
            self._decimal = False
            self._unit = u.degree
            if '.' in value:
                self._precision = len(value) - value.index('.')
                self._fields = 3
            else:
                self._precision = 0
                self._fields = value.count(':') + 1
        elif HMS_RE.match(value) is not None:
            self._decimal = False
            self._unit = u.hourangle
            if '.' in value:
                self._precision = len(value) - value.index('.')
                self._fields = 3
            else:
                self._precision = 0
                self._fields = value.count(':') + 1
        elif DDEC_RE.match(value) is not None:
            self._decimal = True
            self._unit = u.degree
            self._fields = 1
            if '.' in value:
                self._precision = len(value) - value.index('.')
            else:
                self._precision = 0
        else:
            raise ValueError("Invalid format: {0}".format(value))

    @property
    def base_spacing(self):

        if self._fields == 1:
            spacing = 1. * u.degree
        elif self._fields == 2:
            spacing = 1. * u.arcmin
        elif self._fields == 3:
            if self._precision == 0:
                spacing = 1. * u.arcsec
            else:
                spacing = u.arcsec / (10. ** self.precision)

        if self._unit is u.hourangle:
            spacing *= 15

        return spacing

    def locator(self, v1, v2):
        if self.values is not None:
            return self.values
        elif self.spacing is not None:
            return np.arange(np.floor(v1), np.ceil(v2) + 1., self.spacing.degree)
        elif self.number is not None:
            # do what axisartist does
            pass

    def formatter(self, direction, factor, values, **kwargs):
        if len(values) > 0:
            angles = Angle(np.asarray(values) / factor, unit=u.deg)
            string = angles.to_string(unit=self._unit,
                                      precision=self._precision,
                                      decimal=self._decimal,
                                      fields=self._fields).tolist()
            return string
        else:
            return []
