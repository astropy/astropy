import re
import warnings

from astropy import units as u
from astropy.coordinates import Angle
from matplotlib.ticker import Formatter

# Define a formatter based on the Astropy Angle class. Of course, we can add
# formatting options later, this is just a proof of concept.


DMS_RE = re.compile('dd(:mm(:ss(.(s)+)?)?)?')
HMS_RE = re.compile('hh(:mm(:ss(.(s)+)?)?)?')
DDEC_RE = re.compile('d(.(d)+)?')


def re_exact_match(pattern, string):
    m = re.match(pattern, string)
    if m is None:
        return False
    else:
        return m.start() == 0 and m.end() == len(string)


class AngleFormatter(Formatter):

    def __init__(self, precision=1):
        self.precision = precision

    @property
    def format(self):
        return self._format

    @format.setter
    def format(self, value):
        if re_exact_match(DMS_RE, value) or re_exact_match(HMS_RE, value):
            self._decimal = False
            self._unit = u.degree
            if len(value) > 8:
                self._precision = len(value) - 9
            else:
                self._precision = 0
                warnings.warn("This format is not yet supported, falling back to dd:mm:ss")
        elif re_exact_match(HMS_RE, value):
            self._decimal = False
            self._unit = u.hourangle
            if len(value) > 8:
                self._precision = len(value) - 9
            else:
                self._precision = 0
                warnings.warn("This format is not yet supported, falling back to hh:mm:ss")
        elif re_exact_match(DDEC_RE, value):
            self._decimal = True
            self._unit = u.degree
            self._precision = max(0, len(value) - 2)
        else:
            raise ValueError("Invalid format: {0}".format(value))

    def __call__(self, axis, other, value, **kwargs):
        if len(value) > 0:
            angles = Angle(value, unit=u.deg)
            string = angles.to_string(unit=self._unit,
                                      precision=self._precision,
                                      decimal=self._decimal).tolist()
            return string
        else:
            return []
