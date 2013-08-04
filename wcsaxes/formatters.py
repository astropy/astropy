from astropy import units as u
from astropy.coordinates import Angle
from matplotlib.ticker import Formatter

# Define a formatter based on the Astropy Angle class. Of course, we can add
# formatting options later, this is just a proof of concept.


class AngleFormatter(Formatter):

    def __init__(self, precision=1):
        self.precision = precision

    def __call__(self, axis, other, value, **kwargs):
        if len(value) > 0:
            angles = Angle(value, unit=u.deg)
            string = angles.to_string(precision=self.precision).tolist()
            return string
        else:
            return []
