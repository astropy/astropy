from astropy.coords import Angle, Coord

"""
Angles and coordinates alternate API.

Features:

- Exactly one import statement that is the same for all use cases.
- Exactly one way to make an angle and one way to make a coordinate.
- Every output and transformation is done by attribute access.
- Object behavior is controlled by class initialization or object attrs.
- Very similar to Time.
"""

class Angle(angle, format=None, precision=3, seps=None, pad=False):
    """
    Represent one or more angular distance(s).

    ``angle`` can be a single string or number, or an sequence of
    strings or numbers.

    ``format`` specifies the representation of an angle.  This
    is like ``unit`` but more generic.

    angle : input angle(s) in specified format (scalar or sequence)
    format : angle format
    precision : decimal precision for string-output
    seps : sexigesmal separators, can be list (default=':')
    pad : Add leading zeros as needed (default=False)
    """

a = Angle(180.0, 'deg', precision=2)
a.rad  # 3.14159..
a.hour  # 12.0
a.hms  # '12:00:00.00'
a.hms_hour  # 12 (int)
a.hms_min  # 0
a.hms_sec  # 0
a.dms_deg  # 180 (int)
a.dms_min  # 0
a.dms_sec  # 0

a = Angle('12:00:00', 'hms')  # 'hms' required because of ambiguity
a.deg  # 180.0

a = Angle('180:00:00.000')  # unambiguous
a = Angle('180d10m45.25s')

a = Angle(180.0, 'deg', seps=('h', 'm', 's'))
a.hms  # '12h00m00.000s'

# I hate anything but degrees, so let me just declare that
# at the top of my coe.
Angle.set_default_format('deg')

a = Angle([0, 90, 180])
a.hour  # [0.0, 6.0, 12.0]  (really numpy arrays)
a.hms_hour  # [0, 6, 12]

# Formatting control via object initialization or attributes or
# global defaults (as applicable when object is created)
a.pad = True
a.precision = 1
a.hms  # ['00:00:00.0', '06:00:00.0', '12:00:00.0']

# Note: every input or output is a single scalar or one-component array.
# No tuples are ever output because these start getting unwieldy
# for multi-dimensional output.  For the Time class there was
# general consensus that one should be able to input an arbitrarily
# shaped array and get that same shape back out.  If you allow
# for attrs that return a 3-tuple in the case of a scalar, it
# forces people to remember which index the new dimension is in.


class Coord(angle1, angle2, format1=None, format2=None, system=None,
            PLUS OTHERS like equinox, epoch, precision, seps, pad):
    """A coordinate is just two angles with a coordinate system.
    (Is that naive??)

    This stipulates that user have to provide the two coordinates
    as separate values or arrays.  Inputs like "180.25, -23.5" are
    not allowed, but this greatly simplifies the API.  (IMO)
    """

c = Coord(90, 270, 'deg', 'deg', 'equatorial')
c.ra  # 90
c.dec  # 270

c = Coord([90, 180], [270, 271], 'deg', 'deg', 'equatorial' precision=0)
c.ra  # [90, 180]
c.dec  # [270, 271]

c.galactic.b  # whatever that would be
c.b  # Do an implicit transformation to the galactic scale.  Not
     # sure if there are difficulties with this.

# OK, here is where I'm not sure this can work but it would really sweet:

c.ra.hms  # ['6:00:00', '12:00:00']
c.ra.hms_hour  # [6, 12]

# The trick is that c.ra returns an instance of a subclass of
# float or str that has all of the extra format attributes
# available.  OR likewise an ndarray subclass with extra format
# attributes or __getattr__, etc.  In a normal context
# c.ra looks like a float or string or array, but one can then
# get the representation in a different format.

# OTHERWISE, the straightforward solution that requires an extra line

c.format1 = 'hms'
c.ra  # ['6:00:00', '12:00:00']

c2 = c.galactic
c2.l  # galactic "l" values

# Arguments which correspond to a Time can be specified as anything
# that will correctly initialize the Time() class.
c2 = c.precess(equinox=Time(50000.0, format='mjd', scale='utc'))
c2 = c.precess(equinox=Time('J2050.0'))  # Long form
c2 = c.precess(equinox='J2050.0')  # Arg gets passed to Time()

c = Coord(alt, az, 'deg', 'deg', system='alt_az')

# If I'm lucky enough that my source lists are all in decimal degrees...
Coord.set_default_format1('deg')
Coord.set_default_format2('deg')
c = Coord(alt, az, system='alt_az')
