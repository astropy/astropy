import re

import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord


def _sexagesimal(g):
    # convert matched regex groups to sexigesimal array
    sign, h, m, s, frac = g
    sign = [1, -1][sign == '-']
    s = '.'.join((s, frac))
    return sign * np.array([h, m, s], float)


class Jparser(object):
    """
    Class that extracts coordinates from object names containing coordinated
    eg: '2MASS J06495091-0737408'
    """
    basePattern = '([+-]?)(\d{1,2})(\d{2})(\d{2})\.?(\d{0,3})'
    single_parser = re.compile(basePattern)
    fullPattern = '(.*?J)' + basePattern * 2
    parser = re.compile(fullPattern)

    def __call__(self, name):
        """Convert to `name` to `SkyCoords` object"""
        return self.to_skycoord(name)

    def match(self, name, raise_=False):
        """Regex match for coordinates in name"""
        # extract the coordinate data from name
        match = self.parser.search(name)
        if match is None and raise_:
            raise ValueError('No coordinate match found!')
        return match

    def ra_dec(self, name):
        """get RA in hourangle and DEC in degrees by parsing name """
        groups = self.match(name, True).groups()
        prefix, hms, dms = np.split(groups, [1, 6])
        ra = (_sexagesimal(hms) / (1, 60, 60 * 60) * u.hourangle).sum()
        dec = (_sexagesimal(dms) * (u.deg, u.arcmin, u.arcsec)).sum()
        return ra, dec

    def to_skycoord(self, name, frame='icrs'):
        return SkyCoord(*self.ra_dec(name), frame=frame)

    def shorten(self, name):
        """
        Produce a shortened version of the full object name using: the prefix
        (usually the survey name) and RA (hour, minute), DEC (deg, arcmin) parts.
            e.g.: '2MASS J06495091-0737408' --> '2MASS J0649-0737'

        Parameters
        ----------
        name : str
            Full object name with J-coords embedded.

        Returns
        -------
        shortName: str
        """
        match = self.match(name)
        return ''.join(match.group(1, 3, 4, 7, 8, 9))


jparser = Jparser()

if __name__ == '__main__':
    # a few test cases:
    for name in ['CRTS SSS100805 J194428-420209',
                 'MASTER OT J061451.7-272535.5',
                 '2MASS J06495091-0737408',
                 '1RXS J042555.8-194534',
                 'SDSS J132411.57+032050.5',
                 'DENIS-P J203137.5-000511',
                 '2QZ J142438.9-022739',
                 'CXOU J141312.3-652013']:
        print(name, jparser(name))
