# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Planetology constants in SI units.  See :mod:`astropy.constants`
for a complete listing of constants defined in Astropy.
"""

import numpy as np

from .constant import Constant

# PLANETOLOGY CONSTANTS


class Planetology(Constant):
    default_reference = 'Planetology'
    _registry = {}
    _has_incompatible_units = set()


# Planet equatorial radius and mass
# source : https://ssd.jpl.nasa.gov/planets/phys_par.html
R_mercury = Planetology("R_mercury", "Mercury equatorial radius", 2.44053e6, "m",
                        .04e3,
                        "Archinal, B.A. et al. 2018. 'Report of the IAU/IAG Working Group on cartographic coordinates and rotational elements: 2015'",
                        system="si")

M_mercury = Planetology("M_mercury", "Mercury mass", .330103e24, "kg", .000021e24,
        "Anderson, J.D., et al. 1987. 'The mass, gravity field, and ephemeris of Mercury' Icarus 71:337-349.", system="si")

# VENUS
R_venus = Planetology("R_venus", "Venus equatorial radius", 6.0518e6, "m", 1.e3,
                "Archinal, B.A. et al. 2018. 'Report of the IAU/IAG Working Group on cartographic coordinates and rotational elements: 2015'",
                system="si")

M_venus = Planetology("M_venus", "Venus mass", 4.86731e24, "kg", .00023e24,
        "Konopliv, A.S., et al. 1999. 'Venus gravity: 180th degree and order model' Icarus 139:3-18.",
        system="si")

# EARTH
R_earth = Planetology("R_earth", "Earth equatorial radius", 6.3781366e6, "m", .0001e3,
                "Archinal, B.A. et al. 2018. 'Report of the IAU/IAG Working Group on cartographic coordinates and rotational elements: 2015'",
                system="si")

M_earth = Planetology("M_earth", "Earth mass", 5.97217e24, "kg", .00028e24,
        "Folkner, W.M. and Williams, J.G. 2008. 'Mass parameters and uncertainties in planetary ephemeris DE421.' Interoffice Memo. 343R-08-004 (internal document), Jet Propulsion Laboratory, Pasadena, CA",
        system="si")

# MARS
R_mars = Planetology("R_mars", "Mars equatorial radius", 3.396e6, "m", .1e3,
                "Archinal, B.A. et al. 2018. 'Report of the IAU/IAG Working Group on cartographic coordinates and rotational elements: 2015'",
                system="si")

M_mars = Planetology("M_mars", "Mars mass", .641691e24, "kg", .00003e24,
        "Jacobson, R.A. 2008. 'Ephemerides of the Martian Satellites - MAR080' Interoffice Memo. 343R-08-006 (internal document), Jet Propulsion Laboratory, Pasadena, CA",
        system="si")

# JUPYTER
R_jupiter = Planetology("R_jupiter", "Jupyter equatorial radius", 7.1492e7, "m", 4.e3,
                "Archinal, B.A. et al. 2018. 'Report of the IAU/IAG Working Group on cartographic coordinates and rotational elements: 2015'",
                system="si")

M_jupiter = Planetology("M_jupiter", "Jupyter mass", 1898.125e24, "kg", .088e24,
        "Jacobson, R.A. 2013. 'Jovian Satellite ephemeris - JUP310' private communication",
        system="si")

# SATURN
R_saturn = Planetology("R_saturn", "Saturn equatorial radius", 60.268e6, "m", 4.e3,
                "Archinal, B.A. et al. 2018. 'Report of the IAU/IAG Working Group on cartographic coordinates and rotational elements: 2015'",
                system="si")

M_saturn = Planetology("M_saturn", "Saturn mass", 568.317e24, "kg", .026e24,
        "Jacobson, R.A., et al. 2006. 'The gravity field of the Saturnian system from satellite observations and spacecraft tracking data' AJ 132(6):2520-2526",
        system="si")

# URANUS
R_uranus = Planetology("R_uranus", "Uranus equatorial radius", 25.559e6, "m", 4.e3,
                "Archinal, B.A. et al. 2018. 'Report of the IAU/IAG Working Group on cartographic coordinates and rotational elements: 2015'",
                system="si")

M_uranus = Planetology("M_uranus", "Uranus mass", 86.8099e24, "kg", .004e24,
        "Jacobson, R.A. 2014. 'The Orbits of the Uranian Satellites and Rings, the Gravity Field of the Uranian System, and the Orientation of the Pole of Uranus' AJ 148:76-88",
        system="si")

# NEPTUNE
R_neptune = Planetology("R_neptune", "Neptune equatorial radius", 24.764e6, "m", 15.e3,
                "Archinal, B.A. et al. 2018. 'Report of the IAU/IAG Working Group on cartographic coordinates and rotational elements: 2015'",
                system="si")

M_neptune = Planetology("M_neptune", "Neptune mass", 102.4092e24, "kg", .0048e24,
        "Jacobson, R.A. 2009. 'The orbits of the Neptunian satellites and the orientation of the pole of Neptune' AJ 137:4322",
        system="si")



# Dwarf planets

# CERES
R_ceres = Planetology("R_ceres", "Ceres equatorial radius", 482.e3, "m", 0.,
        "R.S. Park, et al. (2016) 'A partially differentiated interior for (1) Ceres deduced from its gravity field and shape' Nature 537:515-517",
                system="si")

M_ceres = Planetology("M_ceres", "Ceres mass", 938.416e18, "kg", .013e18,
        "R.S. Park, et al. (2016) 'A partially differentiated interior for (1) Ceres deduced from its gravity field and shape' Nature 537:515-517",
        system="si")


# PLUTO
R_pluto = Planetology("R_pluto", "Pluto equatorial radius", 1.188e6, "m", 1.6e3,
                "Archinal, B.A. et al. 2018. 'Report of the IAU/IAG Working Group on cartographic coordinates and rotational elements: 2015'",
                system="si")

M_pluto = Planetology("M_pluto", "Pluto mass", 13029.e18, "kg", 27.e18,
        "Brozovic, M. et al. 2015. 'The orbits and masses of satellites of Pluto' Icarus 246:317-329.",
        system="si")


# ERIS
R_eris = Planetology("R_eris", "Eris equatorial radius", 1.200e6, "m", 50.e3,
                "Brown, Michael E.; Schaller, Emily L. (2007) Science 316:1585",
                system="si")

M_eris = Planetology("M_eris", "Eris mass", 16600.e18, "kg", 200.e18,
                "Brown, Michael E.; Schaller, Emily L. (2007) Science 316:1585",
                system="si")

# MAKEMAKE
R_makemake = Planetology("R_makemake", "Makemake equatorial radius", .717e6, "m", 7.e3,
                "Brown, Michael E. (2013) ApJ 767, article id. L7",
                system="si")

M_makemake = Planetology("M_makemake", "Makemake mass", 3100.e18, "kg", 0.e18,
                "Parker, Alex et al. (2018) DPS meeting# 50 id. 509.02",
                system="si")

# HAUMEA
R_haumea = Planetology("R_haumea", "Haumea equatorial radius", .870e6, "m", 0.,
                "Lockwood, Alexandra C.; Brown, Michael E.; Stansberry, John. (2014) Earth, Moon, Planets 111:127-137",
                system="si")

M_haumea = Planetology("M_haumea", "Haumea mass", 4006.e18, "kg", 40.e18,
                "Parker, Alex et al. (2018) DPS meeting# 50 id. 509.02",
                system="si")
