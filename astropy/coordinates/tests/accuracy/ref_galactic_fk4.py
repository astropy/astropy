# Accuracy tests for the ICRS (with no E-terms of aberration) to/from FK5
# conversion, with arbitrary equinoxes and epoch of observation.

import os
import numpy as np
import starlink.Ast as Ast
from astropy.table import Table, Column

np.random.seed(12345)

N = 200

# Sample uniformly on the unit sphere. These will be either the ICRS
# coordinates for the transformation to FK5, or the FK5 coordinates for the
# transformation to ICRS.
lon = np.random.uniform(0., 360., N)
lat = np.degrees(np.arcsin(np.random.uniform(-1., 1., N)))

# Generate random observation epoch and equinoxes
obstime = ["B{0:7.2f}".format(x) for x in np.random.uniform(1950., 2000., N)]
equinox_fk4 = ["J{0:7.2f}".format(x) for x in np.random.uniform(1975., 2025., N)]

lon_gal, lat_gal = [], []
ra_fk4, dec_fk4 = [], []

for i in range(N):

    # Set up frames for AST
    frame_gal = Ast.SkyFrame('System=Galactic,Epoch={epoch}'.format(epoch=obstime[i]))
    frame_fk4 = Ast.SkyFrame('System=FK4,Epoch={epoch},Equinox={equinox_fk4}'.format(epoch=obstime[i], equinox_fk4=equinox_fk4[i]))

    # ICRS to FK5
    frameset = frame_gal.convert(frame_fk4)
    coords = np.degrees(frameset.tran([[np.radians(lon[i])], [np.radians(lat[i])]]))
    ra_fk4.append(coords[0, 0])
    dec_fk4.append(coords[1, 0])

    # FK5 to ICRS
    frameset = frame_fk4.convert(frame_gal)
    coords = np.degrees(frameset.tran([[np.radians(lon[i])], [np.radians(lat[i])]]))
    lon_gal.append(coords[0, 0])
    lat_gal.append(coords[1, 0])

# Write out table to a CSV file
t = Table()
t.add_column(Column('equinox_fk4', equinox_fk4))
t.add_column(Column('obstime', obstime))
t.add_column(Column('lon_in', lon))
t.add_column(Column('lat_in', lat))
t.add_column(Column('ra_fk4', ra_fk4))
t.add_column(Column('dec_fk4', dec_fk4))
t.add_column(Column('lon_gal', lon_gal))
t.add_column(Column('lat_gal', lat_gal))
f = open('galactic_fk4.csv', 'wb')
f.write("# This file was generated with the {0} script, and the reference "
        "values were computed using AST\n".format(os.path.basename(__file__)))
t.write(f, format='ascii', delimiter=',')
