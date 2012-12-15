# Accuracy tests for the FK4 (with no E-terms of aberration) to/from FK5
# conversion, with arbitrary equinoxes and epoch of observation.

import numpy as np
import starlink.Ast as Ast
from astropy.table import Table, Column

np.random.seed(12345)

N = 200

# Sample uniformly on the unit sphere. These will be either the FK4
# coordinates for the transformation to FK5, or the FK5 coordinates for the
# transformation to FK4.
ra = np.random.uniform(0., 360., N)
dec = np.degrees(np.arcsin(np.random.uniform(-1., 1., N)))

# Generate random observation epoch and equinoxes
obstime = ["B{0:7.2f}".format(x) for x in np.random.uniform(1950., 2000., N)]
equinox_in = ["B{0:7.2f}".format(x) for x in np.random.uniform(1925., 1975., N)]
equinox_out = ["J{0:7.2f}".format(x) for x in np.random.uniform(1975., 2025., N)]

ra_fk4, dec_fk4 = [], []
ra_fk5, dec_fk5 = [], []

for i in range(N):

    # Set up frames for AST
    frame_fk4 = Ast.SkyFrame('System=FK4-NO-E,Epoch={epoch},Equinox={equinox_in}'.format(epoch=obstime[i], equinox_in=equinox_in[i]))
    frame_fk5 = Ast.SkyFrame('System=FK5,Equinox={equinox_out}'.format(equinox_out=equinox_out[i]))

    # FK4 to FK5
    frameset = frame_fk4.convert(frame_fk5)
    coords = np.degrees(frameset.tran([[np.radians(ra[i])], [np.radians(dec[i])]]))
    ra_fk5.append(coords[0,0])
    dec_fk5.append(coords[1,0])

    # FK5 to FK4
    frameset = frame_fk5.convert(frame_fk4)
    coords = np.degrees(frameset.tran([[np.radians(ra[i])], [np.radians(dec[i])]]))
    ra_fk4.append(coords[0,0])
    dec_fk4.append(coords[1,0])

# Write out table to a CSV file
t = Table()
t.add_column(Column('equinox_in', equinox_in))
t.add_column(Column('equinox_out', equinox_out))
t.add_column(Column('obstime', obstime))
t.add_column(Column('ra_in', ra))
t.add_column(Column('dec_in', dec))
t.add_column(Column('ra_fk5', ra_fk5))
t.add_column(Column('dec_fk5', dec_fk5))
t.add_column(Column('ra_fk4', ra_fk4))
t.add_column(Column('dec_fk4', dec_fk4))
t.write('fk4_no_e_fk5.csv', format='ascii', delimiter=',')
