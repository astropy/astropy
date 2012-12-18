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
ra = np.random.uniform(0., 360., N)
dec = np.degrees(np.arcsin(np.random.uniform(-1., 1., N)))

# Generate random observation epoch and equinoxes
obstime = ["B{0:7.2f}".format(x) for x in np.random.uniform(1950., 2000., N)]
equinox_fk5 = ["J{0:7.2f}".format(x) for x in np.random.uniform(1975., 2025., N)]

ra_icrs, dec_icrs = [], []
ra_fk5, dec_fk5 = [], []

for i in range(N):

    # Set up frames for AST
    frame_icrs = Ast.SkyFrame('System=ICRS,Epoch={epoch}'.format(epoch=obstime[i]))
    frame_fk5 = Ast.SkyFrame('System=FK5,Epoch={epoch},Equinox={equinox_fk5}'.format(epoch=obstime[i],equinox_fk5=equinox_fk5[i]))

    # ICRS to FK5
    frameset = frame_icrs.convert(frame_fk5)
    coords = np.degrees(frameset.tran([[np.radians(ra[i])], [np.radians(dec[i])]]))
    ra_fk5.append(coords[0, 0])
    dec_fk5.append(coords[1, 0])

    # FK5 to ICRS
    frameset = frame_fk5.convert(frame_icrs)
    coords = np.degrees(frameset.tran([[np.radians(ra[i])], [np.radians(dec[i])]]))
    ra_icrs.append(coords[0, 0])
    dec_icrs.append(coords[1, 0])

# Write out table to a CSV file
t = Table()
t.add_column(Column('equinox_fk5', equinox_fk5))
t.add_column(Column('obstime', obstime))
t.add_column(Column('ra_in', ra))
t.add_column(Column('dec_in', dec))
t.add_column(Column('ra_fk5', ra_fk5))
t.add_column(Column('dec_fk5', dec_fk5))
t.add_column(Column('ra_icrs', ra_icrs))
t.add_column(Column('dec_icrs', dec_icrs))
f = open('icrs_fk5.csv', 'wb')
f.write("# This file was generated with the {0} script, and the reference "
        "values were computed using AST\n".format(os.path.basename(__file__)))
t.write(f, format='ascii', delimiter=',')