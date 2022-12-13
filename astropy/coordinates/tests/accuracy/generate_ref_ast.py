"""
This series of functions are used to generate the reference CSV files
used by the accuracy tests.  Running this as a command-line script will
generate them all.
"""

import os

import numpy as np

from astropy.table import Column, Table


def ref_fk4_no_e_fk4(fnout="fk4_no_e_fk4.csv"):
    """
    Accuracy tests for the FK4 (with no E-terms of aberration) to/from FK4
    conversion, with arbitrary equinoxes and epoch of observation.
    """

    import starlink.Ast as Ast

    np.random.seed(12345)

    N = 200

    # Sample uniformly on the unit sphere. These will be either the FK4
    # coordinates for the transformation to FK5, or the FK5 coordinates for the
    # transformation to FK4.
    ra = np.random.uniform(0.0, 360.0, N)
    dec = np.degrees(np.arcsin(np.random.uniform(-1.0, 1.0, N)))

    # Generate random observation epoch and equinoxes
    obstime = [f"B{x:7.2f}" for x in np.random.uniform(1950.0, 2000.0, N)]

    ra_fk4ne, dec_fk4ne = [], []
    ra_fk4, dec_fk4 = [], []

    for i in range(N):
        # Set up frames for AST
        frame_fk4ne = Ast.SkyFrame(f"System=FK4-NO-E,Epoch={obstime[i]},Equinox=B1950")
        frame_fk4 = Ast.SkyFrame(f"System=FK4,Epoch={obstime[i]},Equinox=B1950")

        # FK4 to FK4 (no E-terms)
        frameset = frame_fk4.convert(frame_fk4ne)
        coords = np.degrees(frameset.tran([[np.radians(ra[i])], [np.radians(dec[i])]]))
        ra_fk4ne.append(coords[0, 0])
        dec_fk4ne.append(coords[1, 0])

        # FK4 (no E-terms) to FK4
        frameset = frame_fk4ne.convert(frame_fk4)
        coords = np.degrees(frameset.tran([[np.radians(ra[i])], [np.radians(dec[i])]]))
        ra_fk4.append(coords[0, 0])
        dec_fk4.append(coords[1, 0])

    # Write out table to a CSV file
    t = Table()
    t.add_column(Column(name="obstime", data=obstime))
    t.add_column(Column(name="ra_in", data=ra))
    t.add_column(Column(name="dec_in", data=dec))
    t.add_column(Column(name="ra_fk4ne", data=ra_fk4ne))
    t.add_column(Column(name="dec_fk4ne", data=dec_fk4ne))
    t.add_column(Column(name="ra_fk4", data=ra_fk4))
    t.add_column(Column(name="dec_fk4", data=dec_fk4))
    f = open(os.path.join("data", fnout), "wb")
    f.write(
        f"# This file was generated with the {os.path.basename(__file__)} script, and"
        " the reference values were computed using AST\n"
    )
    t.write(f, format="ascii", delimiter=",")


def ref_fk4_no_e_fk5(fnout="fk4_no_e_fk5.csv"):
    """
    Accuracy tests for the FK4 (with no E-terms of aberration) to/from FK5
    conversion, with arbitrary equinoxes and epoch of observation.
    """

    import starlink.Ast as Ast

    np.random.seed(12345)

    N = 200

    # Sample uniformly on the unit sphere. These will be either the FK4
    # coordinates for the transformation to FK5, or the FK5 coordinates for the
    # transformation to FK4.
    ra = np.random.uniform(0.0, 360.0, N)
    dec = np.degrees(np.arcsin(np.random.uniform(-1.0, 1.0, N)))

    # Generate random observation epoch and equinoxes
    obstime = [f"B{x:7.2f}" for x in np.random.uniform(1950.0, 2000.0, N)]
    equinox_fk4 = [f"B{x:7.2f}" for x in np.random.uniform(1925.0, 1975.0, N)]
    equinox_fk5 = [f"J{x:7.2f}" for x in np.random.uniform(1975.0, 2025.0, N)]

    ra_fk4, dec_fk4 = [], []
    ra_fk5, dec_fk5 = [], []

    for i in range(N):
        # Set up frames for AST
        frame_fk4 = Ast.SkyFrame(
            f"System=FK4-NO-E,Epoch={obstime[i]},Equinox={equinox_fk4[i]}"
        )
        frame_fk5 = Ast.SkyFrame(
            f"System=FK5,Epoch={obstime[i]},Equinox={equinox_fk5[i]}"
        )

        # FK4 to FK5
        frameset = frame_fk4.convert(frame_fk5)
        coords = np.degrees(frameset.tran([[np.radians(ra[i])], [np.radians(dec[i])]]))
        ra_fk5.append(coords[0, 0])
        dec_fk5.append(coords[1, 0])

        # FK5 to FK4
        frameset = frame_fk5.convert(frame_fk4)
        coords = np.degrees(frameset.tran([[np.radians(ra[i])], [np.radians(dec[i])]]))
        ra_fk4.append(coords[0, 0])
        dec_fk4.append(coords[1, 0])

    # Write out table to a CSV file
    t = Table()
    t.add_column(Column(name="equinox_fk4", data=equinox_fk4))
    t.add_column(Column(name="equinox_fk5", data=equinox_fk5))
    t.add_column(Column(name="obstime", data=obstime))
    t.add_column(Column(name="ra_in", data=ra))
    t.add_column(Column(name="dec_in", data=dec))
    t.add_column(Column(name="ra_fk5", data=ra_fk5))
    t.add_column(Column(name="dec_fk5", data=dec_fk5))
    t.add_column(Column(name="ra_fk4", data=ra_fk4))
    t.add_column(Column(name="dec_fk4", data=dec_fk4))
    f = open(os.path.join("data", fnout), "wb")
    f.write(
        f"# This file was generated with the {os.path.basename(__file__)} script, and"
        " the reference values were computed using AST\n"
    )
    t.write(f, format="ascii", delimiter=",")


def ref_galactic_fk4(fnout="galactic_fk4.csv"):
    """
    Accuracy tests for the ICRS (with no E-terms of aberration) to/from FK5
    conversion, with arbitrary equinoxes and epoch of observation.
    """

    import starlink.Ast as Ast

    np.random.seed(12345)

    N = 200

    # Sample uniformly on the unit sphere. These will be either the ICRS
    # coordinates for the transformation to FK5, or the FK5 coordinates for the
    # transformation to ICRS.
    lon = np.random.uniform(0.0, 360.0, N)
    lat = np.degrees(np.arcsin(np.random.uniform(-1.0, 1.0, N)))

    # Generate random observation epoch and equinoxes
    obstime = [f"B{x:7.2f}" for x in np.random.uniform(1950.0, 2000.0, N)]
    equinox_fk4 = [f"J{x:7.2f}" for x in np.random.uniform(1975.0, 2025.0, N)]

    lon_gal, lat_gal = [], []
    ra_fk4, dec_fk4 = [], []

    for i in range(N):
        # Set up frames for AST
        frame_gal = Ast.SkyFrame(f"System=Galactic,Epoch={obstime[i]}")
        frame_fk4 = Ast.SkyFrame(
            f"System=FK4,Epoch={obstime[i]},Equinox={equinox_fk4[i]}"
        )

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
    t.add_column(Column(name="equinox_fk4", data=equinox_fk4))
    t.add_column(Column(name="obstime", data=obstime))
    t.add_column(Column(name="lon_in", data=lon))
    t.add_column(Column(name="lat_in", data=lat))
    t.add_column(Column(name="ra_fk4", data=ra_fk4))
    t.add_column(Column(name="dec_fk4", data=dec_fk4))
    t.add_column(Column(name="lon_gal", data=lon_gal))
    t.add_column(Column(name="lat_gal", data=lat_gal))
    f = open(os.path.join("data", fnout), "wb")
    f.write(
        f"# This file was generated with the {os.path.basename(__file__)} script, and"
        " the reference values were computed using AST\n"
    )
    t.write(f, format="ascii", delimiter=",")


def ref_icrs_fk5(fnout="icrs_fk5.csv"):
    """
    Accuracy tests for the ICRS (with no E-terms of aberration) to/from FK5
    conversion, with arbitrary equinoxes and epoch of observation.
    """

    import starlink.Ast as Ast

    np.random.seed(12345)

    N = 200

    # Sample uniformly on the unit sphere. These will be either the ICRS
    # coordinates for the transformation to FK5, or the FK5 coordinates for the
    # transformation to ICRS.
    ra = np.random.uniform(0.0, 360.0, N)
    dec = np.degrees(np.arcsin(np.random.uniform(-1.0, 1.0, N)))

    # Generate random observation epoch and equinoxes
    obstime = [f"B{x:7.2f}" for x in np.random.uniform(1950.0, 2000.0, N)]
    equinox_fk5 = [f"J{x:7.2f}" for x in np.random.uniform(1975.0, 2025.0, N)]

    ra_icrs, dec_icrs = [], []
    ra_fk5, dec_fk5 = [], []

    for i in range(N):
        # Set up frames for AST
        frame_icrs = Ast.SkyFrame(f"System=ICRS,Epoch={obstime[i]}")
        frame_fk5 = Ast.SkyFrame(
            f"System=FK5,Epoch={obstime[i]},Equinox={equinox_fk5[i]}"
        )

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
    t.add_column(Column(name="equinox_fk5", data=equinox_fk5))
    t.add_column(Column(name="obstime", data=obstime))
    t.add_column(Column(name="ra_in", data=ra))
    t.add_column(Column(name="dec_in", data=dec))
    t.add_column(Column(name="ra_fk5", data=ra_fk5))
    t.add_column(Column(name="dec_fk5", data=dec_fk5))
    t.add_column(Column(name="ra_icrs", data=ra_icrs))
    t.add_column(Column(name="dec_icrs", data=dec_icrs))
    f = open(os.path.join("data", fnout), "wb")
    f.write(
        f"# This file was generated with the {os.path.basename(__file__)} script, and"
        " the reference values were computed using AST\n"
    )
    t.write(f, format="ascii", delimiter=",")


if __name__ == "__main__":
    ref_fk4_no_e_fk4()
    ref_fk4_no_e_fk5()
    ref_galactic_fk4()
    ref_icrs_fk5()
