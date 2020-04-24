# Script to generate random targets, observatory locations, and times, and
# run these using the Starlink rv command to generate reference values for the
# velocity frame corrections. This requires that Starlink is installed and that
# the rv comamnd is in your PATH. More information about Starlink can be found
# at http://starlink.eao.hawaii.edu/starlink

if __name__ == "__main__":

    from random import choice
    from subprocess import check_output

    import numpy as np

    from astropy.table import QTable
    from astropy.coordinates import SkyCoord, Angle
    from astropy.time import Time
    from astropy import units as u

    np.random.seed(12345)

    N = 100

    tab = QTable()
    target_lon = np.random.uniform(0, 360, N) * u.deg
    target_lat = np.degrees(np.arcsin(np.random.uniform(-1, 1, N))) * u.deg
    tab['target'] = SkyCoord(target_lon, target_lat, frame='fk5')
    tab['obstime'] = Time(np.random.uniform(Time('1997-01-01').mjd, Time('2017-12-31').mjd, N), format='mjd', scale='utc')
    tab['obslon'] = Angle(np.random.uniform(-180, 180, N) * u.deg)
    tab['obslat'] = Angle(np.arcsin(np.random.uniform(-1, 1, N)) * u.deg)
    tab['geocent'] = 0.
    tab['heliocent'] = 0.
    tab['lsrk'] = 0.
    tab['lsrd'] = 0.
    tab['galactoc'] = 0.
    tab['localgrp'] = 0.

    for row in tab:

        # Produce input file for rv command
        with open('rv.input', 'w') as f:
            f.write(row['obslon'].to_string('deg', sep=' ') + ' ' + row['obslat'].to_string('deg', sep=' ') + '\n')
            f.write(f"{row['obstime'].datetime.year} {row['obstime'].datetime.month} {row['obstime'].datetime.day} 1\n")
            f.write(row['target'].to_string('hmsdms', sep=' ') + ' J2000\n')
            f.write('END\n')

        # Run Starlink rv command
        check_output(['rv', 'rv.input'])

        # Parse values from output file
        lis_lines = []
        started = False
        for lis_line in open('rv.lis'):
            if started and lis_line.strip() != '':
                lis_lines.append(lis_line.strip())
            elif 'LOCAL GROUP' in lis_line:
                started = True

        # Some sources are not observable at the specified time and therefore don't
        # have entries in the rv output file
        if len(lis_lines) == 0:
            continue

        # If there are lines, we pick one at random. Note that we can't get rv to
        # run at the exact time we specified in the input, so we will re-parse the
        # actual date/time used and replace it in the table
        lis_line = choice(lis_lines)

        # The column for 'SUN' has an entry also for the light travel time, which
        # we want to ignore. It sometimes includes '(' followed by a space which
        # can cause issues with splitting, hence why we get rid of the space.
        lis_line = lis_line.replace('(  ', '(').replace('( ', '(')
        year, month, day, time, zd, row['geocent'], row['heliocent'], _, \
            row['lsrk'], row['lsrd'], row['galactoc'], row['localgrp'] = lis_line.split()
        row['obstime'] = Time(f'{year}-{month}-{day}T{time}:00', format='isot', scale='utc')

    # We sampled 100 coordinates above since some may not have results - we now
    # truncate to 50 sources since this is sufficient.
    tab[:50].write('reference_rv.ecsv', format='ascii.ecsv')
