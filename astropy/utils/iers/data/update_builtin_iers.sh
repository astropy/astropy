#!/bin/sh

set -euv

# This script should be run every time an astropy release is made.
# It downloads up-to-date versions of the earth rotation and leap
# second tables.

rm Leap_Second.dat
rm eopc04_IAU2000.62-now

# iers.IERS_B_URL
wget http://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now
# iers.IERS_LEAP_SECOND_URL
wget https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat
