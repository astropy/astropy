# -*- coding: utf-8 -*-
"""
=====================================================================
Accessing data stored as a table in a multi-extension FITS (MEF) file
=====================================================================

FITS files can often contain large amount of multi-dimensional data and
tables. This example opens a FITS file with information
from Chandra's HETG-S instrument.

The example uses `astropy.utils.data` to download multi-extension FITS (MEF)
file, `astropy.io.fits` to investigate the header, and
`astropy.table.Table` to explore the data.


*By: Lia Corrales, Adrian Price-Whelan, and Kelle Cruz*

*License: BSD*


"""

##############################################################################
# Use `astropy.utils.data` subpackage to download the FITS file used in this
# example. Also import `~astropy.table.Table` from the `astropy.table` subpackage
# and `astropy.io.fits`

from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table
from astropy.io import fits

##############################################################################
# Download a FITS file

event_filename = get_pkg_data_filename('tutorials/FITS-tables/chandra_events.fits')

##############################################################################
# Display information about the contents of the FITS file.

fits.info(event_filename)

##############################################################################
# Extension 1, EVENTS, is a Table that contains information about each X-ray
# photon that hit Chandra's HETG-S detector.
#
# Use `~astropy.table.Table` to read the table

events = Table.read(event_filename, hdu=1)

##############################################################################
# Print the column names of the Events Table.

print(events.columns)

##############################################################################
# If a column contains unit information, it will have an associated
# `astropy.units` object.

print(events['energy'].unit)

##############################################################################
# Print the data stored in the Energy column.

print(events['energy'])
