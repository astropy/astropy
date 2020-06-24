.. _validation:

Validation and Bounds checking
******************************

Bounds checking is enabled by default, and any computed world
coordinates outside of [-180°, 180°] for longitude and [-90°, 90°] in
latitude are marked as invalid.  To disable this behavior, use
`astropy.wcs.Wcsprm.bounds_check`.