0.2 (unreleased)
----------------

- New `astropy.units` sub-package

- `astropy.cosmology` update to include cosmologies with variable dark
  energy equations of state. (This includes some API incompatibilities
  with the older Cosmology objects)

- Astropy doc themes moved into `astropy.sphinx` to allow affilated
  packages to access them

- Added option to disable building of "legacy" packages (pyfits, vo,
  etc.)

- `astropy.table` I/O infrastructure for custom readers/writers
  implemented

- New `astropy.time` sub-package

- The value of the astronomical unit (au) has been updated to that
  adopted by IAU 2012 Resolution B2, and the values of the pc and kpc
  constants have been updated to reflect this.


0.1.1 (unreleased)
------------------

- `astropy.io.fits`
  - Fixed a bug where opening a file containing compressed image HDUs in
    'update' mode and then immediately closing it without making any changes
    caused the file to be rewritten unncessarily.

  - Fixed two memory leaks that could occur when writing compressed image data,
    or in some cases when opening files containing compressed image HDUs in
    'update' mode.

  - Fixed a bug where ``ImageHDU.scale(option='old')`` wasn't working at
    all--it was not restoring the image to its original BSCALE and BZERO
    values.

  - Fixed a bug when writing out files containing zero-width table columns,
    where the TFIELDS keyword would be updated incorrectly, leaving the table
    largely unreadable.


0.1 (2012-06-19)
----------------

- Inital release.
