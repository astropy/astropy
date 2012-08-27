0.2 (unreleased)
----------------

- `astropy.cosmology` update to include cosmologies with variable dark energy equations of state. (This includes some API incompatibilities with the older Cosmology objects)
- Astropy doc themes moved into `astropy.sphinx` to allow affilated packages to access them
- Added option to disable building of "legacy" packages (pyfits, vo, etc.)
- `astropy.table` I/O infrastructure for custom readers/writers implemented
- Besancon Galaxy models (http://model.obs-besancon.fr) reader implemented in `astropy.ascii`
- new `astropy.time` sub-package


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
