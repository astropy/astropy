0.4 (unreleased)
----------------

- Astropy frame instances to now be passed to ``get_transform`` and ``get_coords_overlay``. [#149]

- Fixed bug that caused plotting to crash if there were no ticks on an axis.

- Fix compatibility with different versions of Matplotlib. [#153]

- Fix bug when plotting overlays on images with 3+ WCS dimensions.

- Fix bug that occurred in some cases where Matplotlib would try and plot a
  grid as soon as the axes were initialized.

- Fix bug with plotting RGB images. [#15]

0.3 (2014-12-07)
----------------

- Fixed a bug that caused axis labels to be incorrectly placed when hiding tick 
  labels using ``set_visible(False)``. [#111]

- Fixed a bug that caused tick labels to not always appear. [#110]

- Fixed a bug that caused the clip path for images to not be updated when
  making subsequent changes after drawing the figure. [#117]

- Fixed a bug that caused the clip path for images to not be updated after
  calling reset_wcs. [#120]

- Added ``CoordinateHelper.get_axislabel()`` method. [#122]

- Improve handling of get/set_xlabel/ylabel. [#94, #126]

- Added new feature to display world coordinates in interactive mode. [#112]

- When converting celestial coordinates, don't compute a new representation if
  the old one was already spherical or unit spherical. [#125]

- WCS instances from ``wcsaxes.wcs_wrapper`` can now be used to instantiate
  Matplotlib plots with the ``projection=`` keyword. [#136]

- Example datasets are now downloaded from http://data.astropy.org. [#144]

- Registering new frame identifiers is now done with the same API as in
  astropy.wcs.utils, and is handled in wcsaxes.wcs_utils. [#145]

0.2 (2014-08-11)
----------------

### New features

- Added option to specify whether overlapping ticks should be displayed. [#85]

- Added minor ticks. [#89]

- Added option to set frame linewidth and color. [#88]

- Added option to set separators for angle coordinates. [#90]

- Added Python ``%`` format to set format for scalar coordinates. [#98]

### Improvements

- Improved performance of grid drawing. [#100]

- Code is now natively compatible with Python 3.

### Bug Fixes

- Fix axis labels overlapping with ticks. [#96]

- Fix drawing grids multiple times for multi-dimensional data. [#99]

- Fixed bug that occurred when providing a custom ``transData``.

0.1 (2014-07-04)
----------------

- Initial Release.
