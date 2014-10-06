0.3 (unreleased)
----------------

- Fixed a bug that caused axis labels to be incorrectly placed when hiding tick 
  labels using ``set_visible(False)``. [#111]

- Fixed a bug that caused tick labels to not always appear. [#110]

- Fixed a bug that caused the clip path for images to not be updated when
  making subsequent changes after drawing the figure. [#117]

- Fixed a bug that caused the clip path for images to not be updated after
  calling reset_wcs. [#120]

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
