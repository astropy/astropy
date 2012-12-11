******************************
Current status of sub-packages
******************************

Astropy has benefited from the addition of widely tested legacy code, as well
as new development, resulting in variations in stability accross
sub-packages. This document summarises the current status of the Astropy
sub-packages, so that users understand where they might expect changes in
future, and which sub-packages they can safely use for production code.

.. |planned| image:: _static/planned.png

.. |dev| image:: _static/dev.png

.. |stable| image:: _static/stable.png

.. |mature| image:: _static/mature.png

The classification is as follows:

.. raw:: html

    <table align='center'>
      <tr>
        <td><img src='_images/planned.png'></td>
        <td>Planned</td>
      </tr>
      <tr>
        <td><img src='_images/dev.png'></td>
        <td>Actively developed, be prepared for API changes</td>
      </tr>
      <tr>
        <td><img src='_images/stable.png'></td>
        <td>Reasonably stable API, no major changes likely</td>
      </tr>
      <tr>
        <td><img src='_images/mature.png'></td>
        <td>Mature</td>
      </tr>
    </table>

The current planned and existing sub-packages are:

+----------------------+--------------+----------------------------------------------------------------------------------+
| Sub-Package          | Status       | Comments                                                                         |
+======================+==============+==================================================================================+
|astropy.constants     |  |dev|       | Constants have recently been changed to ``Quantity`` objects.                    |
+----------------------+--------------+----------------------------------------------------------------------------------+
|astropy.coordinates   |  |dev|       | Recently been added and is under very active development, so API changes are     |
|                      |              | possible.                                                                        |
+----------------------+--------------+----------------------------------------------------------------------------------+
|astropy.cosmology     |  |stable|    | Incremental improvements since v0.1, but mostly stable API.                      |
+----------------------+--------------+----------------------------------------------------------------------------------+
|astropy.io.ascii      |  |mature|    | Originally developed as ``asciitable``, and has maintained a stable API.         |
+----------------------+--------------+----------------------------------------------------------------------------------+
|astropy.io.fits       |  |mature|    | Originally developed as ``pyfits``, and retains an API consistent with the       |
|                      |              | standalone version.                                                              |
+----------------------+--------------+----------------------------------------------------------------------------------+
|astropy.io.votable    |  |mature|    | Originally developed as ``vo.table``, and has a stable API.                      |
+----------------------+--------------+----------------------------------------------------------------------------------+
|astropy.photometry    |  |planned|   |                                                                                  |
+----------------------+--------------+----------------------------------------------------------------------------------+
|astropy.stats         |  |dev|       | Still in development, and does not yet contain much functionality.               |
+----------------------+--------------+----------------------------------------------------------------------------------+
|astropy.table         |  |stable|    | Incremental improvements since v0.1, but mostly stable API.                      |
+----------------------+--------------+----------------------------------------------------------------------------------+
|astropy.time          |  |stable|    | Incremental improvements since v0.1, but mostly stable API.                      |
+----------------------+--------------+----------------------------------------------------------------------------------+
|astropy.units         |  |stable|    | Recently heavily adapted from ``pnbody`` and integrated into Astropy. The API    |
|                      |              | is likely stable.                                                                |
+----------------------+--------------+----------------------------------------------------------------------------------+
|astropy.utils         | |dev|        | This sub-package contains mostly utilities destined for use in other parts of    |
|                      |              | Astropy, and is not yet stable.                                                  |
+----------------------+--------------+----------------------------------------------------------------------------------+
|astropy.vo            |  |planned|   |                                                                                  |
+----------------------+--------------+----------------------------------------------------------------------------------+
|astropy.wcs           |  |stable|    | Originally developed as ``pywcs``, and has a stable API for now. However, there  |
|                      |              | are plans to generalize the WCS interface to accommodate non-FITS WCS            |
|                      |              | transformations, and this may lead to small changes in the user interface.       |
+----------------------+--------------+----------------------------------------------------------------------------------+
