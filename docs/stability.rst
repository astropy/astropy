******************************
Current status of sub-packages
******************************

Astropy has benefited from the addition of widely tested legacy code, as well
as new development, resulting in variations in stability across
sub-packages. This document summarizes the current status of the Astropy
sub-packages, so that users understand where they might expect changes in
future, and which sub-packages they can safely use for production code.

Note that until version 1.0, even sub-packages considered *Mature* could
undergo some user interface changes as we work to integrate the packages
better. Thus, we cannot guarantee complete backward-compatibility between
versions at this stage.

.. |planned| image:: _static/planned.png

.. |dev| image:: _static/dev.png

.. |stable| image:: _static/stable.png

.. |mature| image:: _static/mature.png

The classification is as follows:

.. raw:: html

    <table align='center'>
      <tr>
        <td align='center'><img src='_images/planned.png'></td>
        <td>Planned</td>
      </tr>
      <tr>
        <td align='center'><img src='_images/dev.png'></td>
        <td>Actively developed, be prepared for possible significant changes</td>
      </tr>
      <tr>
        <td align='center'><img src='_images/stable.png'></td>
        <td>Reasonably stable, no major changes likely</td>
      </tr>
      <tr>
        <td align='center'><img src='_images/mature.png'></td>
        <td>Mature</td>
      </tr>
    </table>

The current planned and existing sub-packages are:

.. raw:: html

    <table border="1" class="docutils" align='center'>
        <tr>
            <th class="head">
                Sub-Package
            </th>
            <th class="head">
                &nbsp;
            </th>
            <th class="head">
                Comments
            </th>
        </tr>
        <tr>
            <td>
                astropy.analytic_functions
            </td>
            <td align='center'>
                <img alt="dev" src="_images/dev.png">
            </td>
            <td>
                New in v1.0.
            </td>
        </tr>
        <tr>
            <td>
                astropy.constants
            </td>
            <td align='center'>
                <img alt="dev" src="_images/dev.png">
            </td>
            <td>
                Constants have been changed to <tt class="docutils literal"><span class="pre">Quantity</span></tt> objects in v0.2.
            </td>
        </tr>
        <tr>
            <td>
                astropy.convolution
            </td>
            <td align='center'>
                <img alt="dev" src="_images/dev.png">
            </td>
            <td>
                New top-level package in v0.3 (was previously part of
                <tt class="docutils literal"><span class="pre">astropy.nddata</span></tt>).
                No major changes in v0.4.
            </td>
        </tr>
        <tr>
            <td>
                astropy.coordinates
            </td>
            <td align='center'>
                <img alt="dev" src="_images/stable.png">
            </td>
            <td>
                New in v0.2, major changes in v0.4.  Subsequent versions should
                maintain a stable/backwards-compatible API.
            </td>
        </tr>
        <tr>
            <td>
                astropy.cosmology
            </td>
            <td align='center'>
                <img alt="stable" src="_images/stable.png">
            </td>
            <td>
                Incremental improvements since v0.1, but mostly stable API.
                Pure functional interface deprecated in v0.4.
            </td>
        </tr>
        <tr>
            <td>
                astropy.io.ascii
            </td>
            <td align='center'>
                <img alt="mature" src="_images/mature.png">
            </td>
            <td>
                Originally developed as <tt class="docutils literal"><span class="pre">asciitable</span></tt>, and has maintained a stable API.
            </td>
        </tr>
        <tr>
            <td>
                astropy.io.fits
            </td>
            <td align='center'>
                <img alt="mature" src="_images/mature.png">
            </td>
            <td>
                Originally developed as <tt class="docutils literal"><span class="pre">pyfits</span></tt>, and retains an API consistent with the standalone version.
            </td>
        </tr>
        <tr>
            <td>
                astropy.io.misc
            </td>
            <td align='center'>
                <img alt="mature" src="_images/dev.png">
            </td>
            <td>
                 The functionality that is currently present is stable, but this sub-package will likely see major additions in future.
            </td>
        </tr>
        <tr>
            <td>
                astropy.io.votable
            </td>
            <td align='center'>
                <img alt="mature" src="_images/mature.png">
            </td>
            <td>
                Originally developed as <tt class="docutils literal"><span class="pre">vo.table</span></tt>, and has a stable API.
            </td>
        </tr>
        <tr>
            <td>
                astropy.image
            </td>
            <td align='center'>
                <img alt="dev" src="_images/dev.png">
            </td>
            <td>
                New in v1.0, and in development.
            </td>
        </tr>
        <tr>
            <td>
                astropy.modeling
            </td>
            <td align='center'>
                <img alt="dev" src="_images/dev.png">
            </td>
            <td>
                New in v0.3
            </td>
        </tr>
        <tr>
            <td>
                astropy.nddata
            </td>
            <td align='center'>
                <img alt="dev" src="_images/dev.png">
            </td>
            <td>
                In development, and does not yet contain much functionality apart from a base class for N-dimensional datasets.
            </td>
        </tr>
        <tr>
            <td>
                astropy.photometry
            </td>
            <td align='center'>
                <img alt="planned" src="_images/planned.png">
            </td>
            <td>
                &nbsp;
            </td>
        </tr>
        <tr>
            <td>
                astropy.stats
            </td>
            <td align='center'>
                <img alt="dev" src="_images/dev.png">
            </td>
            <td>
                Still in development, and does not yet contain much functionality.
            </td>
        </tr>
        <tr>
            <td>
                astropy.table
            </td>
            <td align='center'>
                <img alt="stable" src="_images/stable.png">
            </td>
            <td>
                Incremental improvements since v0.1, but mostly stable API.
            </td>
        </tr>
        <tr>
            <td>
                astropy.time
            </td>
            <td align='center'>
                <img alt="stable" src="_images/stable.png">
            </td>
            <td>
                Incremental improvements since v0.1, but mostly stable API.
            </td>
        </tr>
        <tr>
            <td>
                astropy.units
            </td>
            <td align='center'>
                <img alt="stable" src="_images/stable.png">
            </td>
            <td>
                New in v0.2. Adapted from <tt class="docutils literal"><span class="pre">pnbody</span></tt> and integrated into Astropy.
            </td>
        </tr>
        <tr>
            <td>
                astropy.utils
            </td>
            <td align='center'>
                <img alt="dev" src="_images/dev.png">
            </td>
            <td>
                This sub-package contains mostly utilities destined for use in other parts of Astropy, and is not yet stable.
            </td>
        </tr>
        <tr>
            <td>
                astropy.vo
            </td>
            <td align='center'>
                <img alt="dev" src="_images/dev.png">
            </td>
            <td>
                Virtual Observatory service access and validation. Currently, only Simple Cone Search and SAMP are supported.
            </td>
        </tr>
        <tr>
            <td>
                astropy.wcs
            </td>
            <td align='center'>
                <img alt="stable" src="_images/stable.png">
            </td>
            <td>
                Originally developed as <tt class="docutils literal"><span class="pre">pywcs</span></tt>, and has a stable API for now. However, there are plans to generalize the WCS interface to accommodate non-FITS WCS transformations, and this may lead to small changes in the user interface.
            </td>
        </tr>
    </table>

