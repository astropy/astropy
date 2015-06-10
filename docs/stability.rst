******************************
Current status of sub-packages
******************************

Astropy has benefited from the addition of widely tested legacy code, as well
as new development, resulting in variations in stability across
sub-packages. This document summarizes the current status of the Astropy
sub-packages, so that users understand where they might expect changes in
future, and which sub-packages they can safely use for production code.

The classification is as follows:

.. raw:: html

    <style>
         .planned:before {
              color: #cbcbcb;
              content: "⬤";
         }
         .dev:before {
              color: #ffad00;
              content: "⬤";
         }
         .stable:before {
              color: #4e72c3;
              content: "⬤";
         }
         .mature:before {
              color: #03a913;
              content: "⬤";
         }
    </style>

    <table align='center'>
      <tr>
        <td align='center'><span class="planned"></span></td>
        <td>Planned</td>
      </tr>
      <tr>
        <td align='center'><span class="dev"></span></td>
        <td>Actively developed, be prepared for possible significant changes.</td>
      </tr>
      <tr>
        <td align='center'><span class="stable"></span></td>
        <td>Reasonably stable, any significant changes/additions will generally include backwards-compatiblity.</td>
      </tr>
      <tr>
        <td align='center'><span class="mature"></span></td>
        <td>Mature.  Additions/improvements possible, but no major changes planned. </td>
      </tr>
    </table>

The current planned and existing sub-packages are:

.. raw:: html

    <table border="1" class="docutils stability" align='center'>
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
                <span class="dev"></span>
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
                <span class="stable"></span>
            </td>
            <td>
                Constants were changed to <tt class="docutils literal"><span class="pre">Quantity</span></tt> objects in v0.2. Since then on, the package has been stable, with occasional additions of new constants.
            </td>
        </tr>
        <tr>
            <td>
                astropy.convolution
            </td>
            <td align='center'>
                <span class="stable"></span>
            </td>
            <td>
                New top-level package in v0.3 (was previously part of
                <tt class="docutils literal"><span class="pre">astropy.nddata</span></tt>).
                No major changes since, likely will maintain backwards compatibility but possible future additions or improvements.
            </td>
        </tr>
        <tr>
            <td>
                astropy.coordinates
            </td>
            <td align='center'>
                <span class="stable"></span>
            </td>
            <td>
                New in v0.2, major changes in v0.4.  Subsequent versions should
                maintain a stable/backwards-compatible API, following the plan of <a href="https://github.com/astropy/astropy-APEs/blob/master/APE5.rst">APE 5</a>.  Further major additions/enhancements likely, but with basic framework unchanged.
            </td>
        </tr>
        <tr>
            <td>
                astropy.cosmology
            </td>
            <td align='center'>
                <span class="stable"></span>
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
                <span class="mature"></span>
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
                <span class="mature"></span>
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
                <span class="mature"></span>
            </td>
            <td>
                Originally developed as <tt class="docutils literal"><span class="pre">vo.table</span></tt>, and has a stable API.
            </td>
        </tr>
        <tr>
            <td>
                astropy.modeling
            </td>
            <td align='center'>
                <span class="dev"></span>
            </td>
            <td>
                New in v0.3.  Major changes in v1.0, significant additions planned.  Backwards-compatibility likely to be maintained, but not guaranteed.
            </td>
        </tr>
        <tr>
            <td>
                astropy.nddata
            </td>
            <td align='center'>
                <span class="dev"></span>
            </td>
            <td>
                Significantly revised in v1.0 to implement <a href="https://github.com/astropy/astropy-APEs/blob/master/APE7.rst">APE 7</a>. Major changes in the API are not anticipated, broader use may reveal flaws that require API changes.
            </td>
        </tr>
        <tr>
            <td>
                astropy.photometry
            </td>
            <td align='center'>
                <span class="planned"></span>
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
                <span class="dev"></span>
            </td>
            <td>
                Likely to maintain backwards-compatibility, but functionality continually being expanded, so significant additions likely in the future.
            </td>
        </tr>
        <tr>
            <td>
                astropy.table
            </td>
            <td align='center'>
                <span class="stable"></span>
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
                <span class="mature"></span>
            </td>
            <td>
                Incremental improvements since v0.1, API likely to remain stable
                for the foreseeable future.
            </td>
        </tr>
        <tr>
            <td>
                astropy.units
            </td>
            <td align='center'>
                <span class="stable"></span>
            </td>
            <td>
                New in v0.2. Adapted from <tt class="docutils literal"><span class="pre">pnbody</span></tt> and integrated into Astropy. Current functionality stable with intent to maintain backwards compatibility. Significant new functionality is likely to be added in future versions.
            </td>
        </tr>
        <tr>
            <td>
                astropy.utils
            </td>
            <td align='center'>
                <span class="dev"></span>
            </td>
            <td>
                Contains mostly utilities destined for internal use with other parts of Astropy.  Existing functionality generally stable, but regular additions and occasional changes.
            </td>
        </tr>
        <tr>
            <td>
                astropy.visualization
            </td>
            <td align='center'>
                <span class="dev"></span>
            </td>
            <td>
                New in v1.0, and in development.
            </td>
        </tr>
        <tr>
            <td>
                astropy.vo
            </td>
            <td align='center'>
                <span class="stable"></span>
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
                <span class="stable"></span>
            </td>
            <td>
                Originally developed as <tt class="docutils literal"><span class="pre">pywcs</span></tt>, and has a stable API for now. However, there are plans to generalize the WCS interface to accommodate non-FITS WCS transformations, and this may lead to small changes in the user interface.
            </td>
        </tr>
    </table>
