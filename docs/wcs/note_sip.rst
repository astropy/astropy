.. include:: references.rst
.. doctest-skip-all
.. _note_sip: :orphan:


Note about SIP and WCS
**********************

`astropy.wcs` supports the Simple Imaging Polynomial (`SIP`_) convention.
The SIP distortion is defined in FITS headers by the presence of the
SIP specific keywords **and** a ``-SIP`` suffix in ``CTYPE``, for example
``RA---TAN-SIP``, ``DEC--TAN-SIP``.

This has not been a strict convention in the past and the default in
`astropy.wcs` is to always include the SIP distortion if the SIP coefficients
are present, even if ``-SIP`` is not included in CTYPE.
The presence of a ``-SIP`` suffix in CTYPE is not used as a trigger
to initialize the SIP distortion.

It is important that headers implement correctly the SIP convention.
If the intention is to use the SIP distortion, a header should have
the SIP coefficients and the ``-SIP`` suffix in CTYPE.

`astropy.wcs` prints INFO messages when inconsistent headers are detected,
for example when SIP coefficients are present but CTYPE is missing a ``-SIP`` suffix,
see examples below.
`astropy.wcs` will print a message about the inconsistent header
but will create and use the SIP distortion and it will be used in
calls to `~astropy.wcs.wcs.WCS.all_pix2world`. If this was not the intended use
(e.g. it's a drizzled image and has no distortions) it is best to remove the SIP
coefficients from the header. They can be removed temporarily from a WCS object by

>>> wcsobj.sip = None

In addition, if SIP is the only distortion in the header, the two methods,
`~astropy.wcs.wcs.WCS.wcs_pix2world` and `~astropy.wcs.wcs.WCS.wcs_world2pix`,
may be used to transform from pixels to world coordinate system while omitting distortions.

Another consequence of the inconsistent header is that if
`~astropy.wcs.wcs.WCS.to_header()` is called with ``relax=True`` it will return a header
with SIP coefficients and a ``-SIP`` suffix in CTYPE and will not reproduce the original header.

**In conclusion, when astropy.wcs detects inconsistent headers, the recommendation
is that the header is inspected and corrected to match the data.**

Below is an example of a header with SIP coefficients when ``-SIP`` is missing from CTYPE.
The data is drizzled, i.e. distortion free, so the intention is **not** to include the
SIP distortion.

>>> wcsobj = wcs.WCS(header)

INFO::

        Inconsistent SIP distortion information is present in the FITS header and the WCS object:
        SIP coefficients were detected, but CTYPE is missing a "-SIP" suffix.
        astropy.wcs is using the SIP distortion coefficients,
        therefore the coordinates calculated here might be incorrect.

        If you do not want to apply the SIP distortion coefficients,
        please remove the SIP coefficients from the FITS header or the
        WCS object.  As an example, if the image is already distortion-corrected
        (e.g., drizzled) then distortion components should not apply and the SIP
        coefficients should be removed.

        While the SIP distortion coefficients are being applied here, if that was indeed the intent,
        for consistency please append "-SIP" to the CTYPE in the FITS header or the WCS object.


>>> hdr = wcsobj.to_header(relax=True)

INFO::

        Inconsistent SIP distortion information is present in the current WCS:
        SIP coefficients were detected, but CTYPE is missing "-SIP" suffix,
        therefore the current WCS is internally inconsistent.

        Because relax has been set to True, the resulting output WCS will have
        "-SIP" appended to CTYPE in order to make the header internally consistent.

        However, this may produce incorrect astrometry in the output WCS, if
        in fact the current WCS is already distortion-corrected.

        Therefore, if current WCS is already distortion-corrected (eg, drizzled)
        then SIP distortion components should not apply. In that case, for a WCS
        that is already distortion-corrected, please remove the SIP coefficients
        from the header.
