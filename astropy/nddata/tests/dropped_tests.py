def test_to_hdu_long_metadata_item(ccd_data):
    # There is no attempt to try to handle the general problem of
    # a long keyword (that requires HIERARCH) with a long string value
    # (that requires CONTINUE).
    # However, a long-ish keyword with a long value can happen because of
    # auto-logging, and we are supposed to handle that.

    # So, a nice long command:
    from ..core import subtract_dark, _short_names

    dark = CCDData(np.zeros_like(ccd_data.data), unit="adu")
    result = subtract_dark(ccd_data, dark, dark_exposure=30 * u.second,
                           data_exposure=15 * u.second, scale=True)
    assert 'subtract_dark' in result.header
    hdulist = result.to_hdu()
    header = hdulist[0].header
    assert header['subtract_dark'] == _short_names['subtract_dark']
    args_value = header[_short_names['subtract_dark']]
    # Yuck -- have to hand code the ".0" to the numbers to get this to pass...
    assert "dark_exposure={0} {1}".format(30.0, u.second) in args_value
    assert "data_exposure={0} {1}".format(15.0, u.second) in args_value
    assert "scale=True" in args_value


def test_ccddata_header_does_not_corrupt_fits(ccd_data, tmpdir):
    # This test is for the problem described in astropy/ccdproc#165
    # The issue comes up when a long FITS keyword value is in a header
    # that is read in and then converted to a non-fits.Header object
    # that is dict-like, and then you try to write that out again as
    # FITS. Certainly FITS files out to be able to round-trip, and
    # this test checks for that.

    fake_dark = ccd_data.copy()
    # This generates a nice long log entry in the header.
    ccd = subtract_dark(ccd_data, fake_dark, dark_exposure=30*u.second,
                        data_exposure=30*u.second)
    # The write below succeeds...
    long_key = tmpdir.join('long_key.fit').strpath
    ccd.write(long_key)

    # And this read succeeds...
    ccd_read = CCDData.read(long_key, unit="adu")

    # This write failed in astropy/ccdproc#165 but should not:
    rewritten = tmpdir.join('should_work.fit').strpath
    ccd_read.write(rewritten)

    # If all is well then reading the file we just wrote should result in an
    # identical header.
    ccd_reread = CCDData.read(rewritten, unit="adu")
    assert ccd_reread.header == ccd_read.header


def test_ccddata_with_fits_header_as_meta_works_with_autologging(ccd_data,
                                                                 tmpdir):
    tmp_file = tmpdir.join('tmp.fits')
    hdr = fits.Header(ccd_data.header)
    ccd_data.header = hdr
    fake_dark = ccd_data.copy()
    # The combination below will generate a long keyword ('subtract_dark')
    # and a long value (the function signature) in autlogging.
    ccd2 = subtract_dark(ccd_data, fake_dark,
                         dark_exposure=30*u.second,
                         data_exposure=15*u.second,
                         scale=True)
    # This should not fail....
    ccd2.write(tmp_file.strpath)
    # And the header on ccd2 should be a subset of the written header; they
    # do not match exactly because the written header contains information
    # about the array size that is the hdr we created manually.
    ccd2_read = CCDData.read(tmp_file.strpath, unit=u.adu)
    for k, v in six.iteritems(ccd2.header):
        assert ccd2_read.header[k] == v
