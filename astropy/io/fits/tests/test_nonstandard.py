from astropy.io import fits
from . import FitsTestCase


class TestNonstandardHdus(FitsTestCase):
    def test_create_fitshdu(self):
        """
        A round trip test of creating a FitsHDU, adding a FITS file to it,
        writing the FitsHDU out as part of a new FITS file, and then reading
        it and recovering the original FITS file.
        """

        self._test_create_fitshdu(compression=False)

    def test_create_fitshdu_with_compression(self):
        """Same as test_create_fitshdu but with gzip compression enabled."""

        self._test_create_fitshdu(compression=True)

    def _test_create_fitshdu(self, compression=False):
        hdul_orig = fits.open(self.data('test0.fits'),
                              do_not_scale_image_data=True)

        fitshdu = fits.FitsHDU.fromhdulist(hdul_orig, compress=compression)
        # Just to be meta, let's append to the same hdulist that the fitshdu
        # encapuslates
        hdul_orig.append(fitshdu)
        hdul_orig.writeto(self.temp('tmp.fits'), clobber=True)
        del hdul_orig[-1]

        hdul = fits.open(self.temp('tmp.fits'))
        assert isinstance(hdul[-1], fits.FitsHDU)

        wrapped = hdul[-1].hdulist
        assert isinstance(wrapped, fits.HDUList)

        assert hdul_orig.info(output=False) == wrapped.info(output=False)
        assert (hdul[1].data == wrapped[1].data).all()
        assert (hdul[2].data == wrapped[2].data).all()
        assert (hdul[3].data == wrapped[3].data).all()
        assert (hdul[4].data == wrapped[4].data).all()
