from astropy.io.fits.hdu.table import BinTableHDU

__all__ = ["_CompBinTableHDU"]


class _CompBinTableHDU(BinTableHDU):
    _load_variable_length_data = False
    """
    We don't want to always load all the tiles so by setting this option
    we can then access the tiles as needed.
    """

    _manages_own_heap = True

    @classmethod
    def match_header(cls, header):
        return False
