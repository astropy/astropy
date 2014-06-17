from astropy.utils.data import download_file


def msx():
    return download_file("http://astrofrog.github.io/wcsaxes-datasets/msx.fits", cache=True)


def rosat():
    return download_file("http://astrofrog.github.io/wcsaxes-datasets/rosat.fits", cache=True)


def twoMASS_k():
    return download_file("http://astrofrog.github.io/wcsaxes-datasets/2MASS_k.fits", cache=True)


def cube():
    return download_file("http://astrofrog.github.io/wcsaxes-datasets/L1448_13CO.fits", cache=True)


# image1 = download_file("http://astrofrog.github.io/wcsaxes-datasets/msx.fits", cache=True)
# image2 = download_file("http://astrofrog.github.io/wcsaxes-datasets/rosat.fits", cache=True)
# image3 = download_file("http://astrofrog.github.io/wcsaxes-datasets/2MASS_k.fits", cache=True)
# data_cube = download_file("http://astrofrog.github.io/wcsaxes-datasets/L1448_13CO.fits", cache=True)
