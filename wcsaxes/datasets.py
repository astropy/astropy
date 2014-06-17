from astropy.utils.data import download_file

image1 = download_file("http://astrofrog.github.io/wcsaxes-datasets/msx.fits", cache=True)
image2 = download_file("http://astrofrog.github.io/wcsaxes-datasets/rosat.fits", cache=True)
image3 = download_file("http://astrofrog.github.io/wcsaxes-datasets/2MASS_k.fits", cache=True)
data_cube = download_file("http://astrofrog.github.io/wcsaxes-datasets/L1448_13CO.fits", cache=True)
