# Download an example FITS file, create a 2D cutout, and save it to a
# new FITS file, including the updated cutout WCS.
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.utils.data import download_file
from astropy.wcs import WCS


def download_image_save_cutout(url, position, size):
    # Download the image
    filename = download_file(url)

    # Load the image and the WCS
    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)

    # Make the cutout, including the WCS
    cutout = Cutout2D(hdu.data, position=position, size=size, wcs=wcs)

    # Put the cutout image in the FITS HDU
    hdu.data = cutout.data

    # Update the FITS header with the cutout WCS
    hdu.header.update(cutout.wcs.to_header())

    # Write the cutout to a new FITS file
    cutout_filename = 'example_cutout.fits'
    hdu.writeto(cutout_filename, overwrite=True)


if __name__ == '__main__':
    url = 'https://astropy.stsci.edu/data/photometry/spitzer_example_image.fits'

    position = (500, 300)
    size = (400, 400)
    download_image_save_cutout(url, position, size)
