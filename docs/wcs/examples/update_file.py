# Download image from tutorials folder, load the WCS information from the fits header, trim the image using the cutout2D
# function. Finally, replace the original image and WCS information from with the data from the cutout2D object.
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.data import download_file
from astropy.nddata import Cutout2D
from astropy import coordinates, units as u

def download_image_update_from_cutout(file_address, size):

    # Download the image
    image_file = download_file(file_address, cache=False)

    # Load the image and the world coordinates from the header
    hdu = fits.open(image_file)[0]
    wcs = WCS(hdu.header)

    # Get Horsehead Nebula coordinates
    horsehead_coord = coordinates.SkyCoord.from_name('Horsehead Nebula')

    # Make the cutout using the wcs
    cutout = Cutout2D(hdu.data, position=horsehead_coord, size=size, wcs=wcs)

    # Update image data
    hdu.data = cutout.data

    # Update the WCS from the cutout
    hdu.header.update(cutout.wcs.to_header())

    # Replace original image by the new "trimmed" version
    hdu.writeto(image_file, overwrite=True)

    # Load the new image
    hdu = fits.open(image_file)[0]
    wcs = WCS(hdu.header)

    # Plot the image
    fig = plt.figure()
    fig.add_subplot(111, projection=wcs)
    plt.imshow(wcs.data, origin='lower', cmap=plt.cm.viridis)
    plt.xlabel('RA')
    plt.ylabel('Dec')
    plt.show()

if __name__ == '__main__':

    # Horsehead Nebula image from the tutorials folder
    fits_file = 'http://data.astropy.org/tutorials/FITS-images/HorseHead.fits'

    # Defining area around the horsehead nebular "head"
    sizeTrim = (400 * u.pixel, 400 * u.pixel)

    #Run the function
    download_image_update_from_cutout(fits_file, sizeTrim)
