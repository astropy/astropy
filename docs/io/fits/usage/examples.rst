Examples
--------

Converting a 3-color image (JPG) to separate FITS images
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: ../images/Hs-2009-14-a-web.jpg
   :scale: 100 %
   :align: center
   :alt: Starting image

.. container:: figures

    .. figure:: ../images/Red.jpg
       :scale: 50
       :alt: Red color information

       Red color information

    .. figure:: ../images/Green.jpg
       :scale: 50
       :alt: Green color information

       Green color information

    .. figure:: ../images/Blue.jpg
       :scale: 50
       :alt: Blue color information

       Blue color information

::

    #!/usr/bin/env python
    import numpy 
    import Image

    from astropy.io import fits

    #get the image and color information
    image = Image.open('hs-2009-14-a-web.jpg')
    #image.show()
    xsize, ysize = image.size
    r, g, b = image.split()
    rdata = r.getdata() # data is now an array of length ysize\*xsize
    gdata = g.getdata()
    bdata = b.getdata()

    # create numpy arrays
    npr = numpy.reshape(rdata, (ysize, xsize))
    npg = numpy.reshape(gdata, (ysize, xsize))
    npb = numpy.reshape(bdata, (ysize, xsize))

    # write out the fits images, the data numbers are still JUST the RGB
    # scalings; don't use for science
    red = fits.PrimaryHDU(data=npr)
    red.header['LATOBS'] = "32:11:56" # add spurious header info
    red.header['LONGOBS'] = "110:56"
    red.writeto('red.fits')

    green = fits.PrimaryHDU(data=npg)
    green.header['LATOBS'] = "32:11:56"
    green.header['LONGOBS'] = "110:56"
    green.writeto('green.fits')

    blue = fits.PrimaryHDU(data=npb)
    blue.header['LATOBS'] = "32:11:56"
    blue.header['LONGOBS'] = "110:56"
    blue.writeto('blue.fits')
