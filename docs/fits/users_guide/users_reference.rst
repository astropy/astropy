****************
Reference Manual
****************

**Examples**


Converting a 3-color image (JPG) to separate FITS images
========================================================


.. figure:: ../_static/Hs-2009-14-a-web.jpg
   :scale: 100 %
   :align: center
   :alt: Starting image

.. container:: figures

    .. figure:: ../_static/Red.jpg
       :target: ../_static/Red.jpg
       :scale: 50
       :alt: Red color information

       Red color information

    .. figure:: ../_static/Green.jpg
       :target: ../_static/Green.jpg
       :scale: 50
       :alt: Green color information

       Green color information

    .. figure:: ../_static/Blue.jpg
       :target: ../_static/Blue.jpg
       :scale: 50
       :alt: Blue color information

       Blue color information

::

    #!/usr/bin/env python
    import pyfits
    import numpy 
    import Image

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
    red = pyfits.PrimaryHDU()
    red.header.update('LATOBS', "32:11:56") # add spurious header info
    red.header.update('LONGOBS', "110:56")
    red.data = npr
    red.writeto('red.fits')
    green = pyfits.PrimaryHDU()
    green.header.update('LATOBS', "32:11:56")
    green.header.update('LONGOBS', "110:56")
    green.data = npg
    green.writeto('green.fits')
    blue = pyfits.PrimaryHDU()
    blue.header.update('LATOBS', "32:11:56")
    blue.header.update('LONGOBS', "110:56")
    blue.data = npb
    blue.writeto('blue.fits')
