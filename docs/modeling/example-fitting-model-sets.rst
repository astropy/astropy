Fitting Model Sets
==================

Astropy model sets let you fit the same (linear) model to lots of independent
data sets. It solves the linear equations simultaneously, so can avoid looping.
But getting the data into the right shape can be a bit tricky. 

The time savings could be worth the effort. In the example below, if we change
the width*height of the data cube to 500*500 it takes 140 ms on a 2015 MacBook Pro
to fit the models using model sets. Doing the same fit by looping over the 500*500 models 
takes 1.5 minutes, more than 600 times slower.

In the example below, we create a 3D data cube where the first dimension is a ramp --
for example as from non-destructive readouts of an IR detector. So each pixel has a 
depth along a time axis, and flux that results a total number of counts that is 
increasing with time. We will be fitting a 1D polynomial vs. time to estimate the 
flux in counts/second (the slope of the fit).

First, import the necessary libraries:

.. doctest-requires:: matplotlib

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy import stats
    >>> from astropy.modeling import models, fitting

We will use just a small image of 3 rows by 4 columns, with a depth of 10 non-destructive reads::

    >>> depth, width, height = 10, 3, 4  # Time is along the depth axis
    >>> t = np.arange(depth, dtype=np.float64)*10.  # e.g. readouts every 10 seconds

The number of counts in neach pixel is flux*time with the addition of some Gaussian noise::

    >>> fluxes = np.arange(1.*width*height).reshape(height, width)
    >>> image = fluxes[np.newaxis, :, :] * t[:, np.newaxis, np.newaxis]
    >>> image += stats.norm.rvs(0., image*0.05, size=image.shape)  # Add noise
    >>> image.shape
    (10, 3, 4)

Create the models and the fitter. We need N=width*height instances of the same linear,
parametric model (model sets currently only work with linear models and fitters)::

    >>> N = width*height 
    >>> line = models.Polynomial1D(degree=1, n_models=N)
    >>> fit = fitting.LinearLSQFitter()
    >>> print("We created %d models" % len(line))
    We created 12 models

We need to get the data to be fit into the right shape. It's not possible to just feed
the 3D data cube. In this case, the time axis can be one dimensional. 
The fluxes have to be organized into an array that is of shape `width*height,depth` --  in 
other words, we are reshaping to flatten last two axes and transposing to put them first::

    >>> pixels = image.reshape((depth, width*height))
    >>> y = pixels.T
    >>> print("x axis is one dimensional: ",t.shape)
    >>> print("y axis is two dimensional, N by len(x): ", y.shape)
    x axis is one dimensional:  (10,)
    y axis is two dimensional, N by len(x):  (12, 10)

Fit the model. It does the looping over the N models implicitly::

    >>> new_model = fit(line, x=t, y=y)
    >>> print("We fit %d models" % len(new_model))
    We fit 12 models

Fill an array with values computed from the best fit and reshape it to match the original::

    >>> best_fit = new_model(t, model_set_axis=False).T.reshape((depth, height, width))
    print("We reshaped the best fit to dimensions: ", best_fit.shape)
    We reshaped the best fit to dimensions:  (10, 3, 4)

Now inspect the model::

    >>> print(new_model) # doctest: +FLOAT_CMP
    >>> print("The new_model has a param_sets attribute with shape: ",new_model.param_sets.shape)
    >>> print("And values that are the best-fit parameters for each pixel: ")
    >>> print(new_model.param_sets) # doctest: +FLOAT_CMP
    Model: Polynomial1D
    Inputs: ('x',)
    Outputs: ('y',)
    Model set size: 12
    Degree: 1
    Parameters:
                 c0                 c1
        ------------------- ------------------
                        0.0                0.0
         0.3027502259566667  1.011977448587174
        -0.2840478955871922  2.040893880272837
         -2.406760545259093 3.1000703825064604
        0.22342137772278226 3.8945121037774144
          1.651496741874075   4.98888927032699
         3.3670972127381056  6.058033100265409
         1.7403959255727979  7.071466847238193
         2.7731828744191764 7.6671436311988295
         -5.178334595056723  9.116712142530025
         -6.625951349481062 10.128086437447301
        -6.3133581154839264 11.434712892755664
    The new_model has a param_sets attribute with shape:  (2, 12)
    And values that are the best-fit parameters for each pixel:
    [[ 0.          0.30275023 -0.2840479  -2.40676055  0.22342138  1.65149674
       3.36709721  1.74039593  2.77318287 -5.1783346  -6.62595135 -6.31335812]
     [ 0.          1.01197745  2.04089388  3.10007038  3.8945121   4.98888927
       6.0580331   7.07146685  7.66714363  9.11671214 10.12808644 11.43471289]]

Plot the fit along a couple of pixels:

.. doctest-requires:: matplotlib

    >>> def plotramp(t, image, best_fit, row, col):
    >>>     plt.plot(t, image[:, row, col], '.', label='data pixel %d,%d' % (row, col))
    >>>     plt.plot(t, best_fit[:, row, col], '-', label='fit to pixel %d,%d' % (row, col))
    >>>     plt.xlabel('Time')
    >>>     plt.ylabel('Counts')
    >>>     plt.legend(loc='upper left')
    >>> plt.figure(figsize=(10, 5))
    >>> plotramp(t, image, best_fit, 1, 1)
    >>> plotramp(t, image, best_fit, 3, 2)

The data and the best fit model are shown together on one plot.

.. plot::
    
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats
    from astropy.modeling import models, fitting
    
    # Set up the shape of the image and create the time axis
    depth,width,height=10,3,4 # Time is along the depth axis
    t = np.arange(depth, dtype=np.float64)*10.  # e.g. readouts every 10 seconds
    
    # Make up a flux in each pixel 
    fluxes = np.arange(1.*width*height).reshape(height, width)
    # Create the ramps by integrating the fluxes along the time steps
    image = fluxes[np.newaxis, :, :] * t[:, np.newaxis, np.newaxis]
    # Add some Gaussian noise to each sample
    image += stats.norm.rvs(0., image*0.05, size=image.shape)  # Add noise
    
    # Create the models and the fitter
    N = width*height # This is how many instances we need
    line = models.Polynomial1D(degree=1, n_models=N)
    fit = fitting.LinearLSQFitter()
    
    # We need to get the data to be fit into the right shape
    # In this case, the time axis can be one dimensional.
    # The fluxes have to be organized into an array 
    # that is of shape `(width*height, depth)`
    # i.e we are reshaping to flatten last two axes and
    # transposing to put them first.
    pixels = image.reshape((depth, width*height))
    y = pixels.T
    
    # Fit the model. It does the looping over the N models implicitly
    new_model = fit(line, x=t, y=y)
    
    # Fill an array with values computed from the best fit and reshape it to match the original
    best_fit = new_model(t, model_set_axis=False).T.reshape((depth, height, width))
    

    # Plot the fit along a couple of pixels
    def plotramp(t, image, best_fit, row, col):
        plt.plot(t, image[:, row, col], '.', label='data pixel %d,%d' % (row, col))
        plt.plot(t, best_fit[:, row, col], '-', label='fit to pixel %d,%d' % (row, col))
        plt.xlabel('Time')
        plt.ylabel('Counts')
        plt.legend(loc='upper left')    


    plt.figure(figsize=(10, 5))
    plotramp(t, image, best_fit, 1, 1)
    plotramp(t, image, best_fit, 3, 2)
    plt.show()

