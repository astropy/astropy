.. _example-fitting-model-sets:

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
flux in counts/second (the slope of the fit). We will use just a small image
of 3 rows by 4 columns, with a depth of 10 non-destructive reads.

First, import the necessary libraries:

    >>> import numpy as np
    >>> np.random.seed(seed=12345)
    >>> from astropy.modeling import models, fitting

    >>> depth, width, height = 10, 3, 4  # Time is along the depth axis
    >>> t = np.arange(depth, dtype=np.float64)*10.  # e.g. readouts every 10 seconds

The number of counts in neach pixel is flux*time with the addition of some Gaussian noise::

    >>> fluxes = np.arange(1. * width * height).reshape(width, height)
    >>> image = fluxes[np.newaxis, :, :] * t[:, np.newaxis, np.newaxis]
    >>> image += np.random.normal(0., image*0.05, size=image.shape)  # Add noise
    >>> image.shape
    (10, 3, 4)

Create the models and the fitter. We need N=width*height instances of the same linear,
parametric model (model sets currently only work with linear models and fitters)::

    >>> N = width * height
    >>> line = models.Polynomial1D(degree=1, n_models=N)
    >>> fit = fitting.LinearLSQFitter()
    >>> print(f"We created {len(line)} models")
    We created 12 models

We need to get the data to be fit into the right shape. It's not possible to just feed
the 3D data cube. In this case, the time axis can be one dimensional.
The fluxes have to be organized into an array that is of shape ``width*height,depth`` --  in
other words, we are reshaping to flatten last two axes and transposing to put them first::

    >>> pixels = image.reshape((depth, width*height))
    >>> y = pixels.T
    >>> print("x axis is one dimensional: ",t.shape)
    x axis is one dimensional:  (10,)
    >>> print("y axis is two dimensional, N by len(x): ", y.shape)
    y axis is two dimensional, N by len(x):  (12, 10)

Fit the model. It fits the N models simultaneously::

    >>> new_model = fit(line, x=t, y=y)
    >>> print(f"We fit {len(new_model)} models")
    We fit 12 models

Fill an array with values computed from the best fit and reshape it to match the original::

    >>> best_fit = new_model(t, model_set_axis=False).T.reshape((depth, height, width))
    >>> print("We reshaped the best fit to dimensions: ", best_fit.shape)
    We reshaped the best fit to dimensions:  (10, 4, 3)

Now inspect the model::

    >>> print(new_model) # doctest: +FLOAT_CMP
    Model: Polynomial1D
    Inputs: ('x',)
    Outputs: ('y',)
    Model set size: 12
    Degree: 1
    Parameters:
                 c0                 c1
        ------------------- ------------------
	                0.0                0.0
	-0.5206606340901005 1.0463998276552442
         0.6401930368329991 1.9818733492667582
         0.1134712985541639  3.049279878262541
        -3.3556420351251313  4.013810434122983
          6.782223372575449  4.755912707001437
          3.628220497058842  5.841397947835126
        -5.8828309622531565  7.016044775363114
        -11.676538736037775  8.072519832452022
          -6.17932185981594  9.103924115403503
        -4.7258541419613165 10.315295021908833
           4.95631951675311 10.911167956770575

    >>> print("The new_model has a param_sets attribute with shape: ",new_model.param_sets.shape)
    The new_model has a param_sets attribute with shape:  (2, 12)

    >>> print(f"And values that are the best-fit parameters for each pixel:\n{new_model.param_sets}") # doctest: +FLOAT_CMP
    And values that are the best-fit parameters for each pixel:
    [[  0.          -0.52066063   0.64019304   0.1134713   -3.35564204
        6.78222337   3.6282205   -5.88283096 -11.67653874  -6.17932186
       -4.72585414   4.95631952]
     [  0.           1.04639983   1.98187335   3.04927988   4.01381043
        4.75591271   5.84139795   7.01604478   8.07251983   9.10392412
       10.31529502  10.91116796]]

Plot the fit along a couple of pixels:

    >>> def plotramp(t, image, best_fit, row, col):
    ...     plt.plot(t, image[:, row, col], '.', label=f'data pixel {row},{col}')
    ...     plt.plot(t, best_fit[:, row, col], '-', label=f'fit to pixel {row},{col}')
    ...     plt.xlabel('Time')
    ...     plt.ylabel('Counts')
    ...     plt.legend(loc='upper left')
    >>> fig = plt.figure(figsize=(10, 5)) # doctest: +SKIP
    >>> plotramp(t, image, best_fit, 1, 1) # doctest: +SKIP
    >>> plotramp(t, image, best_fit, 2, 1) # doctest: +SKIP

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
    N = width * height # This is how many instances we need
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
        plt.plot(t, image[:, row, col], '.', label=f'data pixel {row},{col}')
        plt.plot(t, best_fit[:, row, col], '-', label=f'fit to pixel {row},{col}')
        plt.xlabel('Time')
        plt.ylabel('Counts')
        plt.legend(loc='upper left')


    plt.figure(figsize=(10, 5))
    plotramp(t, image, best_fit, 1, 1)
    plotramp(t, image, best_fit, 3, 2)
    plt.show()
