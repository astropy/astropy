Fitting Model Sets
==================

Astropy model sets let you fit the same (linear) model to lots of independent
data sets. It does the looping for you, and can be significantly faster than
doing the looping yourself. But getting the data into the right shape can be 
a bit tricky. 

In the example below, we create a 3D data cube where the first dimension is a ramp, 
for example as from non-destructive readouts of an IR detector. So each pixel has a 
depth along a time axis, and flux that results a total number of counts that is 
increasing with time. We will be fitting a 1D polynomial vs. time to estimate the 
flux in counts/second (the slope of the fit).

   .. plot::
       :include-source:

       import numpy as np
       import matplotlib.pyplot as plt
       from scipy import stats
       from astropy.modeling import models, fitting
       
       # Set up the shape of the image and create the time axis
       DEPTH,WIDTH,HEIGHT=10,101,102 # Time is along the DEPTH axis
       t = np.arange(DEPTH,dtype=np.float64)*10. # e.g. readouts every 10 seconds
       
       # Make up a flux in each pixel 
       fluxes = np.arange(1.*WIDTH*HEIGHT).reshape(HEIGHT,WIDTH)
       # Create the ramps by integrating the fluxes along the time steps
       image = np.stack([fluxes*dt for dt in t],axis=0)
       # Add some Gaussian noise to each sample
       image += stats.norm.rvs(0.,image*0.05,size=image.shape) # Add noise
       print("The image shape is: ",image.shape)
       
       # Create the models and the fitter
       N=WIDTH*HEIGHT # This is how many instances we need
       line=models.Polynomial1D(degree=1,n_models=N)
       fit = fitting.LinearLSQFitter()
       print("We created %d models" % len(line))
       
       # We need to get the data to be fit into the right shape
       # In this case, the time axis can be one dimensional.
       # The fluxes have to be organized into an array 
       # that is of shape `(WIDTH*HEIGHT,DEPTH)`
       # i.e we are reshaping to flatten last two axes and 
       # transposing to put them first.
       pixels = image.reshape((DEPTH,WIDTH*HEIGHT))
       y = pixels.T
       print("x axis is one dimensional: ",t.shape)
       print("y axis is two dimensional, N by len(x): ", y.shape)
       
       # Fit the model. It does the looping over the N models implicitly
       new_model = fit(line,x=t,y=y)
       print("We fit %d models: " % len(new_model))
       
       # Fill an array with values computed from the best fit and reshape it to match the original
       best_fit = new_model(t,model_set_axis=False).T.reshape((DEPTH,HEIGHT,WIDTH))
       print("We reshaped the best fit to dimensions: ",best_fit.shape)
       
       # Plot the fit along a couple of pixels
       plt.figure(figsize=(10,5))
       plt.subplot(121)
       X,Y=HEIGHT//5,WIDTH//5
       plt.plot(t,image[:,X,Y],'.')
       plt.plot(t,best_fit[:,X,Y],'-')
       plt.subplot(122)
       X,Y=HEIGHT//2,WIDTH//2
       plt.plot(t,image[:,X,Y],'.')
       plt.plot(t,best_fit[:,X,Y],'-')
       print(new_model.param_sets.shape)
       print(new_model.param_sets)
       print(new_model)
       plt.show()

The printed output of the script is::

    The image shape is:  (10, 102, 101)
    We created 10302 models
    x axis is one dimensional:  (10,)
    y axis is two dimensional, N by len(x):  (10302, 10)
    We fit 10302 models: 
    We reshaped the best fit to dimensions:  (10, 102, 101)

    Model: Polynomial1D
    Inputs: ('x',)
    Outputs: ('y',)
    Model set size: 10302
    Degree: 1
    Parameters:
                 c0                 c1
        ------------------- ------------------
                        0.0                0.0
         1.3093101365940352 0.9725635141702309
        -0.7363862005667433 2.0574413783040253
        -3.2379037283031904  3.103107332439707
         0.1890694350875021 3.9522395588092243
          7.433319082488881  4.739094525515152
          -9.58833113519548  6.205021349646927
         -16.34019588402537  7.718178246509106
         1.3504324911167127  7.850386296532051
           9.66338439344172   8.82219967706662
                        ...                ...
        -11654.336300168248 10703.511450459855
        -1810.7546524482611 10761.276320076036
        -1076.1523247925975 10443.094930390856
        -1657.4021351949136   10300.8517030767
        -7848.9545301854205 10627.101026440047
         6545.5276861968505 10230.805915635534
         -706.4933781219803 10555.927783793251
         -1543.684934784145 10252.098060425496
          1901.922723180854 10365.930306399056
         -6070.755094219937  10249.89171753032
        -11234.448199108616 10559.617582642642
        Length = 10302 rows
    The new_model has a param_sets attribute with shape:  (2, 10302)
    And values that are the best-fit parameters for each pixel:
    [[ 0.00000000e+00  1.30931014e+00 -7.36386201e-01 ...  1.90192272e+03
      -6.07075509e+03 -1.12344482e+04]
     [ 0.00000000e+00  9.72563514e-01  2.05744138e+00 ...  1.03659303e+04
       1.02498917e+04  1.05596176e+04]]
