.. _bounding-boxes:

Efficient Model Rendering with Bounding Boxes
=============================================

.. versionadded:: 1.1

All `astropy.modeling.Model` subclasses have a ``bounding_box`` attribute that
can be used to set the limits over which the model is significant. This greatly
improves the efficiency of evaluation when the input range is much larger than
the characteristic width of the model itself. For example, to create a sky model
image from a large survey catalog, each source should only be evaluated over the
pixels to which it contributes a significant amount of flux. This task can
otherwise be computationally prohibitive on an average CPU.

The `astropy.modeling.render_model` function can be used to evaluate a model on
an input array, or coordinates, limiting the evaluation to the ``bounding_box``
region if it is set. This function will also produce postage stamp images of the
model if no other input array is passed. To instead extract postage
stamps from the data array itself, see :ref:`cutout_images`.

Using the Bounding Box
-----------------------

For basic usage, see `astropy.modeling.Model.bounding_box`.
By default no bounding box is set (``bounding_box`` is `None`), except for
individual model subclasses that have a ``bounding_box_default`` function
defined. ``bounding_box_default`` returns the minimum rectangular region
symmetric about the position that fully contains the model if the model has a
finite extent. If a model does not have a finite extent, the choice for the
``bounding_box_default`` limits is noted in the docstring. For example, see
`astropy.modeling.functional_models.Gaussian2D.bounding_box_default`.

The default function can also be set to any callable. This is particularly
useful for fitting ``custom_model`` or ``CompoundModel`` instances.

    >>> from astropy.modeling import custom_model
    >>> def ellipsoid(x, y, z, x0=0, y0=0, z0=0, a=2, b=3, c=4, amp=1):
    ...     rsq = ((x-x0)/a) ** 2 + ((y-y0)/b) ** 2 + ((z-z0)/c) ** 2
    ...     val = (rsq < 1) * amp
    ...     return val
    ...
    >>> def ellipsoid_bb(self):
    ...     return ((self.z0 - self.c, self.z0 + self.c),
    ...             (self.y0 - self.b, self.y0 + self.b),
    ...             (self.x0 - self.a, self.x0 + self.a))
    ...
    >>> Ellipsoid3D = custom_model(ellipsoid)
    >>> Ellipsoid3D.bounding_box_default = ellipsoid_bb
    >>> model = Ellipsoid3D()
    >>> model.bounding_box = 'auto'
    >>> model.bounding_box
    ((-4.0, 4.0), (-3.0, 3.0), (-2.0, 2.0))

Efficient evaluation with ``render_model``
------------------------------------------

When a model is evaluated over a range much larger than the model itself, it may
be prudent to use `astropy.modeling.render_model` if efficiency is a concern.
The ``render_model`` function can be used to evaluate a model on an array of the
same dimensions. If no array is given, ``render_model`` will return a "postage
stamp" array corresponding to the bounding box region. However, if
``bounding_box`` is `None` an image or coordinates must be passed.

In this example, we generate a 300x400 pixel image of 100 2D
Gaussian sources both with and without using bounding boxes. Using bounding
boxes, the evaluation speed increases by approximately a factor of 10 with
negligible loss of information.

.. plot::
    :include-source:

	import numpy as np
	from time import time
	from astropy.modeling import models, render_model

	import matplotlib as mpl
	import matplotlib.pyplot as plt
	from astropy.visualization import astropy_mpl_style
	astropy_mpl_style['axes.grid'] = False
	astropy_mpl_style['axes.labelcolor'] = 'k'
	mpl.rcParams.update(astropy_mpl_style)

	np.random.seed(0)

	imshape = (300, 400)
	nsrc = 100

	model_params = [
	    dict(amplitude = np.random.uniform(0, 1),
	         x_mean = np.random.uniform(0, imshape[1]),
	         y_mean = np.random.uniform(0, imshape[0]),
	         x_stddev = np.abs(np.random.uniform(3, 6)),
	         y_stddev = np.abs(np.random.uniform(3, 6)),
	         theta = np.random.uniform(0, 2 * np.pi))
	    for i in range(nsrc)]

	model_list = [models.Gaussian2D(**kwargs) for kwargs in model_params]

	#Evaluate all models over their bounded regions and over the full image
	#for comparison.

	def make_image(model_list, shape=imshape, mode='bbox'):
	    image = np.zeros(imshape)
	    t1 = time()
	    for i,model in enumerate(model_list):
	        if mode == 'full': model.bounding_box = None
	        elif mode == 'auto': model.bounding_box = 'auto'
	        image = render_model(model, image)
	    t2 = time()
	    return image, (t2 - t1)

	bb_image, t_bb = make_image(model_list, mode='auto')
	full_image, t_full = make_image(model_list, mode='full')

	flux = full_image.sum()
	diff = (full_image - bb_image)
	max_err = diff.max()

	plt.figure(figsize=(16, 7))
	plt.subplots_adjust(left=.05,right=.97,bottom=.03,top=.97,wspace=0.1)#07)

	plt.subplot(121)
	plt.imshow(full_image, origin='lower')
	plt.axis([0,imshape[1],0,imshape[0]])
	plt.title('Full Models\nTiming: %.2f seconds' % (t_full), fontsize=16)
	plt.xlabel('x', fontsize=14)
	plt.ylabel('y', fontsize=14)

	plt.subplot(122)
	plt.imshow(bb_image, origin='lower')
	for model in model_list:
	    y1,y2,x1,x2 = np.reshape(model.bounding_box_default(),(4,))
	    plt.plot([x1,x2,x2,x1,x1], [y1,y1,y2,y2,y1], 'w-',alpha=.2)

	plt.axis([0,imshape[1],0,imshape[0]])
	plt.title('Bounded Models\nTiming: %.2f seconds' % (t_bb), fontsize=16)
	plt.xlabel('x', fontsize=14)
	plt.ylabel('y', fontsize=14)

	plt.figure(figsize=(16,8))
	plt.imshow(diff, vmin=-max_err, vmax=max_err)
	plt.colorbar(format='%.1e')
	plt.title('Difference Image\nTotal Flux Err = %.0e'
	            %((flux - np.sum(bb_image)) / flux), fontsize=16)
	plt.xlabel('x', fontsize=14)
	plt.ylabel('y', fontsize=14)
	plt.show()

