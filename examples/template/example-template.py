# -*- coding: utf-8 -*-
"""
========================
Title of Example
========================

This example <verb> <active tense> <does something>.

The example uses <packages> to <do something> and <other package> to <do other
thing>. Include links to referenced packages like this: `astropy.io.fits` to
show the astropy.io.fits or like this `~astropy.io.fits`to show just 'fits'

-------------------

*By: <names>*

*License: BSD*

-------------------

"""

##############################################################################
# Make print work the same in all versions of Python, set up numpy,
# matplotlib, and use a nicer set of plot parameters:

import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
# uncomment if including figures:
# import matplotlib.pyplot as plt
# from astropy.visualization import astropy_mpl_style
# plt.style.use(astropy_mpl_style)

##############################################################################
# This code block is executed, although it produces no output. Lines starting
# with a simple hash are code comment and get treated as part of the code
# block. To include this new comment string we started the new block with a
# long line of hashes.
#
# The sphinx-gallery parser will assume everything after this splitter and that
# continues to start with a **comment hash and space** (respecting code style)
# is text that has to be rendered in
# html format. Keep in mind to always keep your comments always together by
# comment hashes. That means to break a paragraph you still need to commend
# that line break.
#
# In this example the next block of code produces some plotable data. Code is
# executed, figure is saved and then code is presented next, followed by the
# inlined figure.

x = np.linspace(-np.pi, np.pi, 300)
xx, yy = np.meshgrid(x, x)
z = np.cos(xx) + np.cos(yy)

plt.figure()
plt.imshow(z)
plt.colorbar()
plt.xlabel('$x$')
plt.ylabel('$y$')

###########################################################################
# Again it is possible to continue the discussion with a new Python string. This
# time to introduce the next code block generates 2 separate figures.

plt.figure()
plt.imshow(z, cmap=plt.cm.get_cmap('hot'))
plt.figure()
plt.imshow(z, cmap=plt.cm.get_cmap('Spectral'), interpolation='none')

##########################################################################
# There's some subtle differences between rendered html rendered comment
# strings and code comment strings which I'll demonstrate below. (Some of this
# only makes sense if you look at the
# :download:`raw Python script <plot_notebook.py>`)
#
# Comments in comment blocks remain nested in the text.


def dummy():
    """Dummy function to make sure docstrings don't get rendered as text"""
    pass

# Code comments not preceded by the hash splitter are left in code blocks.

string = """
Triple-quoted string which tries to break parser but doesn't.
"""

############################################################################
# Output of the script is captured:

print('Some output from Python')

############################################################################
# Finally, I'll call ``show`` at the end just so someone running the Python
# code directly will see the plots; this is not necessary for creating the docs

plt.show()
