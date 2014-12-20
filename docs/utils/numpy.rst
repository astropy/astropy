:orphan:

.. _numpy-compatibility:

NumPy compatibility
===================

NumPy_ forms an essential basis for astropy, and astropy's development has led
to the identification of problems with some of numpy's functionality. Often,
these are corrected in later versions of numpy, but in order for astropy not
to depend on these, work-arounds are made, usually in the code.  If functions
are used in more than one place, however, it can be more convenient to provide
patched routines. Hence, `astropy.utils.compat.numpy`.


Adding a patched routine
------------------------

To ensure that patched code is only used when required, and that it will be
easy to remove it if it is no longer needed for any supported version of
NumPy_, the following procedure should be used to add a patched routine:

* Copy over a correct version of the relevant numpy file to its
  corresponding location below the ``astropy/utils/compat/numpy`` directory.
* In this file, remove everything that does not have to be changed.  If
  necessary, import required pieces from numpy.
* Define a function that tests whether or not a patched version is needed, by
  directly testing whether the desired functionality is present. Suggested
  function names are ``PR####`` with a relevant numpy pull request number,
  or ``GE####`` with a version number.
* Place the redefinition of the relevant piece of code inside an ``if``
  statement that uses the function just defined.  This should ensure that if a
  sufficiently high version of numpy is used, no replacement is made.
* In ``numpy/__init__.py``, import your patched code.
* In ``numpy/tests``, add a new test routine that tests that the patch is used
  when necessary (i.e., test the test function), and that it provides the
  desired functionality. 

For an example, see ``numpy/lib/stride_tricks.py`` and the corresponding
``numpy/tests/test_broadcast_arrays.py``.

Note that patched routines will normally only be considered if they are part 
of NumPy_. Thus, if the patch concerns a new bug discovered in numpy, a `pull
request <https://github.com/numpy/numpy/pulls>`__ should first be made to
NumPy_ (which can of course form the basis of a `pull request 
<https://github.com/astropy/astropy/pulls>`__ to ``astropy``).


Reference/API
=============
.. automodapi:: astropy.utils.compat.numpy
    :no-inheritance-diagram:

.. _Numpy: http://numpy.scipy.org/
