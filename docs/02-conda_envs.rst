An Astropy User's Guide to Managing Conda Environments and Jupyter Notebook
===========================================================================


    Help! My code library works in Python 2 and I don't have time to
    make it work in Python 3

Do not fear! This tutorial will teach you how to use conda environments
to create a *separate* installation of Python 2.7. With conda
environments, you can switch between Python 2 and 3 without having to
worry about version conflicts.

Step 1: Set up a Python 2 environment
-------------------------------------

-  Create a new ``conda`` environment

   ::

       conda create -n python2 python=2.7 anaconda 

**NOTE:** By adding ``anaconda`` at the end, the complete Anaconda Python distribution will be installed. Omitting ``anaconda`` or choosing ``miniconda`` will install a much smaller library of Python packages. You can install any additional packages you need for that environment, after activating it. `Some information about the differences between Anaconda and Miniconda can be found in the conda document pages <https://conda.io/docs/user-guide/install/download.html#anaconda-or-miniconda>`__

-  Activate the Python 2 environment and install any additional packages
   you need to run your favorite code library. Here, we show you how to 
   install the base version of Astropy. 

    ::

        source activate python2
        conda install astropy

**NOTE:** If you want to install Astropy and **all** affiliated packages, you can use `` conda install stsci``

-  When you are ready to return to your default environment:

   ::

       source deactivate

**NOTE:** Some newer versions of Anaconda use ``conda activate`` and
``conda deactivate``. In that case, both ``source`` and ``conda`` will work 
interchangeably when activating or deactivating your chosen environment.

When you want to see all of your available environments:

::

    conda env list

Step 2: Check that your code runs in the new environment
--------------------------------------------------------

-  Now you are ready to work in Python 2! Here's a generic example for
   switching to your Python 2 environment, running your Python 2 script,
   and exiting the environment.

    ::

        cd ~/my-python2-library
        source activate python2
        python my_python2_script.py
        source deactivate

Step 3: Set up a Jupyter Notebook for the new environment
---------------------------------------------------------

-  Activate your custom Python 2 environment:

   ::

       source activate python2

-  Check that you have ipykernel installed

   ::

       conda list | grep ipykernel

If you do not see any output, install it with
``conda install ipykernel``

-  Install that environment for Jupyter notebook. In this case, we are
   choosing a display name, "python2", that matches the environment name, 
   but you may choose something else.

    ::

        python -m ipykernel install --user --name python2 --display-name "python2"`

- Now leave that environement

    ::

        source deactivate

-  Start a Jupyter Notebook session

   ::

       jupyter notebook

-  When you click on *New*, you should see a drop down list of options
   that include "python2", the environment display name we chose above.

-  If you would like to change the environment for an existing Jupyter
   Notebook, click on *Kernel*, then *Change kernel*, and select the
   environment you would like to use.

-  In general, you can view your available Jupyter Notebook kernels by
   running

   ::

       jupyter kernelspec list
