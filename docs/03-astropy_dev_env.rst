How to set up a development version of Astropy
==============================================

Found a bug? Know how to fix it? See a modification you need for your
own project? This tutorial will teach you how to set up a conda
environment for installing and modifying a developer version of Astropy.

Pre-requisites
--------------

-  `Git <https://git-scm.com/>`__ version control system. On a Mac, this
   will come with installing the latest version of the `Xcode
   developer <https://developer.apple.com/xcode/>`__ tools.
-  An account on `GitHub <https://github.com/>`__

This tutorial will show you how to fork, clone, and install the
development version of Astropy from the command line (terminal).

Step 1: Fork the Astropy repo on Github
---------------------------------------

-  Log into your account on Github
-  Go to https://github.com/astropy/astropy and click on **Fork** in the
   upper right-hand corner of the page. This will create a separate repo
   in your Github account under github.com//astropy

Step 2: Clone your forked Astropy repo to your computer
-------------------------------------------------------

-  Go to your forked version of the repository and click on **Clone or
   download** (the green button on the upper right-hand side of the
   page).
-  Highlight and copy the revealed URL. If you prefer using SSH to
   connect to Github, you can toggle between "Use SSH" and "Use HTTPS"
   by clicking on the highlighted text in the upper right-hand corner of
   the box.
-  Open a terminal session
-  Make an astropy folder in the location you would like to develop the
   code. For example:

   ::

       mkdir astropy-dev
       cd astropy-dev
       git clone git@github.com:<user-name>/astropy.git .

**IMPORTANT NOTE:** If it has been a long time since you forked the
Astropy repo for the first time, your fork may out of date. Before
proceeding, follow the instructions for **Keeping your fork up to date**
at the bottom of this page.

Step 3: Create a new conda environment for developing astropy
-------------------------------------------------------------

-  Inside the astropy-dev folder, create a new conda environment called
   *astropy-dev* and activate it.

   ::

       conda create -n astropy-dev --only-deps python=3 astropy
       source activate astropy-dev
       
   The *--only-deps* flag indicates that the new environment will be installed with all of astropy's dependencies without installing astropy itself.

-  Inside the astropy-dev folder, install the Astropy from the source
   code:

   ::

       python setup.py develop
       
   The *develop* command will install the package in way that does not require you to re-install any time you make changes to your astropy-dev library.

-  At this point, you could install any other packages you normally work
   with, *except* for astropy. For example,

   ::

       conda install pandas

Step 4: Make a new branch for editing the source code
-----------------------------------------------------

-  Inside the astropy-dev folder, create a new branch to work in. Try to
   pick a descriptive name for your branch:

   ::

       git branch my-code-update
       git checkout my-code-update

Congratulations! You're ready to edit the astropy source code safely.
Use ``git`` for committing your changes and pushing those changes to
your ``my-code-update`` branch on Github.

Here's an example of creating a new file and pushing it to the new branch:

::

       echo "print('hello world')" > my-new-code.py
       git add my-new-code.py
       git commit -m "first commit"
       git push origin my-code-update

For more help with learning Git and Github, see the `Github
Help <https://help.github.com/>`__ pages, the `Git and Github Learning
Resources <https://help.github.com/articles/git-and-github-learning-resources/>`__
guide, `try.github.io <http://try.github.io/>`_, and the `Github Learning Lab  <https://lab.github.com/>`_

When you are done working on development, don't forget to exit the
astropy-dev environment.

::

    source deactivate

When you want your changes to be incorporated into Astropy, `submit a pull request on Github <https://help.github.com/articles/creating-a-pull-request/>`__. If you're looking for `quick ways to contribute to Astropy <http://www.astropy.org/contribute.html#code>`__, check out the issues page on the main Astropy Github repo.


Keeping your fork up to date
============================

If it has been awhile since you last forked the Astropy repo, your fork
will likely need updating. This tutorial will teach you how to use the
command line to connect to the main Astropy Github repo, update to the
latest development version, and push those changes to your personal
Astropy fork on Github.

Step 1: Add the core Astropy repo to your git config file
---------------------------------------------------------

-  Go to https://github.com/astropy/astropy and click on the green
   "Clone or Download" button in the upper right-hand corner of the
   page.
-  Go to the directory where you have the Astropy repo on your computer,
   for example:

   ::

       cd ~/astropy-dev

-  Add a new remote named *upstream* to your local copy of the Astropy
   repo

   ::

       git remote add upstream git@github.com:astropy/astropy.git

Step 2: Pull any changes from the master branch of the main Astropy repo
------------------------------------------------------------------------

-  Make sure you are in the master branch of your local Astropy repo

   ::

       git checkout master

-  Use ``git pull`` to update your local master branch from the upstream
   Astropy master branch

   ::

       git pull upstream master

Step 3: Push the changes to your fork on Github
-----------------------------------------------

-  Use ``git push`` to update your Github fork of Astropy:

   ::

       git push origin master
   
   If you've already made some changes to your own master branch, you may need to force the push with the `--force` command. This may cause you to lose some changes or issues with your git history. This is why it's good practice to **always develop in a separate branch**.



Congratulations! You are up to date!
