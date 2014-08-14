========================================================
How to create and maintain an Astropy affiliated package
========================================================

If you run into any problems, don't hesitate to ask for help on the
astropy-dev mailing list!

This package provides a template for packages that are affiliated with the
`Astropy`_ project. This package design mirrors the layout of the main
`Astropy`_ repository, as well as reusing much of the helper code used to
organize `Astropy`_.  The instructions below describe how to take this
template and adjust it for your particular affiliated package.

Everywhere below that the text ``yourpkg`` is shown, replace it with the name
of your particular package.

**Note**: The instructions below assume you are using git for version control,
as is used by the Astropy repository. If this is not the case, hopefully it
will be clear from context what to do with your particular VCS.

* Make sure `Astropy`_ is installed, as the template depends in part on
  Astropy to do its setup.

* You may have already done this if you are looking at this file locally, but
  if not, you will need to obtain a copy of the package template.  Assuming
  you have `git`_ installed, just do::

      git clone git://github.com/astropy/package-template.git yourpkg

  This will download the latest version of the template from `github`_ and
  place it in a directory named ``yourpkg``.

* Go into the directory you just created, and open the ``setup.cfg``
  file with your favorite text editor.  Edit the settings in the
  ``metadata`` section.  These values will be used to automatically
  replace special placeholders in the affiliated package template.

  1. Change the ``package_name`` variable to whatever you decide your
     package should be named. By tradition/very strong suggestion,
     python package names should be all lower-case.
  2. Change the ``description`` variable to a short (one or few
     sentence) description of your package.
  3. Add your name and email address by changing the ``author`` and
     ``author_email`` variables.
  4. If your affiliated package has a website, change ``url`` to point
     to that site.  Otherwise, you can leave it pointing to `Astropy`_
     or just delete it.
  5. Exit out of your text editor

* Update the main package docstring in ``packagename/__init__.py``.

* Decide what license you want to use to release your source code. If
  you don't care and/or are fine with the Astropy license, just edit
  the file ``licenses/LICENSE.rst`` with your name (or your
  collaboration's name) at the top as the licensees. Otherwise, make
  sure to replace that file with whatever license you prefer, and
  update the ``license`` variable in ``setup.cfg`` to reflect your
  choice of license. You also may need to update the comment at the
  top of ``packagename/__init__.py`` to reflect your choice of
  license.

* Take a moment to look over the ``packagename/example_mod.py``,
  ``packagename/tests/test_example.py``, ``scripts/script_example``,
  and ``packagename/example_c.pyx`` files, as well as the
  ``packagename/example_subpkg`` directory. These are examples of a
  pure-python module, a test script, an example command-line script, a
  `Cython`_ module, and a sub-package, respectively. (`Cython`_ is a
  way to compile python-like code to C to make it run faster - see the
  project's web site for details). These are provided as examples of
  standard way to lay these out. Once you understand these, though,
  you'll want to delete them (and later replace with your own)::

    git rm packagename/example_mod.py
    git rm scripts/script_example
    git rm packagename/example_c.pyx
    git rm packagename/tests/test_example.py
    git rm -r packagename/example_subpkg
    git commit -m "removed examples from package template"

* Optional: If you're hosting your source code on github, you can
  enable a sphinx extension that will link documentation pages
  directly to github's web site. To do this, set ``edit_on_github`` in
  ``setup.cfg`` to ``True`` and set ``github_project`` to the name of
  your project on github.

* Move the main source directory to reflect the name of your package.
  To tell your DVCS about this move, you should use it, and not ``mv``
  directly, to make the move.  For example, with git::

    git mv packagename yourpkg

* Update the names of the documentation files to match your package's name.
  First open ``docs/index.rst`` in a text editor and change the text
  ``"packagename/index.rst"`` to e.g., ``"yourpkg/index.rst"``.  Then do::

    git add docs/index.rst
    git mv docs/packagename docs/yourpkg

* Edit this file (``README.rst``) and delete all of this content, and replace it
  with a short description of your affiliated package.

*  Open ``docs/yourpkg/index.rst`` and you can start writing the documentation
   for your package, but at least replace ``packagename`` in ``automodapi::``
   with your package name.

* Now tell git to remember the changes you just made::

    git commit -a -m "Adjusted for new project yourpkg"

* (This step assumes your affiliated package is hosted as part of the astropy
  organization on Github.  If it's instead hosted somewhere else, just adjust
  the URL in the instructions below to match wherever your repository lives)
  Now you will want to tell git that it should be pushing and pulling updates
  to the repository of *your* project, rather than the package template::

    git remote rename origin template
    git remote add upstream git@github.com:astropy/yourpkg.git

  Now that it is pointing to the correct master, you should push everything up
  to your project and make sure that your local master is tied to your project
  rather than the template.  You'll only be able to do this if your github
  repository is empty (if not, add the ``-f`` option to the ``push``
  command - that will overwrite whatever is there)::

    git push upstream master
    git branch master --set-upstream upstream/master

* (optional) If you are adopting the standard workflow used by `Astropy`_ with
  github, you will also want to set up a fork of the repo on your own account,
  by going to the Github page https://github.com/astropy/yourpkg and clicking
  the "fork" button on the upper right.  Then run the following commands::

    git remote add origin git@github.com:yourgithubusername/yourpkg.git
    git branch master --set-upstream origin/master

  Now you can push, pull, and branch whatever you want in your local fork
  without affecting the official version, but when you want to push something
  up to the main repository, just switch to the appropriate branch and do
  ``git push upstream master``.

  Additionally, you can set things up to make it easier to pull future
  changes to the package template to your affiliated package.  Add a remote
  for the package template::

    git remote add template git@github.com:astropy/package-template.git

  Then, each time you want to pull in changes to the package template::

    git fetch template
    git fetch upstream

    # Make your master match the upstream master.  This will destroy
    # any unmerged commits on your master (which you shouldn't be doing
    # work on anyway, according to the standard workflow).
    git checkout master
    git reset --hard upstream/master

    # Merge any recent changes from the package-template
    git merge template/master

    # ...possibly resolve any conflicts...

    # Push to upstream master
    git push upstream master

* You should register your package on https://travis-ci.org and modify the
  ``.travis.yml`` file to make the build pass. This will continuously test
  your package for each commit, even pull requests against your main repository
  will be automatically tested, so that you notice when something breaks.
  For further information see
  `here <https://github.com/astropy/astropy/wiki/Continuous-Integration>`__
  and for lot's of example ``.travis.yml`` build configurations see
  `here <https://github.com/astropy/astropy/wiki/travis-ci-test-status>`__.
  Generally you should aim to always have your ``master`` branch work with
  the latest stable as well as the latest development version of astropy
  (i.e. the astropy git master branch) and the same versions of python and
  numpy supported by astropy. The template ``.travis.yml`` covers those
  versions; in some circumstances you may need to limit the versions your
  package covers.

* If you register your package with coveralls.io, then you will need
  to modify the ``coveralls --rcfile`` line in ``.travis.yml`` file to
  replace ``packagename`` with the name of your package.

* If you want the documentation for your project to be hosted by
  `ReadTheDocs <https://readthedocs.org>`_, then you need to setup an
  account there. The following entries in "Advanced Settings" for your
  package on `ReadTheDocs <https://readthedocs.org>`_ should work:

  - activate ``Install your project inside a virtualenv using setup.py install``
  - Requirements file: ``docs/rtd-pip-requirements``
  - activate ``Give the virtual environment access to the global site-packages dir.``

  All other settings can stay on their default value.

* You're now ready to start doing actual work on your affiliated package.  You
  will probably want to read over the developer guidelines of the Astropy
  documentation, and if you are hosting your code in GitHub, you might also
  want to read the `Github help <http://help.github.com/>`_ to ensure you know
  how to push your code to GitHub and some recommended workflows that work for
  the core Astropy project.

* Once you have started work on the affiliated package, you should register
  your package with the Astropy affiliated package registry. Instructions for
  doing this will be provided on the `Astropy`_ website.

* Good luck with your code and your science!

.. _git: http://git-scm.com/
.. _github: http://github.com
.. _Cython: http://cython.org/
