**********************************************************************
How to create and maintain a Python package using the Astropy template
**********************************************************************

If you run into any problems, don't hesitate to ask for help on the
astropy-dev mailing list!

The `package-template`_ repository provides a template for Python
packages. This package design mirrors the layout of the main `Astropy`_
repository, as well as reusing much of the helper code used to organize
`Astropy`_. The instructions below describe how to take this template and
adjust it for your particular package, as well as how to update
your package to the latest version of the package template.

There are two main ways you can use this template layout to create a new
package:

#. simply copy over the files you need manually

#. start from a fork of the package-template repository

There are advantages and disadvantages to both methods. In the first case your
repository history will be clean and you will only include the files you really
need. However, when updating you will need to make sure you update all the
files manually. In the second case, the history of your package will be
cluttered with all the commits from the `package-template`_ repository, but if
done properly this can be easier to update in future since it simply involves
pulling from the latest version of the `package-template`_ repository and then
merging in the changes (and resolving conflicts).

.. note:: The instructions below assume you are using git for version control,
          as is used by the Astropy repository. If this is not the case,
          hopefully it will be clear from context what to do with your
          particular VCS.

Everywhere below that the text ``<packagename>`` is shown, replace it with the
name of your particular package. In fact, you can do this automatically by
entering your package name in the box below:

.. raw:: html

    <form id="myform">
      Package name: <input type="text" name="packagename">
      <button>Update package name</button>
    </form>
    <br>

Managing the template files manually
************************************

Starting a new package
======================

#. Clone the `package-template`_ package so that we can have access to the
   files, but do not go inside it - instead, create another empty repository
   into which we will copy the required files::

    git clone https://github.com/astropy/package-template.git template
    mkdir <packagename>
    cd <packagename>
    git init

#. The `package-template`_ infrastructure relies on the `astropy-helpers`_
   package, and we recommend adding this as a sub-module so as to easily be
   able to bundle it in releases of your packages::

    git submodule add https://github.com/astropy/astropy-helpers.git astropy_helpers

#. Copy over the following files from the package template (these define the
   bare minimum of what is needed) and add them to the repository::

    # .gitignore specifies which types of files to ignore in the repository
    cp ../template/.gitignore .

    # MANIFEST.in specifies which files to include in a tar file
    release. Change all ``'packagename'`` text to be the name of your
    package, <packagename>.
    cp ../template/MANIFEST.in .

    # ah_bootstrap.py is used by setup.py to use the astropy-helpers
    cp ../template/ah_bootstrap.py .

    # ez_setup.py is used by setup.py to be able to get setuptools on-the-fly
    cp ../template/ez_setup.py .

    # setup.cfg contains details about your package - edit it after copying!
    # Note: it is important that the package_name variable matches the name you
    #       are using for your package.
    cp ../template/setup.cfg .

    # edit the VERSION variable and if applicable, the package_data values,
    # the rest can be kept as-is
    # Note: If your package data directory has a sub-directory, you HAVE TO
    #       break 'data/*' to 'data/*.*' and 'data/subdir/*' (two different
    #       appends). If you do not have any package data, you can just comment
    #       out the line that appends 'data/*'.
    cp ../template/setup.py .

   .. important:: Before proceeding, make sure you have edited
                 ``MANIFEST.in``, ``setup.cfg`` and ``setup.py`` as
                 indicated above!

   Once you have edited ``MANIFEST.in``, ``setup.cfg`` and ``setup.py``, you
   can commit the changes::

    git add .gitignore MANIFEST.in ah_bootstrap.py ez_setup.py setup.cfg setup.py

#. Next, you can create a directory for your package's source code, which will
   usually also be called the same name as your package. In this directory
   you can copy over the following files::

    mkdir <packagename>

    cp ../template/packagename/__init__.py <packagename>/
    cp ../template/packagename/_astropy_init.py <packagename>/

    # edit <packagename>/__init__.py to change the docstring and change the
    # example import ``from example_mod import *`` to whatever is needed for
    # your package. If you don't want to try out any specific code yet, just
    # replace the import by ``pass``.

   The main purpose of the ``_astropy_init.py`` file is to set up the
   ``test()`` command at the root of your package so that you can do
   ``<packagename>.test()``. This file is imported into ``__init__``.

   .. important:: Before proceeding, make sure you have edited ``__init__.py`` as
                  indicated above!

   Once you have made the above changes, you can commit the files::

    git add <packagename>/__init__.py
    git add <packagename>/_astropy_init.py

#. In order to benefit from the pytest plugins in Astropy, you should also
   copy over the ``conftest.py`` file to your repository::

    cp ../template/packagename/conftest.py <packagename>/

    git add <packagename>/conftest.py

   You can also uncomment the line ``enable_deprecations_as_exceptions()`` if
   you want deprecation warnings to make tests fail. There are also
   options to customize the information to be printed when running the
   tests. The package template has comments in the ``conftest.py`` file that
   indicate what they are.

#. If you are interested in accurate coverage test results, copy over the
   ``coveragerc`` and the ``setup_package.py`` files to your repository (the
   latter ensures that ``coveragerc`` gets installed with the package::

    mkdir <packagename>/tests/
    cp ../template/packagename/tests/__init__.py <packagename>/tests
    cp ../template/packagename/tests/setup_package.py <packagename>/tests
    cp ../template/packagename/tests/coveragerc <packagename>/tests

    git add <packagename>/tests/__init__.py
    git add <packagename>/tests/setup_package.py
    git add <packagename>/tests/coveragerc

   to your repository. When you run tests with with ``--coverage`` option this
   file will be used to exclude certain files that should not typically be
   included. Note that you don't need to change the ``{packagename}`` string in
   ``coveragerc`` - this gets changed automatically using the package name
   defined in ``setup.cfg``.

   .. note:: the ``python setup.py`` commands will not work until you
             have made your first commit, as shown in the last step of these
             instructions.

#. To set up the infrastructure to build the documentation, copy over the
   following files into a new directory called ``docs``::

    mkdir docs
    cp -r ../template/docs/_templates docs/
    cp ../template/docs/Makefile docs/
    cp ../template/docs/conf.py docs/
    cp ../template/docs/make.bat docs/
    touch docs/index.rst  # creates empty page
    git add docs/_templates docs/Makefile docs/conf.py docs/make.bat docs/index.rst

   you can later start adding content to ``index.rst`` and other documentation
   files.

#. Add a ``README.md`` file to your repository, describing what the package
   does, and for example how to install it and any required dependencies::

    git add README.md

#. Finally, if you plan on using Travis for continuous integration, copy over
   the ``.travis.yml`` file and edit it by replacing "packagname" with your
   package.  (You can do a lot more to customize the travis builds.  More
   detail on that below.)::

    cp ../template/.travis.yml .
    # edit .travis.yml
    git add .travis.yml

   .. important:: Before proceeding, make sure you have edited ``.travis.yml`` as
                  indicated above!

#. Double check that you have changed all the ``'packagename'`` references
   in the configuration files to be <packagename>, the name of your package.

#. Now you are ready to make your first commit::

    git commit -m "Initial layout for package"

#. You can test that your package works correctly by doing e.g.::

    python setup.py build
    python setup.py test --coverage
    python setup.py build_docs

   If you have any issues that you cannot fix, feel free to ask us on the
   `astropy-dev mailing list`_!

Updating to the latest template files
=====================================

From time to time we will make changes to the package-template to fix bugs or
add functionality. Updating to the latest version is simple - simply check
the `TEMPLATE_CHANGES.md`_ file, which provides a changelog of the package
template. You can also re-copy over all the files listed in the above section
and see if any of the changes should be committed (some of the changes will
be reverting some of your edits, so do not include those!). Remember to
update the astropy-helpers sub-module to the latest stable version, and
update the corresponding ``ah_bootstrap.py`` file, for example::

    cd astropy_helpers
    git fetch origin
    git checkout v0.4.3
    cd ..
    cp astropy_helpers/ah_bootstrap.py .
    git add astropy_helpers ah_bootstrap.py
    git commit -m "Updated astropy-helpers to v0.4.3"

You can find out what the latest version of astropy-helpers is by checking the
`astropy-helpers <https://pypi.python.org/pypi/astropy-helpers/>`__ entry on
PyPI.

Customizing the documentation CSS
=================================

As described in the documentation configuration file (`template/docs/conf.py
<https://github.com/astropy/package-template/blob/master/docs/conf.py#L95>`_),
the documentation uses a custom theme based on `bootstrap
<http://getbootstrap.com/css/>`_. You can swap out this theme by editing the
configuration file. You can also tweak aspects of the documentation theme by
creating a custom CSS file in your package documentation.

To do this, create a new CSS file in ``<packagename>/_static/`` -- let's call it
``<packagename>.css``::

    cd <packagename>/_static/
    touch <packagename>.css

We're going to set the HTML style to this new ``<packagename>.css`` stylesheet,
so we need to import the original ``bootstrap-astropy`` style before we start
modifying entries. To the first line of your ``<packagename>.css`` file, import
the default style. We can add any custom CSS below the import. For example, to
hide the Astropy logo and Astropy link from your project's documentation menu
bar:

.. code-block:: css

    @import url("bootstrap-astropy.css");

    div.topbar a.brand {
        background: none;
        background-image: none;
    }

    div.topbar ul li a.homelink {
        background: none;
        background-image: none;
    }

We now have to include the ``<packagename>.css`` in the documentation, and tell
Sphinx to use the new style. To do this, edit your
``<packagename>/docs/conf.py`` file and add the lines::

    # Static files to copy after template files
    html_static_path = ['_static']
    html_style = '<packagename>.css'

Managing the template files via git
***********************************

Starting a new package
======================

Before reading this we recommend reading over the `Managing the template
files manually`_ section since this explains what many of the files do.

#. Make sure `Astropy`_ is installed, as the template depends in part on
   Astropy to do its setup.

#. You may have already done this if you are looking at this file locally, but
   if not, you will need to obtain a copy of the package template.  Assuming
   you have `git`_ installed, just do::

      git clone git://github.com/astropy/package-template.git <packagename>

   This will download the latest version of the template from `github`_ and
   place it in a directory named ``<packagename>``.

#. Go into the directory you just created, and open the ``setup.cfg``
   file with your favorite text editor.  Edit the settings in the
   ``metadata`` section.  These values will be used to automatically
   replace special placeholders in the package template.

   #. Change the ``package_name`` variable to whatever you decide your
      package should be named. By tradition/very strong suggestion,
      python package names should be all lower-case.
   #. Change the ``description`` variable to a short (one or few
      sentence) description of your package.
   #. Add your name and email address by changing the ``author`` and
      ``author_email`` variables.
   #. If your package has a website, change ``url`` to point
      to that site.  Otherwise, you can leave it pointing to `Astropy`_
      or just delete it.
   #. Exit out of your text editor.

#. Move the main source directory to reflect the name of your package.
   To tell your DVCS about this move, you should use it, and not ``mv``
   directly, to make the move.  For example, with git::

    git mv packagename <packagename>

#. Update the main package docstring in ``<packagename>/__init__.py``.

#. Decide what license you want to use to release your source code. If
   you don't care and/or are fine with the Astropy license, just edit
   the file ``licenses/LICENSE.rst`` with your name (or your
   collaboration's name) at the top as the licensees. Otherwise, make
   sure to replace that file with whatever license you prefer, and
   update the ``license`` variable in ``setup.cfg`` to reflect your
   choice of license. You also may need to update the comment at the
   top of ``<packagename>/__init__.py`` to reflect your choice of
   license.

#. Take a moment to look over the ``<packagename>/example_mod.py``,
   ``<packagename>/tests/test_example.py``, and ``<packagename>/example_c.pyx``
   files, as well as the ``<packagename>/example_subpkg`` directory.
   These are examples of a pure-python module, a test script, a
   `Cython`_ module, and a sub-package, respectively. (`Cython`_ is a
   way to compile python-like code to C to make it run faster - see the
   project's web site for details). These are provided as examples of
   standard way to lay these out. Once you understand these, though,
   you'll want to delete them (and later replace with your own)::

      git rm <packagename>/example_c.pyx
      git rm <packagename>/tests/test_example.py
      git rm -r <packagename>/example_subpkg
      git commit -m "removed examples from package template"

#. Optional: If you're hosting your source code on github, you can
   enable a sphinx extension that will link documentation pages
   directly to github's web site. To do this, set ``edit_on_github`` in
   ``setup.cfg`` to ``True`` and set ``github_project`` to the name of
   your project on github.

#. Update the names of the documentation files to match your package's name.
   First open ``docs/index.rst`` in a text editor and change the text
   ``"packagename/index.rst"`` to e.g., ``"<packagename>/index.rst"``.  Then do::

      git add docs/index.rst
      git mv docs/packagename docs/<packagename>

#. Edit the ``README.rst`` file, deleting all of the content and replacing it
   with a short description of your package.

#.  Open ``docs/<packagename>/index.rst`` and you can start writing the documentation
    for your package, but at least replace ``packagename`` in ``automodapi::``
    with your package name.

#. Now tell git to remember the changes you just made::

      git commit -a -m "Adjusted for new project <packagename>"

#. (This step assumes your package is hosted as part of the astropy
   organization on Github.  If it's instead hosted somewhere else, just adjust
   the URL in the instructions below to match wherever your repository lives)
   Now you will want to tell git that it should be pushing and pulling updates
   to the repository of *your* project, rather than the package template::

      git remote rename origin template
      git remote add upstream git@github.com:astropy/<packagename>.git

   Now that it is pointing to the correct master, you should push everything up
   to your project and make sure that your local master is tied to your project
   rather than the template.  You'll only be able to do this if your github
   repository is empty (if not, add the ``-f`` option to the ``push``
   command - that will overwrite whatever is there)::

      git push upstream master
      git branch master --set-upstream upstream/master

#. (optional) If you are adopting the standard workflow used by `Astropy`_ with
   github, you will also want to set up a fork of the repo on your own account,
   by going to the Github page https://github.com/astropy/<packagename> and clicking
   the "fork" button on the upper right.  Then run the following commands::

      git remote add origin git@github.com:yourgithubusername/<packagename>.git
      git branch master --set-upstream origin/master

   Now you can push, pull, and branch whatever you want in your local fork
   without affecting the official version, but when you want to push something
   up to the main repository, just switch to the appropriate branch and do
   ``git push upstream master``.

   Additionally, you can set things up to make it easier to pull future
   changes to the package template to your package.  Add a remote
   for the package template::

      git remote add template git@github.com:astropy/package-template.git

   .. _template-changes-with-gitfu:

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

#. You should register your package on https://travis-ci.org and modify the
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

#. If you register your package with coveralls.io, then you will need
   to modify the ``coveralls --rcfile`` line in ``.travis.yml`` file to
   replace ``packagename`` with the name of your package.

#. If you want the documentation for your project to be hosted by
   `Read the Docs <https://readthedocs.org>`_, then you need to setup an
   account there. The following entries in "Advanced Settings" for your
   package on `Read the Docs <https://readthedocs.org>`_ should work:

   - activate ``Install your project inside a virtualenv using setup.py install``
   - copy these additional files from the package template into the top-level
     directory of your package: ``.rtd-environment.yml`` and ``readthedocs.yml``.
     Edit ``.rtd-environment.yml`` with your package name and requirements.
   - activate ``Give the virtual environment access to the global site-packages dir.``

   All other settings can stay on their default value.

   If you need to mock any Python packages or C libraries that can not be
   installed and built by Read the Docs, you should include the following mocking
   patch before the ``Project information`` section of the ``docs/conf.py`` file::

      class Mock(object):
          def __init__(self, *args, **kwargs):
              pass

          def __call__(self, *args, **kwargs):
              return Mock()

          @classmethod
          def __getattr__(cls, name):
              if name in ('__file__', '__path__'):
                  return '/dev/null'
              elif name[0] == name[0].upper():
                  return type(name, (), {})
              else:
                  return Mock()

      MOCK_MODULES = ['<name of package to mock>', '<name of package to mock>']
      for mod_name in MOCK_MODULES:
          sys.modules[mod_name] = Mock()

#. You're now ready to start doing actual work on your package.  You
   will probably want to read over the developer guidelines of the Astropy
   documentation, and if you are hosting your code in GitHub, you might also
   want to read the `Github help <https://help.github.com/>`_ to ensure you know
   how to push your code to GitHub and some recommended workflows that work for
   the core Astropy Project.

#. Good luck with your code and your science!

Updating to the latest template files
=====================================

See instructions in :ref:`Item 14 above <template-changes-with-gitfu>`.


Releasing a Python package
**************************

You can release a package using the steps given below. In these
instructions, we assume that the release is made from a fresh clone of the
remote "main" repository and not from a forked copy. We also assume that
the changelog file is named ``CHANGES.rst``, like for the astropy core
package. If instead you use Markdown, then you should replace ``CHANGES.rst``
by ``CHANGES.md`` in the instructions.

#. Make sure that Travis and any other continuous integration is passing.

#. Update the ``CHANGES.rst`` file to make sure that all the changes are listed,
   and update the release date, which should currently be set to
   ``unreleased``, to the current date in ``yyyy-mm-dd`` format.

#. Update the version number in ``setup.cfg`` to the version you're about to
   release, without the ``.dev`` suffix (e.g. ``0.1``).

#. Run ``git clean -fxd`` to remove any untracked files (WARNING: this will
   permanently remove any files that have not been previously committed, so
   make sure that you don't need to keep any of these files).

#. Run::

        python setup.py build sdist --format=gztar

   and make sure that generated file is good to
   go by going inside ``dist``, expanding the tar file, going inside the
   expanded directory, and running the tests with::

        python setup.py test

   You may need to add the ``--remote-data`` flag or any other flags that you
   normally add when fully testing your package.

   .. note::

       Running ``python setup.py build sdist`` runs two setup commands in
       succession.  First it runs ``build``, then immediately runs ``sdist``
       to create the source distribution.  The reason to do this is that
       there are several generated source files that must be included in the
       source distribution for it to be valid.  Running ``build`` first
       ensures that those files will be generated and packaged in the source
       distribution.

#. Go back to the root of the directory and remove the generated files with::

        git clean -fxd

#. Add the changes to ``CHANGES.rst`` and ``setup.cfg``::

        git add CHANGES.rst setup.cfg

   and commit with message::

        git commit -m "Preparing release <version>"

#. Tag commit with ``v<version>``, optionally signing with the ``-s`` option::

        git tag v<version>

#. Change ``VERSION`` in ``setup.cfg`` to next version number, but with a
   ``.dev`` suffix at the end (e.g. ``0.2.dev``). Add a new section to
   ``CHANGES.rst`` for next version, with a single entry ``No changes yet``, e.g.::

       0.2 (unreleased)
       ----------------

       - No changes yet

#. Add the changes to ``CHANGES.rst`` and ``setup.cfg``::

        git add CHANGES.rst setup.cfg

   and commit with message::

        git commit -m "Back to development: <next_version>"

#. Check out the release commit with ``git checkout v<version>``.
   Run ``git clean -fxd`` to remove any non-committed files.

#. (optional) Run the tests in an environment that mocks up a "typical user"
   scenario. This is not strictly necessary because you ran the tests above, but
   it can sometimes be useful to catch subtle bugs that might come from you
   using a customized developer environment.  For more on setting up virtual
   environments, see :ref:`virtual_envs`, but for the sake of example we will
   assume you're using `Anaconda <https://conda.io/docs/>`_. Do::

       conda create -n myaffilpkg_rel_test astropy <any more dependencies here>
       source activate myaffilpkg_rel_test
       python setup.py sdist
       cd dist
       pip install myaffilpkg-version.tar.gz
       python -c 'import myaffilpkg; myaffilpkg.test()'
       source deactivate
       cd <back to your source>

   You may want to repeat this for other combinations of dependencies if you think
   your users might have other relevant packages installed.  Assuming the tests
   all pass, you can proceed on.

#. If you did the previous step, do ``git clean -fxd`` again to remove anything
   you made there.  Run ``python setup.py build sdist --format=gztar`` to
   create the files for upload.  Then you can upload to PyPI via ``twine``::

        twine upload dist/*

   as described in `these <https://packaging.python.org/tutorials/distributing-packages/#uploading-your-project-to-pypi>`_
   instructions. Check that the entry on PyPI is correct, and that
   the tarfile is present.

#. Go back to the master branch and push your changes to github::

        git checkout master
        git push --tags origin master

   Once you have done this, if you use Read the Docs, trigger a ``latest`` build
   then go to the project settings, and under **Versions** you should see the
   tag you just pushed. Select the tag to activate it, and save.

#. If your package is available in the ``conda-forge`` conda channel, you
   should also submit a pull request to update the version number in the
   feedstock of your package.

.. note:: The instructions above assume that you do not make use of bug fix
          branches in your workflow. If you do wish to create a bug fix branch,
          we recommend that you read over the more complete astropy
          :doc:`releasing` and adapt these for your package.


Modifications for a beta/release candidate release
==================================================

   Before a new release of your package, you may wish do a "pre-release" of the
   code, for example to allow collaborators to independently test the release.
   If the release you are performing is this kind of pre-release,
   some of the above steps need to be modified.

   The primary modifications to the release procedure is:

   * When entering the new version number, instead of just removing the
     ``.dev``, enter "1.2b1" or "1.2rc1".  It is critical that you follow this
     numbering scheme (``x.yb#`` or ``x.y.zrc#``), as it will ensure the release
     is ordered "before" the main release by various automated tools, and also
     tells PyPI that this is a "pre-release".


.. _git: https://git-scm.com/
.. _github: https://github.com/
.. _Cython: http://cython.org/
.. _package-template: https://github.com/astropy/package-template
.. _astropy-helpers: https://github.com/astropy/astropy-helpers
.. _TEMPLATE_CHANGES.md: https://github.com/astropy/package-template/blob/master/TEMPLATE_CHANGES.md

.. raw:: html

    <script>

    function get_url_vars() {
        var vars = {};
        var parts = window.location.href.replace(/[?&]+([^=&]+)=([^&]*)/gi, function(m,key,value) {
            vars[key] = value;
        });
        return vars;
    }

    packagename = get_url_vars()["packagename"]
    if(packagename) {
      document.body.innerHTML = document.body.innerHTML.replace(/&lt;packagename&gt;/g, packagename);
    }
    </script>
