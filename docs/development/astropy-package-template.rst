**********************************************************************
How to create and maintain a Python package using the Astropy template
**********************************************************************

If you run into any problems, don't hesitate to ask for help on the
astropy-dev mailing list!

The `package-template`_ repository provides a template for Python
packages. This package design mirrors the layout of the main `Astropy`_
repository, as well as reusing much of the helper code used to organize
`Astropy`_. See the
:ref:`package template documentation <packagetemplate:package-template>`
for instructions on using the package template.


.. _simple-release-docs:

Releasing a Python package
**************************

You can release a package using the steps given below. In these
instructions, we assume that the release is made from a fresh clone of the
remote "main" repository and not from a forked copy. We also assume that
the changelog file is named ``CHANGES.rst``, like for the astropy core
package. If instead you use Markdown, then you should replace ``CHANGES.rst``
by ``CHANGES.md`` in the instructions.

.. note:: The instructions below assume that you do not make use of bug fix
          branches in your workflow. If you do wish to create a bug fix branch,
          we recommend that you read over the more complete astropy
          :doc:`releasing` and adapt them for your package.

#. Make sure that continuous integration is passing.

#. Update the ``CHANGES.rst`` file to make sure that all the changes are listed,
   and update the release date, which should currently be set to
   ``unreleased``, to the current date in ``yyyy-mm-dd`` format.

#. Run ``git clean -fxd`` to remove any untracked files (WARNING: this will
   permanently remove any files that have not been previously committed, so
   make sure that you don't need to keep any of these files).

#. At this point, the command to run to build the tar file will depend on
   whether your package has a ``pyproject.toml`` file or not. If it does
   not, then::

        python setup.py build sdist --format=gztar

   If it does, then first make sure the `build <https://pypi.org/project/build/>`_
   package is installed and up-to-date::

        pip install build --upgrade

   then create the source distribution with::

        python -m build --sdist .

   All following instructions will assume you have ``pyproject.toml``.
   If you do not use ``pyproject.toml`` yet, please see
   https://docs.astropy.org/en/v3.2.x/development/astropy-package-template.html
   instead.

   In both cases, make sure that generated file is good to go by going inside
   ``dist``, expanding the tar file, going inside the expanded directory, and
   running the tests with::

        pip install -e .[test]
        pytest

   You may need to add the ``--remote-data`` flag or any other flags that you
   normally add when fully testing your package.

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
       python -m build --sdist .
       cd dist
       pip install myaffilpkg-version.tar.gz
       python -c 'import myaffilpkg; myaffilpkg.test()'
       source deactivate
       cd <back to your source>

   You may want to repeat this for other combinations of dependencies if you think
   your users might have other relevant packages installed.  Assuming the tests
   all pass, you can proceed on.

#. If you did the previous step, do ``git clean -fxd`` again to remove anything
   you made there.  Run ``python -m build --sdist .`` to
   create the files for upload.  Then you can upload to PyPI via ``twine``::

        twine upload dist/*

   as described in `these <https://packaging.python.org/en/latest/tutorials/packaging-projects/#uploading-the-distribution-archives>`_
   instructions. Check that the entry on PyPI is correct, and that
   the tarfile is present.

#. Go back to the main branch and push your changes to github::

        git checkout main
        git push --tags origin main

   Once you have done this, if you use Read the Docs, trigger a ``latest`` build
   then go to the project settings, and under **Versions** you should see the
   tag you just pushed. Select the tag to activate it, and save.

#. If your package is available in the ``conda-forge`` conda channel, you
   should also submit a pull request to update the version number in the
   feedstock of your package.


Modifications for a beta/release candidate release
==================================================

   Before a new release of your package, you may wish do a "pre-release" of the
   code, for example to allow collaborators to independently test the release.
   If the release you are performing is this kind of pre-release,
   some of the above steps need to be modified.

   The primary modifications to the release procedure is:

   * When entering the new version number, instead of just removing the
     ``.dev``, enter "1.2b1" or "1.2rc1".  It is critical that you follow this
     numbering scheme (``X.Yb#`` or ``X.Y.Zrc#``), as it will ensure the release
     is ordered "before" the main release by various automated tools, and also
     tells PyPI that this is a "pre-release".


.. _package-template: https://github.com/astropy/package-template
