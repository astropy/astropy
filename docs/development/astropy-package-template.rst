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


.. _package-template: https://github.com/astropy/package-template
