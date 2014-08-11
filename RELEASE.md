Releasing a new version
=======================

1. Make sure that Travis is passing.

2. Update ``CHANGES.md`` to make sure that all the changes are listed, and
   update ``unreleased`` to the current date.

3. Update the version number in ``setup.py``.

4. Run ``python setup.py clean -fxd`` to remove any non-committed files.

5. Run ``python setup.py sdist`` and make sure that generated file is good to
   go by going inside ``dist``, expanding the tar file, going inside the
   ``wcsaxes-x.x directory, and running the tests with:

        python setup.py test --remote-data

6. Go back to the root of the directory and commit with message:

        Preparing release <version>

7. Tag commit with ``v<version>``, optionally signing with the ``-s`` option.

8. Change ``VERSION`` in ``setup.py`` to next one with ``.dev``. Add a new
   section to ``CHANGES.md`` for next version.

9. Commit with message:

        Back to development: <next_version>

10. Check out the release commit with ``git checkout v<version>``. run ``python setup.py clean -fxd`` to remove any non-committed files, then release with:

        python setup.py register sdist upload
