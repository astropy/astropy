Releasing a new version
=======================

1. Make sure that Travis is passing.

2. Update ``CHANGES.md`` to make sure that all the changes are listed, and
   update ``unreleased`` to the current date.

3. Update the version number in ``setup.py`` to ``v0.x`` (without the ``dev``)

4. Run ``git clean -fxd`` to remove any non-committed files.

5. Run:

        python setup.py sdist --format=gztar
        
   and make sure that generated file is good to
   go by going inside ``dist``, expanding the tar file, going inside the
   ``wcsaxes-x.x directory, and running the tests with:

        python setup.py test --remote-data

6. Go back to the root of the directory and remove the generated files with:

        git clean -fxd

7. Add the changes to ``CHANGES.md`` and ``setup.py``:

        git add CHANGES.md setup.py

   and commit with message:

        git commit -m "Preparing release <version>"

8. Tag commit with ``v<version>``, optionally signing with the ``-s`` option:
 
        git tag v<version>

9. Change ``VERSION`` in ``setup.py`` to next one with ``.dev``. Add a new
   section to ``CHANGES.md`` for next version, with a single entry, ``- No changes yet``.

10. Add the changes to ``CHANGES.md`` and ``setup.py``:

        git add CHANGES.md setup.py


    and commit with message:

        git commit -m "Back to development: <next_version>"

11. Check out the release commit with ``git checkout v<version>``. Run ``git clean -fxd`` to remove any non-committed files, then release with:

        python setup.py register sdist --format=gztar upload
