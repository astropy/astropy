:orphan:

Changelog
=========

This directory contains "news fragments" which are short files that contain a
small **ReST**-formatted text that will be added to the next what's new page.

Make sure to use full sentences with correct case and punctuation.

Each file should be named like ``<PULL REQUEST>.<TYPE>.rst``, where
``<PULL REQUEST>`` is a pull request number, and ``<TYPE>`` is one of:

* ``feature``: New feature.
* ``api``: API change.
* ``bugfix``: Bug fix.
* ``perf``: Performance improvement (this should be significant enough to be measurable using the public API).
* ``other``: Other changes and additions.

If the change concerns a sub-package, the file should go in the sub-directory
relative to this sub-package. Type ``other`` is not allowed in sub-directories.

It is possible to add two files with different types (and text) if both
are relevant. For example a change may add a new feature but introduce an API
change.

So for example: ``123.feature.rst`` would have the content::

    The ``my_new_feature`` option is now available for ``my_favorite_function``.
    To use it, write ``np.my_favorite_function(..., my_new_feature=True)``.

Note the use of double-backticks for code.

If you are unsure what pull request type to use, don't hesitate to ask in your
PR.

You can install ``towncrier`` and run ``towncrier --draft --version 4.3``
if you want to get a preview of how your change will look in the final release
notes.

.. note::

    This README was adapted from the Numpy changelog readme under the terms of
    the MIT licence.
