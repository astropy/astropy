*******************************************
Emacs setup for following coding guidelines
*******************************************

.. _flycheck: https://www.flycheck.org/
.. _flake8: http://flake8.pycqa.org/

The Astropy coding guidelines are listed in :doc:`codeguide`. Here, we describe
how to configure Emacs to help ensure Python code satisfies the guidelines.

For this setup, we add to the standard ``python-mode`` using flycheck_ and the
flake8_ python style checker.  For installation instructions, see their
respective web sites (or install via your distribution; e.g., in Debian/Ubuntu,
the packages are called ``elpa-flycheck`` and ``flake8``).

.. note:: Emacs can be configured in several different ways. So instead of
          providing a drop in configuration file, only the individual
          configurations are presented below.

          The setup below is on purpose minimal.  In principle, it is possible
          to use `Emacs for Python development
          <https://realpython.com/emacs-the-best-python-editor/>`_,
          with, e.g., `elpy <https://elpy.readthedocs.io/>`_.

No tabs
=======

This setting will cause indentation to use spaces rather than tabs for all
files.  For python files, indentation of 4 spaces will be used if the tab key
is pressed.

.. code-block:: scheme

  ;; Don't use TABS for indentations.
  (setq-default indent-tabs-mode nil)

Delete trailing white spaces
============================

One can `delete trailing whitespace
<https://www.emacswiki.org/emacs/DeletingWhitespace#toc3>`_ with ``M-x
delete-trailing-whitespace``. To ensure this is done every time a python file
is saved, use:

.. code-block:: scheme

  ;; Automatically remove trailing whitespace when file is saved.
  (add-hook 'python-mode-hook
  (lambda () (add-to-list 'write-file-functions 'delete-trailing-whitespace)))

If you want to use this for every type of file, you can use
``(add-hook 'before-save-hook 'delete-trailing-whitespace)``.

Flycheck
========

One can make lines that do not satisfy syntax requirements using flycheck_.
When the cursor is on such a line a message is displayed in the mini-buffer.
When mouse pointer is on such a line a "tool tip" message is also shown. By
default, flycheck_ will check if flake8_ is installed and, if so, use that for
its syntax checking. To ensure flycheck_ starts upon opening python files, add:

.. code-block:: scheme

  (add-hook 'python-mode-hook 'flycheck-mode)

Alternatively, you can just use ``(global-flycheck-mode)`` to run flycheck
for all languages it supports.
