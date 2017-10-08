*******************************************
Emacs setup for following coding guidelines
*******************************************

.. _flycheck: http://www.flycheck.org/
.. _flake8: http://flake8.pycqa.org/
.. include:: workflow/known_projects.inc

The Astropy coding guidelines are listed in :doc:`codeguide`. This
document will describe some configuration options for Emacs, that will
help in ensuring that Python code satisfies the guidelines. Emacs can
be configured in several different ways. So instead of providing a drop
in configuration file, only the individual configurations are presented
below.

For this setup, in addition to the standard ``python-mode``, we will
need flycheck_ and the flake8_ python style checker.  For installation
instructions, see their respective web sites (or just install via your
distribution; e.g., in Debian/Ubuntu, the packages are called
``elpa-flycheck`` and ``flake8``).

Global settings
***************

No tabs
=======

This setting will cause all tabs to be replaced with spaces. The number
of spaces to use is set in the :ref:`basic settings` section below.

.. code-block:: scheme

  ;; Don't use TABS for indentations.
  (setq-default indent-tabs-mode nil)

Maximum number of characters in a line
======================================

Emacs will automatically insert a new line after "fill-column" number
of columns. PEP8 specifies a maximum of 79, but this can be set to a
smaller value also, for example 72.

.. code-block:: scheme

  ;; Set the number to the number of columns to use.
  (setq-default fill-column 79)

  ;; Add Autofill mode to mode hooks.
  (add-hook 'text-mode-hook 'turn-on-auto-fill)

  ;; Show line number in the mode line.
  (line-number-mode 1)

  ;; Show column number in the mode line.
  (column-number-mode 1)

Syntax highlighting
===================

Enable syntax highlighting. This will also highlight lines that form a
region.

.. code-block:: scheme

  (global-font-lock-mode 1)

Python specific settings
************************

.. _`basic settings`:

Basic settings
==============

Indentation is automatically added. When a tab is pressed it is
replaced with 4 spaces. When backspace is pressed on an empty line, the
cursor will jump to the previous indentation level.

.. code-block:: scheme

  (load-library "python")

  (autoload 'python-mode "python-mode" "Python Mode." t)
  (add-to-list 'auto-mode-alist '("\\.py\\'" . python-mode))
  (add-to-list 'interpreter-mode-alist '("python" . python-mode))

  (setq interpreter-mode-alist
        (cons '("python" . python-mode)
              interpreter-mode-alist)
        python-mode-hook
        '(lambda () (progn
                      (set-variable 'py-indent-offset 4)
                      (set-variable 'indent-tabs-mode nil))))

Highlight the column where a line must stop
===========================================

The "fill-column" column is highlighted in red. For this to work,
download `column-marker.el
<http://www.emacswiki.org/emacs/column-marker.el>`_ and place it in the
Emacs configuration directory.

.. code-block:: scheme

  ;; Highlight character at "fill-column" position.
  (require 'column-marker)
  (set-face-background 'column-marker-1 "red")
  (add-hook 'python-mode-hook
            (lambda () (interactive)
              (column-marker-1 fill-column)))

Flycheck
========

One can make lines that do not satisfy syntax requirements using
flycheck_. When cursor is on such a line a message is displayed in the
mini-buffer. When mouse pointer is on such a line a "tool tip" message
is also shown. By default, flycheck_ will check if flake8_ is
installed and, if so, use that for its syntax checking. To ensure
flycheck_ starts upon opening python files, add:

.. code-block:: scheme
  (add-hook 'python-mode-hook 'flycheck-mode)

Alternatively, you can just use ``(global-flycheck-mode)`` to run flycheck
for all languages it supports.

Delete trailing white spaces and blank lines
============================================

To manually delete trailing whitespaces, press ``C-t C-w``, which will run
the command "delete-whitespaces`. This command is also run when a file is
saved, and hence all trailing whitespaces will be deleted on saving a Python
file.

To make sure that all "words" are separated by only one space, type
``M-SPC`` (use the ALT key since ``M-SPC`` sometimes brings up a context
menu.).

To collapse a set of blank lines to one blank line, place the cursor on one
of these and press ``C-x C-o``. This is useful for deleting multiple black
lines at the end of a file.

.. code-block:: scheme

  ;; Remove trailing whitespace manually by typing C-t C-w.
  (add-hook 'python-mode-hook
            (lambda ()
              (local-set-key (kbd "C-t C-w")
                             'delete-trailing-whitespace)))

  ;; Automatically remove trailing whitespace when file is saved.
  (add-hook 'python-mode-hook
        (lambda()
          (add-hook 'local-write-file-hooks
                '(lambda()
                   (save-excursion
                     (delete-trailing-whitespace))))))

  ;; Use M-SPC (use ALT key) to make sure that words are separated by
  ;; just one space. Use C-x C-o to collapse a set of empty lines
  ;; around the cursor to one empty line. Useful for deleting all but
  ;; one blank line at end of file. To do this go to end of file (M->)
  ;; and type C-x C-o.

..  LocalWords:  whitespaces
