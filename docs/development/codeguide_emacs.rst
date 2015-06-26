=============================================
 Emacs setup for following coding guidelines
=============================================

.. _flymake: http://www.emacswiki.org/emacs/FlyMake
.. _pyflakes: http://pypi.python.org/pypi/pyflakes
.. _pep8: http://pypi.python.org/pypi/pep8
.. include:: workflow/known_projects.inc

The Astropy coding guidelines are listed in :doc:`codeguide`. This
document will describe some configuration options for Emacs, that will
help in ensuring that Python code satisfies the guidelines. Emacs can
be configured in several different ways. So instead of providing a drop
in configuration file, only the individual configurations are presented
below.

For this setup we will need flymake_, pyflakes_ and the pep8_ Python
script, in addition to ``python-mode``.

Flymake comes with Emacs 23. The rest can be obtained from their websites,
or can be installed using `pip`_.

Global settings
===============

No tabs
-------

This setting will cause all tabs to be replaced with spaces. The number
of spaces to use is set in the :ref:`basic settings` section below.

.. code-block:: scheme

  ;; Don't use TABS for indentations.
  (setq-default indent-tabs-mode nil)

Maximum number of characters in a line
--------------------------------------

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
-------------------

Enable syntax highlighting. This will also highlight lines that form a
region.

.. code-block:: scheme

  (global-font-lock-mode 1)

Python specific settings
========================

.. _`basic settings`:

Basic settings
--------------

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
-------------------------------------------

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

Flymake
-------

Flymake will mark lines that do not satisfy syntax requirements in
red. When cursor is on such a line a message is displayed in the
mini-buffer. When mouse pointer is on such a line a "tool tip" message
is also shown.

For flymake to work with `pep8`_ and `pyflakes`_, create an
executable file named `pychecker`_ with the following contents. This
file must be in the system path.

.. code-block:: sh

  #!/bin/bash

  pyflakes "$1"
  pep8 --ignore=E221,E701,E202 --repeat "$1"
  true

Also download `flymake-cursor.el
<http://www.emacswiki.org/emacs/flymake-cursor.el>`_ and place it in the
Emacs configuration directory.  Then add the following code to the Emacs
configuration:

.. code-block:: scheme

  ;; Setup for Flymake code checking.
  (require 'flymake)
  (load-library "flymake-cursor")

  ;; Script that flymake uses to check code. This script must be
  ;; present in the system path.
  (setq pycodechecker "pychecker")

  (when (load "flymake" t)
    (defun flymake-pycodecheck-init ()
      (let* ((temp-file (flymake-init-create-temp-buffer-copy
                         'flymake-create-temp-inplace))
             (local-file (file-relative-name
                          temp-file
                          (file-name-directory buffer-file-name))))
        (list pycodechecker (list local-file))))
    (add-to-list 'flymake-allowed-file-name-masks
                 '("\\.py\\'" flymake-pycodecheck-init)))

  (add-hook 'python-mode-hook 'flymake-mode)

.. note::

    Flymake will save files with suffix *_flymake* in the current
    directory. If it crashes for some reason, then these files will not
    get deleted.

    Sometimes there is a delay in refreshing the results.

Delete trailing white spaces and blank lines
--------------------------------------------

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
