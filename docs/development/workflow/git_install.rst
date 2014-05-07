:orphan:

.. include:: links.inc
.. _install-git:

==========================
 Install and configure git
==========================


Get git
-------

Installers and instructions for all platforms are available at
http://git-scm.com/downloads

.. _essential_config:

Essential configuration
-----------------------

Though technically not required to install `git`_ and get it running, configure `git`_ so that you get credit for your contributions::

    git config --global user.name "Your Name"
    git config --global user.email you@yourdomain.example.com

.. note::
    Use the same email address here that you used for setting up your GitHub
    account to save yourself a couple of steps later, when you connect your
    git to GitHub.

Check it with::

    $ git config --list
    user.name=Your Name
    user.email=you@yourdomain.example.com
    # ...likely followed by many other configuration values

.. _git_gui_options:

Get a git GUI (optional)
------------------------

There are several good, free graphical interfaces for git.
Even if you are proficient with `git`_ at the command line a GUI can be useful.

Mac and Windows:

+ `SourceTree`_
+ The github client for `Mac`_ or `Windows`_

Linux, Mac and Windows:

+ `git-cola`_

There is a more extensive list of `git GUIs`_, including non-free options, for
all platforms.

.. _git GUIs: http://git-scm.com/downloads/guis
.. _SourceTree: http://www.sourcetreeapp.com/
.. _Mac: http://mac.github.com/
.. _Windows: http://windows.github.com/
.. _git-cola: http://git-cola.github.io/
