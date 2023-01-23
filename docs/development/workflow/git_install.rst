:orphan:

.. include:: links.inc
.. _install-git:

**************************
 Install and configure git
**************************


Get git
=======

Installers and instructions for all platforms are available at
https://git-scm.com/downloads

.. _essential_config:

Essential configuration
=======================

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
