:orphan:

.. _basic-workflow:

****************
Creating patches
****************

Overview
========

If you haven't already configured git::

    git config --global user.name "Your Name"
    git config --global user.email you@yourdomain.example.com

Then, the workflow is the following::

   # Get the repository if you don't have it
   git clone --recursive git://github.com/astropy/astropy.git

   # Make a branch for your patching
   cd astropy
   git branch the-fix-im-thinking-of
   git checkout the-fix-im-thinking-of

   # hack, hack, hack

   # Tell git about any new files you've made
   git add somewhere/tests/test_my_bug.py

   # Commit work in progress as you go
   git commit -am 'BF - added tests for Funny bug'

   # hack hack, hack

   # Commit work
   git commit -am 'BF - added fix for Funny bug'

   # Make the patch files
   git format-patch -M -C master

Then, send the generated patch files to the `astropy-dev mailing list`_ |emdash|
where we will thank you warmly.

In detail
=========

#. Tell git who you are so it can label the commits you've
   made::

    git config --global user.name "Your Name"
    git config --global user.email you@yourdomain.example.com

   This is only necessary if you haven't already done this, and you haven't
   checked to :ref:`check_git_install`.

#. If you don't already have one, clone a copy of the
   Astropy_ repository::

      git clone --recursive git://github.com/astropy/astropy.git
      cd astropy

#. Make a 'feature branch'. This will be where you work on your bug fix. It's
   nice and safe and leaves you with access to an unmodified copy of the code
   in the main branch::

      git branch the-fix-im-thinking-of
      git checkout the-fix-im-thinking-of

#. Do some edits, and commit them as you go::

      # hack, hack, hack

      # Tell git about any new files you've made
      git add somewhere/tests/test_my_bug.py

      # Commit work in progress as you go
      git commit -am 'BF - added tests for Funny bug'

      # hack hack, hack

      # Commit work
      git commit -am 'BF - added fix for Funny bug'

   Note the ``-am`` options to ``commit``. The ``m`` flag just
   signals that you're going to type a message on the command
   line.  The ``a`` flag |emdash| you can just take on faith |emdash|
   or see `why the -a flag?`_.

#. When you have finished, check you have committed all your changes::

      git status

#. Finally, make your commits into patches. You want all the commits since you
   branched from the ``master`` branch::

      git format-patch -M -C master

   You will now have several files named for the commits::

      0001-BF-added-tests-for-Funny-bug.patch
      0002-BF-added-fix-for-Funny-bug.patch

   Send these files to the `astropy-dev mailing list`_.

When you are done, to switch back to the main copy of the
code, just return to the ``master`` branch::

   git checkout master

.. include:: links.inc
