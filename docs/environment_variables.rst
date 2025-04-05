.. currentmodule:: astropy

.. _environment_variables:

*********************
Environment variables
*********************


.. glossary::

    ``XDG_CONFIG_HOME``
        This environment variables control where configuration files are read
        and written. Astropy will look for configuration files in
        ``$XDG_CONFIG_HOME/astropy``.
        If not set, or if ``$XDG_CONFIG_HOME/astropy`` does not exist,
        astropy will default to ``$HOME/.astropy/config``.
        See :ref:`astropy_config` for how to programmatically set or get the
        location of the corresponding directory at runtime.

    ``XDG_CACHE_HOME``
        These environment variables control where data files are cached.
        Astropy will cache files in ``$XDG_CACHE_HOME/astropy``.
        If not set, or if ``$XDG_CACHE_HOME/astropy`` does not exist,
        astropy will default to ``$HOME/.astropy/cache``.
        See :ref:`utils-data` for how to programmatically set or get the
        location of the corresponding directory at runtime.


.. note::
    ``XDG_CONFIG_HOME`` and ``XDG_CACHE_HOME`` come from a Linux-centric
    specification
    (see `here <https://wiki.archlinux.org/index.php/XDG_Base_Directory_support>`_
    for more details), but ``astropy`` will use them on any OS as a more
    general means to know where user-specific configurations should be written.
