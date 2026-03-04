.. Hide the left hand sidebar on the home page because it's empty

.. raw:: html

    <style>

    .bd-sidebar-secondary {
      display: none;
    }

    .bd-header label.secondary-toggle {
      display: none;
    }
    </style>

:tocdepth: 3

.. _astropy-docs-index:

#################################################
astropy: A Community Python Library for Astronomy
#################################################

**Version**: |release| - :ref:`whatsnew-8.0`

**Useful links**:
:ref:`Installation <installing-astropy>` |
`Issues & Ideas <https://github.com/astropy/astropy/issues>`__ |
:ref:`astropy-org-help` |
:ref:`astropy-org-contribute` |
:ref:`astropy-org-about`

The ``astropy`` package contains key functionality and common tools needed for
performing astronomy and astrophysics with Python.  It is at the core of the
:ref:`Astropy Project <astropy-org-about>`, which aims to enable
the community to develop a robust ecosystem of :ref:`astropy-org-affiliated`
covering a broad range of needs for astronomical research, data
processing, and data analysis.

.. Important:: If you use Astropy for work presented in a publication
   or talk please help the project via proper
   :ref:`astropy-org-acknowledge`.  This also applies to use of
   software or :ref:`astropy-org-affiliated` that depend on the astropy core
   package.

.. toctree::
   :maxdepth: 1
   :hidden:

   index_getting_started
   index_user_docs
   index_dev
   index_project_details

.. grid:: 2
    :gutter: 2

    .. grid-item-card:: Getting Started
        :link: index_getting_started
        :link-type: doc
        :text-align: center

        :material-outlined:`directions_run;8em;sd-text-secondary`

        New to Astropy? Check out the getting started guides. They contain an
        introduction to astropy's main concepts and links to additional tutorials.

    .. grid-item-card:: User Guide
        :link: index_user_docs
        :link-type: doc
        :text-align: center

        :material-outlined:`menu_book;8em;sd-text-secondary`

        The user guide provides in-depth information on the key concepts
        of astropy with useful background information and explanation.

    .. grid-item-card:: Learn Astropy
        :link: https://learn.astropy.org
        :link-type: url
        :text-align: center

        :material-outlined:`psychology;8em;sd-text-secondary`

        Learn how to use Python for astronomy through tutorials and guides that cover
        Astropy and other packages in the astronomy Python ecosystem.

    .. grid-item-card:: Astropy Packages
        :link: astropy-org-affiliated
        :link-type: ref
        :text-align: center

        :material-outlined:`inventory_2;8em;sd-text-secondary`

        The Astropy Project ecosystem includes numerous `Coordinated
        <https://www.astropy.org/affiliated/#coordinated-packages>`_ and `Affiliated
        <https://www.astropy.org/affiliated/#affiliated-packages>`_ packages.
        Coordinated packages are maintained by the Project.

    .. grid-item-card:: Contributor's Guide
        :link: index_dev
        :link-type: doc
        :text-align: center

        :material-outlined:`person_add;8em;sd-text-secondary`

        Saw a typo in the documentation? Want to improve
        existing functionalities? The contributing guidelines will show
        you how to improve astropy.

    .. grid-item-card:: Project Details
        :link: index_project_details
        :link-type: doc
        :text-align: center

        :material-outlined:`more_horiz;8em;sd-text-secondary`

        What's new in the latest release, changelog, and other project details.

.. image:: https://github.com/astropy/repo_stats/blob/cache/cache/astropy_user_stats_light.png?raw=true
    :class: only-light
    :target: https://docs.astropy.org/en/latest/impact_health.html
    :alt: Astropy User Statistics

.. image:: https://github.com/astropy/repo_stats/blob/cache/cache/astropy_user_stats_dark.png?raw=true
    :class: only-dark
    :target: https://docs.astropy.org/en/latest/impact_health.html
    :alt: Astropy User Statistics

.. _feedback@astropy.org: mailto:feedback@astropy.org
.. _affiliated packages: https://www.astropy.org/affiliated/
