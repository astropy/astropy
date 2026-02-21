.. _managing-dependencies::

*********************
Managing Dependencies
*********************

Guidelines for defining requirements
====================================

All requirements (including optional and dev-only ones) must be defined
with an explicit lower version bound.

Upper bounds should be reserved as a last resort, only if a real incompatibility has
been identified that cannot be worked around with runtime branching

In order to compensate for missing lower bounds in indirect requirements,
transitive dependencies can also be *constrained* through ``[tool.uv]`` in
``pyproject.toml`` (see ``no-build-packages`` and ``constraint-dependencies``).
These constraints should generally be used as temporary workarounds and upstreamed
as early as possible.


Testing old versions of dependencies (``oldestdeps``)
=====================================================

As part of CI, we verify continuously that our minimal requirements are actually
compatible with the dev state of astropy, and each other.
This is done through the ``oldestdeps`` tox factor. This factor constrains
dependency resolution (the process of discovering a set of *exact* package versions that
satisfies requirements) to prefer oldest available versions of both direct and transitive
dependencies. By construction, the resulting solution is invariant under normal
conditions, which we enumerate here for completion:
- the resolver's algorithm does not change
- the package registry we use as a source (PyPI) does not change (i.e., packages in the
  final solution are neither yanked nor removed)
- our own requirements do not change

This last condition is the most likely reason for the ``oldestdeps`` solution to change.
We'll see how to handle this case in the following section.

Upgrading requirements
======================

In all generality, when a direct requirement is upgraded, as, for instance

.. patch::

    - "matplotlib>=3.5"
    + "matplotlib>=3.6"

then *its* (in)direct requirements (which, to ``astropy``, are *transitive* requirements)
may also change: most often, transitive lower bounds are simply upgraded, but it also
happens that new transitive requirements are added (or removed) as a consequence.
This process affects the set of constraints that dependency managers (like ``pip`` or
``uv``) need to resolve for, both directly and indirectly.
In turn, the ``oldestdeps`` solution, which we discussed in detail in the previous
section, will change, sometimes in ways unexpected by the person performing the upgrade.
In order to mitigate any undesirable side effects of upgrades, a new solution needs to
be compared to the previous one.
This is done using a helper script as follows

.. shell::

    uv run scripts/check-lowest-resolved-tree.py

which will print a diff if anything changed.
Provided the diff is acceptable, one may then upgrade the recorded solution as by
re-running the script with the ``--overwrite`` flag.
This script is also used in CI to prevent stealth upgrades.

Providing guidelines for *all* failure modes in this process is out of scope here.
Maintainers are expected to use their best judgement and collectively
decide what to do on a case by case basis.

However, most common (or obvious) *undesirable* effects include:
- upgrading a direct requirement cascades into upgrading (an)other, user-facing direct
  requirement(s). In this case, consider upgrading the other(s) requirement(s) as well.
- no solution is found. This should be extremely rare, but may be caused by abusive
  upper bounds throughout the dependency graph.
