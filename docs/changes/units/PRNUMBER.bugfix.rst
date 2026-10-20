The ``cds`` format now writes units that share a name with a different unit in
the CDS standard (e.g. gauss, ``u.G``, which is ``gauss`` in CDS since ``G`` is
the gravitational constant there) such that they are read back as the same
unit. Previously such units were silently changed, e.g.
``u.G.to_string("cds")`` gave ``"G"``, the gravitational constant in CDS.
