The unit and angle string parsers are now built with the `lark
<https://lark-parser.readthedocs.io>`_ parsing library, which is a new required
dependency, instead of the bundled copy of PLY, which is no longer maintained
and has been removed together with the pre-generated parser tables. Parsing
behaviour is unchanged, but the parsers now also work when Python is run with
``-OO`` or under tools such as Nuitka.
