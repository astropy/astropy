# Licensed under a 3-clause BSD style license - see LICENSE.rst
import warnings

import numpy as np


def run_schema_example_test(organization, standard, name, version, check_func=None):
    import asdf
    from asdf.exceptions import AsdfDeprecationWarning
    from asdf.schema import load_schema

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            category=AsdfDeprecationWarning,
            message=r"asdf.tests.helpers is deprecated.*",
        )
        from asdf.tests.helpers import yaml_to_asdf

        warnings.filterwarnings(
            "ignore",
            category=AsdfDeprecationWarning,
            message=r"asdf.types.*is deprecated.*",
        )
        from asdf.types import format_tag

        tag = format_tag(organization, standard, version, name)

        warnings.filterwarnings(
            "ignore",
            category=AsdfDeprecationWarning,
            message=r"default_extensions is deprecated.*",
        )
        uri = asdf.extension.default_extensions.extension_list.tag_mapping(tag)

        warnings.filterwarnings(
            "ignore",
            category=AsdfDeprecationWarning,
            message=r"get_default_resolver is deprecated.*",
        )
        r = asdf.extension.get_default_resolver()

    examples = []
    schema = load_schema(uri, resolver=r)
    for node in asdf.treeutil.iter_tree(schema):
        if (
            isinstance(node, dict)
            and "examples" in node
            and isinstance(node["examples"], list)
        ):
            for _, example in node["examples"]:
                examples.append(example)

    for example in examples:
        buff = yaml_to_asdf("example: " + example.strip())
        ff = asdf.AsdfFile(uri=uri)
        # Add some dummy blocks so that the ndarray examples work
        for i in range(3):
            b = asdf.block.Block(np.zeros((1024 * 1024 * 8), dtype=np.uint8))
            b._used = True
            ff.blocks.add(b)
        ff._open_impl(ff, buff, mode="r")
        if check_func:
            check_func(ff)
