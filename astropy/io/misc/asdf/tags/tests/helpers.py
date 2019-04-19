# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import os
import urllib.parse

import yaml

import numpy as np


def run_schema_example_test(organization, standard, name, version, check_func=None):

    import asdf
    from asdf.tests import helpers
    from asdf.types import format_tag
    from asdf.resolver import default_resolver

    tag = format_tag(organization, standard, version, name)
    schema_path = urllib.parse.urlparse(default_resolver(tag)).path

    with open(schema_path, 'rb') as ff:
        schema = yaml.load(ff)

    examples = []
    for node in asdf.treeutil.iter_tree(schema):
        if (isinstance(node, dict) and
            'examples' in node and
            isinstance(node['examples'], list)):
            for desc, example in node['examples']:
                examples.append(example)

    for example in examples:
        buff = helpers.yaml_to_asdf('example: ' + example.strip())
        ff = asdf.AsdfFile(uri=schema_path)
        # Add some dummy blocks so that the ndarray examples work
        for i in range(3):
            b = asdf.block.Block(np.zeros((1024*1024*8), dtype=np.uint8))
            b._used = True
            ff.blocks.add(b)
        ff._open_impl(ff, buff, mode='r')
        if check_func:
            check_func(ff)
