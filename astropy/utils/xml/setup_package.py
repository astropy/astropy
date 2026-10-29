# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import sys
from collections import defaultdict
from os.path import join

from setuptools import Extension

from extension_helpers import pkg_config


def get_extensions(build_type="release"):
    XML_DIR = "astropy/utils/xml/src"

    cfg = defaultdict(list)
    cfg["sources"] = [join(XML_DIR, "iterparse.c")]

    if int(os.environ.get("ASTROPY_USE_SYSTEM_EXPAT", "0")) or int(
        os.environ.get("ASTROPY_USE_SYSTEM_ALL", "0")
    ):
        for k, v in pkg_config(["expat"], ["expat"]).items():
            cfg[k].extend(v)
    else:
        EXPAT_DIR = "cextern/expat/lib"
        expat_sources = [
            "xcs.c",
            "xmlparse.c",
            "xmlrole.c",
            "xmltok.c",
            "xmltok_impl.c",
        ]
        # Source of entropy for the hash salt, matching expat_config.h
        if sys.platform == "win32":
            expat_sources.append("random_rand_s.c")
        else:
            expat_sources.append("random_dev_urandom.c")
        cfg["sources"].extend([join(EXPAT_DIR, fn) for fn in expat_sources])
        cfg["include_dirs"].extend([XML_DIR, EXPAT_DIR])
        if sys.platform.startswith("linux"):
            # This is to ensure we only export the Python entry point
            # symbols and the linker won't try to use the system expat in
            # place of ours.
            cfg["extra_link_args"].extend(
                [f"-Wl,--version-script={join(XML_DIR, 'iterparse.map')}"]
            )
        cfg["define_macros"].append(("HAVE_EXPAT_CONFIG_H", 1))
        if sys.byteorder == "big":
            cfg["define_macros"].append(("BYTEORDER", "4321"))
        else:
            cfg["define_macros"].append(("BYTEORDER", "1234"))
        if sys.platform != "win32":
            cfg["define_macros"].append(("HAVE_UNISTD_H", None))

    return [Extension("astropy.utils.xml._iterparser", **cfg)]
