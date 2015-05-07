# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Command-line script and related functions for reading and installing affiliated
packages.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import sys
import json
import urllib2


DEFAULT_AFFILIATED_REGISTRY = 'http://affiliated.astropy.org/registry.json'


def get_affiliated_json(url=DEFAULT_AFFILIATED_REGISTRY):
    u = urllib2.urlopen(url)
    try:
        return json.loads(u.read())
    finally:
        u.close()


def get_names(regjson=None):
    if regjson is None:
        regjson = get_affiliated_json()

    return [pkg['name'] for pkg in regjson['packages']]


def _do_install(installcmd, pkgname, regjson):
    if regjson is None:
        regjson = get_affiliated_json()

    if pkgname == 'all':
        for pkg in regjson['packages']:
            print("Installing affiliated package "+pkg['name'])

            retcode = os.system('{0} {1}'.format(installcmd, pkg['pypi_name']))
            if retcode != 0:
                print("Failed on installing affiliated package "+pkg['name'])
                return retcode
        return 0
    else:
        for pkg in regjson['packages']:
            if pkg['name'] == pkgname:
                pypiname = pkg['pypi_name']
                break
        else:
            print('Could not find affiliated package "{0}"'.format(pkgname))
            sys.exit(1)

        return os.system('{0} {1}'.format(installcmd, pypiname))


def pip_install(pkgname, regjson=None):
    sys.exit(_do_install('pip install', pkgname, regjson))


def conda_install(pkgname, regjson=None):
    # implicitly this assumes pypi names=conda names... but usually true?
    sys.exit(_do_install('conda install', pkgname, regjson))


def main(args=None):
    from .compat import argparse

    parser = argparse.ArgumentParser(
        description='Operations with astropy affiliated packages.')
    parser.add_argument('--url', default=DEFAULT_AFFILIATED_REGISTRY,
                        help='The URL to use for the affiliated package json '
                             'registry')

    subparsers = parser.add_subparsers(dest='subparser_name')

    parser_list = subparsers.add_parser('list')

    parser_install = subparsers.add_parser('install')
    parser_install.add_argument('pkgname', help='The affiliated package to '
                                                'install or "all" to install '
                                                'all affiliated packages')

    parser_cinstall = subparsers.add_parser('condainstall')
    parser_cinstall.add_argument('pkgname', help='The affiliated package to '
                                                 'install or "all" to install '
                                                 'all affiliated packages')

    parser_pinstall = subparsers.add_parser('pipinstall')
    parser_pinstall.add_argument('pkgname', help='The affiliated package to '
                                                 'install or "all" to install '
                                                 'all affiliated packages')

    args = parser.parse_args(args)

    regjson = get_affiliated_json(args.url)

    if args.subparser_name == 'list':
        names = get_names(regjson=regjson)
        print(names)
    elif args.subparser_name == 'install':
        if 'Anaconda' in sys.version or 'Continuum Analytics' in sys.version:
            conda_install(args.pkgname, regjson=regjson)
        else:
            pip_install(args.pkgname, regjson=regjson)

    elif args.subparser_name == 'pipinstall':
        pip_install(args.pkgname, regjson=regjson)
    elif args.subparser_name == 'condainstall':
        conda_install(args.pkgname, regjson=regjson)
    else:
        print("Encountered unknown subcommand {0}".format(args.subparser_name))
        sys.exit(-1)
