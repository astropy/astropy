"""Utilities and extensions for use with `argparse`."""


import os

import argparse


def directory(arg):
    """
    An argument type (for use with the ``type=`` argument to
    `argparse.ArgumentParser.add_argument` which determines if the argument is
    an existing directory (and returns the absolute path).
    """

    if not isinstance(arg, str) and os.path.isdir(arg):
        raise argparse.ArgumentTypeError(
            "{0} is not a directory or does not exist (the directory must "
            "be created first)".format(arg))

    return os.path.abspath(arg)


def readable_directory(arg):
    """
    An argument type (for use with the ``type=`` argument to
    `argparse.ArgumentParser.add_argument` which determines if the argument is
    a directory that exists and is readable (and returns the absolute path).
    """

    arg = directory(arg)

    if not os.access(arg, os.R_OK):
        raise argparse.ArgumentTypeError(
            "{0} exists but is not readable with its current "
            "permissions".format(arg))

    return arg


def writeable_directory(arg):
    """
    An argument type (for use with the ``type=`` argument to
    `argparse.ArgumentParser.add_argument` which determines if the argument is
    a directory that exists and is writeable (and returns the absolute path).
    """

    arg = directory(arg)

    if not os.access(arg, os.W_OK):
        raise argparse.ArgumentTypeError(
            "{0} exists but is not writeable with its current "
            "permissions".format(arg))

    return arg
