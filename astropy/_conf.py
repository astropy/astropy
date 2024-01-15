import sys
from . import config as _config

__all__ = ["conf"]

class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy`.
    """

    unicode_output = _config.ConfigItem(
        False,
        "When True, use Unicode characters when outputting values, and "
        "displaying widgets at the console.",
    )
    use_color = _config.ConfigItem(
        sys.platform != "win32",
        "When True, use ANSI color escape sequences when writing to the console.",
        aliases=["astropy.utils.console.USE_COLOR", "astropy.logger.USE_COLOR"],
    )
    max_lines = _config.ConfigItem(
        None,
        description=(
            "Maximum number of lines in the display of pretty-printed "
            "objects. If not provided, try to determine automatically from the "
            "terminal size.  Negative numbers mean no limit."
        ),
        cfgtype="integer(default=None)",
        aliases=["astropy.table.pprint.max_lines"],
    )
    max_width = _config.ConfigItem(
        None,
        description=(
            "Maximum number of characters per line in the display of "
            "pretty-printed objects.  If not provided, try to determine "
            "automatically from the terminal size. Negative numbers mean no "
            "limit."
        ),
        cfgtype="integer(default=None)",
        aliases=["astropy.table.pprint.max_width"],
    )


conf = Conf()
