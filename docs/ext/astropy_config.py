import io

from docutils import nodes
from docutils.parsers.rst import Directive

from astropy.config import generate_config


class GenerateConfig(Directive):
    """
    Directive to generate the configuration file for a package and
    include it in the documentation as a literal code block.
    """

    has_content = False
    required_arguments = 1

    def run(self):
        buf = io.StringIO()
        generate_config(pkgname=self.arguments[0], filename=buf)
        text = buf.getvalue()
        node = nodes.literal_block(text, text)
        return [node]


def setup(app):
    app.add_directive("generate_config", GenerateConfig)

    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
