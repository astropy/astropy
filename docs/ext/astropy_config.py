import io

from docutils import nodes
from docutils.parsers.rst import Directive

from astropy.config import generate_config


class AstropyConfig(Directive):

    def run(self):
        buf = io.StringIO()
        generate_config(buf)
        text = buf.getvalue()
        node = nodes.literal_block(text, text)
        # self.add_name(node)
        return [node]


def setup(app):
    app.add_directive("astropy_config", AstropyConfig)

    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
