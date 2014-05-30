# **Please Note**: ``astropy.sphinx`` exists only for backward-compatibility
# purposes - it has now been moved to the separate astropy-helpers package,
# located at https://github.com/astropy/astropy-helpers. Any new development or
# bug fixes should be done there.
from sphinx import addnodes


def fix_toc_entries(app, doctree):
    # Get the docname; I don't know why this isn't just passed in to the
    # callback
    # This seems a bit unreliable as it's undocumented, but it's not "private"
    # either:
    docname = app.builder.env.temp_data['docname']
    if app.builder.env.metadata[docname].get('tocdepth', 0) != 0:
        # We need to reprocess any TOC nodes in the doctree and make sure all
        # the files listed in any TOCs are noted
        for treenode in doctree.traverse(addnodes.toctree):
            app.builder.env.note_toctree(docname, treenode)


def setup(app):
    app.connect('doctree-read', fix_toc_entries)
