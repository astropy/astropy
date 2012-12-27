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
