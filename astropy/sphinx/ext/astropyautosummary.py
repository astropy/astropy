# **Please Note**: ``astropy.sphinx`` exists only for backward-compatibility
# purposes - it has now been moved to the separate astropy-helpers package,
# located at https://github.com/astropy/astropy-helpers. Any new development or
# bug fixes should be done there.
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This sphinx extension builds off of `sphinx.ext.autosummary` to
clean up some issues it presents in the Astropy docs.

The main issue this fixes is the summary tables getting cut off before the
end of the sentence in some cases.

"""
import re

from sphinx.ext.autosummary import Autosummary

# used in AstropyAutosummary.get_items
_itemsummrex = re.compile(r'^([A-Z].*?\.(?:\s|$))')


class AstropyAutosummary(Autosummary):
    def get_items(self, names):
        """Try to import the given names, and return a list of
        ``[(name, signature, summary_string, real_name), ...]``.
        """
        from sphinx.ext.autosummary import (get_import_prefixes_from_env,
            import_by_name, get_documenter, mangle_signature)

        env = self.state.document.settings.env

        prefixes = get_import_prefixes_from_env(env)

        items = []

        max_item_chars = 50

        for name in names:
            display_name = name
            if name.startswith('~'):
                name = name[1:]
                display_name = name.split('.')[-1]

            try:
                real_name, obj, parent = import_by_name(name, prefixes=prefixes)
            except ImportError:
                self.warn('[astropyautosummary] failed to import %s' % name)
                items.append((name, '', '', name))
                continue

            # NB. using real_name here is important, since Documenters
            #     handle module prefixes slightly differently
            documenter = get_documenter(obj, parent)(self, real_name)
            if not documenter.parse_name():
                self.warn('[astropyautosummary] failed to parse name %s' % real_name)
                items.append((display_name, '', '', real_name))
                continue
            if not documenter.import_object():
                self.warn('[astropyautosummary] failed to import object %s' % real_name)
                items.append((display_name, '', '', real_name))
                continue

            # -- Grab the signature

            sig = documenter.format_signature()
            if not sig:
                sig = ''
            else:
                max_chars = max(10, max_item_chars - len(display_name))
                sig = mangle_signature(sig, max_chars=max_chars)
                sig = sig.replace('*', r'\*')

            # -- Grab the summary

            doc = list(documenter.process_doc(documenter.get_doc()))

            while doc and not doc[0].strip():
                doc.pop(0)
            m = _itemsummrex.search(" ".join(doc).strip())
            if m:
                summary = m.group(1).strip()
            elif doc:
                summary = doc[0].strip()
            else:
                summary = ''

            items.append((display_name, sig, summary, real_name))

        return items


def setup(app):
    # need autosummary, of course
    app.setup_extension('sphinx.ext.autosummary')
    # this replaces the default autosummary with the astropy one
    app.add_directive('autosummary', AstropyAutosummary)
