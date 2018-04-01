# Licensed under a 3-clause BSD style license
import os

def get_package_data():
    # Installs the schema files
    schemas = []
    root = os.path.join(os.path.dirname(__file__), 'schemas')
    for node, dirs, files in os.walk(root):
        for fname in files:
            if fname.endswith('.yaml'):
                schemas.append(
                    os.path.relpath(os.path.join(node, fname), root))

    # In the package directory, install to the subdirectory 'schemas'
    schemas = [os.path.join('schemas', s) for s in schemas]

    return {'astropy.io.misc.asdf': schemas}
