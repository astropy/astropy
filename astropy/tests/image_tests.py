try:
    import matplotlib
    from matplotlib import pyplot as plt  # noqa
except ImportError:
    MPL_VERSION = ''
    ROOT = ''
    IMAGE_REFERENCE_DIR = ''
else:
    MPL_VERSION = matplotlib.__version__

    # The developer versions of the form 3.2.x+... contain changes that will only
    # be included in the 3.3.x release, so we update this here.
    if MPL_VERSION[:3] == '3.2' and '+' in MPL_VERSION:
        MPL_VERSION = '3.3'

    ROOT = "http://{server}/testing/astropy/2021-08-25T18:18:36.000000/{mpl_version}/"
    IMAGE_REFERENCE_DIR = (
        ROOT.format(server='data.astropy.org', mpl_version=MPL_VERSION[:3] + '.x') + ',' +
        ROOT.format(server='www.astropy.org/astropy-data', mpl_version=MPL_VERSION[:3] + '.x'))
