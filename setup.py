from distutils.core import setup

long_description = """
Backport of Python 2.7 collections.OrderedDict class which works on 2.6.  See
https://github.com/taldcroft/odict for the (minimal) changes from original code.

Testing has been minimal but it appears to work.  Caveat emptor.
"""

setup(name='odict',
      version='1.0',
      description='Python 2.6 backport of Python 2.7 collections.OrderedDict',
      long_description=long_description,
      author='Tom Aldcroft',
      author_email='aldcroft@head.cfa.harvard.edu',
      license='BSD',
      platforms=['any'],
      classifiers=[
          'Development Status :: 4 - Beta',
          'License :: OSI Approved :: BSD License',
          'Programming Language :: Python :: 2',
          ],
      py_modules=['odict'],
      )
