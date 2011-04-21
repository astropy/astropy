from distutils.core import setup

long_description = """
Backport of Python 2.7 collections.OrderedDict class which works on 2.6.
Only changes from baseline 2.7 code:

Add:
from collections import MutableMapping

Delete:
from _collections import deque, defaultdict
"""

setup(name='odict',
      version=1.0,
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
