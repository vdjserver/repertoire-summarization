#!/usr/bin/env python
"""
Repertoire-summarization setup
"""
# Imports
import os
import sys

try:
    from setuptools import setup
except ImportError:
    sys.exit('Please install setuptools before installing repsum.\n')

try:
    from pip.req import parse_requirements
except ImportError:
    sys.exit('Please install pip before installing repsum.\n')

# Get version, author and license information
info_file = os.path.join('repsum', 'version.py')
__version__, __author__, __license__ = None, None, None
try:
    exec(open(info_file).read())
except:
    sys.exit('Failed to load package information from %s.\n' % info_file)

if __version__ is None:
    sys.exit('Missing version information in %s\n.' % info_file)
if __author__ is None:
    sys.exit('Missing author information in %s\n.' % info_file)
if __license__ is None:
    sys.exit('Missing license information in %s\n.' % info_file)

# Load long package description
#desc_files = ['README.rst', 'INSTALL.rst', 'NEWS.rst']
desc_files = ['README.rst']
long_description = '\n\n'.join([open(f, 'r').read() for f in desc_files])

# TODO: check pip version to avoid problem with parse_requirements(session=False)
# Parse requirements
if os.environ.get('READTHEDOCS', None) == 'True':
    # Set empty install_requires to get install to work on readthedocs
    install_requires = []
else:
    require_file = 'requirements.txt'
    try:
        requirements = parse_requirements(require_file, session=False)
    except TypeError:
        requirements = parse_requirements(require_file)
    install_requires = [str(r.req) for r in requirements]

# Setup
setup(name='repsum',
      version=__version__,
      author=__author__,
      author_email=__email__,
      description='Summarization functions for immune repertoire sequencing data.',
      long_description=long_description,
      zip_safe=False,
      license=__license__,
      url='https://vdjserver.org',
      download_url='https://bitbucket.org/vdjserver/repertoire-summarization/downloads',
      keywords='bioinformatics immunoglobulin lymphocyte sequencing TCR CDR3',
      install_requires=install_requires,
      packages=['repsum'],
      package_dir={'repsum': 'repsum'},
      entry_points={
		'console_scripts': [
        	'repsum=repsum:summary',
        	'repcalc=repsum:main',
    	],
	  },
      classifiers=['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'Natural Language :: English',
                   'Operating System :: OS Independent',
                   'Programming Language :: Python :: 2.7',
                   'Topic :: Scientific/Engineering :: Bio-Informatics'])
