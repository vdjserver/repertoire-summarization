"""
Repertoire calculations
"""
# Imports
import os
import sys
import versioneer

try:
    from setuptools import setup, find_packages
except ImportError:
    sys.exit('Please install setuptools before installing repcalc.\n')

try:
    from pip.req import parse_requirements
except ImportError:
    sys.exit('Please install pip before installing repcalc.\n')

with open('README.rst', 'r') as ip:
    long_description = ip.read()

# Parse requirements
if os.environ.get('READTHEDOCS', None) == 'True':
    # Set empty install_requires to get install to work on readthedocs
    install_requires = []
else:
    with open('requirements.txt') as req:
        install_requires = req.read().splitlines()

# Setup
setup(name='repcalc',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      author='VDJServer Team',
      author_email='vdjserver@utsouthwestern.edu',
      description='Summarization functions for immune repertoire sequencing data.',
      long_description=long_description,
      zip_safe=False,
      license='GNU General Public License v3.0',
      url='https://vdjserver.org',
      download_url='https://bitbucket.org/vdjserver/repertoire-summarization/downloads',
      keywords=['bioinformatics', 'sequencing', 'immunoglobulin', 'antibody', 'lymphocyte',
                'adaptive immunity', 'T cell', 'B cell', 'BCR', 'TCR', 'CDR3'],
      install_requires=install_requires,
      packages=find_packages(),
      entry_points={
                'console_scripts': [
                'repcalc=repcalc.repcalc:main',
                'repcalc_create_config=repcalc.repcalc:create_config',
                'repcalc_group_map=repcalc.repcalc:generate_group_map',
                ]},
      classifiers=['Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'Natural Language :: English',
                   'Operating System :: OS Independent',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3',
                   'Topic :: Scientific/Engineering :: Bio-Informatics'])
