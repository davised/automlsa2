# ################## Maintained by Hatch ####################
# This file is auto-generated by hatch. If you'd like to customize this file
# please add your changes near the bottom marked for 'USER OVERRIDES'.
# EVERYTHING ELSE WILL BE OVERWRITTEN by hatch.
# ###########################################################
from setuptools import find_packages, setup

with open('automlsa2/__version__.py', 'r') as f:
    for line in f:
        if line.startswith('__version__'):
            version = line.strip().split('=')[1].strip(' \'"')
            break
    else:
        version = '0.0.1'

with open('README.rst', 'r', encoding='utf-8') as f:
    readme = f.read()

REQUIRES = ['biopython >= 1.78']

kwargs = {
    'name': 'automlsa2',
    'version': version,
    'description': 'Automated Multi-Locus Sequence Analysis tool',
    'long_description': readme,
    'author': 'Edward Davis',
    'author_email': 'ed@cgrb.oregonstate.edu',
    'maintainer': 'Edward Davis',
    'maintainer_email': 'ed@cgrb.oregonstate.edu',
    'url': 'https://github.com/davised/automlsa2',
    'license': 'Custom',
    'classifiers': [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'License :: Free for non-commercial use',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: Implementation :: CPython',
    ],
    'install_requires': REQUIRES,
    'tests_require': ['coverage', 'pytest'],
    'packages': find_packages(exclude=('tests', 'tests.*')),

}

# ################## BEGIN USER OVERRIDES ####################
# Add your customizations in this section.

# #################### END USER OVERRIDES ####################

setup(**kwargs,
      entry_points={
          "console_scripts": ['automlsa2 = automlsa2.__main__:main']
      }
      )
