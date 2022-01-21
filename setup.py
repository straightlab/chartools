# -*- coding: utf-8 -*-
import os
import os.path
from setuptools import setup, find_packages

classifiers = [
    'Programming Language :: Python :: 3',
]



# Load metadata from __about__.py
base_dir = os.path.dirname(__file__)
about = {}


with open(os.path.join(base_dir, 'chartools', '__about__.py')) as f:
    exec(f.read(), about)

install_requires = [] # ['numpy', 'matplotlib', ]

if __name__ == '__main__':
    setup(
        name = about['__distname__'],
        version = about['__version__'],
        packages = ['chartools'],
        author = about['__author__'],
        classifiers = classifiers,
        install_requires = install_requires,
        entry_points={'console_scripts': ['chartools = chartools.cli.main:main']}
    )