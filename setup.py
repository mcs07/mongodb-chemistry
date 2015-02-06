#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

setup(
    name='mongodb-chemistry',
    version='0.0.1',
    author='Matt Swain',
    author_email='m.swain@me.com',
    license='MIT',
    url='https://github.com/mcs07/mongodb-chemistry',
    packages=['mchem'],
    description='Proof of concept for a MongoDB chemical database',
    keywords='chemistry cheminformatics rdkit',
    zip_safe=False,
    entry_points={'console_scripts': ['mchem = mchem.cli:cli', 'pgchem = mchem.postgres:cli']},
    install_requires=['Click', 'pymongo'],
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
)
