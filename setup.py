#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name='labkeyext-engtools',
    version='0.1.0',
    description='A collection of engineering, database, and miscellaneous dataframe tools.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/dasievers/engtools',
    author='David Sievers',
    packages=['engtools'],
    install_requires=[
                      'numpy',
                      'scipy',
                      'pandas',
                      ],

    classifiers=[
                'Development Status :: 3 - Alpha',
                'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
                'Intended Audience :: Science/Research',
                'Programming Language :: Python :: 3',
                ],
    python_requires=">=3.6",
)
