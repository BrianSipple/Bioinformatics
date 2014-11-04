#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name='BioInformatics Algorithms',
    description='Solutions / Algorithms for solving \
     BioInfomatics problems',
    author='Brian Sipple',
    author_email='Bsipple57@gmail.com',
    url='https://github.com/BrianSipple',
    packages=find_packages(exclude=['tests']),
)