# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

name='rivas_decomposition_py'
version='0.0.19'
short_description="Python packages for Rivas Lab.'s Decomposition project."

########

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name=name,
    version=version,
    description=short_description,
    long_description=readme,
    author='Yosuke Tanigawa',
    author_email='info@yosuketanigawa.com',
    install_requires=requirements,
    url='https://yosuketanigawa.com/software/{}'.format(name),
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)

