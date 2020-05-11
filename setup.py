#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
  name='perf_ssr',
  version='0.4.5',
  description='PERF is an exhaustive repeat finder',
  url='https://github.com/rkmlab/perf',
  keywords='ssr microsatellites',
  author='Divya Tej Sowpati',
  author_email='tej@ccmb.res.in',
  license='MIT',
  packages=find_packages(),
  install_requires=['biopython==1.69', 'tqdm>=4'], # biopython version 1.69 installs numpy
  entry_points={
    'console_scripts': ['PERF=PERF.core:main']
  },
  include_package_data=True, # change path according to package name in MANIFEST.in
)