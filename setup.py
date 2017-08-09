#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
  name='ssr_perf',
  version='0.2.0',
  description='PERF is an exhaustive repeat finder',
  url='http://192.168.14.249/akshay/perf',
  keywords='ssr microsatellites',
  author='Divya Tej Sowpati',
  author_email='tej@ccmb.res.in',
  license='MIT',
  packages=find_packages(),
  install_requires=['biopython>1', 'tqdm>=4'],
  entry_points={
    'console_scripts': ['ssr-perf=ssr_perf.core:main']
  },
  include_package_data=True
)