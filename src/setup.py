from setuptools import setup, Extension

"""
This is a setup for spkmeans C extention
"""


setup(name='spkmeans',
      version='1.0',
      description='spkmeans module for spectral clustering algorithm implemented in C',
      ext_modules=[Extension('spkmeans_module', sources=['spkmeansmodule.c'])])
