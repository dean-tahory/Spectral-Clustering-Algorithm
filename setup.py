from setuptools import setup, Extension

"""
This is a setup for k_means_algorthim C extention
"""


setup(name='mykmeanssp',
      version='1.0',
      description='mykmeanssp module for k-means-algorthim implemented with C lanuage',
      ext_modules=[Extension('spkmeans', sources=['spkmeans.c'])])
