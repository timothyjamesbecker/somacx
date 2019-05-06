#!/usr/env/bin/python
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import pysam
import numpy
import glob
import os

def two_dot(version):
    v = version.split('.')
    return '.'.join(v[0:min(3, len(v))])

def get_version():
    """Extract version number from source file."""
    from ast import literal_eval
    with open('somacx/utils.pyx') as f:
        for line in f:
            if line.startswith('__version__'):
                return literal_eval(line.partition('=')[2].lstrip())
    raise ValueError("__version__ not found")

cythonize('somacx/utils.pyx')
extensions = [Extension('utils',
                        sources=['somacx/utils.pyx'],
                        libraries=['m'],
                        include_dirs=pysam.get_include() + [numpy.get_include()],
                        define_macros=pysam.get_defines(),
                        extra_compile_args=['-ffast-math'])]

setup(
    name='somacx',
    version=get_version(),
    author='Timothy Becker',
    author_email='timothyjamesbecker@gmail.com',
    url='https://github.com/timothyjamesbecker/somacx',
    license='GPL 3 License',
    description='SV calling data fusion framework',
    classifiers=['Intended Audience :: Developers',
                 'License :: GPL 3 License',
                 'Programming Language :: Python :: 2.7',
                 'Programming Language :: Cython',
                 'Programming Language :: C',
                 'Operating System :: POSIX',
                 'Topic :: Software Development :: Libraries :: Python Modules'],
    cmdclass={'build_ext': build_ext},
    ext_modules=extensions,
    packages=['somacx'],
    package_data={'somacx': ['data/*.json', 'data/*.json.gz','data/*.txt', 'data/*.txt.gz',
                              'data/*.bed', 'data/*.bed.gz','data/*.vcf',]},
    scripts=['bin/generator.py','bin/fastq_sample.py'])
