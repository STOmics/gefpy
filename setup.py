#!/usr/bin/env python3
# coding: utf-8
from setuptools import Extension, setup, find_packages
import sys
import os
from pathlib import Path

# Newer packaging standards may recommend removing the current dir from the
# path, add it back if needed.
if '' not in sys.path:
    sys.path.insert(0, '')
import setup_build

if sys.version_info < (3, 7):
    sys.exit('gefpy requires Python >= 3.7')

if os.name == 'nt':
    package_data = {'gefpy': ["*.dll"]}
    exclude_package_data = {'gefpy': ["*.pyx", "*.pxd", "*.cpp"]}
    is_include_package_data = False
else:
    package_data = {'gefpy': [], "gefpy.tests.data_files": ["*.h5"]}
    is_include_package_data = True
    exclude_package_data = {}
# if os.name == 'nt':
#     package_data['gefpy'].append('*.dll')

setup(
    name='gefpy',
    version='0.7.1',
    description='A thin, pythonic wrapper around geftool.',
    long_description=Path('README.md').read_text('utf-8'),
    long_description_content_type="text/markdown",
    url='https://github.com/STOmics/gefpy',
    author='BGIResearch',
    author_email='huangzhibo@genomics.cn',
    python_requires='>=3.7',
    setup_requires=['pkgconfig', 'Cython', 'numpy>=1.20.0,<1.22.0', 'setuptools_scm'],
    install_requires=[
        "h5py <= 3.7.0",
        "numpy >= 1.20.0",
        "matplotlib",
        "seaborn",
        "pandas",
        "geojson",
        "tifffile",
        "opencv-python"
    ],
    extras_require=dict(
        docs=['sphinx>=3.2'],
        test=['pytest>=4.4', 'pytest-nunit'],
    ),
    packages=find_packages(),
    package_data=package_data,
    exclude_package_data=exclude_package_data,
    include_package_data=is_include_package_data,
    ext_modules=[Extension('gefpy.x', ['x.cpp'])],
    classifiers=[
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
    cmdclass={'build_ext': setup_build.gefpy_build_ext}
)
