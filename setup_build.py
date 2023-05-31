#!/usr/bin/env python3
# coding: utf-8
"""
    Implements a custom Distutils build_ext replacement, which handles the
    full extension module build process, from cython code to C compilation and
    linking.
"""

try:
    from setuptools import Extension
except ImportError:
    from distutils.extension import Extension
from distutils.command.build_ext import build_ext

import sys
import os
import os.path as op
from pathlib import Path

from setup_configure import BuildConfig


def localpath(*args):
    return op.abspath(op.join(op.dirname(__file__), *args))


MODULES = ['']

COMPILER_SETTINGS = {
    # 'libraries': ['gef'],
    'include_dirs': [],
    'library_dirs': [],
    'define_macros': [('NPY_NO_DEPRECATED_API', 0)],
    'language': 'c++',
    # 'extra_compile_args': ["-std=c++11", "-Wno-sign-compare"]
}

if sys.platform.startswith('win'):
    COMPILER_SETTINGS['libraries'] = ['geftools', 'zlib', 'opencv_world454', 'hdf5', 'tiff']
    # COMPILER_SETTINGS['include_dirs'].append(localpath('windows'))
    COMPILER_SETTINGS['include_dirs'].append(localpath('win\\include'))
    COMPILER_SETTINGS['library_dirs'].append(localpath('win\\lib'))
    COMPILER_SETTINGS['include_dirs'].append(localpath('win\\include\\hdf5'))
    COMPILER_SETTINGS['include_dirs'].append(localpath('win\\include\\libtiff'))
    COMPILER_SETTINGS['include_dirs'].append(localpath('win\\include\\zlib'))
else:
    COMPILER_SETTINGS['libraries'] = ['gef']
    COMPILER_SETTINGS['extra_compile_args'] = ["-std=c++11", "-Wno-sign-compare"]


class gefpy_build_ext(build_ext):
    """
        Custom distutils command which encapsulates api_gen pre-building,
        Cython building, and C compilation.

        Also handles making the Extension modules, since we can't rely on
        NumPy being present in the main body of the setup script.
    """

    @staticmethod
    def _make_extensions(config):
        """ Produce a list of Extension instances which can be passed to
        cythonize().

        This is the point at which custom directories, MPI options, etc.
        enter the build process.
        """
        import numpy

        settings = COMPILER_SETTINGS.copy()
        
        if not sys.platform.startswith('win'):
            settings['include_dirs'][:0] = config.geftools_includedirs
            settings['library_dirs'][:0] = config.geftools_libdirs
            settings['include_dirs'].extend(config.opencv_includedirs)
            settings['include_dirs'].extend(config.hdf5_includedirs)

        try:
            numpy_includes = numpy.get_include()
        except AttributeError:
            # if numpy is not installed get the headers from the .egg directory
            import numpy.core
            numpy_includes = os.path.join(os.path.dirname(numpy.core.__file__), 'include')

        settings['include_dirs'] += [numpy_includes]

        print(settings)

        # TODO: should this only be done on UNIX?
        if os.name != 'nt':
            settings['runtime_library_dirs'] = settings['library_dirs']

        extensions = [
            Extension('gefpy.gene_exp_cy', [localpath('gefpy', 'gene_exp_cy.pyx')],
                      **settings),
            Extension('gefpy.bgef_reader_cy', [localpath('gefpy', 'bgef_reader_cy.pyx')],
                      **settings),
            Extension('gefpy.bgef_writer_cy', [localpath('gefpy', 'bgef_writer_cy.pyx')],
                      **settings),
            Extension('gefpy.cgef_reader_cy',
                      [localpath('gefpy', 'cgef_reader_cy.pyx')], **settings),
            Extension('gefpy.cgef_writer_cy',
                      [localpath('gefpy', 'cgef_writer_cy.pyx')], **settings),
            Extension('gefpy.cgef_adjust_cy',
                      [localpath('gefpy', 'cgef_adjust_cy.pyx')], **settings),
            Extension('gefpy.gef_to_gem_cy',
                      [localpath('gefpy', 'gef_to_gem_cy.pyx')], **settings),
            Extension('gefpy.bgef_creater_cy',
                      [localpath('gefpy', 'bgef_creater_cy.pyx')], **settings)
        ]

        return extensions

    def run(self):
        """ Distutils calls this method to run the command """

        from Cython import __version__ as cython_version
        from Cython.Build import cythonize
        import numpy

        # This allows ccache to recognise the files when pip builds in a temp
        # directory. It speeds up repeatedly running tests through tox with
        # ccache configured (CC="ccache gcc"). It should have no effect if
        # ccache is not in use.
        os.environ['CCACHE_BASEDIR'] = op.dirname(op.abspath(__file__))
        os.environ['CCACHE_NOHASHDIR'] = '1'

        # Get configuration from environment variables
        config = BuildConfig.from_env()
        config.summarise()

        # Run Cython
        print("Executing cythonize()")
        self.extensions = cythonize(self._make_extensions(config),
                                    force=self.force,
                                    language_level=3)
        print("build_ext")
        # Perform the build
        build_ext.run(self)


def write_if_changed(target_path, s: str):
    """Overwrite target_path unless the contents already match s

    Avoids changing the mtime when we're just writing the same data.
    """
    p = Path(target_path)
    b = s.encode('utf-8')
    try:
        if p.read_bytes() == b:
            return
    except FileNotFoundError:
        pass

    p.write_bytes(b)
