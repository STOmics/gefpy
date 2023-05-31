# coding: utf-8
"""
    Implements a new custom Distutils command for handling library
    configuration.

    The "configure" command here doesn't directly affect things like
    config.pxi; rather, it exists to provide a set of attributes that are
    used by the build_ext replacement in setup_build.py.

    Options from the command line and environment variables are stored
    between invocations in a pickle file.  This allows configuring the library
    once and e.g. calling "build" and "test" without recompiling everything
    or explicitly providing the same options every time.

    This module also contains the auto-detection logic for figuring out
    the currently installed GEFTOOLS version.
"""

import os
import os.path as op
import sys


class BuildConfig:
    def __init__(self, opencv_includedirs, hdf5_includedirs, geftools_includedirs, geftools_libdirs):
        self.opencv_includedirs = opencv_includedirs
        self.hdf5_includedirs = hdf5_includedirs
        self.geftools_includedirs = geftools_includedirs
        self.geftools_libdirs = geftools_libdirs

    @classmethod
    def from_env(cls):
        cv_inc, hdf5_inc, gef_inc, gef_lib = cls._find_geftools_compiler_settings()
        return cls(cv_inc, hdf5_inc, gef_inc, gef_lib)

    @staticmethod
    def _find_geftools_compiler_settings():
        """Get compiler settings from environment or pkgconfig.

        Returns (include_dirs, lib_dirs, define_macros)
        """
        geftools = os.environ.get('GEFTOOLS_DIR')
        opencv = os.environ.get('OpenCV_DIR')
        hdf5 = os.environ.get('HDF5_ROOT')

        cv_inc = []
        hdf5_inc = []
        gef_inc = []
        gef_lib = []

        if not opencv:
            #raise ValueError("Specify OpenCV_DIR")
            print("[WARN] OpenCV_DIR value is empty.")
        else:
            cv_inc = [op.join(opencv, 'include/opencv4')]

        if hdf5:
            hdf5_inc = [op.join(hdf5, 'include')]

        if geftools:
            gef_inc = [op.join(geftools, 'include')]
            gef_lib = [op.join(geftools, 'lib')]
            if sys.platform.startswith('win'):
                gef_lib.append(op.join(geftools, 'bin'))

        return cv_inc, hdf5_inc, gef_inc, gef_lib

    def summarise(self):
        def fmt_dirs(l):
            return '\n'.join((['['] + [f'  {d!r}' for d in l] + [']'])) if l else '[]'

        print('*' * 80)
        print(' ' * 23 + "Summary of the gefpy configuration")
        print('')
        print("OpenCV_DIR include dirs:", fmt_dirs(self.opencv_includedirs))
        print("HDF5_ROOT include dirs:", fmt_dirs(self.hdf5_includedirs))
        print("GEFTOOLS include dirs:", fmt_dirs(self.geftools_includedirs))
        print("GEFTOOLS library dirs:", fmt_dirs(self.geftools_libdirs))
        print('')
        print('*' * 80)
