# -*- coding: utf-8 -*-
# distutils: language=c++
# cython: language_level=3, boundscheck=False
# cython: c_string_type=unicode, c_string_encoding=utf8
# Created by huangzhibo on 2022/01/01
"""
    Provides access to the gef2gem interface.
"""

from .gef_to_gem cimport *

cdef class gefToGem:
    cdef geftogem *c_instance 

    def __cinit__(self, strout, strsn, boutexon=True):
        """
        Init the instance.

        :param strout: set out path
        :param strsn: input SN
        :param boutexon: set exon whether or not
        """

        self.c_instance = new geftogem(strout, strsn, boutexon)

    def __init__(self, strout, strsn, boutexon=True):
        pass

    def __dealloc__(self):
        del self.c_instance

    def bgef2gem(self, strbgef, binsize):
        """
        Create bgem file by bgef.

        :param strbgef: the bgef file path
        :param binsize: set the binsize
        """

        self.c_instance.bgeftogem(strbgef, binsize)

    def cgef2gem(self, strcgef, strbgef):
        """
        Create cgem file by cgef and bgef.

        :param strcgef: the cgef file path
        :param strbgef: the bgef file path
        """

        self.c_instance.cgeftogem(strcgef, strbgef)

    def bgef2cgem(self, strmask, strbgef):
        """
        Create cgem file by mask and bgef.

        :param strcgef: the mask file path
        :param strbgef: the bgef file path
        """

        self.c_instance.bgeftocgem(strmask, strbgef)
