# -*- coding: utf-8 -*-
# distutils: language=c++
# cython: language_level=3, boundscheck=False
# cython: c_string_type=unicode, c_string_encoding=utf8
# Created by zhaozijian on 2022/08/22
"""
    Provides access to the bgef_creater interface.
"""

from .bgef_creater cimport *
import numpy as np
cimport numpy as np
from cython cimport view

cdef class BgefCreater:
    cdef bgefCreater *c_instance 

    def __cinit__(self, thcnt=8):
        self.c_instance = new bgefCreater(thcnt)

    def __init__(self, thcnt=8):
        pass

    def __dealloc__(self):
        del self.c_instance

    def create_bgef(self, strin, bin, strmask, strout):
        self.c_instance.createBgef(strin, bin, strmask, strout)

    def get_stereo_data(self, strin, bin, strmask):

        cdef vector[unsigned int] cell_ind
        cdef vector[unsigned int] gene_ind
        cdef vector[unsigned int] count
        cdef vector[string] gene_names
        cdef vector[unsigned long long] uniq_cell

        self.c_instance.getStereoData(strin, bin, strmask, gene_names, uniq_cell, cell_ind, gene_ind, count)
        return np.asarray(uniq_cell), np.asarray(gene_names), np.asarray(count), np.asarray(cell_ind), np.asarray(gene_ind) 