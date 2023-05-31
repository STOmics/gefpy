# -*- coding: utf-8 -*-
# distutils: language=c++
# cython: language_level=3, boundscheck=False
# cython: c_string_type=unicode, c_string_encoding=utf8

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "main_cgef.h" nogil:
    int generateCgef(const string& cgef_file,
                     const string& bgef_file,
                     const string& mask_file,
                     const int* block_size,
                     int celltype,
                     bool verbose)

    int cgem2cgef(const string &strcgem, const string &strcgef, const int* block_size, int celltype)

    void AddClusterId4Cgef(const string &input_file,
                                    const string &output_file,
                                    const string &cluster_file)

cdef extern from "cgef_writer.h" nogil:
    pass