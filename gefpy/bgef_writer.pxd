# -*- coding: utf-8 -*-
# distutils: language=c++
# cython: language_level=3, boundscheck=False
# cython: c_string_type=unicode, c_string_encoding=utf8

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "main_bgef.h" nogil:
    int generateBgef(const string & input_file,
                     const string & bgef_file,
                     const string &stromics,
                     int n_thread,
                     vector[unsigned int] bin_sizes,
                     vector[int] region)

    int generateBgef(const string & input_file,
                     const string & bgef_file,
                     const string &stromics,
                     int n_thread,
                     vector[unsigned int] bin_sizes)

    void MergeProteinAndRnaMatrices(const string &input_gef, const string &output_gef,
                                    const string &omics_type)

    void Gem2Image(const string &gem_path, const string &tif_path)


cdef extern from "bgef_writer.h" nogil:
    pass