# distutils: language=c++
# cython: language_level=3
import numpy as np
cimport numpy as np

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

from .gef cimport *

cdef extern from "geftogem.h" nogil:

    cdef cppclass geftogem:
        geftogem(const string &strout, const string &strsn, bool boutexon)
        void bgeftogem(const string &strbgef, int binsize)
        void cgeftogem(const string &strcgef, const string &strbgef)
        void bgeftocgem(const string &strmask, const string &strbgef)