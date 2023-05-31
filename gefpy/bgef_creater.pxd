# distutils: language=c++
# cython: language_level=3
import numpy as np
cimport numpy as np

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

from .gef cimport *

cdef extern from "bgefCreater.h" nogil:
    cdef cppclass bgefCreater:
        bgefCreater(int thcnt)
        void createBgef(const string &strin, int bin, const string &strmask, const string &strout)
        void getStereoData(const string &strin, int bin, const string &strmask, vector[string] &vec_gene, vector[unsigned long long] &uniq_cells,
                                 vector[unsigned int] &cell_ind, vector[unsigned int] &gene_ind, vector[unsigned int] &count)