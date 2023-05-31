# -*- coding: utf-8 -*-
# distutils: language=c++
# cython: language_level=3, boundscheck=False
# cython: c_string_type=unicode, c_string_encoding=utf8
# Created by huangzhibo on 2022/01/01
"""
    Provides access to the bgef_reader interface.
    For reading common bin GEF.
"""

from .bgef_reader cimport *
import numpy as np
cimport numpy as np

from cython cimport view


#@Depracted
cdef class GEF:
    cdef BgefReader* c_bgef  # Hold a C++ instance which we're wrapping
    cdef unsigned int exp_len
    cdef unsigned int gene_num

    def __cinit__(self, filepath, bin_size):
        self.c_bgef = new BgefReader(filepath, bin_size, 1, True)
        self.exp_len = self.c_bgef.getExpressionNum()
        self.gene_num = self.c_bgef.getGeneNum()

    def __init__(self, filepath, bin_size):
        """
        A class for reading common bin GEF.

        :param filepath: Input bin GEF filepath.
        :param bin_size: Bin size.
        """
        pass

    def __dealloc__(self):
        del self.c_bgef

    def get_expression_num(self):
        """
        Get the number of expression.
        """
        return self.c_bgef.getExpressionNum()

    def get_cell_num(self):
        """
        Get the number of cell.
        """
        return self.c_bgef.getCellNum()

    def get_gene_num(self):
        """
        Get the number of gene.
        """
        return self.c_bgef.getGeneNum()

    def get_exp_data(self):
        """
        Get sparse matrix indexes of expression data.

        :return: (uniq_cell, cell_index, count)
        """
        cdef unsigned int[::1] cell_index = np.empty(self.exp_len, dtype=np.uint32)
        cdef unsigned int[::1] count = np.empty(self.exp_len, dtype=np.uint32)
        cdef vector[unsigned long long] uniq_cell
        uniq_cell.reserve(self.exp_len)
        self.c_bgef.getSparseMatrixIndicesOfExp(uniq_cell, &cell_index[0], &count[0])
        return np.asarray(uniq_cell), np.asarray(cell_index), np.asarray(count)

    def get_gene_data(self):
        """
        Get gene data.
        :return: (gene_index, gene_names)
        """
        cdef unsigned int[::1] gene_index = np.empty(self.c_bgef.getExpressionNum(), dtype=np.uint32)
        # cdef view.array gene_names = view.array((self.c_bgef.getGeneNum(),),
        #                                    itemsize=32*sizeof(char), format='32s', allocate_buffer=True)

        cdef vector[string] gene_names = self.c_bgef.getSparseMatrixIndicesOfGene(&gene_index[0])
        return np.asarray(gene_index), np.asarray(gene_names)
