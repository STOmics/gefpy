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
import array
from numpy cimport ndarray, import_array, PyArray_DATA

from cython cimport view

# EXP_DTYPE = np.dtype([
#     ('x', np.uint32),
#     ('y', np.uint32),
#     ('count', np.uint32),
# ])

import_array()

cdef class BgefR:
    cdef BgefReader* bgef_instance # Hold a C++ instance which we're wrapping
    cdef unsigned int exp_len
    cdef unsigned int gene_num

    def __cinit__(self, filepath, bin_size, n_thread):
        self.bgef_instance = new BgefReader(filepath, bin_size, n_thread, False)
        self.exp_len = self.bgef_instance.getExpressionNum()
        self.gene_num = self.bgef_instance.getGeneNum()

    def __init__(self, filepath, bin_size, n_thread):
        """
        A class for reading common bin GEF.

        :param filepath: Input bin GEF filepath.
        :param bin_size: Bin size.
        :param n_thread: number of thread
        """
        pass

    def __dealloc__(self):
        del self.bgef_instance

    def get_expression_num(self):
        """
        Get the number of expression.
        """
        return self.bgef_instance.getExpressionNum()

    def get_cell_num(self):
        """
        Get the number of cell.
        """
        return self.bgef_instance.getCellNum()

    def get_gene_num(self):
        """
        Get the number of gene.
        """
        return self.bgef_instance.getGeneNum()

    def get_gene_names(self):
        """
        Get an array of gene names.
        """
        # cdef view.array gene_names = view.array((self.bgef_instance.getGeneNum(),),
        #                                         itemsize=32 * sizeof(char), format='32s', allocate_buffer=True)
        cdef vector[string] gene_names
        gene_names.reserve(self.gene_num)
        self.bgef_instance.getGeneNameList(gene_names)
        return np.asarray(gene_names)

    def get_cell_names(self):
        """
        Get an array of cell ids.
        """
        cdef unsigned long long int[::1] cell_names = np.empty(self.get_cell_num(), dtype=np.uint64)
        # cdef view.array gene_names = view.array((self.bgef_instance.getGeneNum(),),
        #                                    itemsize=32*sizeof(char), format='32s', allocate_buffer=True)
        self.bgef_instance.getCellNameList(&cell_names[0])
        return np.asarray(cell_names)

    def get_cell_names2(self, np.ndarray[np.ulonglong_t, ndim=1] cell_names):
        """
        Get an array of cell ids.

        :param cell_names:    cell names.
        """
        # cdef unsigned long long int[::1] cell_names = np.empty(self.get_cell_num(), dtype=np.uint64)
        # cdef view.array gene_names = view.array((self.bgef_instance.getGeneNum(),),
        #                                    itemsize=32*sizeof(char), format='32s', allocate_buffer=True)
        self.bgef_instance.getCellNameList(<unsigned long long int *>cell_names.data)
        # return np.asarray(cell_names)

    def get_gene_data(self):
        """
        Get gene data.
        :return: (gene_index, gene_names)
        """
        cdef unsigned int[::1] gene_index = np.empty(self.exp_len, dtype=np.uint32)
        cdef vector[string] gene_names = self.bgef_instance.getSparseMatrixIndicesOfGene(&gene_index[0])
        return np.asarray(gene_index), np.asarray(gene_names)

    def get_expression(self):
        """
        Get the all expression from bgef.

        :return: exp
        """
        cdef view.array exp = view.array((self.exp_len,4),
                                         itemsize=4*sizeof(char), format='I', allocate_buffer=False)
        exp.data = <char*>self.bgef_instance.getExpression()
        # arr = np.asarray(exp, dtype=EXP_DTYPE)
        # arr.astype(EXP_DTYPE)
        return np.asarray(exp)

    def get_reduce_expression(self):
        """
        Get the reduce expression
        """
        cdef view.array exp = view.array((self.get_cell_num(),4),
                                         itemsize=4*sizeof(char), format='I', allocate_buffer=False)
        exp.data = <char*>self.bgef_instance.getReduceExpression()
        # arr = np.asarray(exp, dtype=EXP_DTYPE)
        # arr.astype(EXP_DTYPE)
        return np.asarray(exp)

    def get_sparse_matrix_indices(self):
        """
        Gets indices for building csr_matrix.

        Examples:
        from scipy import sparse
        sparse.csr_matrix((count, indices, indptr), shape=(cell_num, gene_num))

        :param indices:  CSR format index array of the matrix. Cell id array, the column indices, is the same size as count.
        :param indptr:   CSR format index pointer array of the matrix. indptr length = gene_num + 1 .
        :param count:    CSR format data array of the matrix. Expression count.
        :return: (indices, indptr, count)
        """
        cdef unsigned int[::1] indices = np.empty(self.exp_len, dtype=np.uint32)
        cdef unsigned int[::1] indptr = np.empty(self.gene_num + 1, dtype=np.uint32)
        cdef unsigned int[::1] count = np.empty(self.exp_len, dtype=np.uint32)
        self.bgef_instance.getSparseMatrixIndices(&indices[0], &indptr[0], &count[0])
        return np.asarray(indices), np.asarray(indptr), np.asarray(count)

    def get_sparse_matrix_indices2(self):
        """
        Gets indices for building csr_matrix.

        Examples:
        from scipy import sparse
        sparse.csr_matrix((count, cell_ind, gene_ind), shape=(cell_num, gene_num))

        :param cell_ind:  CSR format index array of the matrix. same size as count.
        :param gene_ind:  CSR format index array of the matrix. same size as count.
        :param count:     CSR format data array of the matrix. Expression count.
        :return: (cell_ind, gene_ind, count)
        """
        cdef unsigned int[::1] cell_ind = np.empty(self.exp_len, dtype=np.uint32)
        cdef unsigned int[::1] gene_ind = np.empty(self.exp_len, dtype=np.uint32)
        cdef unsigned int[::1] count = np.empty(self.exp_len, dtype=np.uint32)
        self.bgef_instance.getSparseMatrixIndices2(&cell_ind[0], &gene_ind[0], &count[0])
        return np.asarray(cell_ind), np.asarray(gene_ind), np.asarray(count)


    def get_exp_data(self):
        """
        Get sparse matrix indexes of expression data.
    
        :return: (uniq_cell, cell_index, count)
        """
        cdef unsigned int[::1] cell_index = np.empty(self.exp_len, dtype=np.uint32)
        cdef unsigned int[::1] count = np.empty(self.exp_len, dtype=np.uint32)
        cdef vector[unsigned long long] uniq_cell
        uniq_cell.reserve(self.exp_len)
        self.bgef_instance.getSparseMatrixIndicesOfExp(uniq_cell, &cell_index[0], &count[0])
        return np.asarray(uniq_cell), np.asarray(cell_index), np.asarray(count)

    # def get_sparse_matrix_indices_of_gene(self):
    #     """
    #     Get gene data.
    #     :return: (gene_index, gene_names)
    #     """
    #     cdef unsigned int[::1] gene_index = np.empty(self.bgef_instance.getExpressionNum(), dtype=np.uint32)
    #     cdef view.array gene_names = view.array((self.bgef_instance.getGeneNum(),),
    #                                        itemsize=32*sizeof(char), format='32s', allocate_buffer=True)
    #
    #     self.bgef_instance.getSparseMatrixIndicesOfGene(&gene_index[0], <char*>gene_names.data)
    #     return np.asarray(gene_index), np.asarray(gene_names)


    def read_whole_exp_matrix(self,
                              unsigned int offset_x,
                              unsigned int offset_y,
                              unsigned int rows,
                              unsigned int cols,
                              string & key,
                              np.ndarray[np.uint8_t, ndim=1] matrix):
        """
        Get wholeExp matrix.

        :param offset_x: The starting position on the x-axis to be read.
        :param offset_y: The starting position on the y-axis to be read.
        :param rows:     Number of rows to read.
        :param cols:     Number of cols to read.
        :param key: MIDcount or genecount.
        :param matrix: When the value is greater than 255, it will be set to 255.
        :return:
        """
        self.bgef_instance.readWholeExpMatrix(offset_x, offset_y, rows, cols, key, <unsigned char*>matrix.data)


    def read_whole_exp_matrix_all(self, str key):
        """
        Get wholeExp matrix.

        :param key: MIDcount or genecount.
        :return: 2D Matrix: When the value is greater than 255, it will be set to 255.
        """
        cdef np.ndarray[np.uint8_t, ndim = 2] matrix = np.empty(self.get_whole_exp_matrix_shape(), dtype=np.uint8)
        self.bgef_instance.readWholeExpMatrix(key, <unsigned char*>matrix.data)
        return matrix

    def get_whole_exp_matrix_shape(self) -> np.ndarray:
        """
        Get the shape of wholeExp matrix.
        :return: [rows, cols]
        """
        cdef view.array shape = view.array((2,), itemsize=sizeof(unsigned int), format="I", allocate_buffer=False)
        shape.data = <char *>self.bgef_instance.getWholeExpMatrixShape()
        return np.asarray(shape)


    def get_filtered_data(self, region, genelist):
        """
        Get the gene exp matrix.
        
        """
        cdef vector[unsigned int] cell_ind
        cdef vector[unsigned int] gene_ind
        cdef vector[unsigned int] count
        cdef vector[string] gene_names
        cdef vector[unsigned long long] uniq_cell

        self.bgef_instance.getfiltereddata(region, genelist, gene_names, uniq_cell, cell_ind, gene_ind, count)
        return np.asarray(uniq_cell), np.asarray(gene_names), np.asarray(count), np.asarray(cell_ind), np.asarray(gene_ind) 

    def get_filtered_data_exon(self, region, genelist):
        """
        Get the gene exp matrix.
        
        """
        cdef vector[unsigned int] cell_ind
        cdef vector[unsigned int] gene_ind
        cdef vector[unsigned int] count
        cdef vector[unsigned int] exon
        cdef vector[string] gene_names
        cdef vector[unsigned long long] uniq_cell

        self.bgef_instance.getfiltereddata_exon(region, genelist, gene_names, uniq_cell, cell_ind, gene_ind, count, exon)
        return np.asarray(uniq_cell), np.asarray(gene_names), np.asarray(count), np.asarray(cell_ind), np.asarray(gene_ind), np.asarray(exon)

    def get_offset(self):
        """
        Get the offset in bgef.

        :return: [minx, miny]
        """
        cdef int offval[2]
        self.bgef_instance.getOffset(offval)
        return offval[0], offval[1]
	
    def get_exp_attr(self):
        """
        Get the bgef attr.
        
        :return: [minx, miny, maxx, maxy, maxexp, resolution]
        """
        cdef int offval[6]
        self.bgef_instance.getExpAttr(offval)
        return offval[0], offval[1], offval[2], offval[3], offval[4], offval[5]

    def is_Contain_Exon(self):
        return self.bgef_instance.isContainExon()

    def get_dnb_data(self, bfilter, btop, level, data_range):
        cdef unsigned int cols = data_range[1] - data_range[0]
        cdef unsigned int rows = data_range[3] - data_range[2]

        cdef vector[unsigned long long] index

        dt = np.dtype([('x', '<f4'), ('y', '<f4'), ('MIDcount', '<u4'), ('genecount', '<u4'), ('color', '<f4')])
        dnbdata = np.zeros(shape=(rows*cols,), dtype=dt)
        cdef void* data = PyArray_DATA(dnbdata)
        cdef unsigned int cnt = self.bgef_instance.getleveldnb(bfilter, btop, level, data_range[0], data_range[2], rows, cols, data, index)
        return np.resize(dnbdata, cnt), np.asarray(index)

    def get_gene_dnb_data(self, bfilter, btop, level, data_range, gene_list):
        cdef unsigned int cols = data_range[1] - data_range[0]
        cdef unsigned int rows = data_range[3] - data_range[2]

        cdef vector[unsigned long long] index
        # cdef vector[levelgenednb] dnbdata

        self.bgef_instance.GetGenesLevelDnb(bfilter, btop, level, data_range[0], data_range[2], data_range[1], data_range[3], index, gene_list)
        
        data_format="""f:x:
        f:y:
        I:midcnt:
        f:color:
        """
        cdef int lnum = self.bgef_instance.getGeneDnbNum()
        cdef view.array dnbdata
        if lnum != 0:
            dnbdata = view.array((lnum,), itemsize=sizeof(levelgenednb), format=data_format, allocate_buffer=False)
            dnbdata.data = <char*>self.bgef_instance.getGeneDnbData()
            return np.asarray(dnbdata), np.asarray(index)
        else:
            return np.array([]), np.array([])