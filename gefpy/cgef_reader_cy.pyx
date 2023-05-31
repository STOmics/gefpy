# -*- coding: utf-8 -*-
# distutils: language=c++
# cython: language_level=3, boundscheck=False
# cython: c_string_type=unicode, c_string_encoding=utf8
# Created by huangzhibo on 2022/01/01
"""
    Provides access to the cgef_reader interface.
"""

from .cgef_reader cimport *
import numpy as np
cimport numpy as np
from cython cimport view


cdef class CgefR:
    cdef CgefReader* cgef_instance  # Hold a C++ instance which we're wrapping

    def __cinit__(self, filepath, verbose = False):
        self.cgef_instance = new CgefReader(filepath, verbose)

    def __init__(self, filepath, verbose = False):
        """
        A class for reading cell bin GEF.

        :param filepath: Input cell bin GEF filepath.
        """
        pass

    def __dealloc__(self):
        del self.cgef_instance

    def get_expression_num(self):
        """
        Get the number of expression.
        """
        return self.cgef_instance.getExpressionNum()

    def get_cell_num(self):
        """
        Get the number of cell.
        """
        return self.cgef_instance.getCellNum()

    def get_gene_num(self):
        """
        Get the number of gene.
        """
        return self.cgef_instance.getGeneNum()

    def get_gene_names(self):
        """
        Get an array of gene names.
        The type of gene name is bytes.
        """
    #     cdef vector[string] gene_names
    #     gene_names.reserve(self.cgef_instance.getGeneNum())
    #     self.cgef_instance.getGeneNameList(gene_names)
        # below is faster method
        cdef view.array gene_names = view.array((self.cgef_instance.getGeneNum(),),
                                                itemsize=32 * sizeof(char), format='32s', allocate_buffer=True)
        self.cgef_instance.getGeneNames(gene_names.data)
        return np.asarray(gene_names)

    def get_cell_names(self):
        """
        Get an array of cell ids.
        """
        cdef unsigned long long int[::1] cell_names = np.empty(self.get_cell_num(), dtype=np.uint64)
        self.cgef_instance.getCellNameList(&cell_names[0])
        return np.asarray(cell_names)

    def get_cells(self):
        """
        Get cells.

        :return: [cells]
        """
        data_format="""I:id:
        I:x:
        I:y:
        I:offset:
        H:geneCount:
        H:expCount:
        H:dnbCount:
        H:area:
        H:cellTypeID:
        H:clusterID:
        """
        # 可从hdf5 datatset中获取datatype
        # 或许可以参考h5py的方法，自动判断
        cdef int cnum = self.cgef_instance.getCellNum()
        cdef view.array cells
        if cnum != 0:
            cells= view.array((cnum,), itemsize=sizeof(CellData), format=data_format, allocate_buffer=False)
            cells.data = <char*>self.cgef_instance.getCell()
            return np.asarray(cells)
        else:
            return np.array([])

    def get_genes(self):
        """
        Get genes.

        :return: [genes]
        """
        data_format="""32s:geneName:
        I:offset:
        I:cellCount:
        I:expCount:
        H:maxMIDcount:
        """
        cdef view.array genes = view.array((self.cgef_instance.getCellNum(),),
                                           itemsize=sizeof(GeneData), format=data_format, allocate_buffer=False)
        genes.data = <char*>self.cgef_instance.getGene()
        return np.asarray(genes)

    def get_sparse_matrix_indices(self, str order = 'gene'):
        """
        Gets indices for building csr_matrix.

        Examples:
        from scipy import sparse
        sparse.csr_matrix((count, indices, indptr), shape=(cell_num, gene_num))

        :param indices:  CSR format index array of the matrix. Cell id array, the column indices, is the same size as count.
        :param indptr:   CSR format index pointer array of the matrix. indptr length = gene_num + 1 .
        :param count:    CSR format data array of the matrix. Expression count.
        :param order:    Order of count, "gene" or "cell".
        :return: (indices, indptr, count, order)
        """
        indptr_len = self.cgef_instance.getGeneNum() + 1 if order == 'gene' else self.cgef_instance.getCellNum() + 1
        cdef unsigned int[::1] indices = np.empty(self.cgef_instance.getExpressionNum(), dtype=np.uint32)
        cdef unsigned int[::1] indptr = np.empty(indptr_len, dtype=np.uint32)
        cdef unsigned int[::1] count = np.empty(self.cgef_instance.getExpressionNum(), dtype=np.uint32)
        self.cgef_instance.getSparseMatrixIndices(&indices[0], &indptr[0], &count[0], order)
        return np.asarray(indices), np.asarray(indptr), np.asarray(count)

    def get_sparse_matrix_indices2(self):
        """
        Gets indices for building csr_matrix.
        restrict_region and restrict_gene is "invalid" for this function.

        Examples:
        from scipy import sparse
        sparse.csr_matrix((count, cell_ind, gene_ind), shape=(cell_num, gene_num))

        :param cell_ind:  CSR format index array of the matrix. same size as count.
        :param gene_ind:  CSR format index array of the matrix. same size as count.
        :param count:     CSR format data array of the matrix. Expression count.
        :return: (cell_ind, gene_ind, count)
        """
        cdef unsigned int[::1] cell_ind = np.empty(self.cgef_instance.getExpressionNum(), dtype=np.uint32)
        cdef unsigned int[::1] gene_ind = np.empty(self.cgef_instance.getExpressionNum(), dtype=np.uint32)
        cdef unsigned int[::1] count = np.empty(self.cgef_instance.getExpressionNum(), dtype=np.uint32)
        self.cgef_instance.getSparseMatrixIndices2(&cell_ind[0], &gene_ind[0], &count[0])
        return np.asarray(cell_ind), np.asarray(gene_ind), np.asarray(count)

    def restrict_region(self, min_x, max_x, min_y, max_y):
        """
        Restrict to the input region.
        Some member variables (e.g. cell_num) of this class will be updated.

        :param min_x: set the minx
        :param max_x: set the maxx
        :param min_y: set the miny
        :param max_y: set the maxy
        """
        self.cgef_instance.restrictRegion(min_x, max_x, min_y, max_y)

    # TODO fix
    # def restrict_gene(self, vector[string] & gene_list):
    #     """
    #     Restrict to a gene list.
    #
    #     :param gene_list: A list of gene_names
    #     """
    #     self.cgef_instance.restrictGene(gene_list, False)

    def update_gene_info(self):
        """
        Delete genes with zero expression after restrict_region.
        """
        self.cgef_instance.updateGeneInfo()

    def free_restriction(self):
        """
        Free restrict_region and restrict_gene.
        Must call this function before reuse restrict_*.
        """
        self.cgef_instance.freeRestriction()

    def get_cellid_and_count(self):
        """
        Gets cellId and count from the geneExp dataset.

        :return:  (cell_id, count)
        """
        cdef unsigned int[::1] cell_id = np.empty(self.cgef_instance.getExpressionNum(), dtype=np.uint32)
        cdef unsigned short[::1] count = np.empty(self.cgef_instance.getExpressionNum(), dtype=np.uint16)
        self.cgef_instance.getCellIdAndCount(&cell_id[0], &count[0])
        return np.asarray(cell_id), np.asarray(count)

    def get_geneid_and_count(self):
        """
        Gets geneId and count from the cellExp dataset.
        
        :return:  (gene_id, count)
        """
        cdef unsigned int[::1] gene_id = np.empty(self.cgef_instance.getExpressionNum(), dtype=np.uint32)
        cdef unsigned short[::1] count = np.empty(self.cgef_instance.getExpressionNum(), dtype=np.uint16)
        self.cgef_instance.getGeneIdAndCount(&gene_id[0], &count[0])
        return np.asarray(gene_id), np.asarray(count)

    def cgef_close(self):
        """
        Close the cgef.
        """
        self.cgef_instance.closeH5()

    def get_cellborders(self, cellid):
        """
        Gets cell borders.
        
        :return:  [borders]
        """
        
        cdef vector[short] borders
        cdef cnt = self.cgef_instance.getCellBorders(cellid, borders)
        return np.asarray(borders), cnt

    def get_filtered_data(self, region, genelist):

        cdef vector[unsigned int] cell_ind
        cdef vector[unsigned int] gene_ind
        cdef vector[unsigned int] count
        cdef vector[string] gene_names
        cdef vector[unsigned long long] uniq_cell

        self.cgef_instance.getfiltereddata(region, genelist, gene_names, uniq_cell, cell_ind, gene_ind, count)
        return np.asarray(uniq_cell), np.asarray(gene_names), np.asarray(count), np.asarray(cell_ind), np.asarray(gene_ind) 
    
    def get_filtered_data_exon(self, region, genelist):

        cdef vector[unsigned int] cell_ind
        cdef vector[unsigned int] gene_ind
        cdef vector[unsigned int] count
        cdef vector[unsigned int] exon
        cdef vector[string] gene_names
        cdef vector[unsigned long long] uniq_cell

        self.cgef_instance.getfiltereddata_exon(region, genelist, gene_names, uniq_cell, cell_ind, gene_ind, count, exon)
        return np.asarray(uniq_cell), np.asarray(gene_names), np.asarray(count), np.asarray(cell_ind), np.asarray(gene_ind), np.asarray(exon)

    def is_Contain_Exon(self):
        return self.cgef_instance.isContainExon()

