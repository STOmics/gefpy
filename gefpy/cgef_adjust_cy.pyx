# -*- coding: utf-8 -*-
# distutils: language=c++
# cython: language_level=3, boundscheck=False
# cython: c_string_type=unicode, c_string_encoding=utf8
# Created by huangzhibo on 2022/01/01
"""
    Provides access to the cgef_adjust interface.
"""

from .cgef_adjust cimport *
import numpy as np
cimport numpy as np
from cython cimport view

cdef class CgefAdjust:
    cdef cellAdjust *c_instance 

    def __cinit__(self):
        self.c_instance = new cellAdjust()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.c_instance

    def get_cell_data(self, bgef, cgef):
        """
        Get raw cell data from cgef and bgef file.

        :param bgef: the bgef file path
        :param cgef: the cgef file path

        :returns: (genelist, vec_cell)
        """

        self.c_instance.readBgef(bgef)
        self.c_instance.readCgef(cgef)
        cdef vector[cellgem_label] vec_cell
        cdef vector[string] genelist
        self.c_instance.getCellLabelgem(genelist, vec_cell)
        return np.asarray(genelist), np.asarray(vec_cell)

    def write_cgef_adjustdata(self, path, celldata, dnbdata, outline_path=''):
        """
        write the adjust cell data to cgef

        :param path: set the Output path
        :param celldata: input the cell data
        :param dandata: input the dandata
        """
        cdef Cell [:] cell = celldata
        cdef DnbExpression [:] dnb = dnbdata
        self.c_instance.writeCellAdjust(path, outline_path, &cell[0], cell.shape[0], &dnb[0], dnb.shape[0])

    def create_Region_Bgef(self, inpath, outpath, pos):
        cdef vector[vector[int]] vec
        for t in pos:
            vec.push_back(t)

        self.c_instance.readBgef(inpath)
        self.c_instance.getRegionGenedata(vec)
        self.c_instance.createRegionGef(outpath)

    def create_Region_Cgef(self, inpath, outpath, pos):
        cdef vector[vector[int]] vec
        for t in pos:
            vec.push_back(t)
        
        self.c_instance.getRegionCelldata(vec)
        self.c_instance.readRawCgef(inpath)
        self.c_instance.writeToCgef(outpath)

    def get_regiondata_frombgef(self, inpath, bin, thcnt, pos, isIndex = False):
        cdef vector[vector[int]] vec
        for t in pos:
            vec.push_back(t)

        cdef vector[vector[int]] vecdata_idx
        cdef vector[sapBgefData] vecdata
        cdef float region_area = 0.0
        if isIndex:
            self.c_instance.getSapRegionIndex(inpath, bin, thcnt, vec, vecdata_idx)
            return np.asarray(vecdata_idx)
        else:
            self.c_instance.getSapRegion(inpath, bin, thcnt, vec, vecdata, region_area)
            return np.asarray(vecdata), region_area
    
    def get_regiondata_fromcgef(self, input_path, pos):
        cdef vector[vector[int]] vec
        for t in pos:
            vec.push_back(t)
        
        self.c_instance.readRawCgef(input_path)
        self.c_instance.getRegionCelldataSap(vec)
        cdef sapCgefData vecdata
        self.c_instance.getSapCellbinRegion(vecdata)
        return np.asarray(vecdata)

    def get_multilabel_regiondata_bgef(self, inpath, pos, bin=1, thcnt=4):
        cdef vector[vector[int]] vec
        for t in pos:
            vec.push_back(t)

        cdef vector[LabelGeneData] region_data
        cdef unsigned int total_mid = 0
        self.c_instance.getMultiLabelInfoFromBgef(inpath, vec, region_data, total_mid, bin, thcnt)
        return np.asarray(region_data), total_mid

    def get_multilabel_regiondata_cgef(self, inpath, pos):
        cdef vector[vector[int]] vec
        for t in pos:
            vec.push_back(t)
        
        # self.c_instance.readRawCgef(inpath)
        # self.c_instance.getRegionCelldataSap(vec)
        cdef vector[LabelCellData] vecdata
        cdef vector[LabelCellData] total_data
        self.c_instance.getMultiLabelInfoFromCgef(inpath, vec, vecdata, total_data)
        return np.asarray(vecdata), np.asarray(total_data)

    def get_position_by_clusterid(self, inpath, clusterid):
        cdef vector[int] vec
        for t in clusterid:
            vec.push_back(t)

        cdef vector[vector[int]] region_data
        self.c_instance.GetPositionIndexByClusterId(inpath, vec, region_data)
        return np.asarray(region_data)

    def generate_filter_bgef_by_midcnt(self, inpath, outpath, binsize, filter_data, only_filter=False):
        cdef vector[MidCntFilter] filter_genes
        cdef MidCntFilter tmp
        for t in filter_data:
            tmp.gene_name = str(t['Gene'])
            tmp.max_mid = int(t['MaxFilterMID'])
            tmp.min_mid = int(t['MinFilterMID'])
            filter_genes.push_back(tmp)
        ret = self.c_instance.GenerateFilterBgefFileByMidCount(inpath, outpath, binsize, filter_genes, only_filter)
        return ret

    def get_filter_bgef_process_rate(self):
        ret = self.c_instance.GenerateFilterBgefDuration()
        return int(ret)

    def generate_bgef_by_lasso(self, inpath, outpath, pos):
        cdef vector[vector[int]] vec
        for t in pos:
            vec.push_back(t)
        ret = self.c_instance.GenerateBgefByLasso(inpath, outpath, vec)
        return ret

    def get_lasso_bgef_process_rate(self):
        ret = self.c_instance.GenerateLassoBgefDuration()
        return int(ret)