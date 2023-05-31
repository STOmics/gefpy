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

    def get_regiondata_frombgef(self, inpath, bin, thcnt, pos):
        cdef vector[vector[int]] vec
        for t in pos:
            vec.push_back(t)

        cdef vector[sapBgefData] vecdata
        self.c_instance.getSapRegion(inpath, bin, thcnt, vec, vecdata)
        return np.asarray(vecdata)
    
    def get_regiondata_fromcgef(self, input_path, pos):
        cdef vector[vector[int]] vec
        for t in pos:
            vec.push_back(t)

        self.c_instance.getRegionCelldataSap(vec)
        self.c_instance.readRawCgef(input_path)
        cdef sapCgefData vecdata
        self.c_instance.getSapCellbinRegion(vecdata)
        return np.asarray(vecdata)