#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Provides convenient invocation methods for other modules.
"""

import h5py
import numpy as np

def gef_is_square_bin(gef_file):
    """
    Determine if the file is a bgef file

    :param gef_file: input the gef path

    :return: the return True for bgef, False otherwise.
    :rtype: bool
    """

    h5f = h5py.File(gef_file)
    if 'geneExp/bin1/expression' in h5f:
        h5f.close()
        return True
    else:
        h5f.close()
        return False

def gef_is_cell_bin(gef_file):
    """
    Determine if the file is a cgef file

    :param gef_file: input the gef path

    :return: the return True for cgef, False otherwise.
    :rtype: bool
    """

    h5f = h5py.File(gef_file)
    if 'cellBin' in h5f:
        h5f.close()
        return True
    else:
        h5f.close()
        return False

def gef_contain_exon(gef_file):
    h5f = h5py.File(gef_file)
    if 'cellBin' in h5f: #cgef
        if '/cellBin/cellExon' in h5f:
            return True
    else:
        if '/geneExp/bin1/exon' in h5f:
            return True

    return False

def StereoDataToGef(path, bin, expression, gene, attr):
    """
    Write stereodata to cgef file

    :param path: set the output path
    :param bin: set the bin size
    :param expression: input the expression data
    :param gene: input the gene data
    :param attr: input the attr data

    """
    h5f = h5py.File(path,"w")
    geneExp = h5f.create_group("geneExp")
    binsz = "bin"+str(bin)
    bing = geneExp.create_group(binsz)

    #genetype = np.dtype({'names':['gene','offset','count'], 'formats':['S32', 'i','i']})
    #exptype = np.dtype({'names':['x','y','count'], 'formats':['i', 'i','i']})

    geneExp[binsz]["expression"] = expression #np.array([(10,20,2), (20,40,3)], dtype=exptype)
    geneExp[binsz]["gene"] = gene #np.array([("gene1",0,21), ("gene2",21,3)], dtype=genetype)

    bing["expression"].attrs.create("minX", attr[0], dtype=np.uint32)
    bing["expression"].attrs.create("minY", attr[1], dtype=np.uint32)
    bing["expression"].attrs.create("maxX", attr[2], dtype=np.uint32)
    bing["expression"].attrs.create("maxY", attr[3], dtype=np.uint32)
    bing["expression"].attrs.create("maxExp", attr[4], dtype=np.uint32)
    bing["expression"].attrs.create("resolution", attr[5], dtype=np.uint32)

    h5f.attrs.create("version", 2, dtype=np.uint32)
    ver=[0,6,4]
    h5f.attrs.create("geftool_ver", ver, dtype=np.uint32)
    h5f.close()