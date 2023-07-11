# -*- coding: utf-8 -*-
# distutils: language=c++
# cython: language_level=3
# cython: c_string_type=unicode, c_string_encoding=utf8
"""
    Provides access to the cgef_writer interface.
"""
import numpy as np
cimport numpy as np

from .cgef_writer cimport *


def generate_cgef(cgef_file, bgef_file, mask_file, block_size: list):
    """
    Generate cell bin GEF file.

    :param cgef_file: Output CGEF filepath.
    :param bgef_file: Input BGEF filepath.
    :param mask_file: Input make filepath.
    :param block_size: Block size.
    """
    cdef int[::1] bsize = np.array(block_size, dtype=np.int32)
    generateCgef(cgef_file, bgef_file, mask_file, &bsize[0], 0, True)


def cgem_to_cgef(cgem_file, outpath, block_size: list, omics="Transcriptomics"):
    cdef int[::1] bsize = np.array(block_size, dtype=np.int32)
    cgem2cgef(cgem_file, outpath, &bsize[0], 0, omics)

def clusterId_to_cgef(input_file, out_file, cluster_file):
    AddClusterId4Cgef(input_file, out_file, cluster_file)