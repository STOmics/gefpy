# -*- coding: utf-8 -*-
# distutils: language=c++
# cython: language_level=3, boundscheck=False
# cython: c_string_type=unicode, c_string_encoding=utf8


cdef extern from "gef.h":
    ctypedef union Coordinate:
        unsigned int pos[2] # dnb coordinates x, y
        unsigned long long int pos_id

    ctypedef struct Expression:
        int x
        int y
        unsigned int count
        unsigned int exon

    ctypedef struct CellExpData:
        unsigned int gene_id
        unsigned short count

    ctypedef struct GeneExpData:
        unsigned short cell_id
        unsigned short count

    ctypedef struct Gene:
        char gene[32]
        unsigned int offset
        unsigned int count

    ctypedef struct CellData:
        unsigned int id
        int x              # Coordinate X of center point in this cell
        int y              # Coordinate Y of center point in this cell
        unsigned int offset         # Offset of current cell in cellExp, 0-based
        unsigned short gene_count   # The number of gene in this cell
        unsigned short exp_count    #  The total expression count of all genes in this cell
        unsigned short dnb_count    #  Dnb number in this cell
        unsigned short area         #  The polygon area of this cell
        unsigned short cell_type_id #  Cell type ID to index the CellTypeList
        unsigned short cluster_id #  Cell type ID to index the CellTypeList

    ctypedef struct GeneData:
        char gene_name[32]
        unsigned int offset
        unsigned int cell_count
        unsigned int exp_count
        unsigned short max_mid_count

    ctypedef struct Cell:
        unsigned int cellid
        unsigned int offset
        unsigned int count

    ctypedef struct DnbExpression:
        int x
        int y
        unsigned short count
        unsigned int gene_id
