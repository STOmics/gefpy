import sys

import numpy as np

from gefpy.gene_exp_cy import GEF
from gefpy.bgef_reader_cy import BgefR

def main():
    bgef_old = GEF("../test_data/FP200000617TL_B6/stereomics.h5", 50)
    gene_indexes, gene_names = bgef_old.get_gene_data()
    cell_names, cell_index, count = bgef_old.get_exp_data()
    print(gene_names)
    print(cell_names)
    print(bgef_old.get_cell_num())
    print(bgef_old.get_gene_num())

    bgef = BgefR("../test_data/FP200000617TL_B6/FP200000617TL_B6.gefpy.bin1.3.gef", 50, 4)

    # generate gem
    bgef.to_gem("../test_data/FP200000617TL_B6/FP200000617TL_B6.gefpy.bin1.3.gem")

    # exp = np.empty(bgef.get_expression_num()*3, dtype=np.int32)
    exp = bgef.get_expression()
    # exp.astype(EXP_DTYPE)
    # print(bgef.get_gene_num())
    # print(bgef.get_expression_num())
    # print(bgef.get_whole_exp_matrix_shape().shape)
    # print(bgef.get_whole_exp_matrix_shape()[0])
    # print(bgef.get_whole_exp_matrix_shape()[1])
    # gene_index, gene_names = bgef.get_sparse_matrix_indices_of_gene()
    # cell_names, cell_index, count = bgef.get_sparse_matrix_indices_of_exp()
    # print(gene_names)
    # print(cell_names)

    print(bgef.get_gene_names())
    print(bgef.get_cell_names())
    indices, indptr, count = bgef.get_sparse_matrix_indices()
    print(indices)
    print(indptr)
    print(count)
    cell_index, gene_index, count = bgef.get_sparse_matrix_indices2()
    print(cell_index)
    print(gene_index)
    return 0


if __name__ == "__main__":
    sys.exit(main())