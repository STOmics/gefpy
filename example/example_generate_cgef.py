import sys
from gefpy.cgef_writer_cy import generate_cgef

def main():
    mask_file = "../test_data/FP200000617TL_B6/FP200000617TL_B6_mask.tif"
    bgef_file = "../test_data/FP200000617TL_B6/stereomics.h5"
    cgef_file = "../test_data/FP200000617TL_B6/FP200000617TL_B6.gefpy.cgef"
    block_sizes = [256, 256]

    # 从bGEF生成cGEF
    generate_cgef(cgef_file, bgef_file, mask_file, block_sizes)

    return 0


if __name__ == "__main__":
    sys.exit(main())