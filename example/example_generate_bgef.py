import sys
from gefpy.bgef_writer_cy import generate_bgef

def main():
    gem_file = "../test_data/FP200000617TL_B6/FP200000617TL_B6.gem"
    bgef_file = "../test_data/FP200000617TL_B6/FP200000617TL_B6.gefpy.bins.gef"

    # 从gem生成bGEF
    generate_bgef(gem_file, bgef_file)

    # 从gem生成bin1 bGEF
    bin1_bgef_file = "../test_data/FP200000617TL_B6/FP200000617TL_B6.gefpy.bin1.gef"
    bin_sizes = [1]
    generate_bgef(gem_file, bin1_bgef_file, bin_sizes=bin_sizes)

    # 从多bin bGEF提取bin1 bGEF
    bin1_bgef_file2 = "../test_data/FP200000617TL_B6/FP200000617TL_B6.gefpy.bin1.2.gef"
    generate_bgef(bgef_file, bin1_bgef_file2, bin_sizes=bin_sizes)

    # 提取sub bGEF
    region = [1000, 2000, 1000, 2000]
    bin1_bgef_file2 = "../test_data/FP200000617TL_B6/FP200000617TL_B6.gefpy.bins.sub.gef"
    generate_bgef(bgef_file, bin1_bgef_file2, region=region)

    return 0


if __name__ == "__main__":
    sys.exit(main())