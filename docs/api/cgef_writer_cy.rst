gefpy.cgef_writer_cy
===========================

Provides access to the cgef_writer interface.

.. py:module:: cgef_writer_cy

.. py:function:: generate_cgef(cgef_file, bgef_file, mask_file, block_size: list)

    Generate cell bin GEF file from bgef + mask.

    :param cgef_file: Output CGEF filepath.
    :param bgef_file: Input BGEF filepath.
    :param mask_file: Input make filepath.
    :param block_size: Block size list, usually set to [256,256].

.. py:function:: cgem_to_cgef(cgem_file, outpath, block_size: list)

    Generate cell bin GEF file from cgem.

    :param cgem_file: Input cgem path.
    :param outpath: Output cgef path.
    :param block_size: Block size list,  usually set to [256,256].