gefpy.bgef_writer_cy
===========================

Provides access to the bgef_writer interface.

.. py:module:: bgef_writer_cy

.. py:function:: generate_bgef(input_file, bgef_file, stromics="Transcriptomics", n_thread = 8, bin_sizes = None, region = None)

    Function to generate common bin GEF file(.bgef).

    :param input_file:  The input file path of gem file or bin1 bgef.
    :param bgef_file:   Output BGEF filepath.
    :param stromics:    input the omics.
    :param n_thread:    Number of thread, default 8
    :param bin_sizes:   A list of bin sizes, default: 1,10,20,50,100,200,500
    :param region:      A list of region (minX, maxX, minY, maxY)

.. py:function:: gem2tif(gempath, tif_path)

    Function to generate tif file by GEM file(.gem & .gem.gz).

    :param gempath:  The input file path of gem file.
    :param tif_path:   Output tiff filepath.