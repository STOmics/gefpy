gefpy.bgef_creater_cy
===========================

Provides an interface to create bgef format files from .

.. py:module:: bgef_creater_cy

.. py:class:: BgefCreater(thcnt=8)

    .. py:method:: create_bgef(self, strin, bin, strmask, strout)

        Create tisscuecut bgef by bgef/gem and mask.

        :param strin: raw bgef or bgem
        :param bin: mask binsize
        :param strmask: mask path
        :param strout: out path

    .. py:method:: get_stereo_data(self, strin, bin, strmask)

        Get tisscuecut stereo data by bgef/gem and mask. 

        :param strin: raw bgef or bgem
        :param bin: mask binsize
        :param strmask: mask path

        + uniq_cell is list that save all cell, each cell val (exp.x<<32 | exp.y).
        + gene_names is a list of gene names.
        + count is a list that save the midcnt of each expression.
        + cell_index is a list that save the cell idx of each expression.
        + gene_index is a list that records the gene serial number corresponding to each exp.


        :return: (uniq_cell, gene_names, count, cell_index, gene_index)

