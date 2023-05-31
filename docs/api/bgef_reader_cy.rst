gefpy.bgef_reader_cy
===========================

Provides an interface to read bgef format files.

.. py:module:: bgef_reader_cy

.. py:class:: BgefR(filepath, bin_size, n_thread)

    .. py:method:: get_expression_num(self)

        Get the number of expression.

    .. py:method:: get_cell_num(self)

        Get the number of cell.

    .. py:method:: get_gene_num(self)

        Get the number of gene.

    .. py:method:: get_gene_names(self)

        Get a list of gene names.

    .. py:method:: get_cell_names(self)

        Get a list of cell ids, each item is (exp.x<<32 | exp.y)

    .. py:method:: get_gene_data(self)

        Get gene data.

        + gene_index is a list that records the gene serial number corresponding to each exp.
        + gene_names is a list of gene names.
        
        :return: (gene_index, gene_names)

    .. py:method:: get_expression(self)

        Get the all expression from bgef. 

        + explist is a list, each item is (x, y, count, exon).

        :return: explist

    .. py:method:: get_exp_data(self)

        Get sparse matrix indexes of expression data.

        + uniq_cell is list that save all cell, each cell val (exp.x<<32 | exp.y).
        + cell_index is a list that save the cell idx of each expression.
        + count is a list that save the midcnt of each expression.

        :return: (uniq_cell, cell_index, count)

    .. py:method:: get_genedata_in_region(self, min_x, max_x, min_y, max_y, key)

        Get the explist by the specified gene name in the region.

        :param min_x: region minx
        :param max_x: region maxx
        :param min_y: region miny
        :param max_y: region maxy
        :param key: gene name
        :return: explist

    .. py:method:: get_offset(self)

        Get the offset in bgef.

        :return: (minx, miny)

    .. py:method:: get_exp_attr(self)

        Get the bgef attr.

        :return: (minx, miny, maxx, maxy, maxexp, resolution)

    .. py:method:: get_filtered_data(self, region, genelist)

        Get the filtered data from bgef by region or gene.

        :param region: rect region(minx,maxx,miny,maxy)
        :param genelist: gene name list

        + uniq_cell is list that save all cell, each cell val (exp.x<<32 | exp.y).
        + gene_names is a list of gene names.
        + count is a list that save the midcnt of each expression.
        + cell_index is a list that save the cell idx of each expression.
        + gene_index is a list that records the gene serial number corresponding to each exp.

        :return: (uniq_cell, gene_names, count, cell_index, gene_index)

