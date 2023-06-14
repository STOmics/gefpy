gefpy.cgef_adjust_cy
===========================

Provides access to the cgef_adjust interface.

.. py:module:: cgef_adjust_cy

.. py:class:: CgefAdjust()

    .. py:method:: get_cell_data(self, bgef, cgef)

        Get raw cell data from cgef and bgef file.

        :param bgef: the bgef file path
        :param cgef: the cgef file path
        :returns: (genelist, vec_cell)

    .. py:method:: write_cgef_adjustdata(self, path, celldata, dnbdata)

        write the adjust cell data to cgef

        :param path: set the Output path
        :param celldata: input the cell data
        :param dandata: input the dandata
        
    .. py:method:: create_Region_Bgef(self, inpath, outpath, pos)

        generate spatial bin gef file by lasso region datas

        :param inpath: the bgef file path
        :param outpath: set the Output path
        :param pos: lasso region datas

    .. py:method:: create_Region_Cgef(self, inpath, outpath, pos)

        generate cell bin gef file by lasso region datas

        :param inpath: the cgef file path
        :param outpath: set the Output path
        :param pos: lasso region datas

    .. py:method:: get_regiondata_frombgef(self, inpath, bin, thcnt, pos)

        Get gene info from spatial bin gef file by lasso region datas

        :param inpath: the bgef file path
        :param bin: set bin size
        :param thcnt: thread counts
        :param pos: lasso region datas
        :returns vecdata: gene info{genecnt,midcnt,x,y} in region

    .. py:method:: get_regiondata_fromcgef(self, input_path, pos)

        Get cell statistical info from cell bin gef file by lasso region datas

        :param input_path: the cgef file path
        :param pos: lasso region datas
        :returns vecdata: statistical info{cell_count,total_area,average_gene_count,average_exp_count,average_dnb_count,average_area,median_gene_count,median_exp_count,median_dnb_count,median_area} in region
