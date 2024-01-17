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

    .. py:method:: get_multilabel_regiondata_bgef(self, inpath, pos, bin=1, thcnt=4)

       The gene name and MIDcount of multiple labels are returned after the lasso

        :param inpath: the input bgef file path
        :param pos: lasso region datas(contain multi labels)
        :param bin: binsize
        :param thcnt: thread count
        :returns region_data, total_mid: region_data(gene_name, MIDcount), total midcount in region

    .. py:method:: get_multilabel_regiondata_cgef(self, inpath, pos)

        The gene name and MIDcount of multiple labels are returned after the lasso

        :param inpath: the input cgef file path
        :param pos: lasso region datas(contain multi labels)
        :returns vecdata, total_data: vecdata(cluster_id, mid_cnt, area, cell_id, x, y), total_data(cluster_id, mid_cnt, area, cell_id)

    .. py:method:: get_position_by_clusterid(self, inpath, clusterid)

        Get position value(x, y) by cluster id from h5ad file

        :param inpath: the input h5ad file
        :param clusterid: input cluster id need to get position
        :returns region_data: position value(x, y)

    .. py:method:: generate_filter_bgef_by_midcnt(self, inpath, outpath, binsize, filter_data, only_filter=False)

        generate complete bgef file by gene&protein mid count value

        :param inpath: input bgef file
        :param outpath: output bgef file
        :param binsize: current binsize
        :param only_filter: generate bgef only have filter gene&protein
        :param filter_data: filter gene&protein name and mid count
        :returns ret: generate result

    .. py:method:: get_filter_bgef_process_rate(self)

        Get generate process rate, Must be used in conjunction with the generate_filter_bgef_by_midcnt

        :returns ret: current process rate
        
    .. py:method:: generate_bgef_by_lasso(self, inpath, outpath, pos)

        generate complete bgef file by lasso region datas

        :param inpath: the input bgef file path
        :param outpath: the generate bgef file path
        :param pos: lasso region datas
        :returns ret: generate result

    .. py:method:: get_lasso_bgef_process_rate(self)

        Get generate process rate, Must be used in conjunction with the generate_bgef_by_lasso

        :returns ret: current process rate

    .. py:method:: generate_bgef_by_coordinate(self, inpath, outpath, cord, bin_size)

        generate bgef file by input coordinate in lasso region

        :param inpath: the bgef file path
        :param inpath: the bgef file path
        :param cord: lasso region datas
        :param bin_size: input binsize
        :returns ret: generate result

    .. py:method:: generate_cgef_by_coordinate(self, inpath, outpath, cord)

        generate cell bin gef file by coordinate in lasso region

        :param inpath: the input cgef file path
        :param outpath: the output cgef file path
        :param cord: lasso region datas
        :returns ret: generate result
