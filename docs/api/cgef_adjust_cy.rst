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