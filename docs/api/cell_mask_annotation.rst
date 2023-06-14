gefpy.cell_mask_annotation
===========================

Provides a way to generate gem gef tiff file based on contour information .

.. py:module:: cell_mask_annotation

.. py:class:: MaskSegmentation(sampleid, infile, geneFile, outpath, binsize)

    .. py:method:: run_cellMask(self)
 
        generate gem、gef、tiff file based on contour information.
        
        :param sampleid: SampleID for input data
        :param infile: Segmentation mask or geojson
        :param geneFile: Path contains gene expression matrix
        :param outpath: saves the data to the outpath directory
        :param binsize: bin size
