
from gefpy import cell_mask_annotation
import os


sampleid = ""
infile = "mask or geojson file path"
geneFile = "gene expression matrix"
outpath = "output path"
binsize = 1

seg = cell_mask_annotation.MaskSegmentation(sampleid, infile, geneFile, os.path.join(outpath, 'segmentation'), binsize)
seg.run_cellMask()