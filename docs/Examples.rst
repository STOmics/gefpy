Examples
============
This section contains various short examples.

Generate bgef
----------------------
.. code:: python

    from gefpy.bgef_writer_cy import generate_bgef

    gem_file = "FP200000617TL_B6.gem"
    ingef_file = "FP200000617TL_B6.bin1.gef"
    bgef_file = "FP200000617TL_B6.bgef"
    stromics = "Transcriptomics"
    bin_sizes = [1,10,20,50,100,200,500]
    region = [1000, 2000, 1000, 2000]

    #generate bgef by gem
    generate_bgef(input_file=gem_file,
                  bgef_file=bgef_file,
                  stromics=stromics,
                  n_thread=8,
                  bin_sizes=bin_sizes,
                  region=region)

    #generate bgef by bin1 gef
    generate_bgef(input_file=ingef_file,
                  bgef_file=bgef_file,
                  stromics=stromics,
                  n_thread=8,
                  bin_sizes=bin_sizes,
                  region=region)

Read bgef
----------------------
.. code:: python

    from gefpy.bgef_reader_cy import BgefR

    bgef = BgefR("FP200000617TL_B6.gef", 50, 4)

    #get expression num
    expnum = bgef.get_expression_num()

    #get cell num
    cellnum = bgef.get_cell_num()

    #get gene num
    genenum = bgef.get_gene_num()

    #get gene name list
    genelist = bgef.get_gene_names()

    #get cell id list(cellid item is (x<<32|y))
    cellid = bgef.get_cell_names()

    #get gene data
    #gene_index is a list that save the gene idx of each expression 
    #gene_names is a list of gene names
    gene_index, gene_names = bgef.get_gene_data()

    #get the all expression, each item is(x, y, count, exon)
    explist = bgef.get_expression()

    #get sparse matrix indexes of expression data
    #uniq_cell is list that save all cell, each cell val (x<<32|y)
    #cell_index is a list that save the cell idx of each expression
    #count is a list that save the midcnt of each expression
    uniq_cell, cell_index, count = bgef.get_exp_data()

    #get the explist by the specified gene name in the region
    explist = bgef.get_genedata_in_region(minx, maxx, miny, maxy, "xxx")

    #get bgef minx miny
    minx, miny = bgef.get_offset()

    #gef bgef attr
    minx, miny, maxx, maxy, maxexp, resolution = bgef.get_exp_attr()

Generate cgef
----------------------
.. code:: python

    from gefpy.cgef_writer_cy import generate_cgef

    mask_file = "FP200000617TL_B6_mask.tif"
    bgef_file = "FP200000617TL_B6.raw.bgef"
    cgef_file = "FP200000617TL_B6.cgef"
    block_sizes = [256, 256]

    # Generate cgef by bgef and mask
    generate_cgef(cgef_file, bgef_file, mask_file, block_sizes)


Read cgef
----------------------
.. code:: python

    from gefpy.cgef_reader_cy import CgefR

    cgef = CgefR("FP200000617TL_B6.cgef")

    #get the number of expression 
    expnum = cgef.get_expression_num()

    #get cell num
    cellnum = cgef.get_cell_num()

    #get gene num
    genenum = cgef.get_gene_num()

    #get a list of gene names
    genelist = cgef.get_gene_names()

    #get a list of cell ids, each cell id is (cell.x <<32 | cell.y)
    cellidlist = cgef.get_cell_names()

    #get all cell 
    celllist = cgef.get_cells()

    #get all gene
    genelist = cgef.get_genes()

    #get the count of each cell in each gene
    cell_id, count = cgef.get_cellid_and_count()

    #get the count of each gene in each cell
    gene_id, count = cgef.get_geneid_and_count()

    #get the borders
    border = cgef.get_cellborders()


Correct cell
----------------------
.. code:: python

    from gefpy.cgef_adjust_cy import CgefAdjust

    adjust = CgefAdjust()
    #1. get gene and cell info by bgef and cgef
    bgef = "FP200000617TL_B6.raw.bgef"
    cgef = "FP200000617TL_B6.cgef"

    #genelist is a list of gene names
    #cell is a list of cell data, every item include (geneid, x, y, midcnt, cellid)
    genelist, cell = adjust.get_cell_data(bgef, cgef)

    #2. do cell correct in stereopy

    #3. write result to cgef
    path = "FP200000617TL_B6.adjust.cgef"
    celltype = np.dtype({'names':['cellid','offset','count'], 'formats':[np.uint32,np.uint32,np.uint32]})
    dnbtype = np.dtype({'names':['x','y','count','gene_id'], 'formats':[np.int32,np.int32,np.uint16,np.uint16]})
    celldata = np.array([(0, 0, 10),(10,10,5),(13,15,20),...], dtype = celltype)
    dnbdata = np.array([(400,400,6,456),(5000,5000,7,258),...], dtype = dnbtype)
    adjust.write_cgef_adjustdata(path, celldata, dnbdata)


Generate gem by gef
----------------------
.. code:: python

    from gefpy.gef_to_gem_cy import gefToGem

    strout = "FP200000617TL_B6.gem"
    strsn = "FP200000617TL_B6"
    obj = gefToGem(strout, strsn)

    # generate bgem
    strbgef = "FP200000617TL_B6.bgef"
    binsize = 10
    obj.bgef2gem(strbgef, binsize)

    # generate cgem by bgef and cgef
    strcgef = "FP200000617TL_B6.cgef"
    obj.cgef2gem(strcgef, strbgef)

    # generate cgem by bgef and mask
    strmask = "FP200000617TL_B6_mask.tif"
    obj.bgef2cgem(strmask, strbgef)


