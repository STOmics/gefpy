#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    Provides convenient drawing methods for other modules.
"""

import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import h5py
import pandas as pd
import seaborn as sns
from gefpy.bgef_reader_cy import BgefR
from gefpy.bgef_writer_cy import generate_bgef
from gefpy.utils import gef_is_cell_bin
import gzip

def get_exp_from_gemfile(input_file, dot_size):
    df = pd.read_csv(input_file, sep='\t', comment='#')
    fileError = False
    x, y, count = [], [], []
    if 'x' in df.columns and 'y' in df.columns:
        if 'MIDCounts' in df.columns:
            df.rename(columns={"MIDCounts": "values"}, inplace = True)
        elif 'MIDCount' in df.columns:
            df.rename(columns={"MIDCount": "values"}, inplace = True)
        else:
            fileError = True
            print("The file cannot be recognized!")
        if dot_size != 1:
            df['x'] = (df['x']/dot_size).astype('int')*dot_size
            df['y'] = (df['y']/dot_size).astype('int')*dot_size
            dnb_dff_tmp = df['values'].groupby([df['x'], df['y']]).sum().reset_index()
            x = dnb_dff_tmp['x']
            y = dnb_dff_tmp['y']
            count = dnb_dff_tmp['values']
        else:
            x = df['x'][1:].astype(np.int32)
            y = df['y'][1:].astype(np.int32)
            count = df['values'][1:].astype(np.int32)
    else:
        fileError = True
        print("The file cannot be recognized!")
    return fileError, x, y, count

def get_binsize_from_gemfile(input_file):
    isHeaderInFile = False
    binSize = 0
    if input_file.endswith('.gem'):
        with open(input_file,'r') as gemFile:
            for i in range(6):
                headLine = gemFile.readline()
                if '#' not in headLine:
                    break
                else:
                    if 'BinSize' in headLine:
                        binSizes = headLine.split('=')
                        binSize = int(binSizes[1])
                        isHeaderInFile = True
                        break
                    else:
                        continue
    else:
        with gzip.open(input_file,'r') as gemFile:
            for i in range(6):
                headLine = gemFile.readline()
                if '#' not in headLine.decode():
                    break
                else:
                    if 'BinSize' in headLine.decode():
                        binSizes = headLine.decode().split('=')
                        binSize = int(binSizes[1])
                        isHeaderInFile = True
                        break
                    else:
                        continue
    return isHeaderInFile, binSize

def save_exp_heat_map_by_binsize(input_file, output_png, bin_size, scale=2, dpi=72):
    if input_file.endswith('.gem') or input_file.endswith('.gz'):
        isHeaderInFile, binSize = get_binsize_from_gemfile(input_file)
        if isHeaderInFile:
            if binSize == bin_size:
                fileError, d_x, d_y, d_cnt = get_exp_from_gemfile(input_file, int(bin_size))
                if fileError:
                    print("The file cannot be recognized!")
                    return False
            else:
                if binSize == 1:
                    fileError, d_x, d_y, d_cnt = get_exp_from_gemfile(input_file, int(bin_size))
                    if fileError:
                        print("The file cannot be recognized!")
                        return False
                else:
                    print("can not do scatter plot, because binsize mismatch!")
                    return False
        else:
            fileError, d_x, d_y, d_cnt = get_exp_from_gemfile(input_file, 1)
            if fileError:
                print("The file cannot be recognized!")
                return False
            print("No standard header, process as bin1 gene expression matrix!")
    elif input_file.endswith('.gef'):
        h5f = h5py.File(input_file, 'r')
        if ('geneExp/bin1' in h5f) or (f'geneExp/bin{bin_size}' in h5f):
            bgef = BgefR(input_file, bin_size, 4)
            exp = bgef.get_reduce_expression()
            d_x = exp[:, 0]
            d_y = exp[:, 1]
            d_cnt = exp[:, 2]
        else:
            print("The file cannot be recognized!")
            h5f.close()
            return False
        h5f.close()
    else:
        print("The file cannot be recognized!")
        return False

    try:
        cmap = mpl.colors.ListedColormap(['#0C3383', '#0A88BA', '#F2D338', '#F28F38', '#D91E1E'])
        x_range=max(d_x) - min(d_x)
        y_range=max(d_y) - min(d_y)
        #x_num = len(data['x'].drop_duplicates())
        x_num = len(list(set(d_x)))

        plt.figure(figsize=(1*scale,y_range/x_range*scale), facecolor='#262B3D', edgecolor='black') ## 设置图像大小 inch
        ##去掉图像旁边的空白区
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        ##将y的坐标翻转
        plt.gca().xaxis.set_ticks_position('top')
        plt.gca().invert_yaxis()
        plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
        plt.margins(0,0)
        ##添加标题
        r = scale*72/(x_range/bin_size)
        dot_size = r**2
        plt.scatter(d_x, d_y, c=d_cnt, s=dot_size, cmap=cmap)
        plt.axis('off')
        plt.savefig(output_png,facecolor='#262B3D', dpi=dpi, pad_inches = 0)
        return True
    except Exception as e:
        print(e)
        return False

# def find_cutoff(score_list, p):
#     """
#     calculate
#     :param score_list:
#     :param p: expression data for gene
#     :return: the cutoff of E10 or C50
#     """
#     #pdb.set_trace()
#     curve = score_list
#     mu = np.mean(curve)
#     sd = statistics.stdev(curve)
#     cutoff = stats.norm.ppf(p) * sd + mu
#     pkk = stats.norm.cdf(stats.norm.ppf(p))
#     return cutoff

def save_exp_heat_map(input_gef, output_png, scale=2, dpi=72):
    """
    Save bgef expression to heat map.

    :param input_gef: set the input bgef path
    :param output_png: set the out png path
    :param scale: set the scale, default 2
    :param dpi: set the dpi, default 72

    """
    h5f = h5py.File(input_gef, 'r')
    if 'cellBin' in h5f:
        d_x = h5f['cellBin']['cell']['x']
        d_y = h5f['cellBin']['cell']['y']
        d_cnt = h5f['cellBin']['cell']['expCount']
        h5f.close()
    else:
        h5f.close()
        bgef = BgefR(input_gef, 200, 4)
        exp = bgef.get_reduce_expression()
        d_x = exp[:, 0]
        d_y = exp[:, 1]
        d_cnt = exp[:, 2]

    try:
        cmap = mpl.colors.ListedColormap(['#0C3383', '#0A88BA', '#F2D338', '#F28F38', '#D91E1E'])
        x_range=max(d_x) - min(d_x)
        y_range=max(d_y) - min(d_y)
        #x_num = len(data['x'].drop_duplicates())
        x_num = len(list(set(d_x)))

        plt.figure(figsize=(1*scale,y_range/x_range*scale), facecolor='#262B3D', edgecolor='black') ## 设置图像大小 inch
        ##去掉图像旁边的空白区
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        ##将y的坐标翻转
        plt.gca().xaxis.set_ticks_position('top')
        plt.gca().invert_yaxis()
        plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
        plt.margins(0,0)
        ##添加标题
        r = scale*72/(x_range/200)
        dot_size = r**2
        plt.scatter(d_x, d_y, c=d_cnt, s=dot_size, cmap=cmap)
        plt.axis('off')
        plt.savefig(output_png,facecolor='#262B3D', dpi=dpi, pad_inches = 0)
        return True
    except Exception as e:
        print(e)
        return False


def cgef_stat(input_cgef, figpath):
    """
    Save cgef stat to png.

    :param input_cgef: set the input bgef path
    :param figpath: set the out png path

    """
    b = 1
    scapath = os.path.join(figpath, "scatter_{0}x{0}_MID_gene_counts.png".format(b if b != 0 else "cell"))
    violinpath = os.path.join(figpath, "violin_{0}x{0}_gene.png".format(b if b != 0 else "cell"))
    violinpath1 = os.path.join(figpath, "violin_{0}x{0}_MID.png".format(b if b != 0 else "cell"))
    statisticPath1 = os.path.join(figpath, "statistic_{0}x{0}_MID.png".format(b if b != 0 else "cell"))
    statisticPath2 = os.path.join(figpath, "statistic_{0}x{0}_gene.png".format(b if b != 0 else "cell"))
    statisticPath3 = os.path.join(figpath, "statistic_{0}x{0}_DNB.png".format(b if b != 0 else "cell"))
    statisticPath4 = os.path.join(figpath, "statistic_{0}x{0}_cell_area.png".format(b if b != 0 else "cell"))
    plt.figure(figsize=(5, 5))

    cgef = h5py.File(input_cgef)
    df = pd.DataFrame(cgef['cellBin']['cell']['expCount', 'geneCount', 'dnbCount', 'area'])
    df = df.rename(columns={'expCount': 'MID Count', 'geneCount': 'Gene Type', 'dnbCount': 'DNB Number', 'area': 'Cell Area'})
    # write to file
    # df.to_csv(os.path.join(figpath, "cellbin.midAndGeneCount.txt"), sep='\t',index=False, header=True)

    # sns.scatterplot(x=df['n_counts'], y=df['n_genes'], edgecolor="gray", color="gray")
    plt.scatter(df['MID Count'], df['Gene Type'], color="gray", edgecolors="gray", s=0.8)
    plt.grid()
    plt.xlabel("MID Count")
    plt.ylabel("Gene Type")
    plt.savefig(scapath, format="png", bbox_inches="tight")

    plt.figure(figsize=(5, 6))
    # plt.subplot(121)
    sns.violinplot(y=df['MID Count'])
    sns.stripplot(y=df['MID Count'], jitter=0.4, color="black", size=0.8)
    plt.ylabel("")
    plt.title("MID Count")
    plt.savefig(violinpath1, format="png", bbox_inches="tight")
    # plt.subplot(122)
    sns.violinplot(y=df['Gene Type'])
    sns.stripplot(y=df['Gene Type'], jitter=0.4, color="black", size=0.8)
    plt.ylabel("")
    plt.title("Gene Type")
    plt.savefig(violinpath, format="png", bbox_inches="tight")

    g = sns.FacetGrid(pd.melt(df[['MID Count']]), col='variable', hue='variable',
                      sharex=False, sharey=False, height=8, palette='Set1')
    g = (g.map(sns.distplot, "value", hist=False, rug=True, color="red"))
    plt.savefig(statisticPath1)

    g = sns.FacetGrid(pd.melt(df[['Gene Type']]), col='variable', hue='variable',
                      sharex=False, sharey=False, height=8, palette='Set1')
    g = (g.map(sns.distplot, "value", hist=False, rug=True, color="blue"))
    plt.savefig(statisticPath2)

    g = sns.FacetGrid(pd.melt(df[['DNB Number']]), col='variable', hue='variable',
                      sharex=False, sharey=False, height=8, palette='Set1')
    g = (g.map(sns.distplot, "value", hist=False, rug=True, color="green"))
    plt.savefig(statisticPath3)

    g = sns.FacetGrid(pd.melt(df[['Cell Area']]), col='variable', hue='variable',
                      sharex=False, sharey=False, height=8, palette='Set1')
    g = (g.map(sns.distplot, "value", hist=False, rug=True, color="violet"))
    plt.savefig(statisticPath4)

if __name__=='__main__':
    #a = [8,9,10,35,78,6,45,23,11,66,33,24,28,54,32, 26]
    #find_cutoff(a, 0.9)
    rc = save_exp_heat_map(
        # "../test_data/FP200000617TL_B6/FP200000617TL_B6.gefpy.bin1.3.gef",
        # "../../test_data/FP200000617TL_B6/FP200000617TL_B6.bgef.h5.gef",
        "../test_data/FP200000617TL_B6/FP200000617TL_B6.gefpy.cgef",
        "../test_data/FP200000617TL_B6/cellbin.png")
    cgef_stat("../test_data/FP200000617TL_B6/FP200000617TL_B6.gefpy.cgef", "../test_data/FP200000617TL_B6/")
    sys.exit(0 if rc else 2)



