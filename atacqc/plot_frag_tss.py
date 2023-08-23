import os,re
import glob
import argparse, sys
import pandas as pd
import pyranges as pr
import seaborn as sns
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as tkr
def change_axs_kb(x, pos):
    return '{:.1f}'.format(float(x)/1000.0)
parser = argparse.ArgumentParser(description='ATAC-seq Data Quality Control.\n Designed by Tao Zhu.')
parser.add_argument('sample', action = 'store', nargs = '?', type = str, default = sys.stdin, help = 'Sample infomation file, it contains sample name and FASTQ path.')
parser.add_argument('tss', action = 'store', nargs = '?', type = str, default = sys.stdin, help = 'TSS file.')
parser.add_argument('outfile_path', action = 'store', nargs = '?',type = str, default = sys.stdout,help = 'Output path.')
args = parser.parse_args()
def plot_frag_tss(sam, tss, png):
    
            #path, dirs, files = next(os.walk("fragments_size/"))
            #file_count = len(files)
            #fig, ax = plt.subplots(figsize=[5,3])
            #plt.subplots_adjust(wspace =0.5, hspace =0.5)
            #fig, axs = plt.subplots(file_count,2,squeeze=False, figsize=[10,5*file_count])
            #plt.figure(figsize=[10,5*file_count])
            #plt.rcParams['font.sans-serif'] = "Arial"
            #plt.rcParams['font.family'] = "sans-serif"
            #plt.rcParams['axes.labelweight'] = 'bold'
            #plt.rcParams['axes.titleweight'] = 'bold'
            #sns.set(style="whitegrid")
            #plt.rcParams["axes.edgecolor"] = "0.15"
            #plt.rcParams["axes.linewidth"]  = 0.5
            
    fig, ax = plt.subplots(figsize=[5,2])
    plt.subplots_adjust(wspace =0.5, hspace =0.5)
            #fig, axs = plt.subplots(file_count,2,squeeze=False, figsize=[10,5*file_count])
            #plt.figure(figsize=[10,5*file_count])
    plt.rcParams['font.sans-serif'] = "Arial"
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['axes.titleweight'] = 'bold'
    sns.set(style="whitegrid")
    plt.rcParams["axes.edgecolor"] = "0.15"
    plt.rcParams["axes.linewidth"]  = 0.5
                #frag_size=pd.read_table('%s/plot_df/fragments_size/%s*_size.txt'%(bam_path,name),names='chr pos len'.split())
    frag_size=pd.read_table(sam,names='chr pos len'.split())
    name=re.sub('_q30_frag_size.txt','',os.path.basename(sam))
    tss=pd.read_table(tss)
    plt.subplot(1,2,1)
    sns.kdeplot(data=frag_size, x="len",shade=True,common_norm=False,alpha=.8, linewidth=0,color="#585EAA")
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.title('%s'%name,fontsize=10)
                #axs[i,0].set_xticklabels('Fragment size (bp)')
                #axs[i,0].set_yticklabels('Density')
    plt.xlabel('Fragment size (bp)',fontsize=10)
    plt.ylabel('Density',fontsize=10)
    plt.subplot(1,2,2)
    ax=sns.lineplot(x=tss['pos'],y=tss['cut_sum'],linewidth=0.5,color="#438E95")
    plt.xticks(fontsize=8)
                #axs[i,1].set_xticklabels('Distance from TSS (kb)')
                #axs[i,1].set_yticklabels('TSS enrichment score')
    x_format = tkr.FuncFormatter(change_axs_kb)
    ax.xaxis.set_major_formatter(x_format)
    plt.yticks(fontsize=8)
    plt.xlabel('Distance from TSS (kb)',fontsize=10)
    plt.ylabel('TSS enrichment score',fontsize=10)
    plt.title('%s'%name,fontsize=10)
    plt.savefig(png,dpi=300,bbox_inches='tight')
    return()
plot_frag_tss(args.sample, args.tss, args.outfile_path)
