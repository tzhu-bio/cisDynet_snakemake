import pandas as pd
from functools import reduce
import os,re
import argparse, sys
parser = argparse.ArgumentParser(description='ATAC-seq Data Quality Control.\n Designed by Tao Zhu.')
#parser.add_argument('sample', action = 'store', nargs = '?', type = str, default = sys.stdin, help = 'Sample infomation file, it contains sample name and FASTQ path.')
parser.add_argument('outpath', action = 'store', nargs = '?',type = str, default = sys.stdout,help = 'Output path.')
args = parser.parse_args()
dup = pd.read_table("final_result/all_dup.csv",names=['Sample','raw_reads_num','dup_rate'])
map_rate = pd.read_table("final_result/all_mapping.csv",names=['Sample','mapping_rate'])
q30 = pd.read_table("final_result/all_q30.csv",names=['Sample','q30'])
organ = pd.read_table("final_result/all_organ.csv",names=['Sample','organ'])
tss = pd.read_table("final_result/all_tss.csv",names = ['Sample','tss'])
frip = pd.read_table("final_result/all_frip.csv",names = ['Sample','frip'])
lst = [dup,map_rate,q30,organ,tss,frip]
df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['Sample'],how='inner'), lst)
df_merged ['organ_mapping'] = df_merged['organ'] /df_merged['raw_reads_num']
df_final = df_merged[['Sample','raw_reads_num','q30','mapping_rate','organ_mapping','dup_rate','tss','frip']]
df_final['mapping_rate'] = df_final['mapping_rate'].astype(float).map(lambda n: '{:.2%}'.format(n))
df_final['organ_mapping'] = df_final['organ_mapping'].astype(float).map(lambda n: '{:.2%}'.format(n))
df_final['dup_rate'] = df_final['dup_rate'].astype(float).map(lambda n: '{:.2%}'.format(n))
df_final['frip'] = df_final['frip'].astype(float).map(lambda n: '{:.2%}'.format(n))
#df_final = df_final.style.format({'organ_mapping':"{:.2%}","dup_rate":"{:.2%}","raw_reads_num":"{:,.0f}","q30":"{:,.0f}"})
df_final.to_csv(args.outpath,header=True,index=False,sep='\t')
