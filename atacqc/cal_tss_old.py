import argparse, sys
import re,os
import pyranges as pr
import pandas as pd
parser = argparse.ArgumentParser(description='ATAC-seq Data Quality Control.\n Designed by Tao Zhu.')
parser.add_argument('sample', action = 'store', nargs = '?', type = str, default = sys.stdin, help = 'Sample infomation file, it contains sample name and FASTQ path.')
parser.add_argument('tss', action = 'store', nargs = '?', type = str, default = sys.stdin, help = 'TSS position file.')
parser.add_argument('outfile_path', action = 'store', nargs = '?',type = str, default = sys.stdout,help = 'Output path.')
parser.add_argument('tss_score', action = 'store', nargs = '?',type = str, default = sys.stdout,help = 'Output path.')
args = parser.parse_args()
#tss_file=pd.read_table('/public/workspace/zhutao/ref/mm10_tss_pos.csv',usecols=[0,1,2,4],skiprows=1,names='Chromosome Start End Strand'.split())
tss_file=pd.read_table(args.tss,usecols=[0,1,2,4],skiprows=1,names='Chromosome Start End Strand'.split())
tss_gr=pr.PyRanges(tss_file)
cut=pd.read_table(args.sample, usecols=[0,1,2], names='Chromosome Start End'.split())
cut_gr=pr.PyRanges(cut)
dis=cut_gr.nearest(tss_gr).as_df()
dis_fliter_for=dis[dis['Strand']=='+']
dis_fliter_for['new_distance']=dis_fliter_for['Start']-dis_fliter_for['Start_b']
res_for=pd.DataFrame(dis_fliter_for['new_distance'].value_counts()).reset_index().rename(columns={'index':'pos'})
#res_for=pd.DataFrame(dis_fliter_for['new_distance'].value_counts()).reset_index().rename(columns={'new_distance':'pos','count':'new_distance'})
dis_fliter_rev=dis[dis['Strand']=='-']
dis_fliter_rev['new_distance']=dis_fliter_rev['Start_b']-dis_fliter_rev['Start']
# The value_counts() function has been changed in pandas v2.0.3.
res_rev=pd.DataFrame(dis_fliter_rev['new_distance'].value_counts()).reset_index().rename(columns={'index':'pos'})
#res_rev=pd.DataFrame(dis_fliter_rev['new_distance'].value_counts()).reset_index().rename(columns={'new_distance':'pos','count':'new_distance'})
res=res_for.merge(res_rev,on='pos',how='outer')
res=res.fillna(0)
res['cut_sum']=res['new_distance_x']+res['new_distance_y']
mean_val=res[(res['pos']>2500) & (res['pos']<=3000) | (res['pos']>=-3000) &(res['pos']<-2500) ]['cut_sum'].mean()
res['cut_sum']=res['cut_sum'] / mean_val
res_final=res[(res['pos']>=-3000) & (res['pos']<3000)]
res_final.to_csv(args.outfile_path,index=False , sep='\t')
score = open(args.tss_score,'w')
name = re.sub("_q30_cut_sites.bed","",os.path.basename(args.sample))
score.write(name+ '\t' + str(res_final['cut_sum'].max()) + '\n')
score.close()

