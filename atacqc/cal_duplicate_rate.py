import argparse, sys
import re,os
parser = argparse.ArgumentParser(description='ATAC-seq Data Quality Control.\n Designed by Tao Zhu.')
parser.add_argument('sample', action = 'store', nargs = '?', type = str, default = sys.stdin, help = 'Sample infomation file, it contains sample name and FASTQ path.')
parser.add_argument('outfile_path', action = 'store', nargs = '?',type = str, default = sys.stdout,help = 'Output path.')
args = parser.parse_args()
out_df=open(args.outfile_path,'w')
rmdup_txt=open(args.sample,'r')
for line in rmdup_txt:
    if line.startswith('READ'):
        all_reads_num=line.strip().split(' ')[1]
    if line.startswith('DUPLICATE TOTAL'):
        duplicates_num=line.strip().split(' ')[2]
out_df.write(str(re.sub('_rmdup_rate.txt','',os.path.basename(args.sample))) + '\t' + str(round(float(all_reads_num)/2))+ '\t'+ str('%.4f'%(int(duplicates_num)/int(all_reads_num)))+'\n')
out_df.close()
