import os,re
import argparse, sys
parser = argparse.ArgumentParser(description='ATAC-seq Data Quality Control.\n Designed by Tao Zhu.')
parser.add_argument('sample', action = 'store', nargs = '?', type = str, default = sys.stdin, help = 'Sample infomation file, it contains sample name and FASTQ path.')
parser.add_argument('outpath', action = 'store', nargs = '?',type = str, default = sys.stdout,help = 'Output path.')
args = parser.parse_args()
f = open(args.sample,'r')
out = open(args.outpath,'w')
name = re.sub("_flagstat.txt","",os.path.basename(args.sample))
i = 1
for line in f:
    i+=1
    if i==6:
        out.write(name + '\t' + str(round(float(line.strip().split('(')[1].split(':')[0].rstrip(' %'))/100.0,3)) + '\n')
out.close()
