import argparse, sys
import re,os
import pysam
parser = argparse.ArgumentParser(description='ATAC-seq Data Quality Control.\n Designed by Tao Zhu.')
parser.add_argument('sample', action = 'store', nargs = '?', type = str, default = sys.stdin, help = 'Sample infomation file, it contains sample name and FASTQ path.')
parser.add_argument('outfile_path', action = 'store', nargs = '?',type = str, default = sys.stdout,help = 'Output path.')
args = parser.parse_args()
infile = pysam.AlignmentFile(args.sample, "rb")
outfile = open(args.outfile_path,'w')
for sam in infile:
    if sam.is_reverse :
        pos_end = sam.aend - 5
        outfile.write(sam.reference_name + '\t'+ str(pos_end)+'\t'+str(pos_end + 1)+'\t'+'1'+'\n')
    else:
        pos_start = sam.pos + 4
        outfile.write(sam.reference_name +'\t'+str(pos_start)+'\t'+str(pos_start + 1)+'\t'+'1'+'\n')
outfile.close()
