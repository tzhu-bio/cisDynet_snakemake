import argparse, sys
import re,os
parser = argparse.ArgumentParser(description='ATAC-seq Data Quality Control.\n Designed by Tao Zhu.')
parser.add_argument('sample', action = 'store', nargs = '?', type = str, default = sys.stdin, help = 'Sample infomation file, it contains sample name and FASTQ path.')
parser.add_argument('outfile_path', action = 'store', nargs = '?',type = str, default = sys.stdout,help = 'Output path.')
parser.add_argument('--chromosome', type=str, required=True, help='Chromosomes to process (one or two characters, comma-separated).')
args = parser.parse_args()
chromosomes = args.chromosome.split(',')
idx_txt=open(args.sample,'r')
output=open(args.outfile_path,'w')
if len(chromosomes) == 1:
    for line in idx_txt:
        if line.startswith(str(chromosomes[0])):
            organelle_reads=line.strip().split('\t')[1]
    #if line.startswith('chrY'):
    #    pt_reads=line.strip().split('\t')[1]
if len(chromosomes) == 2:
    for line in idx_txt:
        if line.startswith(str(chromosomes[0])):
            organelle_reads1=line.strip().split('\t')[1]
        if line.startswith(str(chromosomes[1])):
            organelle_reads2=line.strip().split('\t')[1]
    organelle_reads = int(organelle_reads1) + int(organelle_reads2)
output.write(str(re.sub('_idxstats.txt','',os.path.basename(args.sample))) + '\t' + str(int(organelle_reads)) + '\n')
output.close()
