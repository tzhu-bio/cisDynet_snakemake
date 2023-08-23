import argparse, sys
import re,os
import subprocess
parser = argparse.ArgumentParser(description='ATAC-seq Data Quality Control.\n Designed by Tao Zhu.')
parser.add_argument('sample', action = 'store', nargs = '?', type = str, default = sys.stdin, help = 'Sample infomation file, it contains sample name and FASTQ path.')
parser.add_argument('outfile_path', action = 'store', nargs = '?',type = str, default = sys.stdout,help = 'Output path.')
args = parser.parse_args()
output=open(args.outfile_path,'w')
name = re.sub('_q30.bam','',os.path.basename(args.sample))
proc = subprocess.Popen("samtools view -c %s" % args.sample,shell=True,stdout=subprocess.PIPE)
number = proc.communicate()
number = number[0].strip().decode('ascii')
f = round(float(number)/2)
output.write(name + '\t'+str(f) + '\n')
output.close()
