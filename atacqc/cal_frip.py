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
parser = argparse.ArgumentParser(description='ATAC-seq Data Quality Control.\n Designed by Tao Zhu.')
parser.add_argument('frip', action = 'store', nargs = '?', type = str, default = sys.stdin, help = 'Sample infomation file, it contains sample name and FASTQ path.')
parser.add_argument('q30', action = 'store', nargs = '?', type = str, default = sys.stdin, help = 'Sample infomation file, it contains sample name and FASTQ path.')
parser.add_argument('outfile', action = 'store', nargs = '?',type = str, default = sys.stdout,help = 'Output file.')
args = parser.parse_args()
frip = open(args.frip)
q30_reads = open(args.q30)
output=open(args.outfile,'w')
num = frip.readline().strip()
name = str(re.sub('_frip.csv','',os.path.basename(args.frip)))
q30 = q30_reads.readline().strip().split("\t")[1]
ratio = round(int(num) / (2 * int(q30)),4)
output.write(name+ '\t' + str(ratio) + '\n')
output.close()
