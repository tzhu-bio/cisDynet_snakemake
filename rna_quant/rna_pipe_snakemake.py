import sys
import yaml
from os.path import join
import subprocess
import random
import shutil
import re, os
import glob
from joblib import Parallel, delayed
import pysam
import pandas as pd
import pyranges as pr
import argparse, sys
from functools import reduce
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as tkr
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
sns.set(style="whitegrid")
def change_axs_kb(x, pos):
    return '{:.1f}'.format(float(x)/1000.0)
def message(mes):
  sys.stderr.write("|---"+ mes + "\n")
configfile: "./rna_sample.yaml"
SAMPLES = config['sample']
base_dir = config['workdir']
workdir:base_dir
message("The current dir is" +base_dir)

rule allout:
        input:
            #directory('canFam3STAR'),
            #expand('{sample}_pass1/SJ.out.tab', sample = config["sample"]),
            #directory('SJ'),
            #expand('SJ/{sample}_pass1SJ.filtered.tab', sample = config["sample"]),
            #expand('{sample}_pass2/Aligned.sortedByCoord.out.bam', sample = config["sample"]),
            #expand('{sample}_HTSeq_union_gff3_no_gene_ID.log', sample = config["sample"]),
            #expand('{sample}_HTSeq_raw.csv', sample = config["sample"]),
            #expand('{sample}_HTSeq_final.csv', sample = config["sample"]),
            #expand('{sample}.rsem.log', sample = config["sample"]),
            'all_sample_name.txt',
            expand('{sample}.genes.results',sample = config["sample"]),
            #expand('quant/{sample}.kallisto.log', sample = config["sample"]),
            'RSEM_TPM_temp_Results.txt',
            'RSEM_TPM_Results.txt',
            'RSEM_FPKM_temp_Results.txt',
            'RSEM_FPKM_Results.txt',
            'RSEM_Counts_temp_Results.txt',
            'RSEM_Counts_Results.txt'
#rule index:
#        input:
#            fa = config['fasta'], # provide your reference FASTA file
#            gtf = config['gtf'] # provide your GTF file
#        output:
#            directory('canFam3STAR') # you can rename the index folder
#        threads: 20 # set the maximum number of available cores
#        shell:
#            'mkdir {output} && '
#            'STAR --runThreadN {threads} '
#            '--runMode genomeGenerate '
#            '--genomeDir {output} '
#            '--genomeFastaFiles {input.fa} '
#            '--sjdbGTFfile {input.gtf} '
#            '--sjdbOverhang 100'

#rule pass1:
#        input:
#            R1 = lambda wildcards: config["sample"][wildcards.sample]["fwd"], # may need adjustment if your fastq file name format is different
#            R2 = lambda wildcards: config["sample"][wildcards.sample]["rev"], # note each sample has 4 fastq files ~ 2 lanes per file
#            refdir = config['star_index']
#        params:
#            outdir = '{sample}_pass1',
#            rmbam = '{sample}_pass1/Aligned.out.bam'
#        output:
#            '{sample}_pass1/SJ.out.tab'
#        message: "----STAR Running----"
#        threads: 20 # set the maximum number of available cores
#        shell:
#            'rm -rf {params.outdir} &&' # be careful with this. I don't know why, but Snakemake had problems without this cleaning.
#            'mkdir {params.outdir} && ' # snakemake had problems finding output files with --outFileNamePrefix, so I used this approach instead
#            'cd {params.outdir} && '
#            'STAR --runThreadN {threads} '
#            '--genomeDir {input.refdir} '
#            '--readFilesIn {input.R1} {input.R2} '
#            '--readFilesCommand zcat '
#           '--outSAMtype BAM Unsorted && rm {params.rmbam} && cd ..'

#rule filter:
#        input:
#            '{sample}_pass1/SJ.out.tab',
#            directory('SJ')
#        output:
#            'SJ/{sample}_pass1SJ.filtered.tab'
#        threads: 1
#        run:
#            shell("awk '{{ if ($7 >= 3) print $0 }}' {input[0]} > {output[0]}")

#rule pass2:
#        input:
#            R1 = lambda wildcards: config["sample"][wildcards.sample]["fwd"],
#            R2 = lambda wildcards: config["sample"][wildcards.sample]["rev"],
#            SJfiles = 'SJ/{sample}_pass1SJ.filtered.tab',
#            refdir = config['star_index']
#        params:
#            outdir = '{sample}_pass2',
#            id = '{sample}'
#        output:
#            '{sample}_pass2/Aligned.sortedByCoord.out.bam'
#        message: "----STAR Running Round 2----"
#        threads: 20 # set the maximum number of available cores
#        shell:
#            'rm -rf {params.outdir} && ' # be careful with this. I don't know why, but Snakemake had problems without this cleaning.
#            'mkdir {params.outdir} && '
#            'cd {params.outdir} && '
#            'STAR --runThreadN {threads} '
#            '--genomeDir {input.refdir} '
#            '--readFilesIn {input.R1} {input.R2}'
#            '--readFilesCommand zcat '
#            '--outSAMtype BAM SortedByCoordinate '
#            '--sjdbFileChrStartEnd ../{input.SJfiles} '
#            '--outSAMattrRGline ID:{params.id} '
#            '--quantMode GeneCounts '

#rule htseq:
#        input:
#            bam = '{sample}_pass2/Aligned.sortedByCoord.out.bam',
#            gff = config['gtf']
#        output:
#            'quant/{sample}_HTSeq_raw.csv',
#            'quant/{sample}_HTSeq_final.csv'
#        message: "----HTseq Counts----"
#        threads: 1
#        run:
#            shell("samtools index {input.bam} && htseq-count -m union -s no -t gene -i ID -r pos -f bam {input.bam} {input.gff} &> {output[0]} ")
#            shell("grep -v 'GFF lines processed.' {output[0]}  > {output[1]}")

rule rsem:
          input:
            fq1 = lambda wildcards: config["sample"][wildcards.sample]["fwd"],
            fq2 = lambda wildcards: config["sample"][wildcards.sample]["rev"]
            #rsem_index = config['rsem_index']
          output:
                #'quant/{sample}.rsem.log',
                '{sample}.genes.results'
                #dirs = directory("quant")
          params:
                name = '{sample}',
                rsem_index = config['rsem_index']
          message: "----RSEM Expression Quantification----"
          threads: 10
          run:
#              shell("rsem-calculate-expression --paired-end  -p 10 <(zcat {input.fq1}) <(zcat {input.fq2}) {input.rsem_index} {params.name} > {output[0]}")
               shell("rsem-calculate-expression --paired-end -p {threads} --bowtie2 --estimate-rspd -output-genome-bam {input.fq1} {input.fq2} {params.rsem_index} {params.name}")

#rule kallisto:
#            input:
#              fq1 = lambda wildcards: config["sample"][wildcards.sample]["fwd"],
#              fq2 = lambda wildcards: config["sample"][wildcards.sample]["rev"],
#              kallisto_index = config['kallisto_index']
#            output:
#                  'quant/{sample}.kallisto.log'
#            params:
#                  name = '{sample}'
#            message: "----Kallisto Expression Quantification----"
#            threads: 10
#            run:
#               shell("kallisto quant -i {input.kallisto_index} -o {params.name}.kallisto {input.fq1} {input.fq2} > {output[0]}")

rule get_sample_name:
            input:
                 quant = expand('{sample}.genes.results',sample = config["sample"])
            output:
                 #dirs = directory("quant"),
                 sample_name = 'all_sample_name.txt'
            run:
               #shell("ls quant/*.genes.results|sed 's|quant/||g' |sed 's/.genes.results//g'|paste -s -d "\t"|awk '{{print "Gene""\t"$0}}' > {output.sample_name}")
                shell("sh {config[software][rna_quant]}/sample_name.sh '*.genes.results' {output.sample_name}")
               # shell("sh /public/workspace/zhutao/pipeline/rna_pipe/sample_name.sh ./*.genes.results {output.sample_name}")

rule merge_data:
           input:
                quant = expand('{sample}.genes.results',sample = config["sample"]),
                name = 'all_sample_name.txt'
           output:
                tpm = 'RSEM_TPM_temp_Results.txt',
                TPM = 'RSEM_TPM_Results.txt',
                fpkm = 'RSEM_FPKM_temp_Results.txt',
                FPKM = 'RSEM_FPKM_Results.txt',
                counts = 'RSEM_Counts_temp_Results.txt',
                COUNTS = 'RSEM_Counts_Results.txt'
           message: "----Merge The TPM Data----"
           run:
              shell("sh {config[software][rna_quant]}/rsem.fpkm.sh './*.genes.results' {output.fpkm}") 
              shell("cat all_sample_name.txt {output.fpkm} > {output.FPKM}")
              shell("sh {config[software][rna_quant]}/rsem.tpm.sh './*.genes.results' {output.tpm}")
              shell("cat all_sample_name.txt {output.tpm} > {output.TPM}")
              shell("sh {config[software][rna_quant]}/rsem.counts.sh './*.genes.results' {output.counts}")
              shell("cat all_sample_name.txt {output.counts} > {output.COUNTS}")
