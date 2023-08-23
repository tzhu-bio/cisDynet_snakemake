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
configfile: "/public/workspace/zhutao/pipeline/atacqc/test_sample.yaml"
SAMPLES = config['sample']
base_dir = "/public/workspace/zhutao/pipeline/"
workdir:base_dir
message("The current dir is" +base_dir)
trima = "/public/home/xbxu/software/Trimmomatic-0.36"
gff = "/public/home/jcli/data/reference/gff3/Rice.gff3"
fasta = "/public/workspace/zhutao/ref/mm10"
#blackbed = "/public/home/yinmliu/test/snakemake/genome_delete_bin_list.bed"
nextera_adapters = "/public/workspace/zhutao/miniconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/NexteraPE-PE.fa"
#nextera_adapters = "/public/home/yinmliu/software/Trimmomatic-0.38/adapters/test.fa"
#q30txt = "/public/home/yinmliu/test/xxb_test/testraw/bwa_out/test.txt"
#organ = "/public/home/yinmliu/test/xxb_test/testraw/Q95_rice_panicle_chr_14_rm_organelles_bin_top10_quantitation.bed"
#black_list = "/public/home/yinmliu/test/xxb_test/testraw/genome_delete_bin_list.bed"
or_fasta = "/public/home/yinmliu/test/rice_organell/Jan_organell.fasta"
#txt = "/public/home/yinmliu/test/xxb_test/testraw/bwa_out_notrim_0929/q30.text"
#for smp in Samples:
 # message("Sample" +smp)
rule all:
   input:
         "multiqc/multiqc_report.html",
         "final_result/all_dup.csv",
         "final_result/all_mapping.csv",
         "final_result/all_q30.csv",
         "final_result/all_organ.csv",
         #"final_result/all_qc.csv",
         "final_result/final_qc.csv",
         "final_result/ATAC_Quality_Control_Result_mqc.txt",
         expand("raw_bam/{smp}_raw.bam",smp = config["sample"]),
         expand("temp_file/{smp}_rmdup.bam",smp = config["sample"]),
         expand("temp_file/{smp}_rmdup_rate.txt",smp = config["sample"]),
         expand("temp_file/{smp}_sort.bam",smp = config["sample"]),
         expand("temp_file/{smp}_dup_rate.csv",smp = config["sample"]),
         expand("temp_file/{smp}_mapping_rate.txt",smp = config["sample"]),
         expand("q30/{smp}_q30.bam",smp = config["sample"]),
         expand("temp_file/{smp}_organell_mapping_rate.csv",smp = config["sample"]),
         expand("temp_file/{smp}_q30_reads_num.csv",smp = config["sample"]),
         expand("temp_file/{smp}_idxstats.txt",smp = config["sample"]),
         expand("fragments_size/{smp}_q30_frag_size.txt",smp = config["sample"]),
         expand("cut_sites/{smp}_q30_cut_sites.bed",smp = config["sample"]),
         expand("tss/{smp}_tss.csv",smp = config["sample"]),
         expand("temp_file/{smp}_flagstat.txt",smp = config["sample"]),
         expand("sort/{smp}_sorted_raw.bam",smp = config["sample"]),
         expand("trimmed/{smp}_R1_001.fastq.paired.gz",smp = config["sample"]),
         expand("trimmed/{smp}_R2_001.fastq.paired.gz",smp = config["sample"]),
         expand("trimmed/{smp}_R1_001.fastq.unpaired.gz",smp = config["sample"]),
         expand("trimmed/{smp}_R2_001.fastq.unpaired.gz",smp = config["sample"]),
         expand("plot/{smp}_Fragment_size_and_TSS_enrichment_mqc.png",smp = config["sample"])
         #expand("plot/{smp}_tss.png",smp = config["sample"]),
         #expand("temp_result/{smp}_dup_rate.txt",smp = config["sample"]),
         #expand("final_result/{smp}_organell_mapping_rate.csv",smp = config["sample"])
         #expand("q30/{smp}_q30_filtered.bam",smp = config["sample"]),
         #expand("temp/{smp}_mapping_rate.txt",smp = config["sample"]),
         #expand("temp/{smp}_idxstats.txt",smp = config["sample"]),
         #expand("fragments_size/{smp}_q30_frag_size.txt",smp = config['sample']),
         #expand("cut_sites/{smp}_q30_cut_sites.bed",smp = config["sample"]),
         #expand("tss/{smp}_tss.csv",smp = config["sample"]),
         #expand("plot/{smp}_tss.png", smp = config["sample"]),
         #expand("plot/{smp}_frag.png", smp = config["sample"])
         #"plot/All_Fragment_size_and_TSS_enrichment_mqc.png",
         #"final_result/ATAC_Quality_Control_Result_mqc.txt"
         #expand("temp/all_q30_reads_num.txt"),
         #expand("final_result/all_mapping_rate.txt"),
         #expand("final_result/organell_mapping_rate.csv")

rule trimming:
     input:
         fwd = lambda wildcards: config["sample"][wildcards.smp]["fwd"],
         rev = lambda wildcards: config["sample"][wildcards.smp]["rev"]
     output:
         fwd = "trimmed/{smp}_R1_001.fastq.paired.gz",
         rev = "trimmed/{smp}_R2_001.fastq.paired.gz",
         fwd1 = "trimmed/{smp}_R1_001.fastq.unpaired.gz",
         rev1 = "trimmed/{smp}_R2_001.fastq.unpaired.gz"
     params: next_adapters = nextera_adapters
     threads:5
     log:     
         "trimmed/{smp}_cut.log"
     shell:
        """
        ~/miniconda3/bin/trimmomatic PE -threads 5 -phred33 -trimlog {log} {input.fwd} {input.rev} {output.fwd} {output.fwd1} {output.rev} {output.rev1} ILLUMINACLIP:{params.next_adapters}:2:30:10:8:TRUE SLIDINGWINDOW:4:15 MINLEN:30
        """



rule bwa:
    input:
        fwd = "trimmed/{smp}_R1_001.fastq.paired.gz",
        rev = "trimmed/{smp}_R2_001.fastq.paired.gz"
    output:
        bam = "raw_bam/{smp}_raw.bam"
    params:
         Fasta = fasta,
         extra = "-R \"@RG\\tID:{smp}\\tSM:{smp}\\tLB:{smp}\\tPL:illumina\""
    threads: 5
    message: """--- BWA mem reads."""
    run:
        shell("~/miniconda3/bin/bwa mem -M -t {threads} -k 32 {params.extra} {params.Fasta} {input.fwd} {input.rev} |samtools view -bhS >{output.bam}")

rule sorted_by_name:
    input:
         bam = "raw_bam/{smp}_raw.bam"
    output:
         sort = "temp_file/{smp}_sort.bam"
    threads: 5
    message: """--- Sort bam by name."""
    run:
        shell("samtools sort -n -@ {threads} {input.bam} -o {output.sort}")


rule remove_duplicates:
    input:
         bam = "temp_file/{smp}_sort.bam"
    output:
         rmdup_rate = "temp_file/{smp}_rmdup_rate.txt",
         rmdup_bam = "temp_file/{smp}_rmdup.bam"
         #out = "final_result/all_dup_rate.csv"
    threads: 5
    message: """--- Remove PCR duplicates."""
    run:
        shell("samtools fixmate -@ {threads} -cm -O bam {input.bam} - |samtools sort -@ {threads} |samtools markdup -r -s - {output.rmdup_bam} 2>{output.rmdup_rate}")
        shell("samtools index {output.rmdup_bam}")



rule cal_duplicate_rate:
    input:
         bam = "temp_file/{smp}_rmdup_rate.txt"
    output:
         out = "temp_file/{smp}_dup_rate.csv"
    threads: 5
    message: """--- calculate PCR duplicates rate."""
    run:
        shell("python ./atacqc/cal_duplicate_rate.py {input.bam} {output.out}")



rule mapping_rate:
    input:
      bam = "raw_bam/{smp}_raw.bam"
    output:
      mapped_rate = "temp_file/{smp}_mapping_rate.txt",
      flagstat = "temp_file/{smp}_flagstat.txt"
      #out_f = "/final_result/all_mapping_rate.csv"
    threads: 5
    message: """--- Calculate mapping rate."""
    run:
         shell("samtools flagstat {input.bam} > {output.flagstat}")
         shell("python ./atacqc/cal_mapping_rate.py {output.flagstat} {output.mapped_rate}")

rule remove_low_quality:
    input:
      rmdup = "temp_file/{smp}_rmdup.bam"
    output:
      bam = "q30/{smp}_q30.bam",
      q30_reads = "temp_file/{smp}_q30_reads_num.csv"
    threads: 5
    message: """--- generate q30 bam and calculate q30 reads number."""
    run:
     shell("samtools view -h -q 30 {input.rmdup}|grep -v 'chr13\|chr14'|samtools view -bF 4 - >{output.bam}")
     shell("python ./atacqc/cal_q30_reads_num.py {output.bam} {output.q30_reads}")


rule organell_mapping_rate:
    input:
      bam = "raw_bam/{smp}_raw.bam"
    output:
      organell = "temp_file/{smp}_organell_mapping_rate.csv",
      idxstat = "temp_file/{smp}_idxstats.txt",
      sort = "sort/{smp}_sorted_raw.bam"
    threads: 5
    message: """--- organell mapping rate."""
    run:
     shell("samtools sort -@ {threads} {input.bam} -o {output.sort}")
     shell("samtools index {output.sort}")
     shell("samtools idxstats -@ {threads} {output.sort} |cut -f1,3 > {output.idxstat}")
     shell("python ./atacqc/cal_mt_pt_reads_num.py {output.idxstat} {output.organell}")

rule fragments_size:
    input:
      bam = "q30/{smp}_q30.bam"
    output:
      txt = "fragments_size/{smp}_q30_frag_size.txt"
    message: """--- fragments size."""
    shell:
      """
      samtools view {input.bam}|awk '$9>0 && $9<800'|cut -f3,4,9 > {output.txt}
      """

rule get_atac_cut_sites:
    input:
      bam = "q30/{smp}_q30.bam"
    output:
      txt = "cut_sites/{smp}_q30_cut_sites.bed"
    message: """--- get tn5 cuts."""
    threads:5
    run:
        shell("python ./atacqc/cal_cut_sites.py {input.bam} {output.txt}")


rule TSS_enrich:
     input:
       bed = "cut_sites/{smp}_q30_cut_sites.bed"
     output:
       tss = "tss/{smp}_tss.csv",
       score = "tss/{smp}_score.csv"
     message: """--- tss enrichment."""
     threads: 5
     run:
        shell("python ./atacqc/cal_tss.py {input.bed} {output.tss} {output.score}")



rule plot_fra_size:
     input:
        frag_size = "fragments_size/{smp}_q30_frag_size.txt"
     output:
        png = "plot/{smp}_Fragment_size_and_TSS_enrichment_mqc.png"
     threads:5
     message: """--- plot fragments size and TSS enrichment."""
     run:
        shell("python ./atacqc/plot_frag_tss.py {input.frag_size} {output.png}")



rule merge_data:
     input:
        dup = expand("temp_file/{smp}_dup_rate.csv",smp = config["sample"]),
        mapping = expand("temp_file/{smp}_mapping_rate.txt",smp = config["sample"]),
        q30 = expand("temp_file/{smp}_q30_reads_num.csv",smp = config["sample"]),
        organ_mapping = expand("temp_file/{smp}_organell_mapping_rate.csv",smp = config["sample"]),
        tss = expand("tss/{smp}_score.csv",smp = config["sample"])
     output:
        all_dup = "final_result/all_dup.csv",
        map_rate = "final_result/all_mapping.csv",
        q30 = "final_result/all_q30.csv",
        organ_map = "final_result/all_organ.csv",
        tss = "final_result/all_tss.csv",
        final = "final_result/final_qc.csv",
        res = "final_result/ATAC_Quality_Control_Result_mqc.txt"
     message: """--- merge dataframe."""
     run:
        shell("cat  temp_file/*_dup_rate.csv > {output.all_dup}")
        shell("cat  temp_file/*_mapping_rate.txt > {output.map_rate}")
        shell("cat  temp_file/*_q30_reads_num.csv > {output.q30}")
        shell("cat  temp_file/*_organell_mapping_rate.csv > {output.organ_map}")
        shell("cat  tss/*_score.csv >{output.tss}")
        shell("python ./atacqc/merge_data.py {output.final}")
        shell("cat ./atacqc/header.txt {output.final} > {output.res}")
 
rule mutiqc:
     input:
          # dup_rate = expand("temp_file/{smp}_rmdup_rate.txt",smp = config["sample"]),
          # mapping_rate = expand("temp_file/{smp}_mapping_rate.txt",smp = config["sample"]),
          # organ_mapping = expand("temp_file/{smp}_organell_mapping_rate.csv",smp = config["sample"]),
          #idx = expand("temp_file/{smp}_flagstat.txt",smp = config["sample"]),
           #tss = expand("plot/{smp}_tss_mqc.png", smp = config["sample"]),
           table = expand("final_result/ATAC_Quality_Control_Result_mqc.txt"),
           frg = expand("plot/{smp}_Fragment_size_and_TSS_enrichment_mqc.png",smp = config["sample"])
     output:
           "multiqc/multiqc_report.html"
     message: """--- MultiQC."""
     run:
         shell("~/miniconda3/bin/multiqc ./ -o multiqc -f -e general_stats -c ./atacqc/multiqc_config.txt")
