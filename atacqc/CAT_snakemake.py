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
configfile: "config.yaml"
base_dir = config['workdir']
message("The current work directory is" + base_dir)
rule all:
   input:
         "multiqc/multiqc_report.html",
         "final_result/all_dup.csv",
         "final_result/all_mapping.csv",
         "final_result/all_q30.csv",
         "final_result/all_organ.csv",
         #"final_result/all_qc.csv",
         "final_result/final_qc.csv",
         #"final_result/chrom_size.txt",
         "final_result/MergedPeaks.bed",
         "final_result/peakheader.txt",
         "BINDetect/run.log",
         "tss/TSS.pos",
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
         expand("plot/{smp}_Fragment_size_and_TSS_enrichment_mqc.png",smp = config["sample"]),
         expand("peaks/{smp}.macs2.log", smp = config["sample"]),
         expand("signal/{smp}.rpkm.bw", smp = config["sample"]),
         expand("signal/{smp}.cpm.bw", smp = config["sample"]),
         #expand("BINDetect/{smp}.log",smp = config["sample"]),
         expand("ATACorrect/{smp}_q30_corrected.bw", smp = config["sample"]),
         expand("ATACorrect/{smp}_footprints.bw", smp = config["sample"]),
         expand("signal/{smp}.bamcoverage.cpm.log", smp = config["sample"]), 
         expand("signal/{smp}.bamcoverage.rpkm.log", smp = config["sample"]),
         expand("peaks/{smp}_peaks.narrowPeak", smp = config["sample"]),
         expand("peaks/{smp}_peaks_unique.narrowPeak.bed", smp = config["sample"]),
         expand("temp_file/{smp}_frip.csv", smp = config["sample"]),
         expand("temp_file/{smp}_frip_final.csv", smp = config["sample"])
         expand("peaks/{smp}_peaks_unique_rm_blacklist.narrowPeak.bed", smp = config["sample"])


rule Trimming:
     input:
         fwd = lambda wildcards: config["sample"][wildcards.smp]["fwd"],
         rev = lambda wildcards: config["sample"][wildcards.smp]["rev"]
     output:
         fwd = "trimmed/{smp}_R1_001.fastq.paired.gz",
         rev = "trimmed/{smp}_R2_001.fastq.paired.gz",
         fwd1 = "trimmed/{smp}_R1_001.fastq.unpaired.gz",
         rev1 = "trimmed/{smp}_R2_001.fastq.unpaired.gz"
     params: next_adapters = config["adapters"]
     threads:5
     log:     
         "trimmed/{smp}_cut.log"
     shell:
        """
        {config[software][trimmomatic]} PE -threads 5 -phred33 -trimlog {log} {input.fwd} {input.rev} {output.fwd} {output.fwd1} {output.rev} {output.rev1} ILLUMINACLIP:{params.next_adapters}:2:30:10:8:TRUE SLIDINGWINDOW:4:15 MINLEN:30
        """

rule BWA:
    input:
        fwd = "trimmed/{smp}_R1_001.fastq.paired.gz",
        rev = "trimmed/{smp}_R2_001.fastq.paired.gz"
    output:
        bam = "raw_bam/{smp}_raw.bam"
    params:
         Fasta = config['fasta'],
         extra = "-R \"@RG\\tID:{smp}\\tSM:{smp}\\tLB:{smp}\\tPL:illumina\""
    threads: 5
    message: """--- BWA mem reads."""
    run:
        shell("{config[software][bwa]} mem -M -t {threads} -k 32 {params.extra} {params.Fasta} {input.fwd} {input.rev} |{config[software][samtools]} view -bhS >{output.bam}")


rule Sorted_by_name:
    input:
         bam = "raw_bam/{smp}_raw.bam"
    output:
         sort = "temp_file/{smp}_sort.bam"
    threads: 5
    message: """--- Sort bam by name."""
    run:
        shell("{config[software][samtools]} sort -n -@ {threads} {input.bam} -o {output.sort}")


rule Remove_duplicates:
    input:
         bam = "temp_file/{smp}_sort.bam"
    output:
         rmdup_rate = "temp_file/{smp}_rmdup_rate.txt",
         rmdup_bam = "temp_file/{smp}_rmdup.bam"
         #out = "final_result/all_dup_rate.csv"
    threads: 5
    message: """--- Remove PCR duplicates."""
    run:
        shell("{config[software][samtools]}  fixmate -@ {threads} -cm -O bam {input.bam} - |{config[software][samtools]}  sort -@ {threads} |{config[software][samtools]}  markdup -r -s - {output.rmdup_bam} 2>{output.rmdup_rate}")
        shell("{config[software][samtools]}  index {output.rmdup_bam}")


rule Cal_duplication_rate:
    input:
         bam = "temp_file/{smp}_rmdup_rate.txt"
    output:
         out = "temp_file/{smp}_dup_rate.csv"
    threads: 5
    message: """--- Calculate PCR duplicates rate."""
    run:
        shell("python {config[software][atacqc]}/cal_duplicate_rate.py {input.bam} {output.out}")


rule Cal_mapping_rate:
    input:
      bam = "raw_bam/{smp}_raw.bam"
    output:
      mapped_rate = "temp_file/{smp}_mapping_rate.txt",
      flagstat = "temp_file/{smp}_flagstat.txt"
      #out_f = "/final_result/all_mapping_rate.csv"
    threads: 5
    message: """--- Calculate mapping rate."""
    run:
         shell("{config[software][samtools]} flagstat {input.bam} > {output.flagstat}")
         shell("python {config[software][atacqc]}/cal_mapping_rate.py {output.flagstat} {output.mapped_rate}")


rule Remove_low_quality:
    input:
      rmdup = "temp_file/{smp}_rmdup.bam"
    output:
      bam = "q30/{smp}_q30.bam",
      q30_reads = "temp_file/{smp}_q30_reads_num.csv"
    params:
      mapp_qual = config["mapping_quality"]
    threads: 5
    message: """--- generate q30 bam and calculate q30 reads number."""
    run:
     shell("{config[software][samtools]} view -h -q {params.mapp_qual} {input.rmdup}|samtools view -bF 4 - >{output.bam}")
     shell("{config[software][samtools]} index {output.bam}")
     shell("python {config[software][atacqc]}/cal_q30_reads_num.py {output.bam} {output.q30_reads}")


rule Stat:
    input:
      bam = "raw_bam/{smp}_raw.bam"
    output:
      sort = "sort/{smp}_sorted_raw.bam",
      idxstat = "temp_file/{smp}_idxstats.txt"
    threads: 5
    message: """--- organell mapping rate."""
    run:
     shell("{config[software][samtools]} sort -@ {threads} {input.bam} -o {output.sort} && {config[software][samtools]} index {output.sort} && {config[software][samtools]} idxstats -@ {threads} {output.sort} |cut -f1,3 > {output.idxstat}")

rule Organell_mapping_rate:
     input: 
          idxstat = "temp_file/{smp}_idxstats.txt"
     output:
          organell = "temp_file/{smp}_organell_mapping_rate.csv"
     params:
      chr_name = config['organell_chr']   
     run:     
       shell("python {config[software][atacqc]}/cal_organelle_reads_num.py {input.idxstat} {output.organell} --chromosome {params.chr_name}")


rule Fragments_size:
    input:
      bam = "q30/{smp}_q30.bam"
    output:
      txt = "fragments_size/{smp}_q30_frag_size.txt"
    message: """--- fragments size."""
    shell:
      """
      {config[software][samtools]} view {input.bam}|awk '$9>0 && $9<800'|cut -f3,4,9 > {output.txt}
      """

rule Get_cut_sites:
    input:
      bam = "q30/{smp}_q30.bam"
    output:
      txt = "cut_sites/{smp}_q30_cut_sites.bed"
    message: """--- get Tn5 cuts."""
    threads:5
    run:
        shell("python {config[software][atacqc]}/cal_cut_sites.py {input.bam} {output.txt}")


rule TSS_pos:
     input:
          gtf = config['gff3']
     output:
          tss_pos = "tss/TSS.pos"
          #dirs = directory("signal"),
          #peak = directory("peaks")
     message: """--- Getting TSS position."""
     shell:
          """cat {input.gtf}|awk '{{if($3=="gene"){{print $0}}}}'|awk -v OFS='\\t' '{{if($7=="+"){{print $1,$4,$4,"1","+","gene"}} else{{print $1,$5,$5,"1","-","gene"}}}}' > {output.tss_pos}"""


rule TSS_Enrichment:
     input:
       bed = "cut_sites/{smp}_q30_cut_sites.bed",
       tss_pos = "tss/TSS.pos"
     output:
       tss = "tss/{smp}_tss.csv",
       score = "tss/{smp}_score.csv"
     message: """--- TSS enrichment."""
     threads: 5
     run:
        shell("python {config[software][atacqc]}/cal_tss.py {input.bed} {input.tss_pos} {output.tss} {output.score}")


rule Plot_fra_size:
     input:
        frag_size = "fragments_size/{smp}_q30_frag_size.txt",
        tss = "tss/{smp}_tss.csv"
     output:
        png = "plot/{smp}_Fragment_size_and_TSS_enrichment_mqc.png"
     threads:5
     message: """--- plot fragments size and TSS enrichment."""
     run:
        shell("python {config[software][atacqc]}/plot_frag_tss.py {input.frag_size} {input.tss} {output.png}")


#rule chrom_size:
      ## pip install pyfaidx
#      input:
#            fa = config["fasta"]
#      output:
#            chromsize = "final_result/chrom_size.txt"
#        run:
#            shell("{config[software][faidx]} {input.fa} -i chromsizes > {output.chromsize}")


rule Callpeak:
     input:
           bam = "q30/{smp}_q30.bam"
     output:
           log = "peaks/{smp}.macs2.log",
           peak = "peaks/{smp}_peaks.narrowPeak"
     params:
         out_dir = "peaks",
         name = "{smp}",
         genome_size = config['effect_genome_size'],
         extsize = config["extsize"],
         shift = config["shift"]
     message: "----Peak calling"
     run:
         shell("{config[software][macs2]} callpeak --call-summits -f BAM --nomodel --mfold 2 20 -t {input.bam} --outdir {params.out_dir} --name {params.name} --gsize {params.genome_size} --qvalue 0.01  --extsize {params.extsize} --shift {params.shift} > {output.log}")


rule Dedup_peak:
      input:
            peak = "peaks/{smp}_peaks.narrowPeak"
      output:
            unique_peak = "peaks/{smp}_peaks_unique.narrowPeak.bed"
      run:
          shell("cat {input.peak}|awk '!seen[$1,$2,$3]++' > {output.unique_peak}")

rule Remove_blacklist:
      input:
            unique_peak = "peaks/{smp}_peaks_unique.narrowPeak.bed" 
      output:
            final_peak = "peaks/{smp}_peaks_unique_rm_blacklist.narrowPeak.bed"
      params:
            blacklist = config["blacklist"]
            overlap_ratio = config["overlap_ratio"]
      run:
          if config["blacklist"] != "NULL":
            shell("{config[software][bedtools]} intersect -v -f {params.overlap_ratio} -r -a {input.unique_peak} -b {params.blacklist} >{output.final_peak}")
          else:
            shell("echo 'No blacklist exists, skip this step.' > {output.final_peak}")

rule FRiP:
     input:
           peak = "peaks/{smp}_peaks_unique.narrowPeak.bed",
           cuts = "cut_sites/{smp}_q30_cut_sites.bed"
     output:
           frip_res = "temp_file/{smp}_frip.csv"
#           out_name = "{smp}"
     #run:
     #   shell("cat {input.peak} {input.cuts} >{output.frip_res}")
     shell:
          """cat {input.cuts} |awk -v OFS='\\t' '{{print $0,1}}'|bedtools sort -i -|bedtools map -a {input.peak} -b - -c 4 -o sum|cut -f11|paste -sd+ - | bc |awk -v OFS='\\t' '{{print $1}}' > {output.frip_res}"""

rule Cal_frip:
     input:
          frip = "temp_file/{smp}_frip.csv",
          q30 = "temp_file/{smp}_q30_reads_num.csv"
     output:
          frip_res = "temp_file/{smp}_frip_final.csv"
     message: "----Calculating the FRiP----"
     run:
        shell("python {config[software][atacqc]}/cal_frip.py {input.frip} {input.q30} {output.frip_res}")

rule Bam2bw:
     input:
           bam = "q30/{smp}_q30.bam"
    # params:
     #      outdir="signal"
     output:
           log = "signal/{smp}.bamcoverage.cpm.log",
           log2 = "signal/{smp}.bamcoverage.rpkm.log",
           bw = "signal/{smp}.cpm.bw",
           bw2 = "signal/{smp}.rpkm.bw"
     params:
           smoothLength = config["smoothLength"],
           binSize = config["binSize"],
           threads = config["threads"]
     message: """--- Bam to bigwig."""
     run:
         shell("{config[software][bamCoverage]} --bam {input.bam} -o {output.bw2} --binSize {params.binSize} --smoothLength {params.smoothLength} --normalizeUsing RPKM -p {params.threads} > {output.log2}")
         shell("{config[software][bamCoverage]} --bam {input.bam} -o {output.bw} --binSize {params.binSize} --smoothLength {params.smoothLength} --normalizeUsing CPM -p {params.threads} > {output.log}")

rule Merge_data:
     input:
        dup = expand("temp_file/{smp}_dup_rate.csv",smp = config["sample"]),
        mapping = expand("temp_file/{smp}_mapping_rate.txt",smp = config["sample"]),
        q30 = expand("temp_file/{smp}_q30_reads_num.csv",smp = config["sample"]),
        organ_mapping = expand("temp_file/{smp}_organell_mapping_rate.csv",smp = config["sample"]),
        tss = expand("tss/{smp}_score.csv",smp = config["sample"]),
        frip = expand("temp_file/{smp}_frip_final.csv",smp = config["sample"])
     output:
        all_dup = "final_result/all_dup.csv",
        map_rate = "final_result/all_mapping.csv",
        q30 = "final_result/all_q30.csv",
        organ_map = "final_result/all_organ.csv",
        tss = "final_result/all_tss.csv",
        final = "final_result/final_qc.csv",
        all_frip = "final_result/all_frip.csv",
        res = "final_result/ATAC_Quality_Control_Result_mqc.txt"
     message: """--- merge dataframe."""
     run:
        shell("cat  temp_file/*_dup_rate.csv > {output.all_dup}")
        shell("cat  temp_file/*_mapping_rate.txt > {output.map_rate}")
        shell("cat  temp_file/*_q30_reads_num.csv > {output.q30}")
        shell("cat  temp_file/*_organell_mapping_rate.csv > {output.organ_map}")
        shell("cat  tss/*_score.csv >{output.tss}")
        shell("cat temp_file/*_frip_final.csv > {output.all_frip}" )
        shell("python {config[software][atacqc]}/merge_data.py {output.final}")
        shell("cat {config[software][atacqc]}/header.txt {output.final} > {output.res}")


rule Mutiqc:
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
         shell("{config[software][multiqc]} ./ -o multiqc -f -e general_stats -c {config[software][atacqc]}/multiqc_config.txt")


rule MergePeaks:
      input:
           peak = expand("peaks/{smp}_peaks_unique.narrowPeak.bed",smp = config["sample"])
      output:
           merge_peaks = "final_result/MergedPeaks.bed",
           header = "final_result/peakheader.txt"
      message: """--- Merge Peaks."""
      shell:
          """
          cat peaks/*_peaks_unique.narrowPeak.bed|{config[software][bedtools]} sort -i -|{config[software][bedtools]} merge -i - >{output.merge_peaks}
          echo -e "peak_chr\tpeak_start\tpeak_end" > {output.header}
          """

rule ATACorrect:
     input:
           bam = "q30/{smp}_q30.bam",
           genome = config['fasta'],
           peaks = "final_result/MergedPeaks.bed"
     output:
           bw = "ATACorrect/{smp}_q30_corrected.bw"
     params:
           out_dir = "ATACorrect"
     message: """--- ATACorrect."""
     threads: 5
     run:
         shell("{config[software][TOBIAS]} ATACorrect --bam {input.bam} --genome {input.genome} --peaks {input.peaks} --outdir {params.out_dir} --cores {threads}")


rule FootprintScores:
     input:
           bw = "ATACorrect/{smp}_q30_corrected.bw",
           peaks = "final_result/MergedPeaks.bed"
     output:
           ftscore = "ATACorrect/{smp}_footprints.bw"
     message: """--- FootprintScores."""
     threads: 5
     run: 
         shell("{config[software][TOBIAS]} FootprintScores --signal {input.bw} --regions {input.peaks} --output {output.ftscore} --cores {threads}")

# rule BINDetect:
#      input:
#            motif = config['motif'],
#            signal = "ATACorrect/{smp}_footprints.bw",
#            genome = config['fasta'],
#            peaks = "final_result/MergedPeaks.bed",
#            header = "final_result/peakheader.txt"
#      output:
#            log = "BINDetect/{smp}/{smp}.log"
#      params:
#            outdir = "BINDetect/{smp}",
#            cond_names = "{smp}"
#      message: """--- BINDetect."""
#      threads: 5
#      run:
#          shell("TOBIAS BINDetect --motifs {input.motif} --signals {input.signal} --genome {input.genome} --peaks {input.peaks} --peak_header {input.header} --outdir {params.outdir} --cond_names {params.cond_names} --cores {threads} >>{output.log}")


rule BINDetect_compare:
     input:
           motif = config['motif'],
           signal = expand("ATACorrect/{smp}_footprints.bw",smp = config["sample"]),
           genome = config['fasta'],
           peaks = "final_result/MergedPeaks.bed",
           header = "final_result/peakheader.txt"
     output:
           log = "BINDetect/run.log"
     params:
            outdir = "BINDetect"
     message: """--- BINDetect Comparing."""
     threads: 10
     run:
         shell("{config[software][TOBIAS]} BINDetect --motifs {input.motif} --signals ATACorrect/*_footprints.bw --genome {input.genome} --peaks {input.peaks} --peak_header {input.header} --outdir {params.outdir} --cores {threads} >> {output.log}")
