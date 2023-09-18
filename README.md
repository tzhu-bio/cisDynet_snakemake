# CAT_snakemake
This is a chromatin accessibility data preprocessing process designed for CAT package.
## 
[Report Example](https://htmlpreview.github.io/?https://github.com/tzhu-bio/CAT_snakemake/blob/main/multiqc_report.html)

## ATAC-seq pipeline
After you have set up your environment and config.yaml file, you can run the process using the following command
```python
snakemake -s CAT_snakemake.py -j 20
```
## RNA-seq quantification pipeline

Before you run the RNA quantification pipeline, you should install the [bowtie2](https://github.com/BenLangmead/bowtie2) and [RSEM](https://github.com/deweylab/RSEM) in your environment.
And RSEM index needs to be built first.
```python
snakemake -s rna_pipe_snakemake.py -j 20
```
