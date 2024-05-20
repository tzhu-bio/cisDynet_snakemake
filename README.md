# cisDynet_snakemake
## Installation
You can use the following command to configure the environment for the pipeline.
```python
conda env create -f cisDynet_env.yaml
```
This is a chromatin accessibility data preprocessing process designed for cisDynet package. After running the pipeline, you will get an HTML report.
## 
[Report Example](https://htmlpreview.github.io/?https://github.com/tzhu-bio/CAT_snakemake/blob/main/multiqc_report.html)

## ATAC-seq pipeline
Before running this pipeline, you need to specify the absolute paths to some necessary software in the config.yaml file.

After you have set up your environment and config.yaml file, you can run the process using the following command
```python
snakemake -s cisDynet_snakemake.py -j 20
```
## RNA-seq quantification pipeline

Before you run the RNA quantification pipeline, you should install the [bowtie2](https://github.com/BenLangmead/bowtie2) and [RSEM](https://github.com/deweylab/RSEM) in your environment.
And RSEM index needs to be built first.
```python
snakemake -s rna_pipe_snakemake.py -j 20
```
## Blacklists

The presence of some anomalous regions on the genome allows for extremely high signals. Most of these regions are due to the presence of repetitive sequences, genome assembly errors, and so on. In mammals, “problematic” regions have been inferred and manually checked and called blacklists, and are widely used in the analysis of genomic data. However, there are no systematic blacklists for plant genomes at present. To fill this gap, we used the [greenscreen](https://academic.oup.com/plcell/article/34/12/4795/6705244) in combination with data collected from our [ChIP-Hub database](https://www.nature.com/articles/s41467-022-30770-1) to obtain “problematic” lists for five plant species: Arabidopsis thaliana, rice, maize, soybean, and tomato. (In the *blacklists* directory)

## Publication
Zhu, Tao, Zhou, Xinkai, You, Yuxin, Wang, Lin, He, Zhaohui, and Chen, Dijun. 2023. “ cisDynet: An Integrated Platform for Modeling Gene-Regulatory Dynamics and Networks.” iMeta e152. https://doi.org/10.1002/imt2.152

<a href="mailto:tzhubio@gmail.com">
  <img src="https://github.com/blackcater/blackcater/raw/main/images/social-gmail.svg" height="40" />
</a>

![Anurag's GitHub stats](https://github-readme-stats.vercel.app/api?username=tzhu-bio&show_icons=true&theme=radical)
