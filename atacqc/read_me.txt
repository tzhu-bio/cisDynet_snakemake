snakemake -s run_snakemake.py --cluster "bsub -R span[hosts=1] -q normal -o logs.txt " -j 5
