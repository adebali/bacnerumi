alias smm='snakemake --snakefile myco.smk -p --use-conda --profile ./config/slurm  --rerun-incomplete'
alias sme='snakemake --snakefile ecoli.smk -p --use-conda --profile ./config/slurm  --rerun-incomplete'
alias smmr='snakemake --snakefile myco.smk --report myco.zip'
alias smer='snakemake --snakefile ecoli.smk --report ecoli.zip'