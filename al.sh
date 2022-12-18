alias smm='snakemake --snakefile myco.smk -p --use-conda --profile ./config/slurm  --rerun-incomplete'
alias rna='snakemake --snakefile myco_rna.smk -p --use-conda --profile ./config/slurm  --rerun-incomplete'
alias sme='snakemake --snakefile ecoli.smk -p --use-conda --profile ./config/slurm  --rerun-incomplete'
alias smmr='snakemake --snakefile myco.smk --report myco.zip && rm -fr ../playground/myco && unzip -o myco.zip -d ../playground'
alias smer='snakemake --snakefile ecoli.smk --report ecoli.zip --report-stylesheet report/custom-stylesheet.css && rm -fr ../playground/ecoli && unzip -o ecoli.zip -d ../playground'
alias publish='cd ../playground && git fetch && git add . && git commit -a -m "update" && git push origin main && cd -'