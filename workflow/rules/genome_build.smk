build = config["build"]

rule genome_build:
    # input: "resources/ref_genomes/build/genome.fa"
    input: f"resources/ref_genomes/{build}/genome.fa"
    # output: "resources/ref_genomes/build/Bowtie2/genome"
    output: f"resources/ref_genomes/{build}/Bowtie2/genome"
    log: f"logs/rule/analysis/{build}/log/bowtie2_build.log"
    benchmark: f"logs/rule/analysis/{build}/log/bowtie2_build.benchmark.txt"
    conda:
        "../envs/align.yaml"
    shell: 
        """
        (echo "`date -R`: Building indexes..." &&
        bowtie2-build  \
        {input} \
        {output} &&
        echo "`date -R`: Success! Indexes are build." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """
