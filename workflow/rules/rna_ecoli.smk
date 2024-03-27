build_ = config["build"]

rule get_fastq_se_gz:
    output:
        "results/rna/se/raw/{accession}.fastq.gz"
    log:
        "logs/se/{accession}.gz.log"
    params:
        extra="--skip-technical"
    threads: 6
    resources:
        memory="64GB",
        cpu=6
    wrapper:
        "v1.23.4/bio/sra-tools/fasterq-dump"

rule trimmomatic:
    input:
        "results/rna/se/raw/{sample}.fastq.gz"  # input and output can be uncompressed or compressed
    output:
        "trimmed/{sample}.fastq.gz"
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["TRAILING:3"],
        # optional parameters
        extra="",
        # optional compression levels from -0 to -9 and -11
        compression_level="-9"
    threads:
        32
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        memory="4G",
        cpu=32
    wrapper:
        "v1.23.4/bio/trimmomatic/se"

rule bowtie2:
    input:
        sample=rules.trimmomatic.output,
        idx=multiext(
            f"resources/ref_genomes/{build_}/Bowtie2/genome",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        "results/rna/se/mapped/{sample}.bam",
    log:
        "logs/bowtie2/rna/{sample}.log",
    params:
        extra="",  # optional parameters
    threads: 8  # Use at least two threads
    resources:
        memory="16GB",
        cpu=8
    wrapper:
        "v1.19.1/bio/bowtie2/align"

rule bam2bed:
    input: rules.bowtie2.output,
    output: "results/rna/se/mapped/{sample}.bed"
    log:
        "logs/bam2bed/{sample}.log",
    resources:
        memory="16GB",
        cpu=1
    conda:
        "../envs/bam2bed.yaml"
    shell:  
        """
        (echo "`date -R`: Processing bam file..." && 
        bedtools bamtobed -i {input} | sort -k1,1 -k2,2n -k3,3n > {output} &&
        echo "`date -R`: Success! Bam file converted to bed format." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1
        """

rule bed2geneCounts_rna:
    input:
        reads=rules.bam2bed.output,
        genes=f"resources/ref_genomes/{build_}/genes.bed",
    output: "results/rna/se/mapped/{sample}.tsv",
    params:
        fields = lambda w: getSampleFields(config["meta"][w.sample])
    log: "logs/bed2geneCounts_rna/{sample}.log",
    resources:
        memory="16GB",
        cpu=1
    conda:
        "../envs/bedtools.yaml"
    shell:
        re.sub(' +', ' ', 
        """
        (echo "`date -R`: Counting reads..." &&
        bedtools intersect -c \
        -a {input.genes} \
        -b {input.reads} \
        -F 0.5 \
        | awk '{{print $0, "\\t{params.fields}"}}' \
        > {output} &&
        echo "`date -R`: Success! {output} is generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        """
        )

rule mergeTSNTScounts:
    input:
        lambda w: input4merge(config["sample"]),
    output:
        out=f"results/rna_ecoli/readCounts.tsv",
    log: f"logs/rna_ecoli/readCounts.log",
    resources:
        memory="16GB",
        cpu=1
    shell:  
        """
        (echo "`date -R`: Merge counts from samples..." &&
        cat {input} > {output.out} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.out}; exit 1; }} ) > {log} 2>&1
        """
