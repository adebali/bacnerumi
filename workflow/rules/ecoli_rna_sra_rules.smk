build_ = config["build"]
project_ = config["project"]

rule sra_se:
    output:
        "resources/samples/{sample}.fastq.gz", 
    params:
        srr="{sample}",
        name="{sample}",
    log:
        "logs/rule/analysis/{sample}/{sample}_se_sra.log",
    # benchmark:
    #     "logs/rule/analysis/{sample}/{sample}_se_sra.benchmark.txt",
    resources:
        memory="16GB",
        cpu=6
    wildcard_constraints:
        samples='|'.join([x for x in config["sample"]])
    conda:
        "../envs/sra_toolkit.yaml"
    # threads:
    #     6
    shell:
        """
        touch resources/samples/{sample}.fastq
        touch {log}

        (echo "`date -R`: Downloading {sample} files..." &&
        prefetch {sample} \
        -O resources/samples/ &&
        vdb-validate resources/samples/{sample} &&
        fastq-dump \
        resources/samples/{sample} \
        --outdir resources/samples/ && \
        echo "`date -R`: Download is successful!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) \
        >> {log} 2>&1 

        (echo "`date -R`: Zipping srr file..." &&
        gzip resources/samples/{sample}.fastq &&
        echo "`date -R`: Zipping is successful!" || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }} ) \
        >> {log} 2>&1
        """

rule trimmomatic:
    input:
        "resources/samples/{sample}.fastq.gz"  # input and output can be uncompressed or compressed
    output:
        "results/{sample}/{sample}_trimmed.fastq.gz"
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
    conda:
        "../envs/trimmomatic.yaml"
    shell:
        """
        (echo "`date -R`: Trim adaptors..." &&
        trimmomatic SE -phred33 {input} {output} ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        echo "`date -R`: Trim is successful!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) \
        >> {log} 2>&1 
        """

rule star_index:
    input:
        fasta=f"resources/ref_genomes/{build_}/genome.fa",
    output: directory(f"resources/ref_genomes/{build_}/star"),
    message:
        "Testing STAR index"
    threads: 1
    params:
        extra="",
    resources:
        memory="16GB",
        cpu=1
    log:
        "logs/star_index.log",
    wrapper:
        "v2.6.0/bio/star/index"

rule star_se:
    input:
        fq1="results/{sample}/{sample}_trimmed.fastq.gz",
        # path to STAR reference genome index
        idx=f"resources/ref_genomes/{build_}/star",
    output:
        # see STAR manual for additional output files
        aln="results/{sample}/{sample}_star.bam",
        log="logs/se/{sample}/Log.out",
        log_final="logs/se/{sample}/Log.final.out",
        unmapped="results/{sample}/{sample}_star/unmapped.fastq",
    log:
        "logs/{sample}.log",
    params:
        # optional parameters
        extra="--outSAMtype BAM Unsorted",
    threads: 8
    resources:
        memory="16GB",
        cpu=8
    wrapper:
        "v3.12.0/bio/star/align"

# rule filterBam:
#     input: 
#         # bam = rules.tophat.output.bam,
#         bam = rules.star_se.output.aln,
#         ref_genome_index=rules.star_index.output,
#     output:
#         bam = "results/{sample}/{sample}_star_filtered.bam",
#     log:
#         "logs/filter/{sample}.log",
#     params:
#         extra="",  # optional parameters
#     threads: 8  # Use at least two threads
#     resources:
#         memory="16GB",
#         cpu=8
#     conda:
#         "../envs/align.yaml"
#     shell:
#         """
#         (echo "`date -R`: Filtering bam file..." && 
#         samtools view -h -t {input.ref_genome_index} {input.bam} |\
#         samtools sort -n - | samtools view -bf 0x2 - -o {output.bam} &&
#         echo "`date -R`: Success! Bam file filtered." || 
#         {{ echo "`date -R`: Process failed..."; rm {output.bam}; exit 1; }}  ) >> {log} 2>&1
#         """

#         # samtools sort -n {output.sam} -o {output.bam} &&
#         # samtools view -bf 0x2 {output.bam}| bedtools bamtobed -i stdin -bedpe > {output.bedpe} &&
#         # samtools index {output.bam} &&

rule bam2bed:
    input: rules.star_se.output.aln,
    # output: "results/{sample}/{sample}_filtered.bed"
    output: "results/{sample}/{sample}.bed"
    log:
        "logs/bam2bedpe/{sample}.log",
    resources:
        memory="16GB",
        cpu=1
    conda:
        "../envs/bam2bed.yaml"
    shell:  
        """
        (echo "`date -R`: Processing bam file..." && 
        bedtools bamtobed -i {input} | sort -k1,1 -k2,2n -k3,3n >  {output} &&
        echo "`date -R`: Success! Bam file converted to bed format." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1
        """

rule bed2geneCounts_rna:
    input:
        reads=rules.bam2bed.output,
        genes=f"resources/ref_genomes/{build_}/genes.bed",
    output: "results/{sample}/{sample}_sense.tsv",
    params:
        fields = lambda w: getSampleFields(config["meta"][w.sample], "sense"),
        extra = "-s"
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
        {params.extra} \
        -a {input.genes} \
        -b {input.reads} \
        -F 0.5 \
        | awk '{{print $0, "\\t{params.fields}"}}' \
        > {output} &&
        echo "`date -R`: Success! {output} is generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """
        )

rule bed2geneCounts_rna_antisense:
    input:
        reads=rules.bam2bed.output,
        genes=f"resources/ref_genomes/{build_}/genes.bed",
    output: "results/{sample}/{sample}_antisense.tsv",
    params:
        fields = lambda w: getSampleFields(config["meta"][w.sample], "antisense"),
        extra = "-S"
    log: "logs/bed2geneCounts_rna_antisense/{sample}.log"
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
        {params.extra} \
        -a {input.genes} \
        -b {input.reads} \
        -F 0.5 \
        | awk '{{print $0, "\\t{params.fields}"}}' \
        > {output} &&
        echo "`date -R`: Success! {output} is generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """
        )

rule mergeTranscriptCounts:
    input:
        lambda w: input4merge(config["sample"]),
    output:
        out=f"results/{project_}/transcriptCounts.tsv",
    log: f"logs/rule/analysis/{project_}/transcriptCounts.log",
    benchmark: f"logs/rule/analysis/{project_}/transcriptCounts.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    shell:  
        """
        (echo "`date -R`: Merge transcript counts..." &&
        cat {input} > {output.out} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.out}; exit 1; }} ) > {log} 2>&1
        """

