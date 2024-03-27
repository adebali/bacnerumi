build_ = config["build"]
project_ = config["project"]
include: "faidx.smk"


rule mergeTSNTScounts:
    input:
        lambda w: input4merge(config["sample"]),
    output: "results/{project}/readCounts.tsv",
    log: "logs/{project}/readCounts.log",
    resources:
        memory="16GB",
        cpu=1
    shell:  
        """
        (echo "`date -R`: Merge counts from samples..." &&
        cat {input} > {output} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }} ) > {log} 2>&1
        """




# rule bowtie2:
#     input:
#         sample= ["resources/RNA-seq/ecoli/{sample}/{sample}_1.fq.gz", "resources/RNA-seq/ecoli/{sample}/{sample}_2.fq.gz"],
#         idx=multiext(
#             f"resources/ref_genomes/{build_}/Bowtie2/genome",
#             ".1.bt2",
#             ".2.bt2",
#             ".3.bt2",
#             ".4.bt2",
#             ".rev.1.bt2",
#             ".rev.2.bt2",
#         ),
#     output:
#         "results/rna/ecoli/{sample}.bam",
#     log:
#         "logs/bowtie2/rna/{sample}.log",
#     params:
#         extra="",  # optional parameters
#     threads: 8  # Use at least two threads
#     resources:
#         memory="16GB",
#         cpu=8
#     wrapper:
#         "v2.6.0/bio/bowtie2/align"

# rule genome_build:
#     input: report(f"resources/ref_genomes/{build_}/genome.fa", category="genome")
#     params: 
#         base = f"resources/ref_genomes/{build_}/Bowtie2/2.2.5.0/genome"
#     output:
#         index_1 = f"resources/ref_genomes/{build_}/Bowtie2/2.2.5.0/genome.1.bt2",
#         index_2 = f"resources/ref_genomes/{build_}/Bowtie2/2.2.5.0/genome.2.bt2",
#         index_3 = f"resources/ref_genomes/{build_}/Bowtie2/2.2.5.0/genome.3.bt2",
#         index_4 = f"resources/ref_genomes/{build_}/Bowtie2/2.2.5.0/genome.4.bt2",
#         index_rev1 = f"resources/ref_genomes/{build_}/Bowtie2/2.2.5.0/genome.rev.1.bt2",
#         index_rev2 = f"resources/ref_genomes/{build_}/Bowtie2/2.2.5.0/genome.rev.2.bt2",
#     log: f"logs/rule/analysis/{build_}/log/bowtie2_2.2.5.0_build.log"
#     benchmark: f"logs/rule/analysis/{build_}/log/bowtie2_2.2.5.0_build.benchmark.txt"
#     conda:
#         "../envs/tophat.yaml"
#     resources:
#         memory="16GB",
#         cpu=1
#     shell: 
#         """
#         (echo "`date -R`: Building indexes..." &&
#         bowtie2-build  \
#         {input} \
#         {params.base} &&
#         echo "`date -R`: Success! Indexes are build." || 
#         {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
#         """

rule tophat:
    input:
        sample= ["resources/RNA-seq/ecoli/{sample}/{sample}_1.fq.gz", "resources/RNA-seq/ecoli/{sample}/{sample}_2.fq.gz"],
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
        bam = "results/rna/ecoli/tophat_out_{sample}/accepted_hits.bam",
    log:
        "logs/bowtie2/rna/{sample}.log",
    params:
        base = f"resources/ref_genomes/{build_}/Bowtie2/genome",
        outputDir = "results/rna/ecoli/tophat_out_{sample}",
        threads = 8, # Use at least two threads
        extra = "--no-discordant --library-type fr-firststrand"
    resources:
        memory="16GB",
        cpu=8
    conda:
        "../envs/tophat.yaml"
    shell:
        """
        (echo "`date -R`: Running tophat aligner..." && 
        tophat2 {params.extra} -o {params.outputDir} -p {params.threads} {params.base} {input.sample} &&
        echo "`date -R`: Success! Bam file created." || 
        {{ echo "`date -R`: Process failed..."; rm {output.bam}; exit 1; }}  ) >> {log} 2>&1
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

rule star_pe_multi:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1=["resources/RNA-seq/ecoli/{sample}/{sample}_1.fq.gz"],
        # paired end reads needs to be ordered so each item in the two lists match
        fq2=["resources/RNA-seq/ecoli/{sample}/{sample}_2.fq.gz"],  #optional
        # path to STAR reference genome index
        idx=rules.star_index.output,
    output:
        # see STAR manual for additional output files
        aln="results/rna/ecoli/star/pe/{sample}/pe_aligned.sam",
        log="results/rna/ecoli/logs/pe/{sample}/Log.out",
        sj="results/rna/ecoli/star/pe/{sample}/SJ.out.tab",
        unmapped=["results/rna/ecoli/star/pe/{sample}/unmapped.1.fastq.gz","star/pe/{sample}/unmapped.2.fastq.gz"],
    log:
        "logs/pe/{sample}.log",
    params:
        # optional parameters
        extra="",
    resources:
        memory="16GB",
        cpu=8
    threads: 8
    wrapper:
        "v2.6.0/bio/star/align"

rule filterBam:
    input: 
        # bam = rules.tophat.output.bam,
        bam = rules.star_pe_multi.output.aln,
        ref_genome_index=rules.samtools_index.output,
    output:
        # bam = "results/rna/ecoli/{sample}_filtered.bam",
        bam = "results/rna/ecoli/{sample}_star_filtered.bam",
        # bedpe = "results/rna/ecoli/{sample}_filtered.bedpe",
        # bai = "results/rna/ecoli/{sample}_filtered.bam.bai",
        # sam = "results/rna/ecoli/{sample}_star_filtered.sam",
    log:
        "logs/filter/rna/{sample}.log",
    params:
        extra="",  # optional parameters
    threads: 8  # Use at least two threads
    resources:
        memory="16GB",
        cpu=8
    conda:
        "../envs/align.yaml"
    shell:
        """
        (echo "`date -R`: Filtering bam file..." && 
        samtools view -h -t {input.ref_genome_index} {input.bam} |\
        samtools sort -n - | samtools view -bf 0x2 - -o {output.bam} &&
        echo "`date -R`: Success! Bam file filtered." || 
        {{ echo "`date -R`: Process failed..."; rm {output.bam}; exit 1; }}  ) >> {log} 2>&1
        """

        # samtools sort -n {output.sam} -o {output.bam} &&
        # samtools view -bf 0x2 {output.bam}| bedtools bamtobed -i stdin -bedpe > {output.bedpe} &&
        # samtools index {output.bam} &&

rule bam2bedpe:
    input: rules.filterBam.output.bam,
    output: "results/rna/ecoli/{sample}_filtered.bedpe"
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
        bedtools bamtobed -i {input} -bedpe | sort -k1,1 -k2,2n -k3,3n >  {output} &&
        echo "`date -R`: Success! Bam file converted to bed format." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1
        """

rule bedpe2bed:
    input: 
        bedpe = rules.bam2bedpe.output,
        sam = rules.star_pe_multi.output.aln
    output: "results/rna/ecoli/{sample}.bed"
    log:
        "logs/bedpe2bed/{sample}.log",
    resources:
        memory="16GB",
        cpu=1
    shell:  
        """
        (echo "`date -R`: Processing bam file..." && 
        python3 workflow/scripts/bedpe2bed.py {input.bedpe} {input.sam} | sort -k1,1 -k2,2n -k3,3n > {output} &&
        echo "`date -R`: Success! Bam file converted to bed format." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1
        """

# rule bam2sam:
#     input: rules.bowtie2.output,
#     output: "results/rna/ecoli/{sample}.sam"
#     log:
#         "logs/bam2sam/{sample}.log",
#     resources:
#         memory="16GB",
#         cpu=1
#     conda:
#         "../envs/bam2bed.yaml"
#     shell:  
#         """
#         (echo "`date -R`: Processing bam file..." && 
#         samtools view {input} -o {output} &&
#         echo "`date -R`: Success! Bam file converted to sam format." || 
#         {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1
#         """

rule sep_strands_original:
    input:
        rules.bedpe2bed.output,
    output:
        plus="results/rna/ecoli/{sample}_plus.bed",
        minus="results/rna/ecoli/{sample}_minus.bed",
    log:
        "logs/sep_strands/{sample}.log",
    resources:
        memory="16GB",
        cpu=1
    shell:  
        """
        (echo "`date -R`: Separating plus stranded reads..." &&
        awk '{{if($6=="+"){{print}}}}' {input} > {output.plus} &&
        echo "`date -R`: Success! Reads are separated." || 
        {{ echo "`date -R`: Process failed..."; rm {output.plus}; exit 1; }}  )  > {log} 2>&1

        (echo "`date -R`: Separating minus stranded reads..." &&
        awk '{{if($6=="-"){{print}}}}' {input} > {output.minus} &&
        echo "`date -R`: Success! Reads are separated." || 
        {{ echo "`date -R`: Process failed..."; rm {output.minus}; exit 1; }}  )  >> {log} 2>&1
        """

rule genomecov:
    input:
        bed="results/rna/ecoli/{sample}_{strand}.bed",
        ref_genome_index=rules.samtools_index.output,
    output:
        bdg="results/rna/ecoli/{sample}_{strand}.bdg",
        bw=report("results/rna/ecoli/{sample}_{strand}.bw", category="bigwig"),
    params:
        read=lambda w, input: mappedReads(str(input['bed'])),
        bwtitle=lambda w: f'results/{project_}/bigwig/' + getTitle(w.sample) + f'_{w.strand}.bw',
    log:
        "logs/rule/analysis/{sample}_{strand}_genomeCov.log",
    benchmark:
        "logs/rule/analysis/{sample}_{strand}_genomeCov.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    conda:
        "../envs/bedtools.yaml"
    shell: 
        re.sub(' +', ' ',
        """
        (echo "`date -R`: Calculating genome coverage of {input.bed}..." &&
        bedtools genomecov \
        -i {input.bed} \
        -g {input.ref_genome_index} \
        -bg \
        -scale $(echo {params.read} | awk '{{print 1000000/$1}}') \
        > {output.bdg} \
        && \
        bedGraphToBigWig {output.bdg} {input.ref_genome_index} {output.bw} && \
        cp {output.bw} {params.bwtitle} &&
        echo "`date -R`: Success! Genome coverage for located damaged is calculated." || 
        echo "`date -R`: Process failed...") > {log} 2>&1 
        """
        )

rule bed2geneCounts_rna:
    input:
        reads=rules.bedpe2bed.output,
        genes=f"resources/ref_genomes/{build_}/genes.bed",
    output: "results/rna/ecoli/{sample}_sense.tsv",
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
        reads=rules.bedpe2bed.output,
        genes=f"resources/ref_genomes/{build_}/genes.bed",
    output: "results/rna/ecoli/{sample}_antisense.tsv",
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
        lambda w: input4mergeTranscriptCounts(config["sample"]),
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


# rule deneme:
#     output:"results/{project}/readCounts.tsv"
#     shell: "touch {output}"