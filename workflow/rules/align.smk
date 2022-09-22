build_ = config["build"]

rule bowtie2_se:
    input:
        sample = rules.cutadapt_se.output.fastq,
        bowtie2=f"resources/ref_genomes/{build_}/Bowtie2/genome.1.bt2",
    output:
        sam=temp("results/{sample}/{sample}.sam"),
        bam="results/{sample}/{sample}.bam",
    params:
        ref_genome=rules.genome_build.output,
        extra="--seed 1",
        sample_id="{sample}"
    # threads: 16  
    log:
        "logs/rule/analysis/{sample}/{sample}_bowtie2.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_bowtie2.benchmark.txt",
    resources:
        memory="64GB",
        cpu=16
    conda:
        "../envs/align.yaml"
    shell:  
        """
        (echo "`date -R`: Aligning fastq file..." &&
        bowtie2 \
        --threads {resources.cpu} \
        {params.extra} \
        -x {params.ref_genome} \
        --rg-id {params.sample_id} \
        --rg SM:{params.sample_id} \
        -U {input.sample} -S {output.sam} &&
        echo "`date -R`: Success! Alignment is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1

        (echo "`date -R`: Converting sam to bam..." &&
        samtools view -Sbh {output.sam} | samtools sort -o {output.bam} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; rm {output.bam}; exit 1; }}  )  >> {log} 2>&1
        """

rule bam2bed_se:
    input:
        bam= rules.bowtie2_se.output.bam,
        fasta="resources/ref_genomes/" + config["build"] + "/singletons.fa"
    output:
        bed="results/{sample}/{sample}.bed",
        bam="results/{sample}/{sample}_sorted.bam",
        sam="results/{sample}/{sample}_sorted.sam",
        idx="results/{sample}/{sample}_sorted.bam.bai",
    params:
        q_trim=config["samtools_q_trim_se"], 
    log:
        "logs/rule/analysis/{sample}/{sample}_bam2bed.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_bam2bed.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    conda:
        "../envs/bedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Sorting (coordinates) bam file..." &&
        samtools view -h {input.bam} | python3 workflow/scripts/filterSam.py -f {input.fasta} -o {output.sam} && 
        samtools view -S {output.sam} -b | samtools sort > {output.bam} &&
        echo "`date -R`: Success! Bam file is filtered and sorted." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Index bam file..." &&
        samtools index {output.bam} {output.idx} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Processing bam file..." && 
        samtools view {params.q_trim} -b {output.bam} |&
        bedtools bamtobed | sort -k1,1 -k2,2n -k3,3n > {output.bed} &&
        echo "`date -R`: Success! Bam file converted to bed format." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1
        """


rule sep_strands:
    input:
        rules.bam2bed_se.output.bed,
    output:
        plus="results/{sample}/{sample}_plus.bed",
        minus="results/{sample}/{sample}_minus.bed",
    log:
        "logs/rule/analysis/{sample}/{sample}_sepStrands.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_sepStrands.benchmark.txt",
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


use rule sep_strands as sep_strands_sim with:
    input: 
        "results/{sample}/simulation/{sample}.bed"
    output:
        plus="results/{sample}/simulation/{sample}_plus.bed",
        minus="results/{sample}/simulation/{sample}_minus.bed",
    log:
        "logs/rule/analysis/{sample}/simulation/{sample}_sepStrands.log",
    benchmark:
        "logs/rule/analysis/{sample}/simulation/{sample}_sepStrands.benchmark.txt",

rule genomecov_plus:
    input:
        bed=rules.sep_strands.output.plus, 
        ref_genome_index=rules.samtools_index.output,
    output:
        bdg=temp("results/{sample}/{sample}_pos.bdg"),
        bw="results/{sample}/{sample}_pos.bw",
    params:
        read=lambda w, input: mappedReads(input['bed']),
    log:
        "logs/rule/analysis/{sample}/{sample}_genomeCovP.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_genomeCovP.benchmark.txt",
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
        bedGraphToBigWig {output.bdg} {input.ref_genome_index} {output.bw} &&
        echo "`date -R`: Success! Genome coverage for located damaged is calculated." || 
        echo "`date -R`: Process failed...") > {log} 2>&1
        """
        )

use rule genomecov_plus as genomecov_minus with:
    input:
        bed=rules.sep_strands.output.minus, 
        ref_genome_index=rules.samtools_index.output,
    output:
        bdg=temp("results/{sample}/{sample}_neg.bdg"),
        bw="results/{sample}/{sample}_neg.bw",
    log:
        "logs/rule/analysis/{sample}/{sample}_genomeCovN.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_genomeCovN.benchmark.txt",

use rule genomecov_plus as genomecov_plus_sim with:
    input:
        bed=rules.sep_strands.output.plus, 
        ref_genome_index=rules.samtools_index.output,
    output:
        bdg=temp("results/{sample}/simulation/{sample}_pos.bdg"),
        bw="results/{sample}/simulation/{sample}_pos.bw",
    params:
        read=lambda w, input: mappedReads(input['bed']),
    log:
        "logs/rule/analysis/{sample}/simulation/{sample}_genomeCovP.log",
    benchmark:
        "logs/rule/analysis/{sample}/simulation/{sample}_genomeCovP.benchmark.txt",


use rule genomecov_minus as genomecov_minus_sim with:
    input:
        bed=rules.sep_strands.output.minus, 
        ref_genome_index=rules.samtools_index.output,
    output:
        bdg=temp("results/{sample}/simulation/{sample}_neg.bdg"),
        bw="results/{sample}/simulation/{sample}_neg.bw",
    params:
        read=lambda w, input: mappedReads(input['bed']),
    log:
        "logs/rule/analysis/{sample}/simulation/{sample}_genomeCovN.log",
    benchmark:
        "logs/rule/analysis/{sample}/simulation/{sample}_genomeCovN.benchmark.txt"

