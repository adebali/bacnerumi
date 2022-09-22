import re
build_ = config["build"]

rule bed2geneCounts_tsnts:
    input:
        reads=rules.bam2bed_se.output.bed,
        genes=f"resources/ref_genomes/{build_}/genes.bed",
    output:
        TS=temp("results/{sample}/{sample}_TS.tsv"),
        NTS=temp("results/{sample}/{sample}_NTS.tsv"),
        mergedTSNTS="results/{sample}/{sample}_TSNTS.tsv",
    log:
        "logs/rule/analysis/{sample}/{sample}_bed2geneTSNTScounts.log",
    benchmark:
        "logs/rule/analysis/{sample}/{sample}_bed2geneTSNTScounts.benchmark.txt",
    resources:
        memory=config["resources"]["memory"],
        cpu=config["resources"]["cpu"],
    params:
        fields = lambda w: getSampleFields(config["meta"][w.sample])
    conda:
        "../envs/bedtools.yaml"
    shell:
        re.sub(' +', ' ', 
        """
        (echo "`date -R`: Counting TS reads {output.TS}..." &&
        bedtools intersect -c \
        -a {input.genes} \
        -b {input.reads} \
        -S \
        > {output.TS} &&
        echo "`date -R`: Success! {output.TS} is generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Counting NTS reads {output.NTS}..." &&
        bedtools intersect -c \
        -a {input.genes} \
        -b {input.reads} \
        -s \
        > {output.NTS} &&
        echo "`date -R`: Success! {output.NTS} is generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Merging TS NTS files..." &&
        paste \
        {output.TS} \
        <(cut -f 7 {output.NTS}) \
        | awk '{{print $0, "\\t{params.fields}"}}' \
        > {output.mergedTSNTS} &&
        echo "`date -R`: Success! {output.mergedTSNTS} is generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        """
        )

rule mergeTSNTScounts:
    input:
        lambda w: input4mergeTSNTScounts(config["sample"], config["srr"]["enabled"], 
            config["srr"]["codes"]),
    output:
        out="results/readCountsTSNTS.tsv",
    log:
        "logs/rule/analysis/readCountsTSNTS.log",
    benchmark:
        "logs/rule/analysis/readCountsTSNTS.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    shell:  
        """
        (echo "`date -R`: Merge TSNTS counts from true sample..." &&
        cat {input} > {output.out} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.out}; exit 1; }} ) > {log} 2>&1
        """

use rule bed2geneCounts_tsnts as bed2geneCounts_tsnts_sim with:
    input:
        reads="results/{sample}/simulation/{sample}.bed",
        genes=f"resources/ref_genomes/{build_}/genes.bed",
    output:
        TS=temp("results/{sample}/simulation/{sample}_TS.tsv"),
        NTS=temp("results/{sample}/simulation/{sample}_NTS.tsv"),
        mergedTSNTS="results/{sample}/simulation/{sample}_TSNTS.tsv",
    log:
        "logs/rule/analysis/{sample}/simulation/{sample}_bed2geneCounts_tsnts.log",
    benchmark:
        "logs/rule/analysis/{sample}/simulation/{sample}_bed2geneCounts_tsnts.benchmark.txt",


use rule mergeTSNTScounts as mergeTSNTScounts_sim with:
    input:
        lambda w: input4mergeTSNTScounts_sim(config["sample"], config["srr"]["enabled"], 
            config["srr"]["codes"]),
    output:
        out="results/simulation/readCountsTSNTS.tsv",
    log:
        "logs/rule/analysis/simulation/readCountsTSNTS.log",
    benchmark:
        "logs/rule/analysis/simulation/readCountsTSNTS.benchmark.txt",