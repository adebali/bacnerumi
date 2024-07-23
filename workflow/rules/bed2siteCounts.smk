import re
build_ = config["build"]
project_ = config["project"]

rule bed2siteCounts:
    input:
        sites=f"results/{project_}/unified.bed",
        bed="results/{sample}/{sample}_shifted.bed",
    output: "results/{sample}/{sample}_siteCounts.tsv"
    log:
        "logs/rule/analysis/{sample}/{sample}bed2siteCounts.log",
    # benchmark:
    #     "logs/rule/analysis/{sample}/{sample}bed2siteCounts.benchmark.txt",
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
        (echo "`date -R`: Mapping reads to sites {output}..." &&
        bedtools intersect -c \
        -a {input.sites} \
        -b {input.bed} \
        -s \
        -f 1 \
        -F 1 \
        | awk '{{print $0, "\\t{params.fields}"}}' \
        > {output} &&
        echo "`date -R`: Success! {output} is generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """
        )

rule mergeSiteCounts:
    input:
        lambda w: input4mergeSiteCounts(config),
    output:
        out=f"results/{project_}/siteCounts.tsv",
    log: f"logs/rule/analysis/{project_}/siteCounts.log",
    # benchmark: f"logs/rule/analysis/{project_}/siteCounts.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    shell:  
        """
        (echo "`date -R`: Merge site counts from samples..." &&
        cat {input} > {output.out} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.out}; exit 1; }} ) > {log} 2>&1
        """