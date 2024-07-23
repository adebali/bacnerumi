import re
build_ = config["build"]
project_ = config["project"]

rule bed2shiftedBed:
    input: 
        bed= "results/{sample}/{sample}_plus.bed",
    output: 
        bed="results/{sample}/{sample}_plus_shifted.bed",
    log:
        "logs/rule/analysis/{sample}/{sample}_bed2shiftedBed.log",
    # benchmark:
    #     "logs/rule/analysis/{sample}/{sample}_bed2shiftedBed.benchmark.txt",
    resources:
        memory=config["resources"]["memory"],
        cpu=config["resources"]["cpu"],
    params:
        title= lambda w: config["meta"][w.sample]["title"]
    shell:
        re.sub(' +', ' ', 
        """
        (echo "`date -R`: Shifting Bed coordinates..." &&
        cat {input.bed} | awk '$3-$2 > 10' | awk '{{print $1"\\t"$2+7"\\t"$2+9"\\t{params.title}\\t.\\t"$6}}' > {output.bed} &&
        echo "`date -R`: Success! Shifted bed - {output.bed} - is generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """
        )

rule bed2shiftedBed_minus:
    input: 
        bed= "results/{sample}/{sample}_minus.bed",
    output: 
        bed="results/{sample}/{sample}_minus_shifted.bed",
    log:
        "logs/rule/analysis/{sample}/{sample}_bed2shiftedBedM.log",
    # benchmark:
    #     "logs/rule/analysis/{sample}/{sample}_bed2shiftedBedM.benchmark.txt",
    resources:
        memory=config["resources"]["memory"],
        cpu=config["resources"]["cpu"],
    params:
        title= lambda w: config["meta"][w.sample]["title"]
    shell:
        re.sub(' +', ' ', 
        """
        (echo "`date -R`: Shifting Bed coordinates..." &&
        cat {input.bed} | awk '$3-$2 > 10' | awk '{{print $1"\\t"$3-9"\\t"$3-7"\\t{params.title}\\t.\\t"$6}}' > {output.bed} &&
        echo "`date -R`: Success! Shifted bed - {output.bed} - is generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """
        )

rule uniqueReads:
    input:
        plus="results/{sample}/{sample}_plus_shifted.bed",
        minus="results/{sample}/{sample}_minus_shifted.bed"
    output:
        plus="results/{sample}/{sample}_plus_shifted_uniq.bed",
        minus="results/{sample}/{sample}_minus_shifted_uniq.bed",
        both="results/{sample}/{sample}_shifted.bed"
    log:
        "logs/rule/analysis/{sample}/{sample}_uniqueReads.log",
    # benchmark:
    #     "logs/rule/analysis/{sample}/{sample}_uniqueReads.benchmark.txt",
    resources:
        memory=config["resources"]["memory"],
        cpu=config["resources"]["cpu"],
    shell:
        re.sub(' +', ' ', 
        """
        (echo "`date -R`: Getting unique reads..." &&
        cat {input.plus} {input.minus} | sort -k1,1 -k2,2n -k3,3n -k6,6 > {output.both} &&
        cat {input.plus} | awk '{{print  $1"\\t"$2"\\t"$3"\\t.\\t.\\t"$6 }}' | uniq > {output.plus} &&
        cat {input.minus} | awk '{{print  $1"\\t"$2"\\t"$3"\\t.\\t.\\t"$6 }}' | uniq  > {output.minus} &&
        echo "`date -R`: Success! Unique reads - {output} - are generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """
        )

rule mergeAndUnifyUniqueReads:
    input: 
        lambda w: input4mergeUnifyReads(config)
    output:
        f"results/{project_}/unified.bed",
    log:
        "logs/rule/analysis/mergeAndUnifyUniqueReads.log",
    # benchmark:
    #     "logs/rule/analysis/mergeAndUnifyUniqueReads.benchmark.txt",
    resources:
        memory=config["resources"]["memory"],
        cpu=config["resources"]["cpu"],
    shell:
        re.sub(' +', ' ', 
        """
        (echo "`date -R`: Getting unique reads..." &&
        cat {input} | sort -k2,2n -k6 -u > {output} &&
        echo "`date -R`: Success! mergeAndUnifyUniqueReads is completed with the output file: {output}" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """
        )