#!/bin/env python

#### configuration file ####
configfile: "config/ecoli_rna.yaml"

def getSampleFields(sample):
    l = []
    for key in sorted(sample.keys()):
        l.append(str(sample[key]))
    return '\t'.join(l)

def input4merge(sampleList):
    inputList = []
    for sample in sampleList:
        inputList.append(f"results/{sample}/{sample}_sense.tsv")
        inputList.append(f"results/{sample}/{sample}_antisense.tsv")
    return inputList

def getSampleFields(sample, additional_value):
    l = []
    for key in sorted(sample.keys()):
        l.append(str(sample[key]))
    l.append(additional_value)
    return '\t'.join(l)

include: "workflow/rules/ecoli_rna_sra_rules.smk"

# output = ["results/rna/ecoli_sra/readCounts.tsv"]
output = []
for sample in config['sample']:
    print(sample)
    # output.append(f"resources/samples/{sample}.fastq.gz")
    # output.append(f"results/{sample}/{sample}_trimmed.fastq.gz")
    # output.append(f"results/{sample}/{sample}_star.bam")
    # output.append(f"results/{sample}/{sample}_star_filtered.bam")
    # output.append(f"results/{sample}/{sample}_filtered.bed")
    # output.append(f"results/{sample}/{sample}_sense.tsv")
    # output.append(f"results/{sample}/{sample}_antisense.tsv")
    output.append(f"results/{project_}/transcriptCounts.tsv")
    

rule all:
    input: output
        
 