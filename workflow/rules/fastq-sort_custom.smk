
rule fastq_sort:
    input:
        "resources/samples/{sample}.fastq",
    output:
        "results/{sample}/{sample}_sorted.fastq"
    log:
        "logs/rule/analysis/{sample}/{sample}_fastqSortCustom.log",
    resources:
        memory="64GB",
        cpu=1
    conda:
        "../envs/fastq-sort-bioawk.yaml"
    shell: "(bioawk -c fastx \'{{print \"@\"$name, $seq, \"+\", $qual}}\' {input} |  sort -k2,2 -T tmp | awk \'{{for(i=1; i<=NF; i++) print $i}}\' > {output}) > {log} 2>&1"

