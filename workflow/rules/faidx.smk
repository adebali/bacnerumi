build_ = config["build"]

rule samtools_index:
    input:
        f"resources/ref_genomes/{build_}/genome.fa",
    output:
        f"resources/ref_genomes/{build_}/genome.fa.fai",
    log:
        f"logs/rule/{build_}/faidx.log",
    benchmark:
        f"logs/rule/{build_}/faidx.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    params:
        extra="",  # optional params string
    wrapper:
        "v1.3.2/bio/samtools/faidx"