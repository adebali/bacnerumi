
rule genome_download:
    output:
        "resources/ref_genomes/{build}/genome.fa"
    params:
        species=config["species"],
        datatype=config["datatype"],
        build=config["build"],
        release=config["release"],
    log:
        "logs/rule/analysis/{build}/log/download.log",
    # benchmark:
        # "logs/rule/analysis/{build}/log/download.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    wrapper:
        "0.69.0/bio/reference/ensembl-sequence"