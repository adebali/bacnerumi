rule prepareSingletons:
    input:
        "resources/ref_genomes/{build}/genome.fa",
    output:
        theoretical="resources/ref_genomes/{build}/theoreticalReads_{damageSite}.fa",
        singletons="resources/ref_genomes/{build}/singletons_{damageSite}.fa",
        bedsingletons=report("resources/ref_genomes/{build}/singletons_{damageSite}.bed", category="genome"),
    log:
        "logs/rule/analysis/{build}/log/prepareSingletons_{damageSite}.log",
    # benchmark:
        # "logs/rule/analysis/{build}/log/prepareSingletons_{damageSite}.benchmark.txt",
    resources:
        memory="16GB",
        cpu=1
    params:
        kmerlist=config['readLength']
    shell:  
        """
        (echo "`date -R`: Preparing theoretical reads..." &&
        python3 workflow/scripts/fa2theoreticalReads.py -i {input} -k {params.kmerlist} -d {wildcards.damageSite} > {output.theoretical} &&
        echo "`date -R`: Success! Theoretical reads file is generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Preparing singletons (reads mapped onto genome only for once)..." &&
        python3 workflow/scripts/fa2singletons.py -i {output.theoretical} -b {output.bedsingletons} > {output.singletons} &&
        echo "`date -R`: Success! Singletons file is generated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """