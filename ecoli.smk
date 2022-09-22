#!/bin/env python

#### configuration file ####
configfile: "config/config_XR_initial.yaml" 
configfile: "config/ecoli.yaml"
include: "workflow/rules/common.smk"

rule all:
    input:
        lambda w: allInput(config["build"], config["sample"], 
            config["srr"]["enabled"], config["srr"]["codes"]),
 