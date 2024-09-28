# Main script to orchestrate QC rules

include: "bwa.smk"

rule reference:
    """Run the reference module"""
    input:
        rules._reference__bwa_create_index.output