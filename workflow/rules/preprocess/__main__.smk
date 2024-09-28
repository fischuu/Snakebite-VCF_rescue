include: "__functions__.smk"
include: "concatenating.smk"
include: "fastp.smk"

rule preprocess:
    """Run the preprocess module"""
    input:
        expand(FASTP / "{sample_id}_R1.fastq.gz", sample_id=samples),