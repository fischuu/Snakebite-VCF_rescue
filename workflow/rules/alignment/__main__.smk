# Main script to orchestrate QC rules

include: "bwa.smk"

rule alignment:
    """Run the alignment module"""
    input:
        expand(ALN / "BWA/{sample_id}.cram", sample_id=samples),