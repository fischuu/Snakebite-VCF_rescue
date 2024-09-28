# Main script to orchestrate QC rules
include: "bcftools.smk"


rule variantcalling:
    """Run the variantcalling module"""
    input:
        expand(VC / "BCFTOOLS/{sample_id}.vcf", sample_id=samples),

