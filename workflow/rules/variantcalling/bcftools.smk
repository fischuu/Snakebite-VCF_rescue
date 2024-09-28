rule _variantcalling__bcftools_mpileup:
    """
    Create the MPILEUP for all provided positions (bcftools).
    """
    input:
        reference=config["reference_genome"] ,
        variants=config["high_quality_vcf"],
        bam=ALN / "BWA/{sample_id}.cram",
    output:
        VC / "BCFTOOLS/{sample_id}.mpileup"
    log:
        VC / "BCFTOOLS/LOGS/{sample_id}.mpileup.log"
    benchmark:
        ALN / "BCFTOOLS/BENCHMARK/{sample_id}.mpileup.tsv"
    threads: resources["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=resources["cpu_per_task"]["multi_thread"],
        mem_per_cpu=resources["mem_per_cpu"]["highmem"] // resources["cpu_per_task"]["multi_thread"],
        time =  resources["time"]["longrun"]
    conda:
        "envs/variants.yaml"
    container:
        docker["variants"]
    shell:"""
        bcftools mpileup -f {input.reference} {input.bam} -T {input.variants} -Ou -o {output} 2>> {log} 1>&2
  	"""
  	
  	
rule _variantcalling__bcftools_calling:
    """
    Create the MPILEUP for all provided positions (bcftools).
    """
    input:
        VC / "BCFTOOLS/{sample_id}.mpileup"
    output:
        VC / "BCFTOOLS/{sample_id}.vcf"
    log:
        VC / "BCFTOOLS/LOGS/{sample_id}.calling.log"
    benchmark:
        ALN / "BCFTOOLS/BENCHMARK/{sample_id}.calling.tsv"
    threads: resources["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=resources["cpu_per_task"]["multi_thread"],
        mem_per_cpu=resources["mem_per_cpu"]["highmem"] // resources["cpu_per_task"]["multi_thread"],
        time =  resources["time"]["longrun"]
    container:
        docker["variants"]
    shell:"""
        bcftools call -mv -A --ploidy 1 -o {input} {output} 2>> {log} 1>&2
  	"""
  