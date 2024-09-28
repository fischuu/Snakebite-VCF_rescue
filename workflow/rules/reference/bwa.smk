rule _reference__bwa_create_index:
    """
    Index the Reference Genome (BWA).
    """
    input:
        config["reference_genome"]
    output:
        config["genome-bwa-index"]
    log:
        BWAL / "bwa_create_index.log"
    benchmark:
        BWAB / "bwa_create_index.tsv"
    threads: resources["cpu_per_task"]["single_thread"]
    resources:
        cpu_per_task = resources["cpu_per_task"]["single_thread"],
        mem_per_cpu = resources["mem_per_cpu"]["highmem"],
        time = resources["time"]["longrun"]
    conda:
        "envs/alignment.yaml"
    container:
        docker["alignment"]
    shell:"""
            bwa index -a bwtsw {input} 2> {log}
            samtools faidx {input} 2> {log}
  	"""