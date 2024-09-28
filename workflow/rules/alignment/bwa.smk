rule _alignment__bwa_alignreads_run:
    """
    Align reads to the Reference Genome (BWA).
    """
    input:
        reference=config["reference_genome"] ,
        refFiles=config["genome-bwa-index"],
        R1=FASTP / "{sample_id}_R1.fastq.gz",
        R2=FASTP / "{sample_id}_R2.fastq.gz",
    output:
        ALN / "BWA/{sample_id}.cram"
    log:
        ALN / "BWA/LOGS/{sample_id}.log"
    benchmark:
        ALN / "BWA/BENCHMARK/{sample_id}.tsv"
    threads: resources["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=resources["cpu_per_task"]["multi_thread"],
        mem_per_cpu=resources["mem_per_cpu"]["highmem"] // resources["cpu_per_task"]["multi_thread"],
        time =  resources["time"]["longrun"]
    conda:
        "envs/alignment.yaml"
    container:
        docker["alignment"]
    shell:"""
        df -h &> {log}
        ( bwa mem -t {threads} -M {input.reference} {input.R1} {input.R2} \
          | samtools sort \
              -l 9 \
              -m 1G \
              -o {output} \
              --reference {input.reference} \
              --threads {threads} \
          ) 2>> {log} 1>&2
  	"""