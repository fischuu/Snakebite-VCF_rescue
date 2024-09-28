rule _preprocess__concatenate_lanes_run:
    """
    Concatenate the demultiplexed fastq files (BASH).
    """
    input:
        R1=get_fastq_for_concatenating_read1,
        R2=get_fastq_for_concatenating_read2
    output:
        R1=CONC / "{sample_id}_R1.concatenated.fastq.gz",
        R1Report=CONC / "{sample_id}_R1.concatenated.fastq.gz.report",
        R1md5=CONC / "{sample_id}_R1.concatenated.fastq.gz.md5",
        R2=CONC / "{sample_id}_R2.concatenated.fastq.gz",
        R2Report=CONC / "{sample_id}_R2.concatenated.fastq.gz.report",
        R2md5=CONC / "{sample_id}_R2.concatenated.fastq.gz.md5",
    log:
        CONCL / "{sample_id}.log",
    benchmark:
        CONCB / "{sample_id}.tsv",
    conda:
        "envs/preprocess.yaml"
    container:
        docker["preprocess"]
    threads: resources["cpu_per_task"]["single_thread"]
    resources:
        cpu_per_task = resources["cpu_per_task"]["single_thread"],
        mem_per_cpu = resources["mem_per_cpu"]["lowmem"],
        time = resources["time"]["shortrun"]
    params:
        outfolder=CONC,
    shell:"""
        mkdir -p {params.outfolder} &> {log}
        cat {input.R1} > {output.R1} 2>> {log}
        ls {input.R1} > {output.R1Report} 2>> {log}
        md5sum {output.R1} > {output.R1md5}
        cat {input.R2} > {output.R2} 2>> {log}
        ls {input.R2} > {output.R2Report} 2>> {log}
        md5sum {output.R2} > {output.R2md5}
  	"""
