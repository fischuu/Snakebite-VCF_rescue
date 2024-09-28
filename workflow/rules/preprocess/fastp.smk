rule _preprocess__fastp__run:
    """Run fastp on one library"""
    input:
        forward_=CONC / "{sample_id}_R1.concatenated.fastq.gz",
        reverse_=CONC / "{sample_id}_R1.concatenated.fastq.gz",
    output:
        forward_=FASTP / "{sample_id}_R1.fastq.gz",
        reverse_=FASTP / "{sample_id}_R2.fastq.gz",
        unpaired1=FASTP / "{sample_id}_u1.fastq.gz",
        unpaired2=FASTP / "{sample_id}_u2.fastq.gz",
        html=FASTP / "{sample_id}.fastp.html",
        json=FASTP / "{sample_id}.fastp.json",
    log:
        FASTPL / "{sample_id}.log",
    benchmark:
        FASTPB / "{sample_id}.tsv"
    conda:
        "envs/preprocess.yaml"
    container:
        docker["preprocess"]
    params:
        adapter_forward=get_forward_adapter,
        adapter_reverse=get_reverse_adapter,
        extra=config["preprocess"]["fastp"]["extra"],
        length_required=config["preprocess"]["fastp"]["length_required"],
        temp_forward_=lambda wildcards: FASTP / f"{wildcards.sample_id}_tmp_1.fq",
        temp_reverse_=lambda wildcards: FASTP / f"{wildcards.sample_id}_tmp_2.fq",
        temp_unpaired1=lambda wildcards: FASTP / f"{wildcards.sample_id}_tmp_u1.fq",
        temp_unpaired2=lambda wildcards: FASTP / f"{wildcards.sample_id}_tmp_u2.fq",
    threads: resources["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=resources["cpu_per_task"]["multi_thread"],
        mem_per_cpu=resources["mem_per_cpu"]["highmem"] // resources["cpu_per_task"]["multi_thread"],
        time =  resources["time"]["longrun"]
    shell:"""
   # Initialize the adapter variables
        adapter_params=""

        # Check if adapters are provided and create the adapter parameters string
        if [ "{params.adapter_forward}" != "" ] && [ "{params.adapter_reverse}" != "" ]; then
            adapter_params="--adapter_sequence {params.adapter_forward} --adapter_sequence_r2 {params.adapter_reverse}"
        fi

        # Run fastp with intermediate files
        fastp \
            --in1 {input.forward_} \
            --in2 {input.reverse_} \
            --out1 {params.temp_forward_} \
            --out2 {params.temp_reverse_} \
            --unpaired1 {params.temp_unpaired1} \
            --unpaired2 {params.temp_unpaired2} \
            --html {output.html} \
            --json {output.json} \
            --verbose \
            $adapter_params \
            --length_required {params.length_required} \
            --thread {threads} \
            {params.extra} \
        2> {log} 1>&2

        # Compress the outputs using bgzip
        bgzip -l 9 -@ {threads} {params.temp_forward_}
        bgzip -l 9 -@ {threads} {params.temp_reverse_}
        bgzip -l 9 -@ {threads} {params.temp_unpaired1}
        bgzip -l 9 -@ {threads} {params.temp_unpaired2}

        # Move the compressed files to the final destination
        mv {params.temp_forward_}.gz {output.forward_}
        mv {params.temp_reverse_}.gz {output.reverse_}
        mv {params.temp_unpaired1}.gz {output.unpaired1}
        mv {params.temp_unpaired2}.gz {output.unpaired2}


        # Check the integrity of the gzipped files and log the output
        echo "Checking integrity of output files" >> {log}
        gzip -t {output.forward_} 2>> {log}
        gzip -t {output.reverse_} 2>> {log}
        gzip -t {output.unpaired1} 2>> {log}
        gzip -t {output.unpaired2} 2>> {log}
        echo "Integrity check completed" >> {log}
        """