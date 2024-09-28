# concatenating
def get_fastq_for_concatenating_read1(wildcards):
    r1 = samplesheet.loc[samplesheet["sample_name"] == wildcards.sample_id]["read1"]
    path = config["rawdata_folder"] + "/"
    output = [path + x for x in r1]
    return output   

def get_fastq_for_concatenating_read2(wildcards):
    r1 = samplesheet.loc[samplesheet["sample_name"] == wildcards.sample_id]["read2"]
    path = config["rawdata_folder"] + "/"
    output = [path + x for x in r1]
    return output  
    
# fastp
def get_forward_adapter(wildcards):
    """Get forward adapter for a sample and library."""
    adapter = samplesheet.loc[samplesheet["sample_name"] == wildcards.sample_id]["forward_adapter"]
    # If the adapter is NaN or empty, return an empty string
    if adapter.empty or pd.isna(adapter.values[0]):
        return ""
    else:
        return adapter.values[0]


def get_reverse_adapter(wildcards):
    """Get reverse adapter for a sample and library."""
    adapter = samplesheet.loc[samplesheet["sample_name"] == wildcards.sample_id]["reverse_adapter"]
    # If the adapter is NaN or empty, return an empty string
    if adapter.empty or pd.isna(adapter.values[0]):
        return ""
    else:
        return adapter.values[0]
