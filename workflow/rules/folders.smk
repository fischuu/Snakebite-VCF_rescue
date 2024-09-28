# Definition of used paths and folders

SCRIPT_FOLDER = os.path.join(config["pipeline_folder"], "workflow", "scripts")

WD = os.getcwd()

# Preprocess folders
PRE = Path("RESULTS/PREPROCESS/")
CONC = PRE / "CONCATENATED/"
CONCL = PRE / "CONCATENATED/LOGS"
CONCB = PRE / "CONCATENATED/BENCHMARK"
FASTP = PRE / "FASTP/"
FASTPL = PRE / "FASTP/LOGS"
FASTPB = PRE / "FASTP/BENCHMARK"

# Reference folders
REF = Path ("RESULTS/REFERENCE")
BWA = REF / "BWA"
BWAL = REF / "BWA/LOGS"
BWAB = REF / "BWA/BENCHMARK"

# Alignment
ALN = Path ("RESULTS/ALIGNMENT")

#variantcalling
VC = Path ("RESULTS/VARIANT_CALLING")