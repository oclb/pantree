#!/bin/bash
#SBATCH -J pantree_all               # Job name
#SBATCH -p medium                     # Partition (medium allows up to 5 days)
#SBATCH -t 72:00:00                   # Time limit (72 hours)
#SBATCH --mem=128G                    # Memory per job
#SBATCH -c 1                          # Cores per job
#SBATCH -o /home/pos149/pantree/logs/pantree_%A_%a.out   # Standard output
#SBATCH -e /home/pos149/pantree/logs/pantree_%A_%a.err   # Standard error
#SBATCH --array=1-96                  # Array for all jobs

# Job layout (96 total jobs):
# Jobs 1-24:   Haplo-contiguous V2.1
# Jobs 25-48:  Weighted DFS V1
# Jobs 49-72:  Weighted DFS V2
# Jobs 73-96:  Weighted DFS V2.1

# Load required modules
module purge
module load gcc/14.2.0
module load python/3.13.1

# Define chromosome list
CHROMOSOMES=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

# Determine job type and parameters based on array task ID
if [ ${SLURM_ARRAY_TASK_ID} -le 24 ]; then
    # Haplo-contiguous V2.1
    METHOD="contiguous"
    VERSION="v2.1"
    CHR_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
    GFA_DIR="/n/data1/hms/dbmi/oconnor/lab/shz311/pangenome/Data/chromosome_gfa_v21_GRCh38"
    VCF_DIR="/n/data1/hms/dbmi/oconnor/lab/pangenome/VCF_V2.1_HaploCont_031026"
    SUMMARY_FILE="/home/pos149/pantree/runtime_memory_haplo_v21.txt"
elif [ ${SLURM_ARRAY_TASK_ID} -le 48 ]; then
    # Weighted DFS V1
    METHOD="max_weight"
    VERSION="v1"
    CHR_INDEX=$((SLURM_ARRAY_TASK_ID - 25))
    GFA_DIR="/n/data1/hms/dbmi/oconnor/lab/shz311/pangenome/Data/chromosome_gfa_v1"
    VCF_DIR="/n/data1/hms/dbmi/oconnor/lab/pangenome/VCF_V1_WeightedDFS_031026"
    SUMMARY_FILE="/home/pos149/pantree/runtime_memory_weighted_v1.txt"
elif [ ${SLURM_ARRAY_TASK_ID} -le 72 ]; then
    # Weighted DFS V2
    METHOD="max_weight"
    VERSION="v2"
    CHR_INDEX=$((SLURM_ARRAY_TASK_ID - 49))
    GFA_DIR="/n/data1/hms/dbmi/oconnor/lab/shz311/pangenome/Data/chromosome_gfa_v2"
    VCF_DIR="/n/data1/hms/dbmi/oconnor/lab/pangenome/VCF_V2_WeightedDFS_031026"
    SUMMARY_FILE="/home/pos149/pantree/runtime_memory_weighted_v2.txt"
else
    # Weighted DFS V2.1
    METHOD="max_weight"
    VERSION="v2.1"
    CHR_INDEX=$((SLURM_ARRAY_TASK_ID - 73))
    GFA_DIR="/n/data1/hms/dbmi/oconnor/lab/shz311/pangenome/Data/chromosome_gfa_v21_GRCh38"
    VCF_DIR="/n/data1/hms/dbmi/oconnor/lab/pangenome/VCF_V2.1_WeightedDFS_031026"
    SUMMARY_FILE="/home/pos149/pantree/runtime_memory_weighted_v21.txt"
fi

CHR=${CHROMOSOMES[$CHR_INDEX]}

# Create output directories if they don't exist
mkdir -p ${VCF_DIR}
mkdir -p /home/pos149/pantree/logs

# Define input and output paths
# V2.1 files are .gfa, V1 and V2 are .gfa.gz
if [ "${VERSION}" == "v2.1" ]; then
    GFA_FILE="${GFA_DIR}/${CHR}.gfa"
else
    GFA_FILE="${GFA_DIR}/${CHR}.gfa.gz"
fi
VCF_FILE="${VCF_DIR}/${CHR}.vcf"
LOG_FILE="${VCF_DIR}/${CHR}.log"

# Print job information
echo "=================================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Method: ${METHOD}"
echo "Version: ${VERSION}"
echo "Chromosome: ${CHR}"
echo "GFA File: ${GFA_FILE}"
echo "VCF File: ${VCF_FILE}"
echo "Start Time: $(date)"
echo "Node: $(hostname)"
echo "=================================================="

PRIORITY_SAMPLES=${PRIORITY_SAMPLES:-"GRCh38,CHM13,HG002"}

export POLARS_SKIP_CPU_CHECK=1

# Record start time
START_TIME=$(date +%s)

# Activate virtual environment
source /home/pos149/pantree/.venv/bin/activate

# Run pantree
pantree gfa2vcf ${GFA_FILE} ${VCF_FILE} \
    --chr-id ${CHR} \
    --ref-name GRCh38 \
    --dfs-method ${METHOD} \
    --priority-samples ${PRIORITY_SAMPLES} \
    --log-path ${LOG_FILE} \
    --verbose

# Record end time
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

# Get peak memory usage
PEAK_MEM_KB=$(sstat -j ${SLURM_JOB_ID}.${SLURM_ARRAY_TASK_ID} --format=MaxRSS -n 2>/dev/null | tail -n 1 | sed 's/K//')
if [ -z "$PEAK_MEM_KB" ]; then
    PEAK_MEM_KB=$(sacct -j ${SLURM_JOB_ID}.${SLURM_ARRAY_TASK_ID} --format=MaxRSS -n | tail -n 1 | sed 's/K//')
fi
PEAK_MEM_GB=$(echo "scale=2; ${PEAK_MEM_KB:-0} / 1024 / 1024" | bc)

# Write runtime and memory usage to summary file
echo "${VERSION} ${METHOD} ${CHR} ${RUNTIME} ${PEAK_MEM_GB}" >> ${SUMMARY_FILE}

echo "=================================================="
echo "End Time: $(date)"
echo "Runtime: ${RUNTIME} seconds ($(echo "scale=2; ${RUNTIME}/3600" | bc) hours)"
echo "Peak Memory: ${PEAK_MEM_GB} GB"
echo "=================================================="
