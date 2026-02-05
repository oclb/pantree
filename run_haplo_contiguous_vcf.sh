#!/bin/bash
#SBATCH -J pantree_haplo_cont        # Job name
#SBATCH -p medium                     # Partition (medium allows up to 5 days, required for 256GB)
#SBATCH -t 48:00:00                   # Time limit (48 hours)
#SBATCH --mem=128G                    # Memory per job (increased for Year 2 graphs)
#SBATCH -c 1                          # Cores per job
#SBATCH -o /home/pos149/pantree/logs/haplo_cont_%A_%a.out   # Standard output
#SBATCH -e /home/pos149/pantree/logs/haplo_cont_%A_%a.err   # Standard error
#SBATCH --array=1-48                  # Array for 24 chromosomes x 2 versions = 48 jobs

# Load required modules (adjust if needed)
module purge
module load gcc/14.2.0
module load python/3.13.1

# Define chromosome list
CHROMOSOMES=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

# Calculate which chromosome and version to process
# Jobs 1-24: Year 1, Jobs 25-48: Year 2
if [ ${SLURM_ARRAY_TASK_ID} -le 24 ]; then
    VERSION="v1"
    CHR_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
    GFA_DIR="/n/data1/hms/dbmi/oconnor/lab/shz311/pangenome/Data/chromosome_gfa_v1"
    VCF_DIR="/n/data1/hms/dbmi/oconnor/lab/pangenome/VCF_V1_HaploCont_020526"
else
    VERSION="v2"
    CHR_INDEX=$((SLURM_ARRAY_TASK_ID - 25))
    GFA_DIR="/n/data1/hms/dbmi/oconnor/lab/shz311/pangenome/Data/chromosome_gfa_v2"
    VCF_DIR="/n/data1/hms/dbmi/oconnor/lab/pangenome/VCF_V2_HaploCont_020526"
fi

CHR=${CHROMOSOMES[$CHR_INDEX]}

# Create output directories if they don't exist
mkdir -p ${VCF_DIR}
mkdir -p /home/pos149/pantree/logs

# Define input and output paths
GFA_FILE="${GFA_DIR}/${CHR}.gfa.gz"
VCF_FILE="${VCF_DIR}/${CHR}.vcf"
LOG_FILE="${VCF_DIR}/${CHR}.log"

# Print job information
echo "=================================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
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

# Activate virtual environment (adjust path to your pantree installation)
source /home/pos149/pantree/.venv/bin/activate

# Run pantree with haplo_contiguous_dfs_tree method
pantree gfa2vcf ${GFA_FILE} ${VCF_FILE} \
    --chr-id ${CHR} \
    --ref-name GRCh38 \
    --dfs-method contiguous \
    --priority-samples ${PRIORITY_SAMPLES} \
    --log-path ${LOG_FILE} \
    --verbose

# Record end time
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

# Get peak memory usage (in KB, convert to GB)
PEAK_MEM_KB=$(sstat -j ${SLURM_JOB_ID}.${SLURM_ARRAY_TASK_ID} --format=MaxRSS -n 2>/dev/null | tail -n 1 | sed 's/K//')
if [ -z "$PEAK_MEM_KB" ]; then
    # If sstat fails, try sacct (for completed jobs)
    PEAK_MEM_KB=$(sacct -j ${SLURM_JOB_ID}.${SLURM_ARRAY_TASK_ID} --format=MaxRSS -n | tail -n 1 | sed 's/K//')
fi
PEAK_MEM_GB=$(echo "scale=2; ${PEAK_MEM_KB:-0} / 1024 / 1024" | bc)

# Write runtime and memory usage to summary file
SUMMARY_FILE="/home/pos149/pantree/runtime_memory_haplo_cont.txt"
echo "${VERSION} ${CHR} ${RUNTIME} ${PEAK_MEM_GB}" >> ${SUMMARY_FILE}

echo "=================================================="
echo "End Time: $(date)"
echo "Runtime: ${RUNTIME} seconds ($(echo "scale=2; ${RUNTIME}/3600" | bc) hours)"
echo "Peak Memory: ${PEAK_MEM_GB} GB"
echo "=================================================="
