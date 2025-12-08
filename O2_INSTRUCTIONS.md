# Instructions for Running Pantree on O2

## Setup Steps

### 1. Transfer the script to O2
```bash
scp run_haplo_contiguous_vcf.sh pos149@transfer.rc.hms.harvard.edu:/home/pos149/pantree/
```

### 2. SSH into O2
```bash
ssh pos149@o2.hms.harvard.edu
```

### 3. Create necessary directories on O2
```bash
mkdir -p /home/pos149/pantree/VCF_V1_HaploCont
mkdir -p /home/pos149/pantree/VCF_V2_HaploCont
mkdir -p /home/pos149/pantree/logs
```

### 4. Make the script executable
```bash
cd /home/pos149/pantree
chmod +x run_haplo_contiguous_vcf.sh
```

### 5. Submit the job array
```bash
sbatch run_haplo_contiguous_vcf.sh
```

## What the script does:

- **Runs 48 parallel jobs**: 24 for Year 1 GFAs + 24 for Year 2 GFAs
- **Array jobs 1-24**: Process Year 1 chromosomes (chr1-chr22, chrX, chrY)
- **Array jobs 25-48**: Process Year 2 chromosomes (chr1-chr22, chrX, chrY)
- **Uses `haplo_contiguous_dfs_tree` method** with GRCh38 as reference
- **Outputs**:
  - VCF files: `/home/pos149/pantree/VCF_V1_HaploCont/` and `/home/pos149/pantree/VCF_V2_HaploCont/`
  - Log files: `/home/pos149/pantree/logs/`
  - Runtime/memory summary: `/home/pos149/pantree/runtime_memory_haplo_cont.txt`

## Resources allocated per job:
- **Memory**: 64GB
- **Time**: 10 hours
- **Cores**: 1
- **Partition**: short

## Monitoring your jobs:

```bash
# Check job status
squeue -u pos149

# Check specific job
squeue -j <JOB_ID>

# View output logs
tail -f /home/pos149/pantree/logs/haplo_cont_<JOB_ID>_<ARRAY_ID>.out

# View error logs
tail -f /home/pos149/pantree/logs/haplo_cont_<JOB_ID>_<ARRAY_ID>.err

# Check runtime summary (after jobs complete)
cat /home/pos149/pantree/runtime_memory_haplo_cont.txt
```

## Cancel jobs if needed:

```bash
# Cancel all array jobs
scancel <JOB_ID>

# Cancel specific array task
scancel <JOB_ID>_<ARRAY_TASK_ID>
```

## Notes:

- The script assumes pantree is installed at `/home/pos149/pantree/.venv/`
- Adjust module loads if O2 has different versions available
- The runtime/memory summary file format: `VERSION CHR RUNTIME_SECONDS PEAK_MEMORY_GB`
