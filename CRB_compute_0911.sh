#!/bin/bash
#
#SBATCH --job-name=Matlab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=15GB
#SBATCH --time=12:00:00
#SBATCH --output=log/slurm_%j.out
#SBATCH --error=log/slurm_%j.err
 
module purge
module load matlab/2016b
 
cd /scratch/kl3141/MRF_matlab/MRF_CRB/oct_hsfp

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_ARRAY_JOB_ID="$SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID="$SLURM_ARRAY_TASK_ID
 
cat<<EOF | srun matlab -nodisplay
parpool('local', $SLURM_CPUS_PER_TASK)
MT_Simulate_Signals_Git_2
exit
EOF
 