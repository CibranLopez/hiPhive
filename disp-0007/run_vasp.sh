#!/usr/bin/env bash
#SBATCH --time=72:00:00
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --job-name=crn-0007
#SBATCH --error=crn.%J.err
#SBATCH --output=crn.%J.out
##SBATCH --exclusive

#SBATCH --constraint=cal

# Amount of RAM needed for this job:
#SBATCH --mem=300gb


# Send e-notification when fails or ends
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=rrurali@icmab.es

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load vasp/5.4.4_2023

#  time mpirun --use-hwthread-cpus -np 64 vasp_std
time mpirun -np $SLURM_NTASKS vasp_std
#time vasp_std
