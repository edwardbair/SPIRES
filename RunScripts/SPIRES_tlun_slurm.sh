#!/bin/sh
## Job name:
#SBATCH --job-name=spires_tlun_slurm
#
## Request one node:
#SBATCH --ntasks=1
#
## Number of processors for single task needed for use case (example):
#SBATCH --cpus-per-task=64
#
## Wall clock limit:
#SBATCH --time=UNLIMITED
## dir to run in
#SBATCH --chdir=/home/snowhydro/nbair/SPIRESScriptsAndInputs/
#
##output file names
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err

## Command(s) to run:

./run_batch_spires_tlun.sh
