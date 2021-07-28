#!/bin/bash
# Job name:
#SBATCH --job-name=ipwhal_nhefs_score
#
# Working directory:
#SBATCH --workdir=/global/home/users/nhejazi/
#
# Partition:
#SBATCH --partition=savio2

# Account:
#SBATCH --account=co_biostat
#
# Processors (1 node = 20 cores):
#SBATCH --nodes=1
#SBATCH --exclusive
#
# Wall clock limit ('0' for unlimited):
#SBATCH --time=48:00:00
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=nhejazi@berkeley.edu
#
# Job output:
#SBATCH --output=slurm.out
#SBATCH --error=slurm.out
#

## Command(s) to run:
export TMPDIR='~/rtmp'
export R_LIBS_USER='/global/scratch/nhejazi/R'
module load gcc/6.3.0 r/3.5.1 r-packages/default
cd ~/ipwhal-meta/application/

#R CMD BATCH --no-save --no-restore \
  #R/00_install_pkgs.R logs/00_install_pkgs.Rout

R CMD BATCH --no-save --no-restore \
  '--args selector_type=score' \
  R/01_ipw_analysis.R logs/01_ipw_analysis.Rout
