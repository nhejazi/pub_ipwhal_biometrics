#!/bin/bash
# Job name:
#SBATCH --job-name=ipwuhal_plateau_test
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
#SBATCH --time=24:00:00
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
export TMPDIR='~/rtmp'  # resolve update issues for compiled packages as per https://github.com/r-lib/devtools/issues/32
export R_LIBS_USER='/global/scratch/nhejazi/R'  # personal package library
module load gcc/6.3.0 r/3.5.1 r-packages/default
cd ~/ipwhal-meta/simulations/
R CMD BATCH --no-save --no-restore \
  '--args dgp=1a Mncv_max=120 ipw=ht Q_reg=Qn_hal selector=plateau' \
  R/mc_ipw_plateau.R logs/plateau_mc_var_test.Rout
