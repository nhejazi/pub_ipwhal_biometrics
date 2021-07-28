#!/bin/bash
# Job name:
#SBATCH --job-name=ipwuhal_dgp1b_ht_nsIPW_Q0
#
# Working directory:
#SBATCH --workdir=/global/home/users/nhejazi/
#
# Account:
#SBATCH --account=co_biostat
#
# Partition:
#SBATCH --partition=savio2
#
# Quality of Service:
#SBATCH --qos=biostat_savio2_normal
#
# Processors (1 node = 20 cores):
#SBATCH --nodes=1
#SBATCH --exclusive
#
# Wall clock limit ('0' for unlimited):
#SBATCH --time=12:00:00
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

## HT with Q0
R CMD BATCH --no-save --no-restore \
  '--args dgp=1b Mncv_max=10 trunc_max=0.2 ipw=ht Q_reg=Q0 g_trunc=profile' \
  R/02_est_npipw_hal.R logs/dgp1b_ht_Q0_profile.Rout
R CMD BATCH --no-save --no-restore \
  '--args dgp=1b Mncv_max=10 trunc_max=0.2 ipw=ht Q_reg=Q0 g_trunc=joint' \
  R/02_est_npipw_hal.R logs/dgp1b_ht_Q0_joint.Rout
