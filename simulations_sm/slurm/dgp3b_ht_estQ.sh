#!/bin/bash
# Job name:
#SBATCH --job-name=ipwuhal_dgp3b_ht_nsIPW_Qn
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
#SBATCH --time=96:00:00
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

## HT based on D_CAR with Qn via HAL
R CMD BATCH --no-save --no-restore \
  '--args dgp=3b Mncv_max=400 trunc_max=0.25 ipw=ht Q_reg=Qn_hal g_trunc=profile' \
  R/02_est_npipw_hal.R logs/dgp3b_ht_Qn_hal_profile.Rout
R CMD BATCH --no-save --no-restore \
  '--args dgp=3b Mncv_max=400 trunc_max=0.25 ipw=ht Q_reg=Qn_hal g_trunc=joint' \
  R/02_est_npipw_hal.R logs/dgp3b_ht_Qn_hal_joint.Rout

## HT based on score with Qn via HAL
R CMD BATCH --no-save --no-restore \
  '--args dgp=3b Mncv_max=400 trunc_max=0.3 ipw=ht Q_reg=Qn_hal selector=score' \
  R/02_est_npipw_hal.R logs/dgp3b_ht_Qn_hal_score.Rout

## HT based on Lepski with Qn via HAL
R CMD BATCH --no-save --no-restore \
  '--args dgp=3b Mncv_max=200 trunc_max=0.3 ipw=ht Q_reg=Qn_hal selector=plateau' \
  R/02_est_npipw_hal.R logs/dgp3b_ht_Qn_hal_plateau.Rout
