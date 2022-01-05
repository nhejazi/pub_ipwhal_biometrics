srun --pty -A co_biostat -p savio3 -t 3:00:00 bash -i
export TMPDIR='/global/scratch/users/nhejazi/rtmp'
module load gcc/6.3.0 r/4.0.3 r-packages/default
cd ~/ipwhal_meta/application
R --vanilla
