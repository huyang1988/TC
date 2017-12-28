#!/bin/sh
# one hour timelimit:
#SBATCH --time 00:30:00

# default queue, 32 processors (two nodes worth)
# -p: which queue

#SBATCH -p allgpu-noecc -N 2 -n 2 -c 6

###########################################
# -n means number of processes
# -N means number of machines 
# salloc -t 100 -p allgpu-noecc: ask for interactive node @ debug queue
# sinfo: shows the queues status
# scancel jid: cancel the job you do not want to continue
# squeue | grep username: shows my allocated node
#################################

# sleep 3600
# You can also put them in ~/.bashrc
# module load openmpi
# module load cuda/toolkit/7.5

srun ./tc /home/huyang/data/twitter
