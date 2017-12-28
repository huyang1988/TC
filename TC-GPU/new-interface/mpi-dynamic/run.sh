#!/bin/sh
# one hour timelimit:
#SBATCH --time 00:30:00

# default queue, 32 processors (two nodes worth)
# -p: which queue

#SBATCH -p allgpu-noecc -N 9 -n 9 

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

# mpirun -n 2 -host localhost ./tc ~/data/toy

mpirun ./tc /home/huyang/data/rmat-28
