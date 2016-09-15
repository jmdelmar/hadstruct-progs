#! /usr/bin/ksh
#@ shell = /usr/bin/ksh
#@ job_type = parallel
#@ initialdir=/home/hpc/pr74yo/di56sof3/hadstruct-runs/cA2.30.24
#@ job_name = LLRUN
#@ class = general 
#@ energy_policy_tag = mg_24c48
#@ minimize_time_to_solution = yes
#@ node_usage = not_shared
#@ wall_clock_limit = 04:00:00
#@ network.MPI = sn_all,,us,,
#@ notification = never
#@ output = LLRUN.out.$(jobid)
#@ error =  LLRUN.err.$(jobid)
#@ island_count=1,1
#@ total_tasks = 5184 
#@ queue
        
. /etc/profile
. /etc/profile.d/modules.sh
        
cd /home/hpc/pr74yo/di56sof3/hadstruct-runs/cA2.30.24

module load hdf5/mpi
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$HOME/loc-install/lib         
export OMP_NUM_THREADS=1
mpiexec ../src/main ini/0420.ini.xml
