#! /usr/bin/ksh
#@ shell = /usr/bin/ksh
#@ job_type = parallel
#@ initialdir=/home/hpc/pr74yo/di56sof3/hadstruct-runs/cA2.0900.64/rms
#@ job_name = LLRUN
#@ class = test
#@ energy_policy_tag = mg_64c128
#@ minimize_time_to_solution = yes
#@ node_usage = not_shared
#@ wall_clock_limit = 00:30:00
#@ network.MPI = sn_all,,us,,
#@ notification = never
#@ output = LLRUN.out.$(jobid)
#@ error =  LLRUN.err.$(jobid)
#@ node = 19
#@ island_count=1,1
#@ total_tasks = 512
#@ queue
        
. /etc/profile
. /etc/profile.d/modules.sh
        
cd /home/hpc/pr74yo/di56sof3/hadstruct-runs/cA2.0900.64/rms

module load hdf5/mpi
export OMP_NUM_THREADS=1
for alpha_gauss in 0.0625 0.0717936471873 0.0824692444233 0.0947322854069 0.108818820412 0.125 0.143587294375 0.164938488847 0.189464570814 0.217637640824 0.25 0.287174588749 0.329876977693 0.378929141628 0.435275281648 0.5 0.574349177499 0.659753955386 0.757858283255 0.870550563296 1.0 1.148698355 1.31950791077 1.51571656651 1.74110112659 2.0 2.29739670999 2.63901582155 3.03143313302 3.48220225318 4.0 4.59479341999 5.27803164309 6.06286626604 6.96440450637 8.0 ; do
	mpiexec ./main /gpfs/work/pr74yo/di49saj/rho_15/conf.0250 ${alpha_gauss} |tee rms-a${alpha_gauss}.0250.log
done 
