#! /usr/bin/ksh
#@ shell = /usr/bin/ksh
#@ output = $(job_name).$(jobid).$(stepid).out
#@ error =  $(job_name).$(jobid).$(stepid).err
#@ initialdir = /home/hpc/pr74yo/di56sof3/hadstruct-progs/cA2.30.24/
#@ job_name = 0420
#@ node_usage = not_shared
##################################################################
#@ step_name = source_0_3
#@ wall_clock_limit = 00:30:00
#@ network.MPI = sn_all,not_shared,us
#@ energy_policy_tag = MG_24c48
#@ minimize_time_to_solution = yes
#@ job_type = parallel
#@ class = test
#@ notification = never
#@ node = 16
#@ island_count = 1,1
#@ total_tasks = 432
#@ queue
##################################################################
#@ step_name = source_4_7
#@ wall_clock_limit = 00:30:00
#@ network.MPI = sn_all,not_shared,us
#@ energy_policy_tag = MG_24c48
#@ minimize_time_to_solution = yes
#@ job_type = parallel
#@ class = test
#@ notification = never
#@ node = 16
#@ island_count = 1,1
#@ total_tasks = 432
#@ queue
##################################################################
#@ step_name = source_8_11
#@ wall_clock_limit = 00:30:00
#@ energy_policy_tag = MG_24c48
#@ minimize_time_to_solution = yes
#@ job_type = parallel
#@ network.MPI = sn_all,not_shared,us
#@ notification = never
#@ class = test
#@ node = 16
#@ island_count = 1,1
#@ total_tasks = 432
#@ queue
##################################################################
#@ step_name = source_12_15
#@ wall_clock_limit = 00:30:00
#@ network.MPI = sn_all,not_shared,us
#@ energy_policy_tag = MG_24c48
#@ minimize_time_to_solution = yes
#@ job_type = parallel
#@ notification = never
#@ class = test
#@ node = 16
#@ island_count = 1,1
#@ total_tasks = 432
#@ queue
        
. /etc/profile
. /etc/profile.d/modules.sh

n=0420
ENSEMBLE=cA2.30.24/
CORR_DIR=$WORK/$ENSEMBLE/Corr        
cd $LOADL_STEP_INITDIR

module load hdf5/mpi
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$HOME/loc-install/lib
export OMP_NUM_THREADS=1

mkdir -p $CORR_DIR/${n}

function run {
    s=$1
    poe ../src/main ./ini/${n}-${s}.ini.xml > ${n}-${s}.log
    retval=$?
    return $retval
}

case $LOADL_STEP_NAME in
    source_0_3 )
	run 0
	;;
    source_4_7 )
	run 1
	;;
    source_8_11 )
	run 2
	;;
    source_12_15 )
	run 3
	;;
esac
