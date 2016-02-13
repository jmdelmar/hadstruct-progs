#!/bin/bash -l
#
# Automates the post processing of the production runs in this
# folder.
# 
# . Takes one argument which is the trajectory number: $traj
# . Assumes python scripts are in $HOME/QuaHoG/utils/
# . Assumes HDF5 correlators are in $WORK/$ENSEMBLE/Corr/$traj
# . Will write the averages in $WORK/$ENSEMBLE/Corr/ave/$traj
# . WORK and ENSEMBLE should be set
#
###

run() {
    eval "$1" &
    proc_counter=$((proc_counter + 1))
    [ $proc_counter -eq $NP ] && wait && proc_counter=0
}

proc_counter=0
module load tools/python/2.7.8
WORK=$HOME/work
ENSEMBLE=cA2.09.48
NP=8

traj=$1
script_dir=$HOME/QuaHoG/utils/
corr_dir=$WORK/$ENSEMBLE/Corr/$traj/
ave_dir=$WORK/$ENSEMBLE/Corr/ave/$traj/

[ ! -d $ave_dir ] && mkdir -p $ave_dir

#
# If we're in a parallel job, prefix commands with aprun
#
pre=$(test -z ${PBS_NUM_NODES+x} || echo aprun -q -n 1 -N 1)

#
# User-specified number of parallel jobs is overriden by number of
# nodes available if we're in a parallel job
#
if [ ! -z ${PBS_NUM_NODES+x} ]; then
    NP=${PBS_NUM_NODES}
fi

# First, fourier-transform the HDF5 files. The subsequent scripts
# assume the HDF5 root group to be the trajectory number, so set it
# here
ft=${script_dir}/corr-ft.py

# three-point functions need inverse FT
#
fnames=($(ls $corr_dir | grep -E thrp_sx..sy..sz..st.._gN50a4p0_aN50a0p5_P[0-9]_dt[0-9]{2}.[ud][pn].h5$))
for fn in ${fnames[*]} ; do
    mf=$(echo $fn | sed 's@_sx@_mom_sx@')
    echo $mf
    run "$pre python $ft -i -r $traj $corr_dir/$fn -o $corr_dir/$mf"
done

# nucleons and mesons
#
fnames=($(ls $corr_dir | grep nucleons_sx..sy..sz..st.._gN50a4p0_aN50a0p5.h5$))
fnames=(${fnames[*]} $(ls $corr_dir | grep mesons_sx..sy..sz..st.._gN50a4p0_aN50a0p5.h5$))
for fn in ${fnames[*]} ; do
    mf=$(echo $fn | sed 's@_sx@_mom_sx@')
    echo $mf
    run "$pre python $ft -r $traj $corr_dir/$fn -o $corr_dir/$mf"
done
wait

# Now average files over source positions. List the files according to
# unique names when ignoring the source position
ave=${script_dir}/corr-ave.py
sets=($(ls $corr_dir | grep -E thrp_mom_sx..sy..sz..st.._gN50a4p0_aN50a0p5_P[0-9]_dt[0-9]{2}.[ud][pn].h5$ | sed 's@_sx.*st..@@'|sort -u))
sets=(${sets[*]} $(ls $corr_dir | grep -E nucleons_mom_sx..sy..sz..st.._gN50a4p0_aN50a0p5.h5$ | sed 's@_sx.*st..@@'|sort -u))
sets=(${sets[*]} $(ls $corr_dir | grep -E mesons_mom_sx..sy..sz..st.._gN50a4p0_aN50a0p5.h5$ | sed 's@_sx.*st..@@'|sort -u))
for s in ${sets[*]} ; do
    fns=$(ls ${corr_dir}/$(echo $s|sed 's@mom_@mom_sx??sy??sz??st??_@')|xargs)
    nsrc=$(echo $fns|wc -w)
    af=$(echo $s|sed "s@_gN50a4p0_aN50a0p5@_gN50a4p0_aN50a0p5_nsrc${nsrc}@")
    echo $af
    run "$pre python $ave -o ${ave_dir}/$af $fns"
done
wait

# Finally, project nucleon two-point function with (1+/-g0)/4.0
proj=${script_dir}/nucl-proj.py
for fn in ${ave_dir}/nucleons_mom_gN50a4p0_aN50a0p5_nsrc??.h5 ; do
    pn=$(echo $fn|sed 's@nucleons@nproj@')
    echo $pn
    run "$pre python $proj -o $pn $fn"
done
wait

