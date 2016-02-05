trajs=($*)
for ((i=0; i<${#trajs[*]}; i++)) ; do
    traj=${trajs[$i]}
    mkdir -p $traj
    cat _main.job_ | sed "s@_TRAJ_@$traj@g" > $traj/main.job
    cd $traj/ && cp -f ../invert.input invert.input && cd ../
    echo $traj/main.job
done
