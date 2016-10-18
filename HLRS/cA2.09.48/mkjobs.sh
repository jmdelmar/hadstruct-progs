trajs=($*)
for ((i=0; i<${#trajs[*]}; i++)) ; do
    traj=${trajs[$i]}
    mkdir -p $traj
    cat _main.job_ | sed "s@_TRAJS_@$traj@g" > $traj/main.job
    for tr in $(echo $traj | sed 's@+@ @g') ; do
	cp ini/${tr}-3.ini.xml $traj/.
    done
    echo $traj/main.job
done
