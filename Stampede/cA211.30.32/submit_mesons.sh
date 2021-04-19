#! /bin/bash

LOG_DIR=${SCRATCH}/cA211.30.32/log

SRC_NUM=4

let SRC_SETS=16/${SRC_NUM}
#SRC_SET=5

while read -r CONF; do

    mkdir -p ${LOG_DIR}/${CONF}

    REPL_TRAJ=$(echo ${CONF} | tr "-" " ")

    SRC=0

    while [ ${SRC} -lt ${SRC_SETS} ]; do
    #while [ ${SRC} -le ${SRC_SET} ]; do

	let SRC_LIST=${SRC}*${SRC_NUM}

	let ISRC=${SRC_LIST}+1

	let SRC_f=${SRC_LIST}+${SRC_NUM}

	while [ ${ISRC} -lt ${SRC_f} ]; do

	    SRC_LIST="${SRC_LIST},${ISRC}"

	    let ISRC=ISRC+1

	done

	INI_FILE=./ini/${CONF}.${SRC}.ini.xml

	python mkinput.py ${REPL_TRAJ} ${SRC_LIST} > ${INI_FILE}

	JOB_FILE=./jobFiles/run_${CONF}.${SRC}.job

	cp _run_job_template_ ${JOB_FILE}

	sed -i s/"_CONF_"/"${CONF}"/g ${JOB_FILE}

	sed -i s/"_SRC_"/"${SRC}"/g ${JOB_FILE}

	WAIT=1

	while [ ${WAIT} -eq 1 ]; do

	    COUNT="$(squeue -u tg838024 | grep -c tg838024)"
	    
	    if [ ${COUNT} -lt 24 ]; then

		sbatch ${JOB_FILE}
		
		echo "${CONF}.${SRC}"
		
		sleep 5

		WAIT=0

	    else

		sleep 180

	    fi

	done

	let SRC=SRC+1

    done

done < ./confs.txt
