#! /bin/bash

LOG_DIR=/scratch/tuf47161/cA211.30.32/log

SRC_NUM=16

let SRC_SETS=16/${SRC_NUM}

while read -r CONF; do

    mkdir -p ${LOG_DIR}/${CONF}

    REPL_TRAJ=$(echo ${CONF} | tr "-" " ")

    SRC=0

    while [ ${SRC} -lt ${SRC_SETS} ]; do

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

	#while [ ${WAIT} -eq 1 ]; do

	COUNT_R="$(qstat -u tuf47161 | grep mesons | grep -c R)"
	    
	COUNT_Q="$(qstat -u tuf47161 | grep mesons | grep -c Q)"

	COUNT=$((COUNT_R + COUNT_Q))

	#if [ ${COUNT} -lt 10 ]; then

	qsub ${JOB_FILE}

	echo "${CONF}.${SRC}"

	sleep 5

	WAIT=0

	#else

	#sleep 180

	#fi

	#done

	let SRC=SRC+1

    done

done < ./confs.txt
