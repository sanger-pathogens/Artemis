#! /bin/bash

###############################################################
# BASH Script to generate EvoSuite JUnit test source code.    #
# It is run from the build_test.xml Ant script using the      #
# create-evotests target. Re-generating all tests can         #
# take several hours.                                         #
# Generated JUnit tests are placed in the evosuite-tests      #
# folder and can be RUN via the ant script.                   #
#                                                             #
# Author: K. Pepper                                           #
#                                                             #
###############################################################

PROGNAME=$0
EVOSUITE_JAR=$1
EVOSUITE_CLASSPATH=$2
EVOSUITE_ARGS=$3
EVOSUITE_REPORT_DIR=$4

if [ "$#" -ne 4 ] 
then
	echo
    echo "Usage: $PROGNAME <EvoSuite jar path> <Classpath for EvoSuite to use> <EvoSuite command line arguments> <Coverage reports folder>"
    echo "NOTE: This script is only intended to be used from an ant build."
    echo
    exit 255
fi

RC=0
LOG_FILE="/tmp/evosuite-run.log"
CLASSPATH=$EVOSUITE_CLASSPATH

#
# Run POOL_SIZE number of EvoSuite invocations is parallel
#
POOL_SIZE=4

CLASS_LIST=`find ../uk/ac/sanger -name *.java -print`


function printLogs {

	FROMIDX=$1
	TOIDX=$2
	
	for ((x=$FROMIDX; x<=$TOIDX; x++))
	do
		if [ -e "${LOG_FILE}.${x}" ]
		then
			echo "$PROGNAME: Replaying log ${x}"
			cat ${LOG_FILE}.${x}
			rm ${LOG_FILE}.${x}
		fi
	done
}

echo "$PROGNAME: Building EvoSuite JUnit test cases with arguments: $EVOSUITE_ARGS"

i=0
BATCH_START_IDX=1
END_OF_BATCH_FLAG=0
for CLASS_FILE in $CLASS_LIST
do

	i=$(($i+1))
	
	if [ $END_OF_BATCH_FLAG -eq 1 ]
	then
		END_OF_BATCH_FLAG=0
		BATCH_START_IDX=$i
	fi
	
	#
	# Determine next class name (including package)
	#
	TMP1=`echo $CLASS_FILE | tr '/' '.'` 
	TMP2=${TMP1%%'.java'}
	CLASS=${TMP2#'...'}
	
	echo "${PROGNAME}: ProcessIndex: ${i} - Executing for class $CLASS"
	
	#
	# Execute the Java EvoSuite jar to create a unit test.
	#
	CMD="java -noverify -jar $EVOSUITE_JAR -class $CLASS $EVOSUITE_ARGS -class $CLASS"
	$CMD 2>&1 > ${LOG_FILE}.${i} &
	RC=$?
	if [ $RC -ne 0 ]
	then
		echo "${PROGNAME} [ERROR]: EvoSuite process returned error status: $RC for class $CLASS"
	fi
	
	CHECK=$(($i % $POOL_SIZE))
	if [ $CHECK -eq 0 ]
	then
		echo "${PROGNAME}: Waiting for $POOL_SIZE background processes to complete..."
		wait
		printLogs $BATCH_START_IDX $i
		END_OF_BATCH_FLAG=1
		echo "${PROGNAME}: Background processes completed"
	fi

done 

if [ $END_OF_BATCH_FLAG -eq 0 ] && [ $i -gt 0 ]
then
	wait
	printLogs $BATCH_START_IDX $i
fi

echo "${PROGNAME}: Finished"

exit 0
