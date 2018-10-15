#! /bin/bash

############################################################
# Run a selection of tests for a Travis build that will 
# fit within the Travis 50 minute build window.
############################################################

usage() {
	echo 
	echo "$0 [-e <Wildcarded Evosuite test exclusions>]"
	echo
}

EVOSUITE_EXCLUSION_ARG=

while getopts "e:h" arg; do
  case $arg in
    e)
      EVOSUITE_EXCLUSION_ARG=${OPTARG}
      ;;
    h)
    	  usage
    	  exit 1
      ;;
  esac
done

echo
echo "EMBOSS_ROOT set to: "$EMBOSS_ROOT
echo "Test exclusion argument: "$EVOSUITE_EXCLUSION_ARG
echo

cd ..

if [[ ! -f ".scaffolding_list.tmp" ]]
then
	echo "ERROR: Cannot find Evosuite .scaffolding_list.tmp file - unit tests will not be sand-boxed! Exiting."
	exit 1
fi

mvn -Devosuite.exclude.filter="${EVOSUITE_EXCLUSION_ARG}" test jacoco:report

