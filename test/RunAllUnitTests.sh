#! /bin/bash

############################################################
# Run all hand-crafted and Evosuite generated unit tests.
# Evosuite unit tests are run in a sandbox by the 
# Evosuite unit test listener which is invoked by Surefire.
#
############################################################
echo
echo "This may take 30 minutes or more..."
echo
echo "EMBOSS_ROOT set to: "$EMBOSS_ROOT
echo

sleep 5
cd ..

#mvn clean compile test-compile

if [[ ! -e ./.scaffolding_list.tmp ]]
then
	echo "ERROR: Cannot find Evosuite .scaffolding_list.tmp file - unit tests will not be sand-boxed! Exiting."
	exit 1
fi 

mvn -Devosuite.exclude.filter='' test jacoco:report

