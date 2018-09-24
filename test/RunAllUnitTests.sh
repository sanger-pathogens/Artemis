#! /bin/bash

echo
echo "This may take 30 minutes or more..."
echo

sleep 5
cd ..

mvn -Dinfo=true -DEMBOSS_ROOT=$EMBOSS_ROOT test-compile evosuite:prepare -P dev-evo

if [[ ! -e ./.scaffolding_list.tmp ]]
then
	echo "ERROR: Cannot find Evosuite .scaffolding_list.tmp file - unit tests will not be sand-boxed! Exiting."
	exit 1
fi 

mvn -Dinfo=true -DEMBOSS_ROOT=$EMBOSS_ROOT test jacoco:report -P dev-evo

