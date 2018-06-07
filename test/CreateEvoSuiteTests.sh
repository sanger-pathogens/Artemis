#! /bin/bash

echo
echo "This may take several hours..."
echo
echo "EMBOSS_ROOT set to: "$EMBOSS_ROOT

mvn -Dconsider_main_methods=false -DEMBOSS_ROOT=$EMBOSS_ROOT evosuite:clean evosuite:generate evosuite:export evosuite:prepare
#mvn -Dcores=1 -DmemoryInMB=5000 -Dconsider_main_methods=false -Duse_separate_classloader=true -DEMBOSS_ROOT=$EMBOSS_ROOT evosuite:clean evosuite:generate evosuite:export evosuite:prepare
#mvn -DmemoryInMB=4000 -Dcores=4 -Dconsider_main_methods=false -DtimeInMinutesPerClass=1 -Duse_separate_classloader=false -DEMBOSS_ROOT=$EMBOSS_ROOT evosuite:clean clean package evosuite:generate evosuite:export evosuite:prepare
#mvn -DmemoryInMB=2000 -Dcores=2 -DEMBOSS_ROOT=$EMBOSS_ROOT evosuite:generate evosuite:export

exit 0
