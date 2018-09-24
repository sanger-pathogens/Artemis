#! /bin/bash

echo
echo "This may take several hours..."
echo
echo "EMBOSS_ROOT set to: "$EMBOSS_ROOT

export JAVA_TOOL_OPTIONS="-Xmx5g"

cd ..

#mvn -X -Dcores=2 -DmemoryInMB=2000 -Dconsider_main_methods=false -Duse_separate_classloader=false -DEMBOSS_ROOT=$EMBOSS_ROOT clean compile evosuite:clean evosuite:generate evosuite:export test-compile evosuite:prepare
#mvn -DEMBOSS_ROOT=$EMBOSS_ROOT test-compile evosuite:prepare
mvn -Dcores=2 -DmemoryInMB=2000 -Dsearch_budget=90 -Dconsider_main_methods=false -Duse_separate_classloader=false -DEMBOSS_ROOT=$EMBOSS_ROOT evosuite:generate evosuite:export test-compile evosuite:prepare


exit 0
