#! /bin/bash

echo
echo "This may take several hours and needs a lot of heap space...."
echo
echo "EMBOSS_ROOT set to: "$EMBOSS_ROOT
echo

export JAVA_TOOL_OPTIONS="-Xmx5g"

cd ..

mvn -Dcores=2 -DmemoryInMB=2000 -Dsearch_budget=90 -Dconsider_main_methods=false -Duse_separate_classloader=false -DEMBOSS_ROOT=$EMBOSS_ROOT evosuite:generate evosuite:export test-compile evosuite:prepare


exit 0
