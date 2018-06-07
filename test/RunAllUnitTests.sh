#! /bin/bash

echo
echo "This may take 30 minutes or more..."
echo

sleep 5
cd ..

# No instrumentation
mvn -Dinfo=true -DEMBOSS_ROOT=$EMBOSS_ROOT clean evosuite:prepare test jacoco:report
# Instrumentation
#mvn -DEMBOSS_ROOT=$EMBOSS_ROOT clean evosuite:prepare test jacoco:restore-instrumented-classes jacoco:report -f pom-jacoco-offline-instrumentation.xml

