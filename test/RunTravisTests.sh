#! /bin/bash

############################################################
# Run a selection of tests for a Travis build that will 
# fit within the Travis 50 minute build window.
############################################################

echo
echo "EMBOSS_ROOT set to: "$EMBOSS_ROOT
echo

cd ..

mvn test jacoco:report

