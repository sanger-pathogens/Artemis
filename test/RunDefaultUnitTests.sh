#! /bin/bash

##########################################################
# Run all non-evosuite unit tests, which is the 
# default option. 
##########################################################

echo "EMBOSS_ROOT set to: "$EMBOSS_ROOT
echo

cd ..

mvn test

