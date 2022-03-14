#! /bin/bash

##########################################################
# Run tests. 
##########################################################

echo "EMBOSS_ROOT set to: "$EMBOSS_ROOT
echo

cd ..

mvn -DEMBOSS_ROOT=${EMBOSS_ROOT} test

