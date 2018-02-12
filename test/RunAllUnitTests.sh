#! /bin/bash

echo
echo "This may take 30 minutes or more..."
echo

sleep 5

SCRIPT_DIR=$(dirname $0)
ant -Dlib=${SCRIPT_DIR}/jacoco-lib -DEMBOSS_ROOT=$EMBOSS_ROOT -buildfile build-test.xml testall

