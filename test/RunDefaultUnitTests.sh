#! /bin/bash

SCRIPT_DIR=$(dirname $0)
ant -Dlib=${SCRIPT_DIR}/jacoco-lib -DEMBOSS_ROOT=$EMBOSS_ROOT -buildfile build-test.xml jacoco-coverage-report

