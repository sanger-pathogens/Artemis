#! /bin/bash

ant -Dlib=jacoco-lib -DEMBOSS_ROOT=$EMBOSS_ROOT -buildfile build-test.xml jacoco-coverage-report

