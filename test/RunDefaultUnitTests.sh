#! /bin/bash

mvn -DEMBOSS_ROOT=$EMBOSS_ROOT evosuite:prepare test -P dev

