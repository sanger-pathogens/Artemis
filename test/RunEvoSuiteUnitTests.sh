#! /bin/bash

echo
echo "This may take 30 minutes or more..."
echo

sleep 5

mvn -DEMBOSS_ROOT=$EMBOSS_ROOT evosuite:prepare test -P dev-evo

