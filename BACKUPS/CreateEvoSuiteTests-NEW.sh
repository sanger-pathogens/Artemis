#! /bin/bash

PACKAGE_CUTS_FILE="/tmp/evosuite-cuts.tmp"
BASEDIR=`dirname "$0"`

cd $BASEDIR/..
BASEDIR="`pwd`"

SCAFFOLD_FILE=$BASEDIR/.scaffolding_list.tmp

echo $BASEDIR

echo
echo "This may take several hours..."
echo
echo "EMBOSS_ROOT set to: "$EMBOSS_ROOT

#MAVEN_CMD="-DmemoryInMB=5000 -Dcores=2 -Dinfo -Ddetail=true -DEMBOSS_ROOT=$EMBOSS_ROOT -Dconsider_main_methods=false -DtimeInMinutesPerClass=4 -DtimeInMinutesPerProject=240 -DcutsFile=$PACKAGE_CUTS_FILE evosuite:generate evosuite:export"
#MAVEN_CMD="-DmemoryInMB=5000 -Dcores=2 -Dinfo -Ddetail=true -Dheadless_mode=true -Dconsider_main_methods=false -DtimeInMinutesPerClass=4 -DEMBOSS_ROOT=$EMBOSS_ROOT -DcutsFile=$PACKAGE_CUTS_FILE -Dp_reflection_on_private=0.5 -Dreflection_start_percent=0.8 -Dp_functional_mocking=0.8 -Dfunctional_mocking_percent=0.5 -Pdev-evo evosuite:generate evosuite:export "
MAVEN_CMD="-DmemoryInMB=5000 -Dcores=2 -Dinfo -Ddetail=true -Dheadless_mode=true -Dconsider_main_methods=false -DtimeInMinutesPerClass=4 -DEMBOSS_ROOT=$EMBOSS_ROOT -DcutsFile=$PACKAGE_CUTS_FILE -Pdev-evo evosuite:generate evosuite:export "

write_cuts_package_file() 
{
	search_path=$1
	output_file=$2
	
	rm -f $PACKAGE_CUTS_FILE
	
	ls -1 ${BASEDIR}/src/main/java/$search_path/*.java | egrep -oh "uk.*" | awk -F. '{print $1}' | tr '/' '.' | tr '\n' ',' > $output_file

}

log_package()
{
	package=$1
	echo "Creating tests for package ${package}..."
}


#
# Create tests for application Java package in turn...
# Evosuite gobbles memory so we do it one package at a time.
#

PACKAGE="uk/ac/sanger/artemis"
log_package $PACKAGE
write_cuts_package_file $PACKAGE $PACKAGE_CUTS_FILE
mvn $MAVEN_CMD

PACKAGE="uk/ac/sanger/artemis/sequence"
log_package $PACKAGE
write_cuts_package_file $PACKAGE $PACKAGE_CUTS_FILE
mvn $MAVEN_CMD

PACKAGE="uk/ac/sanger/artemis/circular"
log_package $PACKAGE
write_cuts_package_file $PACKAGE $PACKAGE_CUTS_FILE
mvn $MAVEN_CMD

PACKAGE="uk/ac/sanger/artemis/circular/digest"
log_package $PACKAGE
write_cuts_package_file $PACKAGE $PACKAGE_CUTS_FILE
mvn $MAVEN_CMD

PACKAGE="uk/ac/sanger/artemis/components"
log_package $PACKAGE
write_cuts_package_file $PACKAGE $PACKAGE_CUTS_FILE
mvn $MAVEN_CMD

PACKAGE="uk/ac/sanger/artemis/components/alignment"
log_package $PACKAGE
write_cuts_package_file $PACKAGE $PACKAGE_CUTS_FILE
mvn $MAVEN_CMD

PACKAGE="uk/ac/sanger/artemis/components/filetree"
log_package $PACKAGE
write_cuts_package_file $PACKAGE $PACKAGE_CUTS_FILE
mvn $MAVEN_CMD

PACKAGE="uk/ac/sanger/artemis/components/genebuilder"
log_package $PACKAGE
write_cuts_package_file $PACKAGE $PACKAGE_CUTS_FILE
mvn $MAVEN_CMD

PACKAGE="uk/ac/sanger/artemis/components/genebuilder/cv"
log_package $PACKAGE
write_cuts_package_file $PACKAGE $PACKAGE_CUTS_FILE
mvn $MAVEN_CMD

PACKAGE="uk/ac/sanger/artemis/components/genebuilder/gff"
log_package $PACKAGE
write_cuts_package_file $PACKAGE $PACKAGE_CUTS_FILE
mvn $MAVEN_CMD

PACKAGE="uk/ac/sanger/artemis/components/genebuilder/ortholog"
log_package $PACKAGE
write_cuts_package_file $PACKAGE $PACKAGE_CUTS_FILE
mvn $MAVEN_CMD

PACKAGE="uk/ac/sanger/artemis/components/variant"
log_package $PACKAGE
write_cuts_package_file $PACKAGE $PACKAGE_CUTS_FILE
mvn $MAVEN_CMD

PACKAGE="uk/ac/sanger/artemis/editor"
log_package $PACKAGE
write_cuts_package_file $PACKAGE $PACKAGE_CUTS_FILE
mvn $MAVEN_CMD

PACKAGE="uk/ac/sanger/artemis/io"
log_package $PACKAGE
write_cuts_package_file $PACKAGE $PACKAGE_CUTS_FILE
mvn $MAVEN_CMD

PACKAGE="uk/ac/sanger/artemis/j2ssh"
log_package $PACKAGE
write_cuts_package_file $PACKAGE $PACKAGE_CUTS_FILE
mvn $MAVEN_CMD

PACKAGE="uk/ac/sanger/artemis/plot"
log_package $PACKAGE
write_cuts_package_file $PACKAGE $PACKAGE_CUTS_FILE
mvn $MAVEN_CMD

PACKAGE="uk/ac/sanger/artemis/util"
log_package $PACKAGE
write_cuts_package_file $PACKAGE $PACKAGE_CUTS_FILE
mvn $MAVEN_CMD

echo "Setting separate classloader flag to false for generated tests"
find $BASEDIR/src/test/evosuite -type f -exec sed -i "" 's/separateClassLoader \= true/separateClassLoader \= false/g' {} +
    
# Manually create the scaffold list file [due to the test partitioning]
#
rm -f $SCAFFOLD_FILE
cd $BASEDIR/src/test/evosuite
find . -name '*_scaffolding.java' -type f -print | tr '/' '.' | sed 's/\.\.//g' | sed 's/\.java//g' > $SCAFFOLD_FILE
cd $BASEDIR

echo "Finished"
exit 0
