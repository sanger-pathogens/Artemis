---
runme:
  document:
    relativePath: README.md
  session:
    id: 01JNFN0AQWHADNKNGP97SGWKRE
    updated: 2025-03-03 23:33:53-05:00
---

# The Artemis Software

The Artemis Software is a set of software tools for genome browsing and annotation.

[![Build Status](ht************************************************************er)](ht******************************************is)  
[![License: GPL v3](ht*********************************************************vg)](ht***********************************************************SE)  
[![install with bioconda](ht**********************************************************************************re)](ht*************************************************ml)  
[![st**us](ht********************************************************************************************vg)](ht********************************************44)  
[![st**us](ht*****************************************************************************************vg)](ht*****************************************03)  
[![st**us](ht*****************************************************************************************vg)](ht*****************************************29)  
[![st**us](ht********************************************************************vg)](ht*******************************24)  
[![st**us](ht*****************************************************************************************vg)](ht*****************************************53)  
[![st**us](ht*******************************************************************vg)](ht******************************73)  
[![st**us](ht*****************************************************************************************vg)](ht*****************************************10)  
[![st**us](ht*****************************************************************************************vg)](ht*****************************************78)

![GitHub Releases (by Release)](ht**************************************************************************al)

## Contents

* [Introduction](#introduction)
* [Installation](#installation)
* [Usage](#usage)
* [Building](#building)
* [License](#license)
* [Documentation](#documentation)

## Introduction

The Artemis Software includes:

* Artemis
* ACT
* BamView
* DNAPlotter

Artemis is a free genome browser and annotation tool that allows visualisation of sequence features, next generation data and the results of analyses within the context of the sequence, and also its six-frame translation. Artemis is written in Java, and is available for UNIX, Macintosh and Windows systems. It can read EMBL and GENBANK database entries or sequence in FASTA, indexed FASTA or raw format. Other sequence features can be in EMBL, GENBANK or GFF format.

ACT is a free tool for displaying pairwise comparisons between two or more DNA sequences. It can be used to identify and analyse regions of similarity and difference between genomes and to explore conservation of synteny, in the context of the entire sequences and their annotation.

DNAPlotter generates images of circular and linear DNA maps.

BamView is a standalone BAM/CRAM file viewer.

For further information and details on how to download the Artemis software, please see the [Artemis GitHub page](ht*************************************is/)

## Prerequisites

Java version 9 or later must be installed. Ideally Java 11 should be used. The easiest way to install a non-commercial open source Java version is from [AdoptOpenJDK](ht**********************************ml) - just select the OpenJDK 11 and Hotspot options for the relevant platform. See the [GitHub pages and manuals](#documentation) for further options. A Java installation is not required if installing via Bioconda or Docker.

## Installation

Please refer to the relevant documentation in the [Documentation](#documentation) section below.

## Usage

Please refer to the relevant documentation in the [Documentation](#documentation) section below.

## Building

If you would prefer to build the applications from scratch rather than use the pre-built releases, then you will firstly need the source code from GitHub. This can be obtained by cloning the repository:

```sh
git clone ht*******************************************it
```

or by downloading the source zip file for a particular release. The latest version of Apache Maven (Java 11 compatible, e.g. v3**.0) will need to be installed beforehand, in order to build. The application jars can be built by issuing the following commands in the top-level folder:

When building for the first time:

```sh
mvn validate
```

This is required to install some legacy libraries to the Maven repository and will be phased out in future releases.

And then to build/rebuild the applications:

```sh
mvn clean package
```

Clean is only required if you wish to do a complete rebuild.
If you do not wish to run tests then add the <i>-Dskip.tests=true</i> flag to the above command.

Note that if you are running a build from behind a proxy you will need to add the proxy parameters to the mvn command line, e.g.

```sh
mvn -Dhttps.proxyHost=myproxyhost -Dhttps.proxyPort=myproxyport -DproxySet=true package

[**FO] Scanning for projects...
[**FO] 
[**FO] ------------------------< r:ar***is >------------------------
[**FO] ******ng Artemis 18.2.0
[**FO]   from pom.xml
[**FO] --------------------------------[ jar ]---------------------------------
[*****NG] *******er 'skipAssembly' is unknown for plugin 'ma**************in:3.1.0:shade (package-win-artemis-jar)'
[*****NG] *******er 'skipAssembly' is unknown for plugin 'ma**************in:3.1.0:shade (package-default-artemis-jar)'
[*****NG] *******er 'skipAssembly' is unknown for plugin 'ma**************in:3.1.0:shade (package-win-act-jar)'
[*****NG] *******er 'skipAssembly' is unknown for plugin 'ma**************in:3.1.0:shade (package-default-act-jar)'
[*****NG] *******er 'skipAssembly' is unknown for plugin 'ma**************in:3.1.0:shade (package-win-dnaplotter-jar)'
[*****NG] *******er 'skipAssembly' is unknown for plugin 'ma**************in:3.1.0:shade (package-default-dnaplotter-jar)'
[*****NG] *******er 'skipAssembly' is unknown for plugin 'ma**************in:3.1.0:shade (package-win-bamview-jar)'
[*****NG] *******er 'skipAssembly' is unknown for plugin 'ma**************in:3.1.0:shade (package-default-bamview-jar)'
[**FO] 
[**FO] --- r:1.4.1:en***ce **************-3) @ s ---
[**FO] 
[**FO] --- l:2.5.2:in********le ******************ar) @ s ---
[**FO] pom.xml not found in biojava.jar
[**FO] Installing /Users/nilarmyint/Documents/GitHub/Awesome/LogicAppsUX/Artemis/lib/biojava.jar to /Users/nilarmyint/.m2/repository/org/biojava/biojava/1.6/bi***********ar
[**FO] Installing /var/folders/g7/8d7h7xrj0mg58zvpbcdlhcc40000gn/T/mv***************************1.pom to /Users/nilarmyint/.m2/repository/org/biojava/biojava/1.6/bi***********om
[**FO] 
[**FO] --- l:2.5.2:in********le *******************ar) @ s ---
[**FO] pom.xml not found in jemAlign.jar
[**FO] Installing /Users/nilarmyint/Documents/GitHub/Awesome/LogicAppsUX/Artemis/lib/jemAlign.jar to /Users/nilarmyint/.m2/repository/org/emboss/jemboss/1.0/je***********ar
[**FO] Installing /var/folders/g7/8d7h7xrj0mg58zvpbcdlhcc40000gn/T/mv**************************0.pom to /Users/nilarmyint/.m2/repository/org/emboss/jemboss/1.0/je***********om
[**FO] 
[**FO] --- l:2.5.2:in********le *********************ar) @ s ---
[**FO] pom.xml not found in j2**********ar
[**FO] Installing /Users/nilarmyint/Documents/GitHub/Awesome/LogicAppsUX/Artemis/lib/j2**********ar to /Users/nilarmyint/.m2/repository/com/sshtools/j2******re/0.2.9/j2****************ar
[**FO] Installing /var/folders/g7/8d7h7xrj0mg58zvpbcdlhcc40000gn/T/mv***************************5.pom to /Users/nilarmyint/.m2/repository/com/sshtools/j2******re/0.2.9/j2****************om
[**FO] 
[**FO] --- l:2.5.2:in********le *****************ar) @ s ---
[**FO] pom.xml not found in ib****************ar
[**FO] Installing /Users/nilarmyint/Documents/GitHub/Awesome/LogicAppsUX/Artemis/lib/ib****************ar to /Users/nilarmyint/.m2/repository/com/ibatis/ibatis/2.*****26/ib****************ar
[**FO] Installing /var/folders/g7/8d7h7xrj0mg58zvpbcdlhcc40000gn/T/mv***************************6.pom to /Users/nilarmyint/.m2/repository/com/ibatis/ibatis/2.*****26/ib****************om
[**FO] 
[**FO] --- s:1.0.0:re*******************es ***********************es) @ s ---
[**FO] Quiet processing - ignoring properties cannot be loaded from File: /Users/nilarmyint/Documents/GitHub/Awesome/LogicAppsUX/Artemis/maven.properties
[**FO] 
[**FO] --- n:1.8:run *************************n-nu**er) @ s ---
[**FO] Executing tasks

main:
[**FO] Executed tasks
[**FO] 
[**FO] --- o:0.8.2:pr*********nt ************st) @ s ---
[**FO] Skipping JaCoCo execution because property jacoco.skip is set.
[**FO] 
[**FO] --- s:3.0.2:re*****es ****************es) @ s ---
[**FO] Using 'UTF-8' encoding to copy filtered resources.
[**FO] Copying 35 resources
[**FO] 
[**FO] --- s:3.0.2:co**********es ******************es) @ s ---
[**FO] Using 'UTF-8' encoding to copy filtered resources.
[**FO] Copying 37 resources
[**FO] 
[**FO] --- s:3.0.2:co**********es ****************er) @ s ---
[**FO] Using 'UTF-8' encoding to copy filtered resources.
[**FO] Copying 5 resources
[**FO] 
[**FO] --- r:3.8.0:co***le **************le) @ s ---
[**FO] Nothing to compile - all classes are up to date
[**FO] 
[**FO] --- s:3.0.2:te*********es ********************es) @ s ---
[**FO] Using 'UTF-8' encoding to copy filtered resources.
[**FO] Copying 60 resources
[**FO] 
[**FO] --- r:3.8.0:te*******le ******************le) @ s ---
[**FO] Nothing to compile - all classes are up to date
[**FO] 
[**FO] --- e:2.**.1:test ***********st) @ s ---
[**FO] 
[**FO] -------------------------------------------------------
[**FO]  T E S T S
[**FO] -------------------------------------------------------
[**FO] Running uk*********************ot.**********st
[**FO] ***ts run: 3, Failures: 0, Errors: 0, Skipped: 0, Time elapsed: 1.715 s - in uk*********************ot.**********st
[**FO] Running uk*********************il.*************st
[**FO] ***ts run: 4, Failures: 0, Errors: 0, Skipped: 0, Time elapsed: 1.649 s - in uk*********************il.*************st
[**FO] Running uk****************is.****************st
[**FO] ***ts run: 6, Failures: 0, Errors: 0, Skipped: 0, Time elapsed: 1.339 s - in uk****************is.****************st
[**FO] Running uk*******************io.*************************st

```

This will build the application jars and place them in target/jars and win-jars folders. The win-jars jars have the etc folder files bundled in, for Windows.

## Building Release Artifacts

To build .zip or .gz installables with unsigned jars, use the following command:

```sh
mvn -Djarsigner.skip=true clean package -P release
```

If building on a Mac then .app packages and a .dmg image can additionally be produced using:

```sh
mvn -Djarsigner.skip=true clean package verify -P release
```

The release artifacts can then be found in the target/release-artifacts directory.

Jar signing can be included in the build cycle by <b>not</b> supplying a -Djarsigner.skip flag. A <i>maven.properties</i> file needs to be placed in the top-level Artemis folder containing the following keystore properties:

```ini
signer-keystore-path=<The path to a .jks keystore>
signer-keystore-alias=<certificate alias name>
signer-keystore-password=<keystore password>
signer-keystore-type=JKS
```

## Searching and Using Local Sequence Databases

Artemis and ACT can optionally be used with local sequence databases on non-Windows systems.

[Note: this may require a lot of disk space]

If you wish to do this then carry out the following steps:

1. To download and install BLAST+ locally, check here:

```sh
ft******************************************************ST/
```

to get the latest Mac OS X release:

```sh
ft******************************************************************************************gz
```

Uncompress, if you have downloaded a .tar.gz, e.g.

```sh
tar -zxvf nc****************************ar.gz
```

For a Mac OS X install you can mount the resulting .dmg image by double-clicking it and then selecting the .pkg installer file to actually install the BLAST+ executables into the default location (/usr/local/ncbi/blast). Comprehensive installation instructions for all operating systems are provided  [here](ht****************************************90/).

2. Ensure that the BLAST+ executables are available in your PATH and that local
   databases can be found by Artemis:

On the command line type:

```sh
blastp -help
```

If the executables are not found then their installation bin folder will need to be added to your .profile PATH variable. This can be achieved by uncommenting and setting the relevant environment variables in the **Variables for locally installed Blast databases** section of the Artemis **setenv** script. The location of the databases is assumed by default to be &lt;home directory&gt;/blast-data but this can be overridden by setting the BLASTDB variable (See Blast+ documentation). The setenv script should then be called from your login **.profile** file. Alternatively copy the relevant variables into your profile.

Note that database locations can also be specified directly within the Artemis/ACT applications using the run menu - set options functionality.

4. Download and format the uniprot database(s):

You can do this by using the Artemis etc/setup_uniprot_dbs.sh script. This script can be changed as desired, to for instance, include additional databases.

```sh
mkdir $HOME/blast-data
cd blast-data
cp <Artemis folder>/etc/setup_blast_dbs.sh .
./setup_blast_dbs.sh
```

(Note: the environment variable http_proxy may need to be set for the ftp downloads to work).

## License

Artemis is free software, licensed under [GPLv3](ht***********************************************************SE).

For more information on how to download Artemis, please see the [Artemis GitHub page](ht*************************************is/)

## Documentation

Artemis GitHub pages:
ht*************************************is/

The Artemis user manual is at:
ht*********************************************is/

The ACT user manual is at:
ht*****************************************CT/

The DNA plotter documentation is at:
ht************************************************er/

The BAM View documentation is at:
ht*********************************************ew/
