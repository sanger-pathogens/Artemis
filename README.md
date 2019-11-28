# The Artemis Software
The Artemis Software is a set of software tools for genome browsing and annotation.

[![Build Status](https://travis-ci.org/sanger-pathogens/Artemis.svg?branch=master)](https://travis-ci.org/sanger-pathogens/Artemis)  
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/Artemis/blob/master/LICENSE)  
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/artemis/README.html)  
![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/sangerpathogens/artemis)  
[![status](https://img.shields.io/badge/Bioinformatics-10.1093%2Fbioinformatics%2F16.10.944-brightgreen.svg)](https://doi.org/10.1093/bioinformatics/16.10.944)  
[![status](https://img.shields.io/badge/Bioinformatics-10.1093%2Fbioinformatics%2Fbtr703-brightgreen.svg)](https://doi.org/10.1093/bioinformatics/btr703)   
[![status](https://img.shields.io/badge/Bioinformatics-10.1093%2Fbioinformatics%2Fbtn529-brightgreen.svg)](https://doi.org/10.1093/bioinformatics/btn529)   
[![status](https://img.shields.io/badge/BIB-10.1093%2Fbib%2F4.2.124-brightgreen.svg)](https://doi.org/10.1093/bib/4.2.124)   
[![status](https://img.shields.io/badge/Bioinformatics-10.1093%2Fbioinformatics%2Fbti553-brightgreen.svg)](https://doi.org/10.1093/bioinformatics/bti553)   
[![status](https://img.shields.io/badge/BIB-10.1093%2Fbib%2Fbbr073-brightgreen.svg)](https://doi.org/10.1093/bib/bbr073)   
[![status](https://img.shields.io/badge/Bioinformatics-10.1093%2Fbioinformatics%2Fbtq010-brightgreen.svg)](https://doi.org/10.1093/bioinformatics/btq010)   
[![status](https://img.shields.io/badge/Bioinformatics-10.1093%2Fbioinformatics%2Fbtn578-brightgreen.svg)](https://doi.org/10.1093/bioinformatics/btn578)  

![GitHub Releases (by Release)](https://img.shields.io/github/downloads/sanger-pathogens/artemis/v18.1.0/total)

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

For further information and details on how to download the Artemis software, please see the [Artemis GitHub page](http://sanger-pathogens.github.io/Artemis/)

## Prerequisites

Java version 9 or later must be installed. Ideally Java 11 should be used. The easiest way to install a non-commercial open source Java version is from [AdoptOpenJDK](https://adoptopenjdk.net/releases.html) - just select the OpenJDK 11 and Hotspot options for the relevant platform. See the [GitHub pages and manuals](#documentation) for further options. A Java installation is not required if installing via Bioconda or Docker.

## Installation

Please refer to the relevant documentation in the [Documentation](#documentation) section below.

## Usage

Please refer to the relevant documentation in the [Documentation](#documentation) section below.

## Building

If you would prefer to build the applications from scratch rather than use the pre-built releases, then you will firstly need the source code from GitHub. This can be obtained by cloning the repository:
```
git clone https://github.com/sanger-pathogens/Artemis.git
```
or by downloading the source zip file for a particular release. The latest version of Apache Maven (Java 11 compatible, e.g. v3.6.0) will need to be installed beforehand, in order to build. The application jars can be built by issuing the following commands in the top-level folder:

When building for the first time:
```
mvn validate
```
This is required to install some legacy libraries to the Maven repository and will be phased out in future releases.

And then to build/rebuild the applications:
```
mvn clean package
```
Clean is only required if you wish to do a complete rebuild.
If you do not wish to run tests then add the <i>-Dskip.tests=true</i> flag to the above command.

Note that if you are running a build from behind a proxy you will need to add the proxy parameters to the mvn command line, e.g.
```
mvn -Dhttps.proxyHost=myproxyhost -Dhttps.proxyPort=myproxyport -DproxySet=true package
```
This will build the application jars and place them in target/jars and win-jars folders. The win-jars jars have the etc folder files bundled in, for Windows.

## Building Release Artifacts
To build .zip or .gz installables with unsigned jars, use the following command:
```
mvn -Djarsigner.skip=true clean package -P release
```
If building on a Mac then .app packages and a .dmg image can additionally be produced using:
```
mvn -Djarsigner.skip=true clean package verify -P release
```
The release artifacts can then be found in the target/release-artifacts directory.

Jar signing can be included in the build cycle by <b>not</b> supplying a -Djarsigner.skip flag. A <i>maven.properties</i> file needs to be placed in the top-level Artemis folder containing the following keystore properties:
```
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

  ```
  ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
  ```
  to get the latest Mac OS X release:

  ```
  ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.7.1+-x64-macosx.tar.gz
  ```
  Uncompress, if you have downloaded a .tar.gz, e.g.

  ```
  tar -zxvf ncbi-blast-2.7.1+-x64-macosx.tar.gz
  ```
  For a Mac OS X install you can mount the resulting .dmg image by double-clicking it and then selecting the .pkg installer file to actually install the BLAST+ executables into the default location (/usr/local/ncbi/blast). Comprehensive installation instructions for all operating systems are provided  [here](https://www.ncbi.nlm.nih.gov/books/NBK279690/).

2. Ensure that the BLAST+ executables are available in your PATH and that local
   databases can be found by Artemis:

  On the command line type:

  ```
  blastp -help
  ```
  If the executables are not found then their installation bin folder will need to be added to your .profile PATH variable. This can be achieved by uncommenting and setting the relevant environment variables in the **Variables for locally installed Blast databases** section of the Artemis **setenv** script. The location of the databases is assumed by default to be &lt;home directory&gt;/blast-data but this can be overridden by setting the BLASTDB variable (See Blast+ documentation). The setenv script should then be called from your login **.profile** file. Alternatively copy the relevant variables into your profile.

  Note that database locations can also be specified directly within the Artemis/ACT applications using the run menu - set options functionality.

4. Download and format the uniprot database(s):

  You can do this by using the Artemis etc/setup_uniprot_dbs.sh script. This script can be changed as desired, to for instance, include additional databases.

  ```
  mkdir $HOME/blast-data
  cd blast-data
  cp <Artemis folder>/etc/setup_blast_dbs.sh .
  ./setup_blast_dbs.sh
  ```

  (Note: the environment variable http_proxy may need to be set for the ftp downloads to work).

## License
Artemis is free software, licensed under [GPLv3](https://github.com/sanger-pathogens/artemis/blob/master/LICENSE).

For more information on how to download Artemis, please see the [Artemis GitHub page](http://sanger-pathogens.github.io/Artemis/)

## Documentation

Artemis GitHub pages:
  http://sanger-pathogens.github.io/Artemis/

The Artemis user manual is at:
  http://sanger-pathogens.github.io/Artemis/Artemis/

The ACT user manual is at:
  http://sanger-pathogens.github.io/Artemis/ACT/

The DNA plotter documentation is at:
  http://sanger-pathogens.github.io/Artemis/Artemis/DNAPlotter/

The BAM View documentation is at:
  http://sanger-pathogens.github.io/Artemis/Artemis/BamView/
