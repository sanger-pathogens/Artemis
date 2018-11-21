# The Artemis Software
The Artemis Software is a set of software tools for genome browsing and annotation.

[![Build Status](https://travis-ci.org/sanger-pathogens/Artemis.svg?branch=master)](https://travis-ci.org/sanger-pathogens/Artemis)  
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/Artemis/blob/master/LICENSE)  
[![status](https://img.shields.io/badge/Bioinformatics-10.1093%2Fbioinformatics%2F16.10.944-brightgreen.svg)](https://doi.org/10.1093/bioinformatics/16.10.944)  
[![status](https://img.shields.io/badge/Bioinformatics-10.1093%2Fbioinformatics%2Fbtr703-brightgreen.svg)](https://doi.org/10.1093/bioinformatics/btr703)   
[![status](https://img.shields.io/badge/Bioinformatics-10.1093%2Fbioinformatics%2Fbtn529-brightgreen.svg)](https://doi.org/10.1093/bioinformatics/btn529)   
[![status](https://img.shields.io/badge/BIB-10.1093%2Fbib%2F4.2.124-brightgreen.svg)](https://doi.org/10.1093/bib/4.2.124)   
[![status](https://img.shields.io/badge/Bioinformatics-10.1093%2Fbioinformatics%2Fbti553-brightgreen.svg)](https://doi.org/10.1093/bioinformatics/bti553)   
[![status](https://img.shields.io/badge/BIB-10.1093%2Fbib%2Fbbr073-brightgreen.svg)](https://doi.org/10.1093/bib/bbr073)   
[![status](https://img.shields.io/badge/Bioinformatics-10.1093%2Fbioinformatics%2Fbtq010-brightgreen.svg)](https://doi.org/10.1093/bioinformatics/btq010)   
[![status](https://img.shields.io/badge/Bioinformatics-10.1093%2Fbioinformatics%2Fbtn578-brightgreen.svg)](https://doi.org/10.1093/bioinformatics/btn578)   

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

For further information, please see the [Artemis GitHub page](http://sanger-pathogens.github.io/Artemis/)

## Prerequisites

Java 9+ must be installed. Ideally Java 11.

## Installation

Please refer to the relevant documentation in the Documentation section below.

## Usage

Please refer to the relevant documentation in the Documentation section below.

## Building

If you would prefer to build the applications from scratch rather than use the pre-built releases in GitHub, then you will firstly need the source code from GitHub. This can be obtained by cloning the repository:
```
git clone https://github.com/sanger-pathogens/Artemis.git
```
or by downloading the source zip file for a particular release. The latest version of Apache Maven (Java 11 compatible, e.g. v3.6.0) will need to be installed beforehand, in order to build. The application jars can be built by issuing the following commands in the top-level folder:

When building for the first time:
```
mvn validate
```
This is required to install some legacy libraries to the Maven repository and will be phased out in future releases.

And then to build/rebuild the applications :
```
mvn clean package
```
Clean is only required if you wish to do a complete rebuild.
If you do not wish to run tests then add the <i>-Dskip.tests=true</i> flag to the above command.

Note that if you are running a build from behind a proxy you will need to add the proxy parameters to the mvn command line, e.g.
```
mvn -Dhttps.proxyHost=myproxyhost -Dhttps.proxyPort=myproxyport -DproxySet=true package
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

2. To download and install Fasta locally, check here:

  ```
  ftp://ftp.ebi.ac.uk/pub/software/unix/fasta/
  ```

  and download the latest release, e.g.

  ```
  ftp://ftp.ebi.ac.uk/pub/software/unix/fasta/fasta3/fasta-36.3.8g.tar.gz
  ```
  Uncompress the tar:

  ```
  tar -zxvf fasta-36.3.8g.tar.gz
  ```

  Build the executables, e.g. for a Mac:

  ```
  cd fasta/src
  make -f ../make/Makefile.os_x86_64 all
  ```
  Move the fasta folder to a convenient area, preferably /usr/local/ for Macs:

3. Ensure that the BLAST+ and Fasta executables are available in your PATH and that local databases can be found by Artemis:

  On the command line type:

  ```
  blastp -help
  fasta36 -help
  ```
  If the executables are not found then their respective installation bin folders will need to be added to your .profile PATH variable. This can be achieved by uncommenting and setting the relevant environment variables in the **Variables for locally installed Blast databases** section of the Artemis **setenv** script. The location of the databases is assumed by default to be &lt;home directory&gt;/blast-data but this can be overridden by setting the BLASTDB and FASTLIBS variables. The setenv script should then be called from your login **.profile** file. Alternatively copy the relevant variables into your profile.

  Note that database locations can also be specified directly within the Artemis/ACT applications.

4. Download and format the uniprot database(s):

  You can do this by using the Artemis etc/setup_uniprot_dbs.sh script. This script can be changed as desired, to for instance, include additional databases or create a FASTLIBS file.

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

The Artemis user manual is at:
  http://sanger-pathogens.github.io/Artemis/Artemis/

The ACT user manual is at:
  http://sanger-pathogens.github.io/Artemis/ACT/

The DNA plotter documentation is at:
  http://sanger-pathogens.github.io/Artemis/Artemis/DNAPlotter/

The BAM View documentation is at:
  http://sanger-pathogens.github.io/Artemis/Artemis/BamView/
