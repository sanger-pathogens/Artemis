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

For further information, please see the [Artemis GitHub page](http://sanger-pathogens.github.io/Artemis/)

## Installation

Please refer to the relevant documentation in the Documentation section below.

## Usage

Please refer to the relevant documentation in the Documentation section below.

## Building

If you would prefer to build the applications from scratch rather than use the pre-built releases in GitHub, then you will firstly need the source code from GitHub. This can be obtained by cloning the repository:

git clone http://github.com/sanger-pathogens/Artemis.git

or by downloading the source zip file for a particular release. Apache Maven and Java 8 will need to be installed beforehand, in order to build. The applications can be built by issuing the following command in the top-level folder:

mvn validate clean test package

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
