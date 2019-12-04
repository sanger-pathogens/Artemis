# The Artemis Software
## A set of software tools for genome browsing and annotation

## Overview
The Artemis Software is a set of software tools for genome browsing and annotation. It includes:

* [Artemis](Artemis/)
* [Artemis Comparison Tool (ACT)](ACT/)
* [BamView](BamView/)
* [DNAPlotter](DNAPlotter)

## Software Availability
The Artemis Software is available under GPL3. The source code can be found on [GitHub](https://github.com/sanger-pathogens/Artemis).

### Download

The latest release of Artemis can be downloaded by clicking on the relevant link below:

* [UNIX](https://github.com/sanger-pathogens/Artemis/releases/download/v18.1.0/artemis-unix-release-18.1.0.tar.gz)
* [MacOS](https://github.com/sanger-pathogens/Artemis/releases/download/v18.1.0/artemis-macosx-release-18.1.0.dmg.gz)
* [MacOS (CHADO) - for out of the box CHADO database connectivity](https://github.com/sanger-pathogens/Artemis/releases/download/v18.1.0/artemis-macosx-chado-release-18.1.0.dmg.gz)
* [Windows](https://github.com/sanger-pathogens/Artemis/releases/download/v18.1.0/artemis-windows-release-18.1.0.zip)

__Note on Java versions__: The old v17.0.1 version of the Artemis software required Java version 1.8 to run. All recent releases from v18.0.0 onwards require a minimum of Java 9 and ideally Java 11. This must be installed first. The easiest way to install a non-commercial open source Java version is from [AdoptOpenJDK](https://adoptopenjdk.net/releases.html) - just select the OpenJDK version and Hotspot options for the relevant platform. See the [manuals](#documentation) for further options. A Java installation is not required for the Bioconda or Docker options.

__Note for MacOSX__: occasionally a browser decides to display the contents of the .dmg.gz archive file rather than downloading it. If this happens hold down the <control> key and click on the download link. A popup menu should appear, containing several choices. One of the choices should be something like "Save Link As" (or perhaps "Download Link...", "Save Link to Desktop", or a variation on this theme). Select that option, and the archive file should be download correctly.

For older versions of the Artemis software please see the [Artemis FTP site](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/)

### Bioconda

Via [Bioconda](https://bioconda.github.io). Simply use:
```
conda config --add channels bioconda     (to add the bioconda channel)
conda config --add channels conda-forge  (to add the conda-forge channel)
conda install artemis
```
from the command line, and you're ready to go (no Java installation required).

### Docker (experimental)

Via [Docker](https://hub.docker.com/r/sangerpathogens/artemis):

Firstly, an X server must be installed and running on the host machine. XQuartz can be used on Mac OSX. Attention should be paid to authentication and any timeout settings.

To pull the Docker image from [Docker Hub](https://hub.docker.com/r/sangerpathogens/artemis), use:
```
docker pull sangerpathogens/artemis
```
To run a Docker container use the following command:
```
docker run -d -e DISPLAY="<your display name>:0" \
  -v <your data folder>:/artemis -v <your user home directory>:/home/artuser \
  --user $(id -u):$(id -g) \
  -e ARTEMIS_JVM_FLAGS="-Duser.home=/home/artuser -Djava.io.tmpdir=/tmp" \
  --tmpfs /tmp \
  --rm artemis <program name> \
  <Any additonal program arguments>
```
Where **\<program name\>** is one of art, act, bamview or dnaplotter and **\<your data folder\>** represents the path of a folder that contains the data files that you wish to use in the Artemis applications. Memory allocation can be changed by passing the -mx and -ms arguments in the ARTEMIS_JVM_FLAGS variable above, for example: -mx1g for a maximum 1Gb memory allocation.

Current Docker image limitations:
1. Printing does not work
2. FTP download of bam files does not work

Files will be available from the /artemis folder when running an application via Docker.

Firefox is invoked for display of some results such as rfam/pfam. The Firefox instance
is not meant for general web browsing, only for viewing analysis results.

## Documentation

Specific information is provided on the individual application pages (see [overview](#overview) section above). Application manuals can be found here:

* Artemis manual - [online](https://sanger-pathogens.github.io/Artemis/Artemis/artemis-manual.html) or [PDF](https://sanger-pathogens.github.io/Artemis/Artemis/artemis-manual.pdf).

* ACT manual - [online](https://sanger-pathogens.github.io/Artemis/ACT/act-manual.html) or [PDF](https://sanger-pathogens.github.io/Artemis/ACT/act-manual.pdf).

The GitHub project [README](https://github.com/sanger-pathogens/Artemis/blob/master/README.md) also has some useful information.

## Courses

The Wellcome Genome Campus currently offers online Bacterial Genome Artemis & ACT courses via FutureLearn:

* [Bacterial Genomes: From DNA to Protein Function Using Bioinformatics](https://www.futurelearn.com/courses/bacterial-genomes-bioinformatics)
* [Bacterial Genomes: Accessing and Analysing Microbial Genome Data](https://www.futurelearn.com/courses/bacterial-genomes-access-and-analysis)

A [Bacterial Genomics Introductory Course](https://www.futurelearn.com/courses/introduction-to-bacterial-genomics) is also available.

A complete list of courses can be found on the [Wellcome Genome Campus advanced courses web site](https://coursesandconferences.wellcomegenomecampus.org/event-type/courses/).
