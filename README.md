PLEASE NOTE THAT THIS SITE IS UNDER CONSTRUCTION
# The Artemis Software - a set of software tools for genome browsing and annotation

## Overview
The Artemis Software is a set of software tools for genome browsing and annotation. It includes:

* Artemis
* [ACT](ACT/)
* BamView
* DNAPlotter

# Artemis
Artemis is a free genome browser and annotation tool that allows visualisation of sequence features, next generation data and the results of analyses within the context of the sequence, and also its six-frame translation.
Artemis is written in Java, and is available for UNIX, Macintosh and Windows systems. It can read EMBL and GENBANK database entries or sequence in FASTA, indexed FASTA or raw format. Other sequence features can be in EMBL, GENBANK or GFF format.
Full information about the latest release of Artemis can be found in the [Artemis manual](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/artemis.pdf) and the [current release notes](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/release_notes.txt).

## Citation

Artemis is a DNA sequence visualization and annotation tool that allows the results of any analysis or sets of analyses to be viewed in the context of the sequence and its six-frame translation. Artemis is especially useful in analysing the compact genomes of bacteria, archaea and lower eukaryotes, and will cope with sequences of any size from small genes to whole genomes. It is implemented in Java, and can be run on any suitable platform. Sequences and annotation can be read and written directly in EMBL, GenBank and GFF format.
   
If you make use of this software in your research, please cite as:
   
__Artemis: sequence visualization and annotation.__
Rutherford K, Parkhill J, Crook J, Horsnell T, Rice P, Rajandream MA and Barrell B
Bioinformatics (Oxford, England) 2000;16;10;944-5 PUBMED: [11120685](http://europepmc.org/abstract/MED/11120685)

## Software Availability

Artemis is available under GPL3. The source code can be found on [GitHub](https://github.com/sanger-pathogens/Artemis).
  
The latest release of Artemis can be downloaded from our ftp site:
  
* [UNIX](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/artemis.tar.gz)
* [MacOS](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/artemis.dmg.gz)
* [Windows](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/artemis.jar)
  
For older versions of the Artemis software please see the [Artemis FTP site](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/)
  
For further instructions please see the [Artemis manual](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/artemis.pdf) and our [GitHub page](https://github.com/sanger-pathogens/Artemis/).
   
__Notes for Download__
   
__Note for all operating systems__: The Artemis tool suite requires Java version 1.8 to run. This must be installed first.
   
__Note for MacOSX__: occasionally a browser decides to display the contents of the .dmg.gz archive file rather than downloading it. If this happens hold down the <control> key and click on the download link. A popup menu should appear, containing several choices. One of the choices should be something like "Save Link As" (or perhaps "Download Link...", "Save Link to Desktop", or a variation on this theme). Select that option, and the archive file should be download correctly.
  
__Note for Windows__: you may need to shift-click on this link to download the file. Also please ensure that after downloading this file is called "act.jar" with no other hidden extensions. Some windows systems change the name so it looks like a zip file, if this happens then do not unzip it just rename it to "act.jar".

## Installation 
For installation instructions please see the [Artemis manual](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/artemis.pdf).

## Usage
For instructions please see the [Artemis manual](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/artemis.pdf)

## Contact
For issues encountered with installing the software please contact your local system administrator. For all other issues, please report them to our [Github issues page](https://github.com/sanger-pathogens/Artemis/issues) or email artemis@sanger.ac.uk.

## FAQ
### Why does Artemis run out of memory on UNIX or GNU/Linux even though the machine has lots of memory?

The Java Virtual Machine (JVM) on UNIX has a fixed upper limit on the amount of memory that is available for applications, but this limit can be changed at runtime. As shipped Artemis will use a maximum of 2GB of memory.

You will need to change the art script to fix this problem. Find the line that reads: FLAGS="-mx2g -ms100m -noverify" and change the 2g (2 gigabytes) to a bigger number dependent on your machine memory (try 3g).

To change the memory allocated to Artemis on MacOSX, set the value in the file Info.plist in the directory Artemis.app/Contents. Towards the bottom of the file there are a couple of lines that look like this:

```
<key>VMOptions</key> <string>-Xmx800m</string>
```
Changing the value after -Xmx will change the memory used by Artemis.

### Why does Artemis run out of memory on Windows even though the machine has lots of memory?

Normally the Java virtual machine artificially limits the amount of memory that Artemis can use. The fix is as follows:

Create a shortcut to the Artemis JAR file. Edit the properties of the shortcut and add java -mx250m -jar to the start of the Target: field. -mx250m sets the maximum memory Java will allocate to Artemis. We recommend choosing a number that is about 50 megabytes less than the total amount of memory in the machine (to allow for the overhead of windows and the Java virtual machine).

You will need to use the shortcut to run Artemis from then on.

## Publications
__Artemis: an integrated platform for visualization and analysis of high-throughput sequence-based experimental data.__
Carver T, Harris SR, Berriman M, Parkhill J and McQuillan JA
Bioinformatics (Oxford, England) 2012;28;4;464-9
PUBMED: [22199388](http://ukpmc.ac.uk/abstract/MED/22199388); PMC: [3278759](http://ukpmc.ac.uk/articles/PMC3278759); DOI: [10.1093/bioinformatics/btr703](http://dx.doi.org/10.1093/bioinformatics/btr703)

__Artemis and ACT: viewing, annotating and comparing sequences stored in a relational database.__
Carver T, Berriman M, Tivey A, Patel C, Böhme U et al.
Bioinformatics (Oxford, England) 2008;24;23;2672-6
PUBMED: [18845581](http://ukpmc.ac.uk/abstract/MED/18845581); PMC: [2606163](http://ukpmc.ac.uk/articles/PMC2606163); DOI: [10.1093/bioinformatics/btn529](http://dx.doi.org/10.1093/bioinformatics/btn529)

__Viewing and annotating sequence data with Artemis.__
Berriman M and Rutherford K
Briefings in bioinformatics 2003;4;2;124-32
PUBMED: [12846394](http://ukpmc.ac.uk/abstract/MED/12846394)

__Artemis: sequence visualization and annotation.__
Rutherford K, Parkhill J, Crook J, Horsnell T, Rice P et al.
Bioinformatics (Oxford, England) 2000;16;10;944-5
PUBMED: [11120685](http://ukpmc.ac.uk/abstract/MED/11120685)