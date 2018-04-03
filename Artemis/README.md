# Artemis
Artemis is a free genome browser and annotation tool that allows visualisation of sequence features, next generation data and the results of analyses within the context of the sequence, and also its six-frame translation.
Artemis is written in Java, and is available for UNIX, Macintosh and Windows systems. It can read EMBL and GENBANK database entries or sequence in FASTA, indexed FASTA or raw format. Other sequence features can be in EMBL, GENBANK or GFF format.
Full information about the latest release of Artemis can be found in the [Artemis manual](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/artemis.pdf) and the [current release notes](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/release_notes.txt).

## Publication

> High-throughput sequencing (HTS) technologies have made low-cost sequencing of large numbers of samples commonplace. An explosion in the type, not just number, of sequencing experiments has also taken place including genome re-sequencing, population-scale variation detection, whole transcriptome sequencing and genome-wide analysis of protein-bound nucleic acids.We present Artemis as a tool for integrated visualization and computational analysis of different types of HTS datasets in the context of a reference genome and its corresponding annotation.

If you make use of this software in your research, please cite as:

__Artemis: an integrated platform for visualization and analysis of high-throughput sequence-based experimental data.__  
Carver T, Harris SR, Berriman M, Parkhill J and McQuillan JA  
_Bioinformatics (Oxford, England) 2012;28;4;464-9_  
PUBMED: [22199388](http://ukpmc.ac.uk/abstract/MED/22199388); PMC: [3278759](http://ukpmc.ac.uk/articles/PMC3278759); DOI: [10.1093/bioinformatics/btr703](http://dx.doi.org/10.1093/bioinformatics/btr703)

__Artemis and ACT: viewing, annotating and comparing sequences stored in a relational database.__  
Carver T, Berriman M, Tivey A, Patel C, Böhme U et al.  
_Bioinformatics (Oxford, England) 2008;24;23;2672-6_  
PUBMED: [18845581](http://ukpmc.ac.uk/abstract/MED/18845581); PMC: [2606163](http://ukpmc.ac.uk/articles/PMC2606163); DOI: [10.1093/bioinformatics/btn529](http://dx.doi.org/10.1093/bioinformatics/btn529)

__Viewing and annotating sequence data with Artemis.__  
Berriman M and Rutherford K  
_Briefings in bioinformatics 2003;4;2;124-32_  
PUBMED: [12846394](http://ukpmc.ac.uk/abstract/MED/12846394)   

__Artemis: sequence visualization and annotation.__  
_Rutherford K, Parkhill J, Crook J, Horsnell T, Rice P, Rajandream MA and Barrell B_  
Bioinformatics (Oxford, England) 2000;16;10;944-5 PUBMED: [11120685](http://europepmc.org/abstract/MED/11120685)

## Software Availability
The Artemis Software is available under GPL3. The source code can be found on [GitHub](https://github.com/sanger-pathogens/Artemis).

The latest release of Artemis can be downloaded from our ftp site:

* [UNIX](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/artemis.tar.gz)
* [MacOS](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/artemis.dmg.gz)
* [Windows](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/artemis-windows.zip)

For older versions of the Artemis software please see the [Artemis FTP site](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/)


### Notes for Download

__Note for UNIX/Windows operating systems__: The Artemis tool suite requires Java version 1.8 to run. This must be installed first.

__Note for MacOSX__: occasionally a browser decides to display the contents of the .dmg.gz archive file rather than downloading it. If this happens hold down the <control> key and click on the download link. A popup menu should appear, containing several choices. One of the choices should be something like "Save Link As" (or perhaps "Download Link...", "Save Link to Desktop", or a variation on this theme). Select that option, and the archive file should be download correctly.

## Installation
### For UNIX/Linux
Copy the artemis.tar.gz file to the directory that you wish to install to. Change to that directory and uncompress and untar the file. On UNIX, the command is:
```
tar zxf artemis.tar.gz
```
### For MacOSX
For MacOSX users an archive artemis.dmg.gz disk image is provided. Click on this file in your Downloads folder. This should automatically unzip and mount the disk image, opening a Finder window that contains the Artemis, ACT, BamView and DNAPlotter applications. These apps can then be dragged to any desired location, for example, your desktop.

### For Windows
Copy the artemis-windows.zip file to the directory that you wish to install to and then unzip using an application such as WinZip.
This should unpack the artemis.jar, act.jar, bamview.jar and dnaplotter.jar application files.

## Usage
### Running Artemis on UNIX/Linux Systems
The easiest way to run the program is to run the script called art in the Artemis installation directory, like this:
```
artemis/art
```
Alternatively you can start Artemis with the name of a sequence file or embl file, eg:
```
artemis/art artemis/etc/c1215.embl
```
For other available options use:
```
artemis/art -help
```
### Running Artemis on Macintosh Systems
On MacOSX machines, Artemis can be started by double clicking on the Artemis icon.

### Running Artemis on Windows Systems
Artemis can be started by double clicking on the artemis.jar icon.

__Note__: For additional information please see the [Artemis manual](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/artemis.pdf) and our [GitHub page](https://github.com/sanger-pathogens/Artemis/).

## Contact
For issues encountered with installing the software please contact your local system administrator. For all other issues, please report them to our [Github issues page](https://github.com/sanger-pathogens/Artemis/issues) or email <artemis-help@sanger.ac.uk>.

## FAQ
### Why does Artemis run out of memory on UNIX or GNU/Linux even though the machine has lots of memory?

The Java Virtual Machine (JVM) on UNIX has a fixed upper limit on the amount of memory that is available for applications, but this limit can be changed at runtime. As shipped Artemis will use a maximum of 2GB of memory.

There are two ways of fixing this problem:
1. Change the art script. Find the line that reads: FLAGS="-mx2g -ms100m -noverify" and change the 2g (2 gigabytes) to a bigger number dependent on your machine memory (try 3g), or
2. Create an ARTEMIS_JVM_FLAGS environment variable set to "-mx2g -ms100m", adjusting the mx value as required. No script change is required for this, but it would need to be added to your environment.

### Why does Artemis run out of memory on MacOSX even though the machine has lots of memory?
To change the memory allocated to Artemis on MacOSX, set the value in the file Artemis.cfg in the directory Artemis.app/Contents/Java. There are a couple of lines that look like this:

```
[JVMOptions]
-Xmx2g
```
Changing the value after -Xmx will change the memory used by Artemis.

### Why does Artemis run out of memory on Windows even though the machine has lots of memory?

Normally the Java virtual machine artificially limits the amount of memory that Artemis can use. The fix is as follows:

Create a shortcut to the Artemis JAR file. Edit the properties of the shortcut and add java -mx2g -jar to the start of the Target: field. -mx2g sets the maximum memory Java will allocate to Artemis (2 gigabytes in this case). We recommend choosing a number that is about 50 megabytes less than the total amount of memory in the machine (to allow for the overhead of windows and the Java virtual machine).

You will need to use the shortcut to run Artemis from then on.
