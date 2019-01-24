# Artemis Comparison Tool (ACT)
ACT is a Java application for displaying pairwise comparisons between two or more DNA sequences. It can be used to identify and analyse regions of similarity and difference between genomes and to explore conservation of synteny, in the context of the entire sequences and their annotation. It can read complete EMBL, GENBANK and GFF entries or sequences in FASTA or raw format.

### Related Software
Two separate ACT-related web applications have been developed.

* [DoubleACT](http://www.hpa-bioinfotools.org.uk/pise/double_actv2.html) was developed by Anthony Underwood and Jonathan Green at the [Public Public Health England](http://www.hpa.org.uk/) and allows you to paste or upload sequences to generate ACT comparison files.
* [WebACT](http://www.webact.org/WebACT/home) was developed at Imperial College, and allows users to generate comparisons from their sequence files or use a set of pre-computed comparisons. All comparison files can be downloaded for local use or launched in ACT. As of Dec 2018 the web site was under redevelopment.

There are also several other types of comparison file that ACT can read, e.g. tblastx hit table files from the NCBI Blast web site.

## Publication
> The Artemis Comparison Tool (ACT) allows an interactive visualisation of comparisons between complete genome sequences and associated annotations. The comparison data can be generated with several different programs; BLASTN, TBLASTX or Mummer comparisons between genomic DNA sequences, or orthologue tables generated by reciprocal FASTA comparison between protein sets. It is possible to identify regions of similarity, insertions and rearrangements at any level from the whole genome to base-pair differences. ACT uses Artemis components to display the sequences and so inherits powerful searching and analysis tools. ACT is part of the Artemis distribution and is similarly open source, written in Java and can run on any Java enabled platform, including UNIX, Macintosh and Windows.

If you make use of this software in your research, please cite as:

__ACT: the Artemis Comparison Tool.__  
Carver TJ, Rutherford KM, Berriman M, Rajandream MA, Barrell BG and Parkhill  
_Bioinformatics (Oxford, England) 2005;21;16;3422-3_   
PUBMED: [15976072](http://ukpmc.ac.uk/abstract/MED/15976072); DOI: [10.1093/bioinformatics/bti553](http://dx.doi.org/10.1093/bioinformatics/bti553)

## Software Availability
ACT is packaged as part of the Artemis Software. The Artemis Software is available under GPL3. The source code can be found on [GitHub](https://github.com/sanger-pathogens/Artemis).

The latest release of Artemis can be downloaded by clicking on the relevant link below:

* [UNIX](https://github.com/sanger-pathogens/Artemis/releases/download/v18.0.1/artemis-unix-release-18.0.1.tar.gz)
* [MacOS](https://github.com/sanger-pathogens/Artemis/releases/download/v18.0.1/artemis-macosx-release-18.0.1.dmg.gz)
* [MacOS (CHADO) - for out of the box CHADO database connectivity](https://github.com/sanger-pathogens/Artemis/releases/download/v18.0.1/artemis-macosx-chado-release-18.0.1.dmg.gz)
* [Windows](https://github.com/sanger-pathogens/Artemis/releases/download/v18.0.1/artemis-windows-release-18.0.1.zip)

Or via [Bioconda](https://bioconda.github.io). Simply use:
```
conda install artemis -c bioconda
```
from the command line, and you're ready to go.

For older versions of the software please see the [Artemis FTP site](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/)

The previous v17.0.1 version of the Artemis software required Java version 1.8 to run. All subsequent versions from v18.0.0 onwards require a minimum of Java 9 and preferably Java 11. This must be installed first. A Java installation is not required for the Bioconda route.

## Installation
### For UNIX/Linux
Change directory to the directory you wish to install the Artemis software in. We will use ~/ in this example and in the next chapter.

Uncompress and untar the artemis-unix-release-{version}.tar.gz or artemis-unix-release-{version}.zip file. On UNIX the command is:
```
tar zxf artemis-unix-release-{version}.tar.gz
```
This will create a directory called ~/artemis which will contain all the files necessary for running ACT and the other Artemis tools.

### For MacOSX
For MacOSX users, an artemis-macosx-release-{version}.dmg.gz disk image is provided. Double-click on this file in your Downloads folder to unzip it. Then double-click the unzipped artemis-macosx-release-{version}.dmg to mount the "Artemis_Tools" image and display its contents - the Artemis, ACT, BamView and DNAPlotter applications. These apps can then be dragged to any desired location, for example, your dock or desktop. The download file can be unzipped from the command line using gunzip, if necessary:
```
gunzip artemis-macosx-release-{version}.dmg.gz
```
There’s also an artemis-macosx-chado-release-{version}.dmg disk image that will start up Artemis with a Chado connection window displayed, if you wish to work connected to a Chado database in Artemis or ACT. This is installed in exactly the same way.

### For Windows
Copy the artemis-windows-release-{version}.zip file to the directory that you wish to install to and then unzip using an application such as WinZip.
This should unpack the artemis.jar, act.jar, bamview.jar and dnaplotter.jar application files.

## Usage
### Running ACT on UNIX/Linux Systems
The easiest way to run the program is to run the script called act in the Artemis installation directory, like this:
```
artemis/act
```
Alternatively you can start ACT by passing the comparison filenames as arguments, eg:
```
artemis/act sequence1 blast_output.crunch sequence2
```
For other available options use:
```
artemis/act -help
```
### Running ACT on Macintosh Systems
On MacOSX machines, ACT can be started by double clicking on the ACT icon.

You may find that when trying to run the ACT app for the first time, that you get a security warning window displayed stating that "the application cannot be opened because it is from an unidentified developer" (because the apps are not obtained from the app store). If that's the case, then just okay the window. Go into your System Preferences via the apple symbol at top left of screen. Then select "Security and Privacy". You should then see a button called "Open Anyway", next to some text saying "ACT was blocked from opening because it is not from an identified developer". Click on the "Open Anyway" button, which will then display the security warning window again - click the "Open" button on it and ACT should then start. The application will then open straight away after this, without any further security warnings.

### Running ACT on Windows Systems
ACT can be started by double clicking on the act.jar icon.

### Running ACT via Bioconda
The ACT start script is available in the path, so in a terminal window just use:
```
act
```

## The User Manual
For additional information please see the [ACT manual](https://sanger-pathogens.github.io/Artemis/ACT/act-manual.html). A PDF version of the manual is also available for download [here](https://sanger-pathogens.github.io/Artemis/ACT/act-manual.pdf).

## Contact
For issues encountered with installing the software please contact your local system administrator. For all other issues, please report them to our [Github issues page](https://github.com/sanger-pathogens/Artemis/issues) or email <artemis-help@sanger.ac.uk>.

## FAQ
### Why does ACT run out of memory on UNIX or GNU/Linux even though the machine has lots of memory?

The Java Virtual Machine (JVM) on UNIX has a fixed upper limit on the amount of memory that is available for applications, but this limit can be changed at runtime. As shipped ACT will use a maximum of 2GB of memory.

There are two ways of fixing this problem:
1. Change the act script. Find the line that reads: FLAGS="-mx2g -ms100m -noverify" and change the 2g (2 gigabytes) to a bigger number dependent on your machine memory (try 3g), or
2. Create an ARTEMIS_JVM_FLAGS environment variable set to "-mx2g -ms100m", adjusting the mx value as required. No script change is required for this, but it would need to be added to your environment.

### Why does ACT run out of memory on MacOSX even though the machine has lots of memory?
To change the memory allocated to ACT on MacOSX, set the value in the file Info.plist in the directory ACT.app/Contents. Towards the bottom of the file you will see these lines:
```
<key>JVMOptions</key>
<array>
<string>-mx2g</string>
```
Changing the value after -mx will change the max memory used by ACT. The default is 2Gb.

### Why does ACT run out of memory on Windows even though the machine has lots of memory?

Normally the Java virtual machine artificially limits the amount of memory that ACT can use. The fix is as follows:

Create a shortcut to the act.jar JAR file. Edit the properties of the shortcut and add java -mx2g -jar to the start of the Target: field. -mx2g sets the maximum memory Java will allocate to ACT (2 gigabytes in this case). We recommend choosing a number that is about 50 megabytes less than the total amount of memory in the machine (to allow for the overhead of windows and the Java virtual machine).

You will need to use the shortcut to run ACT from then on.

## References
__WebACT: an online genome comparison suite.__   
Abbott JC, Aanensen DM and Bentley SD  
_Methods in molecular biology (Clifton, N.J.) 2007;395;57-74_  
PUBMED: [17993667](http://ukpmc.ac.uk/abstract/MED/17993667)