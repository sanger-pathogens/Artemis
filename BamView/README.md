# BamView
BamView is an interactive Java application for visualising read-alignment data stored in BAM or CRAM files.

## Publication
> So-called next-generation sequencing (NGS) has provided the ability to sequence on a massive scale at low cost, enabling biologists to perform powerful experiments and gain insight into biological processes. BamView has been developed to visualize and analyse sequence reads from NGS platforms, which have been aligned to a reference sequence. It is a desktop application for browsing the aligned or mapped reads [Ruffalo, M, LaFramboise, T, Koyutürk, M. Comparative analysis of algorithms for next-generation sequencing read alignment. Bioinformatics 2011;27:2790-6] at different levels of magnification, from nucleotide level, where the base qualities can be seen, to genome or chromosome level where overall coverage is shown. To enable in-depth investigation of NGS data, various views are provided that can be configured to highlight interesting aspects of the data. Multiple read alignment files can be overlaid to compare results from different experiments, and filters can be applied to facilitate the interpretation of the aligned reads. As well as being a standalone application it can be used as an integrated part of the Artemis genome browser, BamView allows the user to study NGS data in the context of the sequence and annotation of the reference genome. Single nucleotide polymorphism (SNP) density and candidate SNP sites can be highlighted and investigated, and read-pair information can be used to discover large structural insertions and deletions. The application will also calculate simple analyses of the read mapping, including reporting the read counts and reads per kilobase per million mapped reads (RPKM) for genes selected by the user.

If you make use of this software in your research, please cite as:

__BamView: visualizing and interpretation of next-generation sequencing read alignments.__   
Carver T, Harris SR, Otto TD, Berriman M, Parkhill J and McQuillan JA   
_Briefings in bioinformatics 2013;14;2;203-12_   
PUBMED: [22253280](http://ukpmc.ac.uk/abstract/MED/22253280); PMC: [3603209](http://ukpmc.ac.uk/articles/PMC3603209); DOI: [10.1093/bib/bbr073](http://dx.doi.org/10.1093/bib/bbr073)

__BamView: viewing mapped read alignment data in the context of the reference sequence.__  
Carver T, Böhme U, Otto TD, Parkhill J and Berriman M  
_Bioinformatics (Oxford, England) 2010;26;5;676-7_  
PUBMED: [20071372](http://ukpmc.ac.uk/abstract/MED/20071372); PMC: [2828118](http://ukpmc.ac.uk/articles/PMC2828118); DOI: [10.1093/bioinformatics/btq010](http://dx.doi.org/10.1093/bioinformatics/btq010)

## Software Availability
Bamview is packaged as part of the Artemis Software. The Artemis Software is available under GPL3. The source code can be found on [GitHub](https://github.com/sanger-pathogens/Artemis).

The latest release of Artemis can be downloaded from our ftp site:

* [UNIX](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/artemis.tar.gz)
* [MacOS](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/artemis.dmg.gz)
* [MacOS (CHADO)](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/artemis.dmg.gz)
* [Windows](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/artemis-windows.zip)

For older versions of the software please see the [Artemis FTP site](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/)

The previous v17.0.1 version of the Artemis software required Java version 1.8 to run. All subsequent versions from v18.0.0 onwards require a minimum of Java 9 and preferably Java 11. This must be installed first.

### Notes for Download

__Note for MacOSX__: We recommend using Firefox or Chrome to download from the links above.

## Installation
### For UNIX/Linux
Change directory to the directory you wish to install the Artemis software in. We will use ~/ in this example and in the next chapter.

Uncompress and untar the artemis-unix-release-{version}.tar.gz or artemis-unix-release-{version}.zip file. On UNIX the command is:
```
tar zxf artemis-unix-release-{version}.tar.gz
```
This will create a directory called ~/artemis which will contain all the files necessary for running Bamview and the other Artemis tools.

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
### Running BamView on UNIX/Linux Systems
The easiest way to run the program is to run the script called bamview in the Artemis installation directory, like this:
```
artemis/bamview
```
Alternatively you can start BamView with the name of a BAM or CRAM, eg:
```
artemis/bamview -a bam_or_cram_file
artemis/bamview -a 'http://127.0.0.1/bam_or_cram_file'
artemis/bamview -a 'ftp://ftp.myftpsite/bam_or_cram_file'
artemis/bamview -a cram_file -r reference_file
```
For other available options use:
```
artemis/bamview -help
```
### Running BamView on Macintosh Systems
On MacOSX machines, BamView can be started by double clicking on the BamView icon.

### Running BamView on Windows Systems
BamView can be started by double clicking on the bamview.jar icon.

## The User Manual
For additional information on installation and viewing of alignment files please see the [Artemis manual](https://sanger-pathogens.github.io/Artemis/Artemis/artemis-manual.html) and our [GitHub page](https://github.com/sanger-pathogens/Artemis/). A PDF version of the manual is also available for download [here](https://sanger-pathogens.github.io/Artemis/Artemis/artemis-manual.pdf).

## Contact
For issues encountered with installing the software please contact your local system administrator. For all other issues, please report them to our [Github issues page](https://github.com/sanger-pathogens/Artemis/issues) or email <artemis-help@sanger.ac.uk>.

## FAQ
### CRAM Reference Lookup
BamView will attempt to download CRAM reference sequences from EBI if none is specified on startup or via CRAM headers. This behaviour can be further governed using two environment variables REF_PATH and REF_CACHE. The usage of these variables is the same as for samtools. Please refer to [this article](http://www.htslib.org/workflow/) for details.

### Troubleshooting
The alignment file needs to be sorted and indexed. Samtools may be used for this:

to sort:
samtools sort <in.bam> <out.prefix>

and then index:
samtools index <in.bam> [<out.index>]

Check that you have placed the index file in the same directory. BamView assumes that the name of the index file is the same as the alignment file but with the added .bai [BAM] or .crai [CRAM] extension.

If the BamView window is actually opening but just blank then try changing to a different view by right clicking on the BamView window. You may need to select the correct reference from the top left drop down list.

BamView in Artemis:
The asynchronous option (when selected from the menu) means that when you scroll along the sequence the BamView window only updates when scrolling stops. This makes scrolling faster when the coverage is high in the region being viewed. However this can be turned off to see the reads as you scroll along.

### Colour schemes
Stack view and Strand stack view: paired reads are blue; single reads or reads with an unmapped pair are black; duplicate reads are green.

Inferred size view: paired reads are blue and those with an inversion are red. Reads that do not have a mapped mate are black and are optionally shown in the inferred insert size view.

Paired stack view: paired reads are blue and those with an inversion are red.

Zoomed in to the nucleotide level: the bases can be coloured by their mapping quality score:
blue <10; green <20; orange <30; black ≥30.

### Why does BamView run out of memory on UNIX or GNU/Linux even though the machine has lots of memory?

The Java Virtual Machine (JVM) on UNIX has a fixed upper limit on the amount of memory that is available for applications, but this limit can be changed at runtime. As shipped BamView will use a maximum of 2GB of memory.

There are two ways of fixing this problem:
1. Change the bamview script. Find the line that reads: FLAGS="-mx2g -ms100m -noverify" and change the 2g (2 gigabytes) to a bigger number dependent on your machine memory (try 3g), or
2. Create an ARTEMIS_JVM_FLAGS environment variable set to "-mx2g -ms100m", adjusting the mx value as required. No script change is required for this, but it would need to be added to your environment.

### Why does BamView run out of memory on MacOSX even though the machine has lots of memory?
To change the memory allocated to BamView on MacOSX, set the value in the file Info.plist in the directory BamView.app/Contents. Towards the bottom of the file you will see these lines:
```
<key>JVMOptions</key>
<array>
<string>-mx2g</string>
```
Changing the value after -mx will change the max memory used by BamView. The default is 2Gb.

### Why does BamView run out of memory on Windows even though the machine has lots of memory?

Normally the Java virtual machine artificially limits the amount of memory that BamView can use. The fix is as follows:

Create a shortcut to the bamview.jar JAR file. Edit the properties of the shortcut and add java -mx2g -jar to the start of the Target: field. -mx2g sets the maximum memory Java will allocate to BamView (2 gigabytes in this case). We recommend choosing a number that is about 50 megabytes less than the total amount of memory in the machine (to allow for the overhead of windows and the Java virtual machine).

You will need to use the shortcut to run BamView from then on.
