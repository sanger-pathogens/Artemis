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

__Note__: For additional information on viewing alignment files please see the [Artemis manual](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/artemis.pdf) and our [GitHub page](https://github.com/sanger-pathogens/Artemis/).

## Contact
For issues encountered with installing the software please contact your local system administrator. For all other issues, please report them to our [Github issues page](https://github.com/sanger-pathogens/Artemis/issues) or email <artemis-help@sanger.ac.uk>.

## FAQ
### CRAM Reference Lookup
BamView will attempt to download CRAM reference sequences from EBI if none is specified on startup or via CRAM headers. This behaviour can be further governed using two environment variables REF_PATH and REF_CACHE. The usage of these variables is the same as for samtools. Please refer to [this article](http://www.htslib.org/workflow/) for details.

### Troubleshooting
The alignment file needs to be indexed and sorted. Check that you have the index file in the same directory. BamView assumes that the name of the index file is the same as the alignment file but with the added .bai [BAM] or .crai [CRAM] extension.
If the BamView window is actually opening but just blank then try changing to a different view by right clicking on the BamView window. You may need to select the correct reference from the top left drop down list.

### Memory usage
If you need to increase the maximum memory available to the application then follow the same procedure detailed in the Artemis page FAQ, for the BamView application files.
