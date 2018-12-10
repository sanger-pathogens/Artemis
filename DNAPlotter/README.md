# DNAPlotter
DNAPlotter can be used to generate images of circular and linear DNA maps to display regions and features of interest. The images can be inserted into a document or printed out directly. As this uses [Artemis](/Artemis/Artemis/) it can read in the common file formats EMBL, GenBank and GFF3.

## Citation
> DNAPlotter is an interactive Java application for generating circular and linear representations of genomes. Making use of the Artemis libraries to provide a user-friendly method of loading in sequence files (EMBL, GenBank, GFF) as well as data from relational databases, it filters features of interest to display on separate user-definable tracks. It can be used to produce publication quality images for papers or web pages.

If you make use of this software in your research, please cite as:

__DNAPlotter: circular and linear interactive genome visualization.__  
Carver T, Thomson N, Bleasby A, Berriman M and Parkhill J  
_Bioinformatics (Oxford, England) 2009;25;1;119-20_  
PUBMED: [18990721](http://ukpmc.ac.uk/abstract/MED/18990721); PMC: [2612626](http://ukpmc.ac.uk/articles/PMC2612626); DOI: [10.1093/bioinformatics/btn578](http://dx.doi.org/10.1093/bioinformatics/btn578)

## Software Availability
DNAPlotter is packaged as part of the Artemis Software. The Artemis Software is available under GPL3. The source code can be found on [GitHub](https://github.com/sanger-pathogens/Artemis).

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
This will create a directory called ~/artemis which will contain all the files necessary for running DNAPlotter and the other Artemis tools.

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
### Running DNAPlotter on UNIX/Linux Systems
The easiest way to run the program is to run the script called dnaplotter in the Artemis installation directory, like this:
```
artemis/dnaplotter
```
Alternatively you can start DNAPlotter with the name of a template file, eg:
```
artemis/dnaplotter -t templatefile
```
For other available options use:
```
artemis/dnaplotter -help
```
### Running DNAPlotter on Macintosh Systems
On MacOSX machines, DNAPlotter can be started by double clicking on the DNAPlotter icon.

### Running DNAPlotter on Windows Systems
DNAPlotter can be started by double clicking on the dnaplotter.jar icon.

## The User Manual
For additional information please see the [Artemis manual](https://sanger-pathogens.github.io/Artemis/Artemis/artemis-manual.html) and our [GitHub page](https://github.com/sanger-pathogens/Artemis/). A PDF version of the manual is also available for download [here](https://sanger-pathogens.github.io/Artemis/Artemis/artemis-manual.pdf).

## Contact
For issues encountered with installing the software please contact your local system administrator. For all other issues, please report them to our [Github issues page](https://github.com/sanger-pathogens/Artemis/issues) or email <artemis-help@sanger.ac.uk>.

## FAQ
### Why does DNAPlotter run out of memory on UNIX or GNU/Linux even though the machine has lots of memory?

The Java Virtual Machine (JVM) on UNIX has a fixed upper limit on the amount of memory that is available for applications, but this limit can be changed at runtime. As shipped DNAPlotter will use a maximum of 2GB of memory.

There are two ways of fixing this problem:
1. Change the dnaplotter script. Find the line that reads: FLAGS="-mx2g -ms100m -noverify" and change the 2g (2 gigabytes) to a bigger number dependent on your machine memory (try 3g), or
2. Create an ARTEMIS_JVM_FLAGS environment variable set to "-mx2g -ms100m", adjusting the mx value as required. No script change is required for this, but it would need to be added to your environment.

### Why does DNAPlotter run out of memory on MacOSX even though the machine has lots of memory?
To change the memory allocated to DNAPlotter on MacOSX, set the value in the file Info.plist in the directory DNAPlotter.app/Contents. Towards the bottom of the file you will see these lines:
```
<key>JVMOptions</key>
<array>
<string>-mx2g</string>
```
Changing the value after -mx will change the max memory used by DNAPlotter. The default is 2Gb.

### Why does DNAPlotter run out of memory on Windows even though the machine has lots of memory?

Normally the Java virtual machine artificially limits the amount of memory that DNAPlotter can use. The fix is as follows:

Create a shortcut to the dnaplotter.jar JAR file. Edit the properties of the shortcut and add java -mx2g -jar to the start of the Target: field. -mx2g sets the maximum memory Java will allocate to DNAPlotter (2 gigabytes in this case). We recommend choosing a number that is about 50 megabytes less than the total amount of memory in the machine (to allow for the overhead of windows and the Java virtual machine).

You will need to use the shortcut to run DNAPlotter from then on.
