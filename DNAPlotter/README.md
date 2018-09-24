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

__Note__: or additional information please see the [Artemis manual](https://kpepper.github.io/Artemis/Artemis/artemis-manual.html) and our [GitHub page](https://github.com/sanger-pathogens/Artemis/). A PDF version of the manual is also available for download [here](https://kpepper.github.io/Artemis/Artemis/artemis-manual.pdf).

## Contact
For issues encountered with installing the software please contact your local system administrator. For all other issues, please report them to our [Github issues page](https://github.com/sanger-pathogens/Artemis/issues) or email <artemis-help@sanger.ac.uk>.

## FAQ
### Why does DNAPlotter run out of memory on UNIX or GNU/Linux even though the machine has lots of memory?

The Java Virtual Machine (JVM) on UNIX has a fixed upper limit on the amount of memory that is available for applications, but this limit can be changed at runtime. As shipped DNAPlotter will use a maximum of 2GB of memory.

There are two ways of fixing this problem:
1. Change the dnaplotter script. Find the line that reads: FLAGS="-mx2g -ms100m -noverify" and change the 2g (2 gigabytes) to a bigger number dependent on your machine memory (try 3g), or
2. Create an ARTEMIS_JVM_FLAGS environment variable set to "-mx2g -ms100m", adjusting the mx value as required. No script change is required for this, but it would need to be added to your environment.

### Why does DNAPlotter run out of memory on MacOSX even though the machine has lots of memory?
To change the memory allocated to DNAPlotter on MacOSX, set the value in the file DNAPlotter.cfg in the directory DNAPlotter.app/Contents/Java. There are a couple of lines that look like this:

```
[JVMOptions]
-Xmx2g
```
Changing the value after -Xmx will change the memory used by DNAPlotter.

### Why does DNAPlotter run out of memory on Windows even though the machine has lots of memory?

Normally the Java virtual machine artificially limits the amount of memory that DNAPlotter can use. The fix is as follows:

Create a shortcut to the dnaplotter.jar JAR file. Edit the properties of the shortcut and add java -mx2g -jar to the start of the Target: field. -mx2g sets the maximum memory Java will allocate to DNAPlotter (2 gigabytes in this case). We recommend choosing a number that is about 50 megabytes less than the total amount of memory in the machine (to allow for the overhead of windows and the Java virtual machine).

You will need to use the shortcut to run DNAPlotter from then on.
