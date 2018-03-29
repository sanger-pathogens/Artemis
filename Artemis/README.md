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

## Installation 
For installation instructions please see the [Artemis manual](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/artemis.pdf).

## Usage
For instructions please see the [Artemis manual](ftp://ftp.sanger.ac.uk/pub/resources/software/artemis/artemis.pdf)

## Contact
For issues encountered with installing the software please contact your local system administrator. For all other issues, please report them to our [Github issues page](https://github.com/sanger-pathogens/Artemis/issues) or email <artemis-help@sanger.ac.uk>.

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
