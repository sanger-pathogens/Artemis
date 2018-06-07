
Artemis, ACT & DNAPlotter
=========================

Artemis is a free genome viewer and annotation tool that allows visualization of sequence features and the results of analyses within the context of the sequence, and its six-frame translation. 

ACT (Artemis Comparison Tool) is a DNA sequence comparison viewer based on Artemis. 

DNAPlotter generates images of circular and linear DNA maps.

BamView is a standalone BAM/CRAM file viewer.

Installation on MacOSX
======================

1. Drag and drop the Artemis, ACT, DNAPlotter and BamView icons to an Applications directory, for example an 'Applications' directory in your home directory. 

2. Then eject the disk image (right click and select 'Eject'). 

3. Drag the icons from where you have installed them, into your dock.

4. To run 'art' and 'act' from the command line add the following to your UNIX PATH:

$ARTEMIS/Artemis.app/Contents

where $ARTEMIS is the directory you have installed to. This can be done by adding a line to ~/.bashrc (if you are using bash shell) for example:

PATH=~/Applications/Artemis.app/Contents; export PATH

Searching and Using Local Sequence Databases (Optional)
=======================================================

1. Download and install BLAST:

ftp://ftp.ncbi.nih.gov/blast/executables/release/

get the latest release, e.g.:

ftp://ftp.ncbi.nih.gov/blast/executables/release/2.2.18/blast-2.2.18-universal-macosx.tar.gz

Untar and add the contents of the bin directory into Artemis.app/Contents/blast.

2. Download and install Fasta
ftp://ftp.ebi.ac.uk/pub/software/unix/fasta

Get the latest, e.g.:

ftp://ftp.ebi.ac.uk/pub/software/unix/fasta/fasta3/fasta-34.26.5.shar.Z

Put this in Artemis.app/Contents/fasta and unwrap and compile:
uncompress fasta-34.26.5.shar.Z
cat fasta-34.26.5.shar | sh
make -f Makefile.os_x all

3. Download and format uniprot database

You can do this by using the setup_dbs.sh script in Artemis.app/Contents/blast-data

cd Artemis.app/Contents/blast-data
./setup_dbs.sh

(Note the environment variable http_proxy may need to be set for the download to work).

Additional databases can be added here.
To select to run searches locally ensure that under the 'Options' menu 'Send Searches via SSH' is not toggled.

BamView CRAM reference lookup
=============================

BamView now has the ability to download CRAM reference sequences if none is specified on startup.
This behaviour can be further governed using two environment variables REF_PATH and REF_CACHE.
The usage of these variables is the same as for samtools. Please refer to this article for details:

http://www.htslib.org/workflow/


References:
===========

"K. Rutherford, J. Parkhill, J. Crook, T. Horsnell, P. Rice, M-A. Rajandream and B. Barrell (2000) Artemis: sequence visualisation and annotation." Bioinformatics 16 (10) 944-945.

T.J. Carver, K.M. Rutherford, M. Berriman, M-A. Rajandream, B.G. Barrell and J. Parkhill, "ACT: the Artemis Comparison Tool", Bioinformatics, 2005.

