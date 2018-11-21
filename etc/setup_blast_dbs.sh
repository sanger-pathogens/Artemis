#!/bin/sh

###############################################################
#                                                             #
# Download uniprot_sprot + uniprot_trembl to create uniprot   #
# and format for blast.                                       #
#                                                             #
###############################################################

BLAST=`which makeblastdb`
if [ ! -x "$BLAST" ]
then
	BLAST=/usr/local/ncbi/blast/bin/makeblastdb
fi

if [ ! -x "$BLAST" ]
then
	echo "ERROR: Unable to locate makeblastdb program. Exiting."
	echo "Please ensure that this is available from /usr/local/ncbi/blast/bin"
	echo "OR has been added to your PATH environment variable." 
	exit 1
fi

export BLAST

echo "Installing blast databases to $HOME/blast-data..."

INSTALL_AREA=$HOME/blast-data
cd $HOME
mkdir -p $INSTALL_AREA
cd $INSTALL_AREA

# Uncomment the following section to create a FASTLIBS file.
# 
#CWD=`pwd`
#
#cat > $CWD/pubseqgbs <<EOF
#UNIPROT\$0U@$CWD/uniprot.nam 12
#UNIPROT_EUKARYOTA\$0E@$CWD/uniprot_eukaryota.nam 12
#UNIPROT_BACTERIA\$0B@$CWD/uniprot_bacteria.nam 12
#FALC_CHAB_BERB_YOEL_PROTEINS\$0F@$CWD/falc_chab_berb_yoel_proteins.nam 12
#EOF
#
#cat > $CWD/uniprot.nam <<EOF
#<$CWD
#uniprot 0
#EOF
#
#cat > $CWD/uniprot_bacteria.nam <<EOF
#<$CWD
#uniprot_bacteria 0
#EOF
#
#cat > $CWD/uniprot_eukaryota.nam <<EOF
#<$CWD
#uniprot_eukaryota 0
#EOF
#
#cat > $CWD/falc_chab_berb_yoel_proteins.nam <<EOF
#<$CWD
#falc_chab_berb_yoel_proteins 0
#EOF
#

echo "Downloading and formatting databases..."

if test ! -f uniprot; then
  curl -O ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz
  gunzip -f uniprot_sprot.fasta.gz

  curl -O ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.fasta.gz
  gunzip -f uniprot_trembl.fasta.gz

  cat uniprot_sprot.fasta uniprot_trembl.fasta > uniprot
  rm -f uniprot_sprot.fasta uniprot_trembl.fasta

  $BLAST -in uniprot -dbtype prot -parse_seqids -title uniprot &
fi

# format uniprot_bacteria if we have it
if test -f uniprot_bacteria.gz; then
  gunzip -f uniprot_bacteria.gz
  $BLAST -in uniprot_bacteria -dbtype prot -parse_seqids -title uniprot_bacteria &
fi

# format uniprot_eukaryota if we have it
if test -f uniprot_eukaryota.gz; then
  gunzip -f uniprot_eukaryota.gz
  $BLAST -in uniprot_eukaryota -dbtype prot -parse_seqids -title uniprot_eukaryota &
fi

# format falc_chab_berb_yoel_proteins (aeb) if we have it
#if test -f falc_chab_berb_yoel_proteins.gz; then
#  gunzip -f falc_chab_berb_yoel_proteins.gz
#  $BLAST -i falc_chab_berb_yoel_proteins -o T -p T -t falc_chab_berb_yoel_proteins &
#fi

wait

echo "Completed"

exit 0

