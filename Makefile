# This is a GNU Makefile for Artemis

# $Header: //tmp/pathsoft/artemis/Makefile,v 1.47 2009-09-21 15:32:03 tjc Exp $

SHELL=/bin/sh

#OPT_FLAGS = -g -deprecation

JAVAC := javac -source 1.5 -target 1.5 $(OPT_FLAGS) $(EXTRA_FLAGS)

REAL_CLASSPATH := CLASSPATH=lib/biojava.jar:lib/jemAlign.jar:lib/j2ssh/j2ssh-core.jar:lib/ibatis/ibatis-2.3.4.726.jar:lib/ibatis/log4j-1.2.14.jar:lib/postgresql-8.4-701.jdbc3.jar:lib/picard/picard.jar:lib/picard/sam.jar:lib/commons-net-2.2.jar:lib/batik/batik-awt-util.jar:lib/batik/batik-dom.jar:lib/batik/batik-ext.jar:lib/batik/batik-svggen.jar:lib/batik/batik-util.jar:lib/batik/batik-xml.jar:.

# NAMES:= \
# 	uk/ac/sanger/artemis/OptionChangeListener \
# 	uk/ac/sanger/artemis/OptionChangeEvent \
# 	uk/ac/sanger/artemis/Options \
# 	uk/ac/sanger/artemis/Selection \
# 	uk/ac/sanger/artemis/components/ArtemisMain \
#   uk/ac/sanger/artemis/components/ActMain \
# 	uk/ac/sanger/artemis/components/Splash \
# 	uk/ac/sanger/artemis/components/ListDialog \
# 	uk/ac/sanger/artemis/components/FeatureListFrame \
# 	uk/ac/sanger/artemis/components/FeatureList \
# 	uk/ac/sanger/artemis/components/FeatureDisplay \
# 	uk/ac/sanger/artemis/components/CanvasPanel \
# 	uk/ac/sanger/artemis/components/EntryGroupDisplay \
# 	uk/ac/sanger/artemis/components/EntryHeaderEdit \
# 	uk/ac/sanger/artemis/components/BasePlotGroup \
# 	uk/ac/sanger/artemis/components/EntryEditVector \
# 	uk/ac/sanger/artemis/components/EntryEdit \
# 	uk/ac/sanger/artemis/components/WritableEMBLCorbaEntrySource \
# 	uk/ac/sanger/artemis/components/DbfetchEntrySource \
# 	uk/ac/sanger/artemis/components/EMBLCorbaEntrySource \
# 	uk/ac/sanger/artemis/components/CorbaEntrySource \
# 	uk/ac/sanger/artemis/components/BioJavaEntrySource \
# 	uk/ac/sanger/artemis/components/FileDialogEntrySource \
# 	uk/ac/sanger/artemis/components/EntryActionListener \
# 	uk/ac/sanger/artemis/components/DisplayComponent \
# 	uk/ac/sanger/artemis/components/EntryFileDialog \
# 	uk/ac/sanger/artemis/components/StickyFileChooser \
# 	uk/ac/sanger/artemis/components/LogReadListener \
# 	uk/ac/sanger/artemis/components/Selector \
# 	uk/ac/sanger/artemis/components/Navigator \
# 	uk/ac/sanger/artemis/components/EntryGroupPanel \
# 	uk/ac/sanger/artemis/components/FeaturePopup \
# 	uk/ac/sanger/artemis/components/WriteMenu \
# 	uk/ac/sanger/artemis/components/SelectMenu \
# 	uk/ac/sanger/artemis/components/GraphMenu \
# 	uk/ac/sanger/artemis/components/RunMenu \
# 	uk/ac/sanger/artemis/components/AddMenu \
# 	uk/ac/sanger/artemis/components/ViewMenu \
# 	uk/ac/sanger/artemis/components/GotoMenu \
# 	uk/ac/sanger/artemis/components/EditMenu \
# 	uk/ac/sanger/artemis/components/SelectionMenu \
# 	uk/ac/sanger/artemis/components/SelectionViewer \
# 	uk/ac/sanger/artemis/components/FeaturePlotGroup \
# 	uk/ac/sanger/artemis/components/FeaturePlot \
# 	uk/ac/sanger/artemis/components/FeatureEdit \
# 	uk/ac/sanger/artemis/components/FeatureViewer \
# 	uk/ac/sanger/artemis/components/SearchResultViewer \
# 	uk/ac/sanger/artemis/components/LogViewer \
# 	uk/ac/sanger/artemis/components/FileViewer \
# 	uk/ac/sanger/artemis/components/MessageFrame \
# 	uk/ac/sanger/artemis/components/MessageDialog \
# 	uk/ac/sanger/artemis/components/YesNoDialog \
# 	uk/ac/sanger/artemis/components/KeyChooser \
# 	uk/ac/sanger/artemis/components/KeyChoice \
# 	uk/ac/sanger/artemis/components/QualifierEditor \
# 	uk/ac/sanger/artemis/components/QualifierTextArea \
# 	uk/ac/sanger/artemis/components/QualifierChoice \
# 	uk/ac/sanger/artemis/components/ChoiceFrame \
# 	uk/ac/sanger/artemis/components/SelectionInfoDisplay \
# 	uk/ac/sanger/artemis/components/EntryGroupMenu \
# 	uk/ac/sanger/artemis/components/EntryGroupInfoDisplay \
# 	uk/ac/sanger/artemis/components/Utilities \
# 	uk/ac/sanger/artemis/components/MarkerRangeRequester \
# 	uk/ac/sanger/artemis/components/MarkerRangeRequesterListener \
# 	uk/ac/sanger/artemis/components/MarkerRangeRequesterEvent \
# 	uk/ac/sanger/artemis/components/TextRequester \
# 	uk/ac/sanger/artemis/components/TextRequesterListener \
# 	uk/ac/sanger/artemis/components/TextRequesterEvent \
# 	uk/ac/sanger/artemis/components/TextDialog \
# 	uk/ac/sanger/artemis/components/BasePlot \
# 	uk/ac/sanger/artemis/components/ProcessWatcher \
# 	uk/ac/sanger/artemis/components/ProcessWatcherEvent \
# 	uk/ac/sanger/artemis/components/ProcessWatcherListener \
# 	uk/ac/sanger/artemis/components/ExternalProgramOptions \
# 	uk/ac/sanger/artemis/components/Plot \
# 	uk/ac/sanger/artemis/components/PlotMouseListener \
# 	uk/ac/sanger/artemis/components/FeatureAminoAcidViewer \
# 	uk/ac/sanger/artemis/components/FeatureBaseViewer \
# 	uk/ac/sanger/artemis/components/SequenceViewer \
# 	uk/ac/sanger/artemis/components/FeatureInfo \
# 	uk/ac/sanger/artemis/components/DisplayAdjustmentListener \
# 	uk/ac/sanger/artemis/components/DisplayAdjustmentEvent \
# 	uk/ac/sanger/artemis/components/ScoreChanger \
# 	uk/ac/sanger/artemis/components/ScoreScrollbar \
# 	uk/ac/sanger/artemis/components/ScoreChangeListener \
# 	uk/ac/sanger/artemis/components/ScoreChangeEvent \
# 	uk/ac/sanger/artemis/components/InputStreamProgressDialog \
# 	uk/ac/sanger/artemis/plot/KarlinSigAlgorithm \
# 	uk/ac/sanger/artemis/plot/UserDataAlgorithm \
# 	uk/ac/sanger/artemis/plot/Codon12CorrelationAlgorithm \
# 	uk/ac/sanger/artemis/plot/ATDeviationAlgorithm \
# 	uk/ac/sanger/artemis/plot/GCDeviationAlgorithm \
# 	uk/ac/sanger/artemis/plot/GCFrameAlgorithm \
# 	uk/ac/sanger/artemis/plot/CodonUsageAlgorithm \
# 	uk/ac/sanger/artemis/plot/CodonUsageFormatException \
# 	uk/ac/sanger/artemis/plot/CodonUsageWeight \
# 	uk/ac/sanger/artemis/plot/CodonWeight \
# 	uk/ac/sanger/artemis/plot/AGWindowAlgorithm \
# 	uk/ac/sanger/artemis/plot/GCSDWindowAlgorithm \
# 	uk/ac/sanger/artemis/plot/GCWindowAlgorithm \
# 	uk/ac/sanger/artemis/plot/HydrophilicityAlgorithm \
# 	uk/ac/sanger/artemis/plot/HydroAlgorithm \
# 	uk/ac/sanger/artemis/plot/HydrophobicityAlgorithm \
# 	uk/ac/sanger/artemis/plot/CoilFeatureAlgorithm \
# 	uk/ac/sanger/artemis/plot/FeatureAlgorithm \
# 	uk/ac/sanger/artemis/plot/BaseAlgorithm \
# 	uk/ac/sanger/artemis/plot/Algorithm \
# 	uk/ac/sanger/artemis/Logger \
# 	uk/ac/sanger/artemis/ExternalProgramListener \
# 	uk/ac/sanger/artemis/ExternalProgramException \
# 	uk/ac/sanger/artemis/ExternalProgramVector \
# 	uk/ac/sanger/artemis/SimpleExternalProgramMonitor \
# 	uk/ac/sanger/artemis/ProcessMonitor \
# 	uk/ac/sanger/artemis/ProcessMonitor \
# 	uk/ac/sanger/artemis/ExternalProgramMonitor \
# 	uk/ac/sanger/artemis/ExternalProgram \
# 	uk/ac/sanger/artemis/EntryGroupChangeListener \
# 	uk/ac/sanger/artemis/EntryGroupChangeEvent \
# 	uk/ac/sanger/artemis/EntryChangeListener \
# 	uk/ac/sanger/artemis/EntryChangeEvent \
# 	uk/ac/sanger/artemis/FilteredEntryGroup \
# 	uk/ac/sanger/artemis/SimpleEntryGroup \
# 	uk/ac/sanger/artemis/EntryGroup \
# 	uk/ac/sanger/artemis/EntryVector \
# 	uk/ac/sanger/artemis/EntrySourceVector \
# 	uk/ac/sanger/artemis/EntrySource \
# 	uk/ac/sanger/artemis/Entry \
# 	uk/ac/sanger/artemis/LastSegmentException \
# 	uk/ac/sanger/artemis/FeatureFromVectorPredicate \
# 	uk/ac/sanger/artemis/FeatureKeyQualifierPredicate \
# 	uk/ac/sanger/artemis/FeatureKeyPredicate \
# 	uk/ac/sanger/artemis/FeaturePatternPredicate \
# 	uk/ac/sanger/artemis/FeaturePredicateConjunction \
# 	uk/ac/sanger/artemis/FeaturePredicate \
# 	uk/ac/sanger/artemis/FeaturePredicateVector \
# 	uk/ac/sanger/artemis/FeatureEnumeration \
# 	uk/ac/sanger/artemis/FeatureSegmentVector \
# 	uk/ac/sanger/artemis/FeatureChangeListener \
# 	uk/ac/sanger/artemis/FeatureChangeEvent \
# 	uk/ac/sanger/artemis/FeatureVector \
# 	uk/ac/sanger/artemis/Feature \
# 	uk/ac/sanger/artemis/FeatureSegment \
# 	uk/ac/sanger/artemis/ActionController \
# 	uk/ac/sanger/artemis/ActionVector \
# 	uk/ac/sanger/artemis/Action \
# 	uk/ac/sanger/artemis/ChangeListener \
# 	uk/ac/sanger/artemis/ChangeEventVector \
# 	uk/ac/sanger/artemis/ChangeEvent \
# 	uk/ac/sanger/artemis/GotoListener \
# 	uk/ac/sanger/artemis/GotoEvent \
# 	uk/ac/sanger/artemis/GotoEventSource \
# 	uk/ac/sanger/artemis/SimpleGotoEventSource \
# 	uk/ac/sanger/artemis/SelectionChangeListener \
# 	uk/ac/sanger/artemis/SelectionChangeEvent \
# 	uk/ac/sanger/artemis/Selectable \
# 	uk/ac/sanger/artemis/sequence/NoSequenceException \
# 	uk/ac/sanger/artemis/sequence/MarkerRangeVector \
# 	uk/ac/sanger/artemis/sequence/MarkerRange \
# 	uk/ac/sanger/artemis/sequence/MarkerChangeListener \
# 	uk/ac/sanger/artemis/sequence/MarkerChangeEvent \
# 	uk/ac/sanger/artemis/sequence/Marker \
# 	uk/ac/sanger/artemis/sequence/Strand \
# 	uk/ac/sanger/artemis/sequence/BasePatternFormatException \
# 	uk/ac/sanger/artemis/sequence/BasePattern \
# 	uk/ac/sanger/artemis/sequence/Bases \
# 	uk/ac/sanger/artemis/sequence/AminoAcidSequence \
# 	uk/ac/sanger/artemis/sequence/SequenceChangeListener \
# 	uk/ac/sanger/artemis/sequence/SequenceChangeEvent \
# 	uk/ac/sanger/artemis/io/BioJavaFeature \
# 	uk/ac/sanger/artemis/io/BioJavaSequence \
# 	uk/ac/sanger/artemis/io/BioJavaEntry \
# 	uk/ac/sanger/artemis/io/GenbankStreamSequence \
# 	uk/ac/sanger/artemis/io/RWCorbaEntry \
# 	uk/ac/sanger/artemis/io/CorbaEntry \
# 	uk/ac/sanger/artemis/io/DocumentEntryFactory \
# 	uk/ac/sanger/artemis/io/BlastDocumentEntry \
# 	uk/ac/sanger/artemis/io/MSPcrunchDocumentEntry \
# 	uk/ac/sanger/artemis/io/GFFDocumentEntry \
# 	uk/ac/sanger/artemis/io/GenbankDocumentEntry \
# 	uk/ac/sanger/artemis/io/EmblDocumentEntry \
# 	uk/ac/sanger/artemis/io/PublicDBDocumentEntry \
# 	uk/ac/sanger/artemis/io/SimpleDocumentEntry \
# 	uk/ac/sanger/artemis/io/DocumentEntry \
# 	uk/ac/sanger/artemis/io/Entry \
# 	uk/ac/sanger/artemis/io/ReadOnlyEntry \
# 	uk/ac/sanger/artemis/io/QualifierInfoException \
# 	uk/ac/sanger/artemis/io/QualifierInfo \
# 	uk/ac/sanger/artemis/io/QualifierInfoVector \
# 	uk/ac/sanger/artemis/io/QualifierInfoHash \
# 	uk/ac/sanger/artemis/chado/DbSqlConfig \
# 	uk/ac/sanger/artemis/circular/DNADraw \
# 	uk/ac/sanger/artemis/circular/digest/Utils \
# 	uk/ac/sanger/artemis/circular/digest/CircularGenomeController \
# 	uk/ac/sanger/artemis/io/GffToEMBL

ARTEMIS_DIRS = uk/ac/sanger/artemis \
uk/ac/sanger/artemis/chado \
uk/ac/sanger/artemis/circular \
uk/ac/sanger/artemis/circular/digest \
uk/ac/sanger/artemis/components \
uk/ac/sanger/artemis/components/alignment \
uk/ac/sanger/artemis/components/database \
uk/ac/sanger/artemis/components/filetree \
uk/ac/sanger/artemis/components/genebuilder \
uk/ac/sanger/artemis/components/genebuilder/cv \
uk/ac/sanger/artemis/components/genebuilder/gff \
uk/ac/sanger/artemis/components/genebuilder/ortholog \
uk/ac/sanger/artemis/components/variant \
uk/ac/sanger/artemis/editor \
uk/ac/sanger/artemis/io \
uk/ac/sanger/artemis/j2ssh \
uk/ac/sanger/artemis/plot \
uk/ac/sanger/artemis/sequence \
uk/ac/sanger/artemis/util

#CLASSES := $(NAMES:%=%.class)
SOURCES := $(foreach DIR,$(ARTEMIS_DIRS),$(wildcard $(DIR)/*.java))
CLASSES := $(SOURCES:%.java=%.class)

all: idl code

# Utils needs to be built before controller
uk/ac/sanger/artemis/circular/digest/CircularGenomeController.class:uk/ac/sanger/artemis/circular/digest/Utils.class
	$(REAL_CLASSPATH) $(JAVAC) $(@:%.class=%.java)

code: $(CLASSES)

topdown: idl
	$(REAL_CLASSPATH) $(JAVAC) uk/ac/sanger/artemis/components/ArtemisMain.java

%.class : %.java
	$(REAL_CLASSPATH) $(JAVAC) $<

idl : type/*.java nsdb/*.java seqdb/*.java

IDL = idlj
IDLCMD = $(IDL) -Icorba

type/*.java : corba/types.idl
	$(IDLCMD) corba/types.idl

nsdb/*.java : corba/nsdb.idl corba/nsdb_write.idl
	$(IDLCMD) corba/nsdb.idl
	$(IDLCMD) corba/nsdb_write.idl

seqdb/*.java : corba/seqdb.idl
	$(IDLCMD) corba/seqdb.idl

doc :
	$(REAL_CLASSPATH) javadoc -J-mx200m -version \
		AppGlobal.java \
		uk.ac.sanger.artemis uk.ac.sanger.artemis.components \
		uk.ac.sanger.artemis.sequence uk.ac.sanger.artemis.plot \
		uk.ac.sanger.artemis.util uk.ac.sanger.artemis.io

manual :
	(cd docs; make)

CLASS_FILES := `find org uk nsdb type seqdb -name '*.class' -print`

OTHER_FILES := `find images/PSUlogo.gif images/icon.gif COPYING README`

dist :
	rm -rf artemis_compiled.tar.gz tar_build
	mkdir tar_build
	mkdir tar_build/artemis
	rm -f artemis_compiled_latest.tar.gz
	tar cf - $(OTHER_FILES) act art Makefile corba etc | (cd tar_build/artemis; tar xf -)
	tar cf - artemis_sqlmap dnaplotter uk org nsdb type seqdb lib | (cd tar_build/artemis; tar xf -)
	(cd tar_build; find . -name 'CVS' -print | xargs rm -rf; find . -name '.svn' -print | xargs rm -rf; tar cvf ../artemis_compiled.tar artemis)

jar : all artemis.jar

artemis.jar : $(CLASSES)
	mkdir jar_build
	rm -f artemis.jar
	cd jar_build; \
	if [ ! -d org ]; then \
	  for fileJar in ../lib/*.jar; do \
	    jar xvf $$fileJar; \
	    rm -rf META-INF/MANIFEST.MF; \
	  done; \
          for fileJar in ../lib/j2ssh/*.jar; do \
            jar xvf $$fileJar; \
            rm -rf META-INF/MANIFEST.MF; \
          done; \
          for fileJar in ../lib/ibatis/*.jar; do \
            jar xvf $$fileJar; \
            rm -rf META-INF/MANIFEST.MF; \
          done; \
          for fileJar in ../lib/batik/*.jar; do \
            jar xvf $$fileJar; \
            rm -rf META-INF/MANIFEST.MF; \
          done; \
          for fileJar in ../lib/picard/*.jar; do \
            jar xvf $$fileJar; \
            rm -rf META-INF/MANIFEST.MF; \
          done; \
        fi; \
	cp -R ../lib/LICENSE.Apache ../uk ../org ../nsdb ../type ../seqdb ../etc ../images ../lib/j2ssh/j2ssh.properties \
	      ../images/PSUlogo.gif ../images/icon.gif ../README ../artemis_sqlmap .
	find jar_build -name '*.java' -print | xargs rm -f
	find jar_build -name '.svn' -print | xargs rm -rf
	cd jar_build; \
	rm -rf META-INF/MANIFEST.MF; \
	echo "Main-Class: uk.ac.sanger.artemis.components.ArtemisMain\nPermissions: all-permissions" > manifest-art; \
	jar cmf manifest-art artemis.jar META-INF/services images/PSUlogo.gif images/icon.gif README etc \
	                     artemis_sqlmap org uk com net nsdb type seqdb LICENSE.Apache j2ssh.properties; \
        echo "Main-Class: uk.ac.sanger.artemis.circular.DNADraw\nPermissions: all-permissions" > manifest-circular; \
        jar cmf manifest-circular DNAPlotter.jar images/PSUlogo.gif README etc \
                             uk org/gmod org/w3c org/apache org/biojava/bio/ com/ibatis/common/jdbc/ net/sf/samtools/ LICENSE.Apache j2ssh.properties; \
	echo "Main-Class: uk.ac.sanger.artemis.components.alignment.BamView\nPermissions: all-permissions" > manifest-bamview; \
	jar cmf manifest-bamview BamView.jar META-INF/services etc uk org/apache org/biojava org/biojavax org/gmod org/w3c net/sf com/ibatis; \
	echo "Main-Class: uk.ac.sanger.artemis.components.ActMain\nPermissions: all-permissions" > manifest-act; \
	jar cmf manifest-act act.jar META-INF/services images/PSUlogo.gif images/icon.gif README etc \
	                     artemis_sqlmap org uk com net nsdb type seqdb LICENSE.Apache j2ssh.properties; \
	rm -f etc/log4j.properties; \
	jar cmf manifest-art artemis_mac.jar images/PSUlogo.gif images/icon.gif README \
	        uk org/gmod nsdb type seqdb LICENSE.Apache artemis_sqlmap

clean :
	-rm -rf *.html artemis.jar seqdb nsdb type resources uk/ac/sanger/jcon/ jar_build tar_build  artemis_compiled.tar
	-rm -rf TAGS* *.o
	-find . -name '*.class' -print | xargs rm -f
