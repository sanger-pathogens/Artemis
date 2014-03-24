# This is a GNU Makefile for Artemis

# $Header: //tmp/pathsoft/artemis/Makefile,v 1.47 2009-09-21 15:32:03 tjc Exp $

SHELL=/bin/sh

#OPT_FLAGS = -g -deprecation

JAVAC := javac -source 1.5 -target 1.5 $(OPT_FLAGS) $(EXTRA_FLAGS)

REAL_CLASSPATH := CLASSPATH=lib/commons-lang-2.6.jar:lib/biojava.jar:lib/jemAlign.jar:lib/j2ssh/j2ssh-core.jar:lib/ibatis/ibatis-2.3.4.726.jar:lib/ibatis/log4j-1.2.14.jar:lib/postgresql-8.4-701.jdbc3.jar:lib/picard/picard.jar:lib/picard/sam.jar:lib/commons-net-2.2.jar:lib/batik/batik-awt-util.jar:lib/batik/batik-dom.jar:lib/batik/batik-ext.jar:lib/batik/batik-svggen.jar:lib/batik/batik-util.jar:lib/batik/batik-xml.jar:.

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
