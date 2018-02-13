# This is a GNU Makefile for Artemis

SHELL=/bin/sh

#OPT_FLAGS = -g -deprecation

JAVAC := javac -source 1.8 -target 1.8 $(OPT_FLAGS) $(EXTRA_FLAGS)

REAL_CLASSPATH := CLASSPATH=lib/commons-lang-2.6.jar:lib/biojava.jar:lib/jemAlign.jar:lib/j2ssh/j2ssh-core.jar:lib/ibatis/ibatis-2.3.4.726.jar:lib/ibatis/log4j-1.2.14.jar:lib/postgresql-8.4-701.jdbc3.jar:lib/picard/picard.jar:lib/commons-net-3.6.jar:lib/batik/batik-awt-util.jar:lib/batik/batik-dom.jar:lib/batik/batik-ext.jar:lib/batik/batik-svggen.jar:lib/batik/batik-util.jar:lib/batik/batik-xml.jar:.

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
uk/ac/sanger/artemis/components/ref \
uk/ac/sanger/artemis/editor \
uk/ac/sanger/artemis/io \
uk/ac/sanger/artemis/j2ssh \
uk/ac/sanger/artemis/plot \
uk/ac/sanger/artemis/sequence \
uk/ac/sanger/artemis/util

SOURCES := $(foreach DIR,$(ARTEMIS_DIRS),$(wildcard $(DIR)/*.java))
CLASSES := $(SOURCES:%.java=%.class)

all: code

# Utils needs to be built before controller
uk/ac/sanger/artemis/circular/digest/CircularGenomeController.class:uk/ac/sanger/artemis/circular/digest/Utils.class
	$(REAL_CLASSPATH) $(JAVAC) $(@:%.class=%.java)

code: $(CLASSES)

topdown:
	$(REAL_CLASSPATH) $(JAVAC) uk/ac/sanger/artemis/components/ArtemisMain.java

%.class : %.java
	$(REAL_CLASSPATH) $(JAVAC) $<

doc :
	$(REAL_CLASSPATH) javadoc -J-mx200m -version \
		AppGlobal.java \
		uk.ac.sanger.artemis uk.ac.sanger.artemis.components \
		uk.ac.sanger.artemis.sequence uk.ac.sanger.artemis.plot \
		uk.ac.sanger.artemis.util uk.ac.sanger.artemis.io

manual :
	(cd docs; make)

CLASS_FILES := `find org uk -name '*.class' -print`

OTHER_FILES := `find images/Wellcome_Sanger_Institute_Logo_Landscape_Digital_RGB_Full_Colour.jpg images/icon.gif README.md`

dist :
	rm -rf artemis_compiled.tar.gz tar_build
	mkdir tar_build
	mkdir tar_build/artemis
	tar cf - $(OTHER_FILES) act art dnaplotter bamview setenv Makefile etc | (cd tar_build/artemis; tar xf -)
	tar cf - artemis_sqlmap dnaplotter uk org lib | (cd tar_build/artemis; tar xf -)
	(cd tar_build; find . -name 'CVS' -print | xargs rm -rf; find . -name '.svn' -print | xargs rm -rf; tar cvf ../artemis_compiled.tar artemis)

jar : all artemis.jar

artemis.jar : $(CLASSES)
	mkdir -p jar_build/build
	rm -f *.jar
	cd jar_build/build; \
    for fileJar in ../../lib/*.jar; do \
      jar xvf $$fileJar; \
      rm -rf META-INF/MANIFEST.MF; \
    done; \
    for fileJar in ../../lib/j2ssh/*.jar; do \
      jar xvf $$fileJar; \
      rm -rf META-INF/MANIFEST.MF; \
    done; \
    for fileJar in ../../lib/ibatis/*.jar; do \
      jar xvf $$fileJar; \
      rm -rf META-INF/MANIFEST.MF; \
    done; \
    for fileJar in ../../lib/batik/*.jar; do \
      jar xvf $$fileJar; \
      rm -rf META-INF/MANIFEST.MF; \
    done; \
    for fileJar in ../../lib/picard/*.jar; do \
      jar xvf $$fileJar; \
      rm -rf META-INF/MANIFEST.MF; \
    done; \
	cp -R ../../lib/LICENSE.* ../../uk ../../org ../../etc ../../images ../../lib/j2ssh/j2ssh.properties \
	      ../../images/Wellcome_Sanger_Institute_Logo_Landscape_Digital_RGB_Full_Colour.jpg ../../images/icon.gif ../../README.md ../../artemis_sqlmap .
	find jar_build/build -name '*.java' -print | xargs rm -f
	find jar_build/build -name '.svn' -print | xargs rm -rf
	find jar_build/build -name '*.DS_Store' -print | xargs rm -rf
	cd jar_build/build; \
	rm -rf META-INF/MANIFEST.MF; \
	echo "Main-Class: uk.ac.sanger.artemis.components.ArtemisMain\nPermissions: all-permissions" > manifest-art; \
	jar cmf manifest-art ../artemis.jar META-INF/services images/Wellcome_Sanger_Institute_Logo_Landscape_Digital_RGB_Full_Colour.jpg images/icon.gif README.md etc \
	                     artemis_sqlmap org uk com net htsjdk picard gov joptsimple ngs freemarker LICENSE.Apache LICENSE.Picard LICENSE.JDBC LICENSE j2ssh.properties; \
    echo "Main-Class: uk.ac.sanger.artemis.circular.DNADraw\nPermissions: all-permissions" > manifest-circular; \
    jar cmf manifest-circular ../DNAPlotter.jar images/Wellcome_Sanger_Institute_Logo_Landscape_Digital_RGB_Full_Colour.jpg README.md etc \
                         org uk com net htsjdk picard gov joptsimple ngs freemarker LICENSE.Apache LICENSE.Picard LICENSE.JDBC LICENSE j2ssh.properties; \
	echo "Main-Class: uk.ac.sanger.artemis.components.alignment.BamView\nPermissions: all-permissions" > manifest-bamview; \
	jar cmf manifest-bamview ../bamview.jar META-INF/services etc org uk com net htsjdk picard gov joptsimple ngs freemarker LICENSE.Apache LICENSE.Picard LICENSE.JDBC LICENSE; \
	echo "Main-Class: uk.ac.sanger.artemis.components.ActMain\nPermissions: all-permissions" > manifest-act; \
	jar cmf manifest-act ../act.jar META-INF/services images/Wellcome_Sanger_Institute_Logo_Landscape_Digital_RGB_Full_Colour.jpg images/icon.gif README.md etc \
	                     artemis_sqlmap org uk com net htsjdk picard gov joptsimple ngs freemarker LICENSE.Apache LICENSE.Picard LICENSE.JDBC LICENSE j2ssh.properties; \
	
	rm -rf jar_build/build;

clean :
	-rm -rf *.html artemis.jar resources uk/ac/sanger/jcon/ jar_build tar_build  artemis_compiled.tar
	-rm -rf TAGS* *.o
	-find . -name '*.class' -print | xargs rm -f
