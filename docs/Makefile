# $Header: //tmp/pathsoft/artemis/docs/Makefile,v 1.2 2009-01-07 09:25:57 tjc Exp $

SHELL=/bin/sh

SANGER_CAT = /software/pathogen/projects/artemis/docbook-sgml/dtd_4.1/CATALOG

JADE = jade
JADE_PARAM = -t sgml -c $(SANGER_CAT) -d ../artemis.dsl\#html
JADETEX_PARAM = -t tex -c $(SANGER_CAT) -d ../artemis.dsl\#tex

html : art_html act_html 

art_html : art_html_build
act_html : act_html_build

tex : art_tex act_tex

art_tex : art_tex_build
act_tex : act_tex_build

art_html_build : *.sgml *.dsl *.css
	(if [ -d art_html_build ] ; \
	 then \
	   touch art_html_build; \
	 else \
	   mkdir art_html_build; \
	 fi; \
	 cd art_html_build; ln -sf ../*.sgml ../*.gif ../*.png ../*.eps ../*.css .; \
	 $(JADE) $(JADE_PARAM) manual.sgml;\
	 $(JADE) -V nochunks $(JADE_PARAM) manual.sgml > artemis_manual_complete.html)

act_html_build : *.sgml *.dsl *.css
	(if [ -d act_html_build ] ; \
	 then \
	   touch act_html_build; \
	 else \
	   mkdir act_html_build; \
	 fi; \
	 cd act_html_build; ln -sf ../*.sgml ../*.gif ../*.png ../*.eps ../*.css .; \
	 $(JADE) $(JADE_PARAM) act_manual.sgml;\
	 $(JADE) -V nochunks $(JADE_PARAM) act_manual.sgml > act_manual_complete.html)

art_tex_build : *.sgml *.dsl *.css
	(if [ -d art_tex_build ] ; \
	 then \
	   touch art_tex_build; \
	 else \
	   mkdir art_tex_build; \
	 fi; \
	 cd art_tex_build; ln -sf ../*.sgml ../*.gif ../*.ps ../*.eps .; \
	 $(JADE) $(JADETEX_PARAM) manual.sgml)

act_tex_build : *.sgml *.dsl *.css
	(if [ -d act_tex_build ] ; \
	 then \
	   touch act_tex_build; \
	 else \
	   mkdir act_tex_build; \
	 fi; \
	 cd act_tex_build; ln -sf ../*.sgml ../*.gif ../*.ps ../*.eps .; \
	 $(JADE) $(JADETEX_PARAM) act_manual.sgml)

tex_format : art_tex_format act_tex_format

art_tex_format : 
	(cd art_tex_build; latex -fmt jadetex manual.tex; \
	 latex -fmt jadetex manual.tex)

act_tex_format :
	(cd act_tex_build; latex -fmt jadetex act_manual.tex; \
	 latex -fmt jadetex act_manual.tex)

all : clean html tex tex_format

clean :
	rm -rf art_html_build act_html_build art_tex_build act_tex_build 
