<!DOCTYPE style-sheet PUBLIC "-//James Clark//DTD DSSSL Style Sheet//EN" [
<!ENTITY html-docbook-ss SYSTEM "/nfs/disk222/yeastpub/bio-soft/docbook-sgml/stylesheets_1.64/html/docbook.dsl" CDATA DSSSL>
<!ENTITY tex-docbook-ss SYSTEM "/nfs/disk222/yeastpub/bio-soft/docbook-sgml/stylesheets_1.64/print/docbook.dsl" CDATA DSSSL>
]>

<!-- 
$Header: //tmp/pathsoft/artemis/docs/artemis.dsl,v 1.1 2004-06-10 09:20:11 tjc Exp $
-->

<style-specification id="html" use="html-spec">
  <style-specification-body>
(define %body-attr%
  ;; REFENTRY body-attr
  ;; PURP What attributes should be hung off of BODY?
  ;; DESC
  ;; A list of the the BODY attributes that should be generated.
  ;; The format is a list of lists, each interior list contains the
  ;; name and value of a BODY attribute.
  ;; /DESC
  ;; AUTHOR N/A
  ;; /REFENTRY
  (list
   (list "BGCOLOR" "#FFFFFF")
   (list "TEXT" "#000000")))

<!-- (define %html-header-tags%
  ;; REFENTRY html-header-tags
  ;; PURP What additional HEAD tags should be generated?
  ;; DESC
  ;; A list of the the HTML HEAD tags that should be generated.
  ;; The format is a list of lists, each interior list consists
  ;; of a tag name and a set of attribute/value pairs:
  ;; '(("META" ("NAME" "name") ("CONTENT" "content")))
  ;; /DESC
  ;; AUTHOR N/A
  ;; /REFENTRY
  '())
-->

(define html-index
  ;; REFENTRY html-index
  ;; PURP HTML indexing?
  ;; DESC
  ;; Turns on HTML indexing.  If true, then index data will be written
  ;; to the file defined by 'html-index-filename'.  This data can be
  ;; collated and turned into a DocBook index with bin/collateindex.pl.
  ;; /DESC
  ;; AUTHOR N/A
  ;; /REFENTRY
  #t)

(define html-manifest
  ;; REFENTRY html-manifest
  ;; PURP Write a manifest?
  ;; DESC
  ;; If not '#f' then the list of HTML files created by the stylesheet
  ;; will be written to the file named by 'html-manifest-filename'.
  ;; /DESC
  ;; AUTHOR N/A
  ;; /REFENTRY
  #t)

(define nochunks
  ;; REFENTRY nochunks
  ;; PURP Suppress chunking of output pages
  ;; DESC
  ;; If true, the entire source document is formatted as a single HTML
  ;; document and output on stdout.
  ;; (This option can conveniently be set with '-V nochunks' on the
  ;; Jade command line).
  ;; /DESC
  ;; AUTHOR N/A
  ;; /REFENTRY
  #f)

(define rootchunk
  ;; REFENTRY rootchunk
  ;; PURP Make a chunk for the root element when nochunks is used
  ;; DESC
  ;; If true, a chunk will be created for the root element, even though
  ;; nochunks is specified. This option has no effect if nochunks is not
  ;; true.
  ;; (This option can conveniently be set with '-V rootchunk' on the
  ;; Jade command line).
  ;; /DESC
  ;; AUTHOR N/A
  ;; /REFENTRY
  #f)

(define %link-mailto-url%
  ;; REFENTRY link-mailto-url
  ;; PURP Mailto URL for LINK REL=made
  ;; DESC
  ;; If not '#f', the '%link-mailto-url%' address will be used in a
  ;; LINK REL=made element in the HTML HEAD.
  ;; /DESC
  ;; AUTHOR N/A
  ;; /REFENTRY
  "artemis@sanger.ac.uk")

(define %show-comments%
  ;; REFENTRY show-comments
  ;; PURP Display Comment elements?
  ;; DESC
  ;; If true, comments will be displayed, otherwise they are suppressed.
  ;; Comments here refers to the 'Comment' element, which will be renamed
  ;; 'Remark' in DocBook V4.0, not SGML/XML comments which are unavailable.
  ;; /DESC
  ;; AUTHOR N/A
  ;; /REFENTRY
  #f)

(define %stylesheet%
  ;; REFENTRY stylesheet
  ;; PURP Name of the stylesheet to use
  ;; DESC
  ;; The name of the stylesheet to place in the HTML LINK TAG, or '#f' to
  ;; suppress the stylesheet LINK.
  ;; /DESC
  ;; AUTHOR N/A
  ;; /REFENTRY
  "psu.css")

(define %html-ext% ".html")

; Includes a summary at the beginning of an item. 
(define %generate-article-toc% #t)

(define %use-id-as-filename% #t)

  </style-specification-body>
</style-specification>


<style-specification id="tex" use="tex-spec">

  <style-specification-body>
; Includes a summary at the beginning of an item.    
(define %generate-article-toc% #t)

  </style-specification-body>
</style-specification>

<external-specification id="html-spec" document="html-docbook-ss">
<external-specification id="tex-spec" document="tex-docbook-ss">
