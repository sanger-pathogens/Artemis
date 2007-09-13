/* QualifierTextArea.java
 *
 * created: Tue Oct 23 2001
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2000  Genome Research Limited
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/QualifierTextArea.java,v 1.5 2007-09-13 10:15:28 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.editor.BrowserControl;
import uk.ac.sanger.artemis.editor.DataCollectionPane;
import uk.ac.sanger.artemis.io.QualifierParseException;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.EmblStreamFeature;
import uk.ac.sanger.artemis.util.StringVector;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Point;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.IOException;
import java.io.StringReader;

import javax.swing.JTextPane;
import javax.swing.text.BadLocationException;
import javax.swing.text.Style;
import javax.swing.text.StyleConstants;
import javax.swing.text.StyledDocument;

/**
 *  This component is a JTextPane that understands qualifiers.
 *  It provides hyperlinks to databases using a StyledDocument
 *  (rather than using HyperlinkListener) so it can remain editable.
 **/
public class QualifierTextArea extends JTextPane
{
  private static String[] DATABASES = 
          { "SWALL", "EMBL", "UniProt", "PMID", "PubMed", "InterPro" };

  /**
   *  Create a new QualifierTextArea containing no text.
   **/
  public QualifierTextArea () 
  {
    super();
    initStyles();
    int nrows = (Options.getOptions ().getPropertyTruthValue ("alicat_mode") ||
            Options.getOptions ().getPropertyTruthValue ("val_mode") ?
            40 :
            18);
    int ncolumns = 81;
    
    setPreferredSize(new Dimension( ncolumns*getColumnWidth(), nrows*getRowHeight() ));
    setBackground (Color.white);

    try  // no such method in java1.3
    {
      setDragEnabled(true);
    }
    catch(java.lang.NoSuchMethodError err){}
    
    addMouseListener(new MouseAdapter()
    {
      public void mouseClicked(MouseEvent e)
      {
        if(e.getClickCount() == 2) 
          handleMouseDoubleClick(e);
      }
    });
  }

  
  public void append(final String s)
  {
    StyledDocument doc = super.getStyledDocument();
    try
    {
      doc.insertString(doc.getLength(), s, getLogicalStyle());
      
      for(int i=0; i<DATABASES.length; i++)
        setStyleForHyperLinks(s, DATABASES[i]);
    }
    catch(BadLocationException e)
    {
      e.printStackTrace();
    } 
  }
  
  /**
   * Add hyperlink style
   */
  private void initStyles()
  {
    // Makes text Blue
    Style style = addStyle("Blue", null);
    StyleConstants.setForeground(style, Color.blue);
    
    // Inherits from "Blue"; makes text red and underlined
    style = addStyle("Blue Underline", style);
    StyleConstants.setUnderline(style, true);
  }
  
  /**
   * Override to ensure hyperlinks are set
   */
  public void setText(String text)
  {
    super.setText(text);
    for(int i=0; i<DATABASES.length; i++)
      setStyleForHyperLinks(text, DATABASES[i]);
  }

  /**
   *  Parse and return the qualifiers in this TextArea in a QualifierVector.
   **/
  public QualifierVector
    getParsedQualifiers (final EntryInformation entry_information)
      throws QualifierParseException 
  {
    final String qualifier_string = getText ();
    return getQualifiersFromString (qualifier_string,
                                    entry_information);
  }

  /**
   *  Return a QualifierVector containing the qualifiers from a String.
   *  @param qual_string contains the qualifiers to parse
   */
  private static QualifierVector
    getQualifiersFromString (final String qual_string,
                             final EntryInformation entry_information)
      throws QualifierParseException 
  {

    try 
    {
      final StringReader string_reader = new StringReader (qual_string);
      final QualifierVector embl_qualifiers =
        EmblStreamFeature.readQualifiers (string_reader,
                                          entry_information);

      string_reader.close();
      return embl_qualifiers;
    } 
    catch (IOException exception) 
    {
      throw (new QualifierParseException (exception.getMessage ()));
    }
  }
  
  /**
   * Analogous to JTextArea column
   * @return
   */
  private int getColumnWidth() 
  {
    FontMetrics metrics = getFontMetrics(getFont());
    return metrics.charWidth('m');
  }
  
  /**
   * Analogous to JTextArea row
   * @return
   */
  private int getRowHeight() 
  {
    FontMetrics metrics = getFontMetrics(getFont());
    return metrics.getHeight();
  }
  
  /**
   * This routine sets the hyerlink style for a given database
   * @param s
   * @param db
   */
  private void setStyleForHyperLinks(final String s, 
                                     final String db)
  {
    int ind = 0;
    while((ind = s.indexOf(db, ind)) > -1)
    {
      int ind2 = getEndOfLink(s,ind);
      getStyledDocument().setCharacterAttributes(ind, ind2-ind, 
                             getStyle("Blue Underline"), true);
      ind = ind+1;
    }
  }
  
  private int getEndOfLink(String s, int ind)
  {
    int ind2 = s.indexOf(" ",ind);
    if(ind2 == -1)
      ind2 = s.indexOf(";",ind);

    int ind3 = s.indexOf(")",ind);
    if(ind3>-1 && (ind3<ind2 || ind2 == -1))
      ind2 = ind3;
    ind3 = s.indexOf(";",ind);
    if(ind3>-1 && (ind3<ind2 || ind2 == -1))
      ind2 = ind3;
    ind3 = s.indexOf("\"",ind);
    if(ind3>-1 && (ind3<ind2 || ind2 == -1))
      ind2 = ind3;
    ind3 = s.indexOf("\n",ind);
    if(ind3>-1 && (ind3<ind2 || ind2 == -1))
      ind2 = ind3;
    return ind2;
  }
  
  private int getStartOfLink(final String s)
  {
    int lastIndexLink = -1;
    for(int i=0; i<DATABASES.length; i++)
    {
      int index = s.lastIndexOf(DATABASES[i]);
      if(index > lastIndexLink)
        lastIndexLink = index;
    }
    return lastIndexLink;
  }
  
  public static String getPubMedSite()
  {
    StringVector pubmed = Options.getOptions().getOptionValues("pubmed_url");
    if(pubmed != null)
      return (String)pubmed.elementAt(0);
    return "http://www.ncbi.nlm.nih.gov/sites/entrez?Db=pubmed&Cmd=ShowDetailView&TermToSearch=";
  }
  
  public static String getInterProSite()
  {
	return "http://www.ebi.ac.uk/interpro/ISearch?query=";
  }
  
  /**
   * Process double click event.
   * @param e
   */
  private void handleMouseDoubleClick(final MouseEvent e)
  {
    Point pt = new Point(e.getX(), e.getY());
    int pos = viewToModel(pt);
    int start;
    if(pos < 15)
      start = 0;
    else
      start = pos - 15;
    
    int length = 30;
    
    if( (start+30) >  getStyledDocument().getLength() )
      length = getStyledDocument().getLength()-start;
    
    try
    {
      String textAtPosition = getStyledDocument().getText(start, length);

      int indEnd = getEndOfLink(textAtPosition, (pos-start));
      textAtPosition = textAtPosition.substring(0, indEnd);
      
      int intStart = getStartOfLink(textAtPosition);
      if(intStart < 0)
        return;
      
      textAtPosition = textAtPosition.substring(intStart);
      
      String cmd;
      if(textAtPosition.indexOf("PMID")   > -1 ||
         textAtPosition.indexOf("PubMed") > -1)
      {
        String id[] = textAtPosition.split(":");
        if(id.length < 2)
          return;
        cmd = getPubMedSite()+id[1];
      }
      else if(textAtPosition.indexOf("InterPro") > -1)
      {
        String id[] = textAtPosition.split(":");
        if(id.length < 2)
          return;
        cmd = getInterProSite()+id[1];
      }
      else
      {
        cmd = DataCollectionPane.getSrsSite()+"/wgetz?-e+["+textAtPosition+"]";

        // link to uniprot accession
        if( (intStart = cmd.indexOf("UniProt:")) > -1)
          cmd = cmd.substring(0,intStart+7)+"-acc:"+
                cmd.substring(intStart+8);
        if(cmd.indexOf("ebi.ac.uk") > -1)
          cmd = cmd + "+-vn+2";
      }
      BrowserControl.displayURL(cmd);
    }
    catch(BadLocationException e1)
    {
      e1.printStackTrace();
    }
  }
}
