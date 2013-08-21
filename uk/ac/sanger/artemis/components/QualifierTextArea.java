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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/QualifierTextArea.java,v 1.17 2008-08-04 08:33:45 tjc Exp $
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
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Point;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.io.IOException;
import java.io.StringReader;
import java.util.Vector;

import javax.swing.JComponent;
import javax.swing.JTextPane;
import javax.swing.text.BadLocationException;
import javax.swing.text.Style;
import javax.swing.text.StyleConstants;
import javax.swing.text.StyleContext;
import javax.swing.text.StyledDocument;

/**
 *  This component is a JTextPane that understands qualifiers.
 *  It provides hyperlinks to databases using a StyledDocument
 *  (rather than using HyperlinkListener) so it can remain editable.
 **/
public class QualifierTextArea extends JTextPane
    implements MouseMotionListener
{
  private static final long serialVersionUID = 1L;
  private static Cursor cbusy = new Cursor(Cursor.WAIT_CURSOR);
  private static Cursor cdone = new Cursor(Cursor.DEFAULT_CURSOR);
  private static Cursor chand = new Cursor(Cursor.HAND_CURSOR);
  private static Style DEFAULT_STYLE = 
    StyleContext.getDefaultStyleContext().getStyle(StyleContext.DEFAULT_STYLE);
  private static StringVector dbsLinks;
  public static Vector<String> DATABASES = new Vector<String>();
  
  private boolean useHyperlinks = true;
  
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
    
    if(dbsLinks == null)
    {
      dbsLinks = Options.getOptions().getOptionValues("hyperlinks");
      for(int i=0; i<dbsLinks.size(); i+=2)
      {
        String dbs[] = ((String)dbsLinks.get(i)).split("\\+");
        for(int j=0; j<dbs.length; j++)
          DATABASES.add(dbs[j]);
      }
    }
    
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
        if(e.getClickCount() == 1)
          handleMouseSingleClick(getHyperlinkTextAtMouseEvent(e),
              QualifierTextArea.this);
      }
    });
    
    super.addMouseMotionListener(this);
  }

  
  public void append(final String s)
  {
    StyledDocument doc = super.getStyledDocument();
    try
    {
      doc.insertString(doc.getLength(), s, getLogicalStyle());
      
      for(int i=0; i<DATABASES.size(); i++)
        setStyleForHyperLinks(s, (String)DATABASES.get(i));
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
    // ensure we have the default style set
    getStyledDocument().setCharacterAttributes(0, text.length(), 
                 DEFAULT_STYLE, true);
    for(int i=0; i<DATABASES.size(); i++)
      setStyleForHyperLinks(text, (String)DATABASES.get(i));
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
    if(!isUseHyperlinks())
      return;
    
    int ind = 0;
    while((ind = indexOfIgnoreCase(s, db+":", ind)) > -1)
    {
      if(ind == 0 || s.charAt(ind-1) != '%')
      {
        int ind2 = getEndOfLink(s,ind);
        int ind3 = s.indexOf('/', ind);
        
        if(ind3 < 0 || ind3 > ind2)
          getStyledDocument().setCharacterAttributes(ind, ind2-ind, 
                                 getStyle("Blue Underline"), true);
      }

      ind = ind+1;
    }
  }
  
  public static int getEndOfLink(String s, int ind)
  {
    final char endOfLinkChar[] = { 
        ' ',
        ';',
        ')',
        ']',
        ';',
        ',',
        '\"',
        '|',
        '\n' };
    
    
    int ind2 = s.indexOf(endOfLinkChar[0],ind);
    int ind3;
    
    for(int i=1; i<endOfLinkChar.length; i++)
    {
      ind3 = s.indexOf(endOfLinkChar[i],ind);
      if(ind3>-1 && (ind3<ind2 || ind2 == -1))
        ind2 = ind3;
    }
    
    return ind2;
  }
  
  private int getStartOfLink(final String s)
  {
    int lastIndexLink = -1;
    for(int i=0; i<DATABASES.size(); i++)
    {
      int index = lastIndexOfIgnoreCase(s, (String)DATABASES.get(i)+":");
      if(index > lastIndexLink)
        lastIndexLink = index;
    }
    return lastIndexLink;
  }

  
  /**
   * Process double click event.
   * @param hyperlinkText
   * @param c
   */
  public static void handleMouseSingleClick(
      final String hyperlinkText,
      final JComponent c)
  {
    if(hyperlinkText == null)
      return;
    
    for(int i=0; i<dbsLinks.size(); i+=2)
    {
      String names[] = ((String)dbsLinks.get(i)).split("\\+");
      
      for(int j=0; j<names.length; j++)
      {
        if(indexOfIgnoreCase(hyperlinkText, names[j], 0) > -1)
        {
          String id[] = hyperlinkText.split(":");
          if(id.length < 2)
            return;
          
          String link = (String)dbsLinks.get(i+1);
          
          if(link.equals("srs_url"))
            sendToBrowser( getSrsLink(hyperlinkText), c );
          else
            sendToBrowser(link + id[1], c);
          return;
        }
      }
    }
  }

  /**
   * Get SRS hyperlink
   * @param hyperlinkText
   * @return
   */
  private static String getSrsLink(String hyperlinkText)
  {
    String cmd = DataCollectionPane.getSrsSite() + 
               "/wgetz?-e+[" + hyperlinkText + "]";

    int ind = cmd.indexOf("UniProt:");
    // link to uniprot accession
    if(ind > -1)
      cmd = cmd.substring(0, ind + 7) + "-acc:"
                       + cmd.substring(ind + 8);

    ind = cmd.indexOf("UniProtKB:");
    // link to uniprotkb accession
    if(ind > -1)
      cmd = cmd.substring(0, ind + 7) + "-acc:"
                      + cmd.substring(ind + 10);

    if(cmd.indexOf("ebi.ac.uk") > -1)
      cmd = cmd + "+-vn+2";
    return cmd;
  }

  /**
   * Send hyperlink to browser
   * @param cmd
   */
  private static void sendToBrowser(final String cmd, final JComponent c)
  {
    SwingWorker browserLaunch = new SwingWorker()
    {
      public Object construct()
      {
        c.setCursor(cbusy);
        BrowserControl.displayURL(cmd);
        c.setCursor(cdone);
        return null;
      }
    };
    browserLaunch.start();
  }
  
  
  private int lastIndexOfIgnoreCase(String a, String b)
  {
    return a.toLowerCase().lastIndexOf(b.toLowerCase());
  }
  
  private static int indexOfIgnoreCase(String a, String b, int fromPos)
  {
    return a.toLowerCase().indexOf(b.toLowerCase(), fromPos);
  }
  
  public void mouseDragged(MouseEvent e){}
  public void mouseMoved(MouseEvent e)
  {
    String hyperlinkText = getHyperlinkTextAtMouseEvent(e);
    if(hyperlinkText == null)
      setCursor(cdone);
    else
      setCursor(chand);
  }
  
  /**
   * Get hyperlink text from a MouseEvent position. Return null if
   * no hyperlink found.
   * @param e
   * @return
   */
  private String getHyperlinkTextAtMouseEvent(final MouseEvent e)
  {
    final Point pt = new Point(e.getX(), e.getY());
    final int pos = viewToModel(pt);
    
    try
    {
      final int viewPos = modelToView(pos).x;
      if(Math.abs(viewPos-e.getX())>5)
        return null;
    }
    catch(BadLocationException e2){}
    
    
    final int start;
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
      if(indEnd < 0)
        return null;
        
      textAtPosition = textAtPosition.substring(0, indEnd);
         
      int indStart = getStartOfLink(textAtPosition);
      if( indStart < 0 || 
         (indStart > 0 && textAtPosition.charAt(indStart-1) == '%') )
        return null;
      return textAtPosition.substring(indStart);
    }
    catch(BadLocationException e1)
    {
      //e1.printStackTrace();
    }
    return null;
  }
  
  private boolean isUseHyperlinks()
  {
    return useHyperlinks;
  }


  public void setUseHyperlinks(boolean useHyperlinks)
  {
    this.useHyperlinks = useHyperlinks;
  }
}
