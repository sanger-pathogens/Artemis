/* 
 *
 * created: Wed Aug 3 2004
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2000  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or(at your option) any later version.
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
 */

package uk.ac.sanger.artemis.editor;

import java.util.StringTokenizer;
import javax.swing.event.*;
import javax.swing.text.html.*;
import javax.swing.text.DefaultStyledDocument;
import javax.swing.text.StyledDocument;
import javax.swing.text.Document;
import javax.swing.*;
import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.Cursor;
import java.awt.Dimension;
import java.io.BufferedReader;
import java.io.StringReader;
import java.io.IOException;
import javax.swing.text.BadLocationException;
import java.net.URL;

public class Annotation extends JEditorPane
                        implements HyperlinkListener
{
  private int startRange;
  private int endRange;
  /** busy cursor */
  private Cursor cbusy = new Cursor(Cursor.WAIT_CURSOR);
  /** done cursor */
  private Cursor cdone = new Cursor(Cursor.DEFAULT_CURSOR);
  /** desktop pane */
  private JDesktopPane desktop = null;

  public Annotation(JDesktopPane desktop)
  {
    super();

    this.desktop = desktop;
    setEditable(false);
    setContentType("text/html");
    Font font = new Font("Monospaced",Font.PLAIN,12);
    setFont(font);
    addHyperlinkListener(this);
  }

  
  public Annotation(URL url) throws IOException
  {
    super(url);

    setEditable(false);
    addHyperlinkListener(this);
  }


  protected void setAnnotation(String text)
  {
    setText("<html><body>"+text+"</html></body>");
    reportHTML();
    startRange = getDocument().getLength();
  }

  protected void reportHTML()
  {
    try
    {
      String txt = ((HTMLDocument)getDocument()).getText(0,getDocument().getLength());
      System.out.println("TXT:\n"+txt);
      System.out.println("\nHTML:\n"+getText()+"\n\n");
    }
    catch(BadLocationException ble)
    {
      ble.printStackTrace();
    }
  }

  protected void insert(String s, boolean ortholog)
  {
    s = getDatabaseHTML(s,"SWALL:");
    s = getDatabaseHTML(s,"UNIPROT:");
    s = getDatabaseHTML(s,"EMBL:");

    int ind = s.indexOf("/gene");
    
    Document doc = getDocument();
    int offset = doc.getLength();
    if(ortholog)
      offset = startRange;

    insert(s,offset);
    reportHTML();
  }


  protected void insert(String s, int offset)
  {
    try
    {
      HTMLEditorKit edKit = (HTMLEditorKit)getEditorKit();
//    ((HTMLDocument)getDocument()).insertString(offset,"\n",null);
      edKit.insertHTML((HTMLDocument)getDocument(),offset,"<BR>\n"+s,0,0,HTML.Tag.BR);
//    HTMLEditorKit.InsertHTMLTextAction("similarity",s,HTML.Tag.BODY,HTML.Tag.P);
    }
    catch(BadLocationException ble)
    {
      ble.printStackTrace();
    }
    catch(Exception exp)
    {
      exp.printStackTrace();
    }
  }


  private String getDatabaseHTML(String s, String db)
  {
    int ind = s.indexOf(db);
    if(ind>-1)
    {
      String startStr = s.substring(0,ind);
      int ind2 = s.indexOf(" ",ind);
      int ind3 = s.indexOf(")",ind);     
      if(ind3>-1 && ind3<ind2)
        ind2 = ind3;
      ind3 = s.indexOf(";",ind);
      if(ind3>-1 && ind3<ind2)
        ind2 = ind3;

      String midStr = s.substring(ind,ind2);
      String endStr = s.substring(ind2);

      String srscmd = "http://srs.sanger.ac.uk/srsbin/cgi-bin/wgetz?-e+" +
                      "["+midStr+"]";

      s = startStr + "<a href=\""+srscmd+"\">" +
          midStr   + "</a>" + endStr;
    }
    return s;
  }

  private void replaceRange(String newStr,int start,int end)
  {
    HTMLDocument doc = (HTMLDocument)getDocument();

    try
    {
      doc.remove(start,(end-start));
      insert(newStr,start);
    }
    catch(BadLocationException ble)
    {
      ble.printStackTrace();
    }
    setDocument(doc);
  }

  /**
  *
  * Deletes the annotation line that contains an ID.
  *
  */
  protected void delete(String id, boolean ortholog)
  {
//  try
//  {
//    reportHTML();

//    String txt = ((HTMLDocument)getDocument()).getText(0,getDocument().getLength());
//    String line = null;
//    int eol = 0;
//    int len = 0;
//    BufferedReader buffRead = new BufferedReader(new StringReader(txt));
//    while((line = buffRead.readLine()) != null)
//    {
//      len = line.length()+1;
//      if(line.indexOf("SWALL:"+id) > -1)
//      {
//        len += eol;
//       
//        if(ortholog)
//        {
//          line = buffRead.readLine();
//          if(line != null && line.startsWith("/gene="))
//            len += line.length();
//        }

//        if(len > txt.length())
//          len = txt.length();

//        replaceRange("",eol-1,len);
//        reportHTML();
//        return;
//      }
//      eol += len;
//    }
      String txt = getText();
      int indID = txt.indexOf("SWALL:"+id);
      if(indID == -1)
        indID = txt.indexOf("UNIPROT:"+id);

      int ind1 = 0;
      int ind2 = 0;
      
      while((ind2 = txt.indexOf("<br>",ind1)) > -1)
      {
        if(ind2 < indID)
          ind1 = ind2+1; 
        else
          break;
      }
      
      ind2 = txt.indexOf("<br>",indID);

      // if ortholog then delete gene line as well
      if(ortholog)
        ind2 = txt.indexOf("<br>",ind2+4);
      
      if(ind2 == -1)
        ind2 = txt.length();

      setText(txt.substring(0,ind1-1)+txt.substring(ind2));
//  }
//  catch(BadLocationException ble) { ble.printStackTrace(); }
//  catch(IOException ioe) { ioe.printStackTrace(); }
    
  }


  /**
  *
  * Method to handle hyper link events.
  * @param event        hyper link event
  *
  */
  public void hyperlinkUpdate(HyperlinkEvent event)
  {
    if (event.getEventType() == HyperlinkEvent.EventType.ACTIVATED)
    {
      setCursor(cbusy);
      try
      {
        URL url = event.getURL();
        
        int ind1 = event.getDescription().indexOf("[");
        int ind2 = event.getDescription().lastIndexOf("]");

        String search = "";
        if(ind1 > -1 && ind2 > -1)
          search = event.getDescription().substring(ind1+1,ind2);
        
        if(desktop != null)
        {
          if(BigPane.srsTabPane.isSelected())
          {
            if(BigPane.srsFrame == null)
              BigPane.setUpSRSFrame((2*desktop.getHeight())/3,desktop);
            Annotation edPane = new Annotation(url);
            JScrollPane jsp = new JScrollPane(edPane);
            JTabbedPane jtab = (JTabbedPane)BigPane.srsFrame.getContentPane().getComponent(0);
            jtab.insertTab(search, null,jsp,null,0);
            BigPane.srsFrame.setVisible(true);
          }

          if(BigPane.srsWin.isSelected())
          {
            int hgt = (2*desktop.getHeight())/3;
            Annotation edPane = new Annotation(url);
            JScrollPane jsp = new JScrollPane(edPane);
            JInternalFrame jif = new JInternalFrame("SRS "+search,
                                                   true, //resizable
                                                   true, //closable
                                                   true, //maximizable
                                                   true);//iconifiable);
            JMenuBar menuBar = new JMenuBar();
            menuBar.add(new CommonMenu(jif));
            jif.setJMenuBar(menuBar);
            jif.getContentPane().add(jsp);
            jif.setLocation(0,0);
            jif.setSize(800,hgt);
            jif.setVisible(true);
            desktop.add(jif);
          }
 
          if(BigPane.srsBrowser.isSelected())
            BrowserControl.displayURL(event.getDescription());
        }
        else
        {
          setPage(url); 
        }
      }
      catch(IOException ioe)
      {
        setCursor(cdone);
        ioe.printStackTrace();
//      ("Can't follow link to " +
//                event.getURL().toExternalForm() );
      }

      setCursor(cdone);
    }
  }


}

