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
import java.util.Vector;
import javax.swing.border.*;
import javax.swing.event.*;
import javax.swing.text.html.*;
import javax.swing.text.DefaultStyledDocument;
import javax.swing.text.StyledDocument;
import javax.swing.text.Document;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
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
  /** back option */
  private Vector back = new Vector();
  /** popup menu */
  private JPopupMenu popup;
  /** known qualifiers */
  private Vector qualifier = new Vector();

  public Annotation(JDesktopPane desktop)
  {
    super();

    this.desktop = desktop;
    setEditable(false);
    setContentType("text/html");
    setFont(BigPane.font);
    addHyperlinkListener(this);
  }

  
  public Annotation(URL url) throws IOException
  {
    super(url);

    setEditable(false);
    addHyperlinkListener(this);

// popup 
    addMouseListener(new PopupListener());
    popup = new JPopupMenu();
    JMenuItem backMenu = new JMenuItem("Back");
    popup.add(backMenu);
    backMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        goBack();
      }
    });
    back.add(url);
  }
  

  protected void setAnnotation(String text)
  {
//  setText("<html><body>"+text+"</html></body>");
//  reportHTML();
//  startRange = getDocument().getLength();
    String line = null;

// record qualifiers used
    try
    {
      BufferedReader buffRead = new BufferedReader(new StringReader(text));
      while((line = buffRead.readLine()) != null)
      {
        int ind = line.indexOf("=");
        if(ind > -1)
          qualifier.add(line.substring(0,ind+1).toLowerCase());
      }
    }
    catch(IOException ioe){}
    qualifier.add("/similarity=");
    qualifier.add("/gene=");

    text = getDatabaseHTML(text,"SWALL:");
    text = getDatabaseHTML(text,"UNIPROT:");
    text = getDatabaseHTML(text,"EMBL:");
    setText("<html><body>"+text+"</html></body>");
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

  
  protected String getFeatureText()
  {
    String txt = "";
    try
    {
      txt = ((HTMLDocument)getDocument()).getText(0,getDocument().getLength()).trim();
      StringBuffer buff = new StringBuffer();
      StringTokenizer tok = new StringTokenizer(txt,"/");
      int ntok = 0;

      while(tok.hasMoreTokens())
      {
        String tokTxt = "/"+tok.nextToken();

        int ind = tokTxt.indexOf("=");
        if(ntok != 0 && ind > -1 && qualifier.contains(tokTxt.substring(0,ind+1)))
          buff.append("\n"+tokTxt);
        else
          buff.append(tokTxt);

        ntok++;
      }

      txt = buff.toString();
    }
    catch(BadLocationException ble)
    {
      ble.printStackTrace();
    }

    return txt;
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
//  reportHTML();
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


  protected void setUpSRSFrame(URL url, String search)
                 throws IOException
  {
    if(BigPane.srsFrame == null)
    {
      BigPane.setUpSRSFrame((2*desktop.getHeight())/3,desktop);
      Border loweredbevel = BorderFactory.createLoweredBevelBorder();
      Border raisedbevel = BorderFactory.createRaisedBevelBorder();
      Border compound = BorderFactory.createCompoundBorder(raisedbevel,loweredbevel);

      JTextField statusField = new JTextField();
      statusField.setBorder(compound);
      statusField.setEditable(false);
      BigPane.srsFrame.getContentPane().add(statusField, BorderLayout.SOUTH);
    }

    Annotation edPane = new Annotation(url);
    JScrollPane jsp = new JScrollPane(edPane);
    JTabbedPane jtab = (JTabbedPane)BigPane.srsFrame.getContentPane().getComponent(0);
    jtab.insertTab(search, null,jsp,null,0);
    BigPane.srsFrame.setVisible(true);
  }

  protected void goBack()
  {
    if(back.size() < 2)
      return;

    try
    {
      URL url = (URL)back.get(back.size()-2);
      back.remove(back.size()-1);
      setPage(url);
    }
    catch(IOException ioe) 
    {
      ioe.printStackTrace();
    }
  }

  /**
  *
  * Method to handle hyper link events.
  * @param event        hyper link event
  *
  */
  public void hyperlinkUpdate(HyperlinkEvent event)
  {
    if(event.getEventType() == HyperlinkEvent.EventType.ACTIVATED)
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
            setUpSRSFrame(url,search);

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
          back.add(url);
        }
      }
      catch(IOException ioe)
      {
        String msg = event.getDescription();
        if(msg.length() > 50)
          msg = msg.substring(0,50)+"....";

        JOptionPane.showMessageDialog(this,
                       "Cannot reach URL:\n"+msg,
                       "Cannot Connect",
                       JOptionPane.INFORMATION_MESSAGE);
//      ioe.printStackTrace();
//      ("Can't follow link to " +
//                event.getURL().toExternalForm() );
      }

      setCursor(cdone);
    }
    else if(event.getEventType() == HyperlinkEvent.EventType.ENTERED)
    {
      try
      {
        JTextField statusField = (JTextField)BigPane.srsFrame.getContentPane().getComponent(1);
        statusField.setText(event.getDescription());
      }
      catch(Exception exp){}
      
    }
    else if(event.getEventType() == HyperlinkEvent.EventType.EXITED)
    {
      try
      {
        JTextField statusField = (JTextField)BigPane.srsFrame.getContentPane().getComponent(1);
        statusField.setText("");
      }
      catch(Exception exp){}

    }

  }

  /**
  *
  * Popup menu listener
  *
  */
  class PopupListener extends MouseAdapter
  {
    public void mousePressed(MouseEvent e)
    {
      maybeShowPopup(e);
    }

    public void mouseReleased(MouseEvent e)
    {
      maybeShowPopup(e);
    }

    private void maybeShowPopup(MouseEvent e)
    {
      if(e.isPopupTrigger() && back.size() > 1)
        popup.show(e.getComponent(),
                e.getX(), e.getY());
    }
  }

}

