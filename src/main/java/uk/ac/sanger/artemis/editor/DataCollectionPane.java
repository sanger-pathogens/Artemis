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

import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.util.StringVector;

import javax.swing.*;
import javax.swing.border.Border;
import java.util.StringTokenizer;

import java.io.StringReader;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.io.InputStreamReader;

import java.awt.*;
import java.net.URL;
import java.net.MalformedURLException;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.util.Vector;
import java.util.Enumeration;
import java.util.Hashtable;

public class DataCollectionPane extends JScrollPane
                                implements FastaListener
{

  /**  */
  private static final long serialVersionUID = 1L;
  private FastaTextPane fastaTextPane;
  private Annotation ann;
  private JDesktopPane desktop;
  private DataViewInternalFrame dataView;
  protected static String srs_url = getSrsSite();
  private static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(FastaTextPane.class);
  
  /**
  *
  * @param fastaTextPane   fasta/blast display.
  * @param ann		   annotation display.
  * @param desktop	   desktop pane.
  *
  */
  public DataCollectionPane(FastaTextPane fastaTextPane,
                            Annotation ann, JDesktopPane desktop,
                            DataViewInternalFrame dataView)
  {
    super();
    this.fastaTextPane = fastaTextPane;
    this.ann = ann;
    this.desktop  = desktop;
    this.dataView = dataView;

//  getSrsSite();
    Box bdown = Box.createVerticalBox();
    ScrollPanel scrollPanel = new ScrollPanel();
    scrollPanel.add(bdown);

    Vector hitInfoCollection = fastaTextPane.getHitCollection();
    Hashtable goHash = getIDHash(hitInfoCollection);
    getGoHash("/nfs/disk222/yeastpub/analysis/pathogen/GO/go.flat",goHash);

    setResultLines(bdown,hitInfoCollection,goHash);
    setViewportView(scrollPanel);
    setPreferredSize(new Dimension(500,300));
  }

  public static String getSrsSite()
  {
    StringVector srs = Options.getOptions().getOptionValues("srs_url");
    if(srs != null)
      srs_url = (String)srs.elementAt(0);
    else
      srs_url = "http://srs.ebi.ac.uk/srsbin/cgi-bin/";
    return srs_url;
  }

  /**
  *
  * fasta/blast results changed.
  *
  */
  public void update()
  {
    Box bdown = Box.createVerticalBox();
    ScrollPanel scrollPanel = new ScrollPanel();
    scrollPanel.add(bdown);

    Vector hitInfoCollection = fastaTextPane.getHitCollection();
    Hashtable goHash = getIDHash(hitInfoCollection);
    getGoHash("/nfs/disk222/yeastpub/analysis/pathogen/GO/go.flat",goHash);

    setResultLines(bdown,hitInfoCollection,goHash);
    setViewportView(scrollPanel);
  }


  /**
  *
  * @param bdown   		vertical Box containing all one line results.
  * @param hitInfoCollection	Collection of HitInfo from results.
  * @param goHash 		Hash of GO ID's with description.
  *
  */
  private void setResultLines(Box bdown, Vector hitInfoCollection, 
	                      Hashtable goHash)
  {
    final Vector orthoCheckBox = new Vector();
    Enumeration hitEnum = hitInfoCollection.elements();

    while(hitEnum.hasMoreElements())
    {
      Box bacross = Box.createHorizontalBox();
      final HitInfo hit = (HitInfo)hitEnum.nextElement();

// ortholog / paralog
      final JCheckBox orthoBox = new JCheckBox("ORTH");
      orthoBox.setFont(BigPane.font_sm);
      final JCheckBox paraBox  = new JCheckBox("PARA");
      paraBox.setFont(BigPane.font_sm);

      orthoBox.setMargin(new Insets(1,1,1,1));
      paraBox.setMargin(new Insets(1,1,1,1));
  
      orthoBox.setActionCommand(hit.getAcc());
      orthoCheckBox.add(orthoBox);

      bacross.add(orthoBox);
      orthoBox.setFont(BigPane.font);
      orthoBox.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          Enumeration checkEnum = orthoCheckBox.elements();
          while(checkEnum.hasMoreElements())
          {
            JCheckBox cb = (JCheckBox)checkEnum.nextElement();
            if( cb.isSelected() &&
               !cb.getActionCommand().equals(orthoBox.getActionCommand()))
            {
              cb.setSelected(false);
              ann.delete(cb.getActionCommand(),true);
            }
          }

          if(orthoBox.isSelected())
          {
            if(paraBox.isSelected())
            {
              paraBox.setSelected(false);
              ann.delete(hit.getAcc(),false);
            }

            try
            {
              setAnnotation(hit,ann,fastaTextPane.getFormat(),true);
            }
            catch(NullPointerException npe)
            {
              JOptionPane.showMessageDialog(DataCollectionPane.this,
                           "There may be a probelem retrieving "+hit.getAcc()+
                           "\nfrom SRS",
                           "Connection Error to SRS?",
                           JOptionPane.WARNING_MESSAGE);
            }
          }
          else
            ann.delete(hit.getAcc(),true);

          ann.deleteNote();
          if(BigPane.addNote.isSelected())
            dataView.updateNote();  
        }
      });

      bacross.add(paraBox);
      paraBox.setFont(BigPane.font);
      paraBox.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          if(paraBox.isSelected())
          {
            if(orthoBox.isSelected())
            {
              orthoBox.setSelected(false);
              ann.delete(hit.getAcc(),true);
            }

            try
            {
              setAnnotation(hit,ann,fastaTextPane.getFormat(),false);
            }
            catch(NullPointerException npe)
            {
              JOptionPane.showMessageDialog(DataCollectionPane.this,
                           "There may be a probelem retrieving "+hit.getAcc()+
                           "\nfrom SRS",
                           "Connection Error to SRS?",
                           JOptionPane.WARNING_MESSAGE);
            }
          }
          else
            ann.delete(hit.getAcc(),false);

          ann.deleteNote();
          if(BigPane.addNote.isSelected())
            dataView.updateNote();
        }
      });


// go to alignment
      MouseOverButton butt = new MouseOverButton();
      butt.setPreferredSize(new Dimension(12,12));
      butt.setBorderPainted(false);
      butt.setToolTipText("Go to alignment");
      butt.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          fastaTextPane.show(hit);
        }
      });
      bacross.add(butt);
      bacross.add(Box.createHorizontalStrut(2));

// SRS entry
      MouseOverButton srsRetrieve = new MouseOverButton(hit);
      srsRetrieve.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          getSRSEntry(hit,desktop);
        }
      });
     srsRetrieve.setForeground(Color.blue);
     srsRetrieve.setToolTipText("");
     srsRetrieve.setFont(BigPane.font);
     srsRetrieve.setMargin(new Insets(1,1,1,1));
     srsRetrieve.setBorderPainted(false);
     bacross.add(srsRetrieve);

//
//    String org = hit.getOrganism();
//    if(org != null)
//    {
//      JLabel orgLabel = new JLabel(org);
//      orgLabel.setFont(BigPane.font);
//      orgLabel.setForeground(Color.green);
//      bacross.add(orgLabel);
//    }

// heading
      String head = hit.getHeader();
      if(head.startsWith(hit.getID()))
        head = head.substring(hit.getID().length());        

      JLabel hiLabel = new JLabel(head);
      hiLabel.setFont(BigPane.font);
      hiLabel.setForeground(Color.black);
//    hiLabel.setForeground(Color.red);

// align button
//    final JButton selectButt = new JButton("ALIGN");
//    selectButt.setMargin(new Insets(1,1,1,1));
//    selectButt.addActionListener(new ActionListener()
//    {
//      public void actionPerformed(ActionEvent e)
//      {
//        fastaTextPane.show(hit);
//      }
//    });
//    selectButt.setFont(BigPane.font);

// retrieve srs entry
//    JButton srsButt = new JButton("->SRS");
//    srsButt.setMargin(new Insets(1,1,1,1));

//    srsButt.addActionListener(new ActionListener()
//    {
//      public void actionPerformed(ActionEvent e)
//      {
//        getSRSEntry(hit,desktop);
//      }
//    });
//    srsButt.setFont(BigPane.font);

      bacross.add(hiLabel);
//    bacross.add(selectButt);
//    bacross.add(srsButt);
     
      bacross.add(Box.createHorizontalGlue());
      bdown.add(bacross);

// go entry
      Vector gov = hit.getGO();
      if(gov != null)
      {
        Enumeration gov_enum = gov.elements();
        while(gov_enum.hasMoreElements())
        {
          final String go_id = ((String)gov_enum.nextElement()).trim();
          bacross = Box.createHorizontalBox();

          final String goLine = (String)goHash.get(go_id);
          final JCheckBox goBox = new JCheckBox();
          goBox.setSelected(false);
          goBox.addActionListener(new ActionListener()
          {
            public void actionPerformed(ActionEvent e)
            {
              if(goBox.isSelected())
              {
                String go_term = hit.getGoAssociation(go_id);
                if(go_term == null)
                { 
                  go_term = setGoAnnotation(ann,hit,go_id,goLine);
 //               hit.setGoAssociation(go_id,go_term);
                }
                else
                  ann.insert(go_term,false);
              }
              else
              {
                // possibly 2 line to delete
                ann.deleteGo(hit.getAcc(),go_id);
                ann.deleteGo(hit.getAcc(),go_id);
              }
            }
          });
          goBox.setMargin(new Insets(0,1,0,1));

          MouseOverButton goButton = new MouseOverButton("GO:"+go_id);
          goButton.setFont(BigPane.font);
          goButton.addActionListener(new ActionListener()
          {
            public void actionPerformed(ActionEvent e)
            {
              String go_cmd = "http://www.godatabase.org/cgi-bin/amigo/go.cgi?query=GO%3A"+go_id;
              BrowserControl.displayURL(go_cmd);

              if(BigPane.srsTabPane.isSelected())
              {
                try
                {
                  setUpSRSFrame(new URL(go_cmd),go_id,desktop);
                }
                catch(java.net.ConnectException connect)
                {
                  JOptionPane.showMessageDialog(DataCollectionPane.this,
                           "Cannot retrieve "+go_id+
                           "\nConnection failed to:\n"+go_cmd,
                           "Connection Error",
                           JOptionPane.WARNING_MESSAGE);
                }
                catch(Exception exp)
                {
                  exp.printStackTrace();
                }
              }

            }
          });
          goButton.setForeground(Color.blue);
          goButton.setMargin(new Insets(1,1,1,1));
          goButton.setBorderPainted(false);

          bacross.add(Box.createHorizontalStrut(10));
          bacross.add(goBox);
          bacross.add(goButton);
          goButton.setMargin(new Insets(0,1,0,1));

          JLabel goLabel = new JLabel(goLine);
          goLabel.setFont(BigPane.font);
          bacross.add(goLabel);
          bdown.add(bacross);
          bacross.add(Box.createHorizontalGlue());
        }
      }
    }
  }


  /**
  *
  * @param hit		HitInfo for a single hit.
  * @param desktop	desktop pane.
  *
  */
  private void getSRSEntry(HitInfo hit, JDesktopPane desktop)
  {
    String srscmd = "/wgetz?-e+";

    if(hit.getID() != null)
    {
      String search = hit.getID();
      srscmd = srscmd.concat("[{uniprot}-ID:"+search+"*]");
      if(hit.getAcc() != null)
        srscmd = srscmd.concat("|[{uniprot}-AccNumber:"+hit.getAcc()+"*]");

      if(srs_url.indexOf("ebi.ac.uk") > -1)
        srscmd = srscmd + "+-vn+2";
      
      if(BigPane.srsBrowser.isSelected())
        BrowserControl.displayURL(srs_url+srscmd);

      try
      {
        URL url = new URL(srs_url+srscmd);

        if(BigPane.srsTabPane.isSelected())
          setUpSRSFrame(url,search,desktop);

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
      }
      catch(java.net.ConnectException connect)
      {
        JOptionPane.showMessageDialog(DataCollectionPane.this,
                 "Cannot retrieve "+search+
                 "\nConnection failed to:\n"+srs_url+srscmd,
                 "Connection Error",
                 JOptionPane.WARNING_MESSAGE);
      }
      catch(Exception exp)
      {
        exp.printStackTrace();
      }
    }
    else
      JOptionPane.showMessageDialog(DataCollectionPane.this,
               "No ID to retrieve SRS entry!", "Missing ID",
                            JOptionPane.INFORMATION_MESSAGE);
  }
  

  /**
  *
  * @param url		URL to be displayed.
  * @param name		URL name.
  * @param desktop	desktop pane.
  *
  */
  protected static void setUpSRSFrame(URL url, String name, JDesktopPane desktop)
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
    jtab.insertTab(name, null,jsp,null,0);
    BigPane.srsFrame.setVisible(true);
  }


  /*private HitInfo findHitInfo(String acc)
  {
    Enumeration hitEnum = fastaTextPane.getHitCollection().elements();
    while(hitEnum.hasMoreElements())
    {
      HitInfo hit = (HitInfo)hitEnum.nextElement();
      if(acc.equals(hit.getID()))
        return hit;
    }

    return null;
  }*/


  /**
  *
  * @param hit		HitInfo for a single hit.
  * @param ortholog	true if ortholog is selected.
  *
  */
  private static void getzCall(HitInfo hit, boolean ortholog)
  {
    final String env[] = { "PATH=/usr/local/pubseq/bin/" };

    if(hit.getOrganism() == null ||
       hit.getDescription() == null)
    {
      String res = null;
      
      if(FastaTextPane.isRemoteMfetch())
      {
        String cmd   = 
          "mfetch -f \"org des gen\" -d uniprot -i \"acc:"+hit.getID()+"\"" ;
        uk.ac.sanger.artemis.j2ssh.SshPSUClient ssh =
          new uk.ac.sanger.artemis.j2ssh.SshPSUClient(cmd);
        ssh.run();
        res = ssh.getStdOut();
      }
      else if(!FastaTextPane.getMfetchExecutable().exists())
      {
        try
        {
          final String queryURL = DataCollectionPane.srs_url+
                                  "/wgetz?-f+acc%20org%20description%20gen+"+
                                  "[uniprot:"+hit.getAcc()+"]|[uniprot:"+hit.getID()+"]";
          logger4j.debug(queryURL);
          URL wgetz = new URL(queryURL);
          
          InputStream in = wgetz.openStream();

          BufferedReader strbuff = new BufferedReader(new InputStreamReader(in));
          StringBuffer resBuff = new StringBuffer();
          String line;
          while((line = strbuff.readLine()) != null)
            resBuff.append(line);
          strbuff.close();
          in.close();

          res = resBuff.toString();
          res = FastaTextPane.insertNewline(res, "OS ");
          res = FastaTextPane.insertNewline(res, "DE ");
          res = FastaTextPane.insertNewline(res, "GN ");
          res = FastaTextPane.insertNewline(res, "AC ");
        }
        catch(MalformedURLException e) {System.err.println(e);}
        catch(IOException e) {System.err.println(e);}
      }
      else
      {
        String cmd[]   = { "mfetch", "-f", "org des gen",
            "-d", "uniprot", "-i", "acc:"+hit.getID() };

        ExternalApplication app = new ExternalApplication(cmd,
                                                    env,null);
        res = app.getProcessStdout();
      }

      StringTokenizer tok = new StringTokenizer(res,"\n");
      while(tok.hasMoreTokens())
      {
        String token = tok.nextToken();
        String tokenline = token.substring(3).trim();
        if(tokenline.equals(""))
          continue;

        if(token.startsWith("OS "))
          hit.setOrganism(tokenline);
        else if(token.startsWith("DE "))
          hit.appendDescription(tokenline);
        else if(token.startsWith("GN "))
        {
          StringTokenizer tokGN = new StringTokenizer(tokenline,";");
          while(tokGN.hasMoreTokens())
          { 
            token = tokGN.nextToken();
            if(token.startsWith("Name="))          
              hit.setGeneName(tokenline.substring(5));
//          else
//            hit.appendDescription(token);
          }
        }
      }
    }

    if(hit.getEMBL() == null)
    {
      File fgetz = new File("/usr/local/pubseq/bin/getz");
      String res = FastaTextPane.getUniprotLinkToDatabase(fgetz, 
                   FastaTextPane.getMfetchExecutable().exists(), hit, env, "embl");
  
      int ind1 = res.indexOf("ID ");
      if(ind1 > -1)
      {
        StringTokenizer tok = new StringTokenizer(res);
        tok.nextToken();
        hit.setEMBL(tok.nextToken());
      }
    }

  }

  /**
  *
  * Sets the GO annotation for a given ID, first by querying
  * Amigo and if not successful querying SRS.
  * @param ann		annotation display
  * @param id		database id
  * @param go_id	GO id
  * @return annotation 
  *
  */
  private String setGoAnnotation(Annotation ann, HitInfo hit, 
                                 String go_id, String goLine)
  {
    String go_ann = new String("/GO=\"GOid=GO:"+go_id+
                               "; with="+hit.getDB()+":"+hit.getAcc()+
                               "; "+goLine+"\"");
  
    String prog = DataCollectionPane.class.getResource("/etc/go_associations.pl").getPath(); 
    String cmd[]   = { prog, "-assoc", hit.getAcc() };
    ExternalApplication app = new ExternalApplication(cmd,
                                                null,null);
    String res = app.getProcessStdout();
    boolean found = false;

    try
    {
      String line;
      BufferedReader buffRead = new BufferedReader(new StringReader(res));
      while((line = buffRead.readLine()) != null)
      {

        int ind1 = line.indexOf("GO:");
        int ind2 = line.indexOf(" ",ind1);
        if(ind1 > -1 && ind2 > -1)
        {
          String this_go_id = line.substring(ind1+3,ind2).trim();
          
// see http://intweb.sanger.ac.uk/help/wiki/html/Intweb/PSUEukaryoticQualifiers.html
          int ref = line.indexOf("db_xref=");
          String db_xref = null;
          if(ref > -1)
          {
            db_xref = line.substring(ref);
            line = line.substring(0,ref);
          }

          // build GO line
          StringBuffer goBuff = new StringBuffer();
          if(line.startsWith("/GO_component"))
            goBuff.append("/GO=\"aspect=component; ");
          else if(line.startsWith("/GO_process"))
            goBuff.append("/GO=\"aspect=process; ");
          else
            goBuff.append("/GO=\"aspect=function; ");
          goBuff.append("GOid=GO:"+this_go_id+"; ");

          ind1 = line.indexOf("(");
          ind2 = line.indexOf(")")+1; 
          if(ind1 > -1 && ind2 > -1)
            goBuff.append("term="+line.substring(ind1,ind2)+"; ");

          ind1 = ind2+1;
          ind2 = line.indexOf(";",ind1);
          if(ind1 > -1 && ind2 > -1)
            goBuff.append("evidence="+line.substring(ind1,ind2)+"; ");

          if(db_xref != null)
            goBuff.append(db_xref+" ");

          goBuff.append("with="+hit.getDB()+":"+hit.getAcc()+"; "); 

          if(go_id.equals(this_go_id))
            go_ann = line + "\"<br>" + goBuff.toString() +"\"";

          hit.setGoAssociation(this_go_id,line + "\"<br>" + goBuff.toString() +"\"");
          found = true;
        }
      }
    }
    catch(IOException ioe) { ioe.printStackTrace(); }

    if(!found)  // try SRS
    {
      String env[]  = { "PATH=/usr/local/pubseq/bin/" };
      String cmd2[] = { "getz", "-f", "dbxref", "[uniprot:"+hit.getAcc()+"]" };
      app = new ExternalApplication(cmd2,env,null);
      res = app.getProcessStdout();

      try
      {
        String line;
        BufferedReader buffRead = new BufferedReader(new StringReader(res));
        while((line = buffRead.readLine()) != null)
        {

          int ind1 = line.indexOf("GO:");
          int ind2 = line.indexOf(";",ind1);
          if(ind1 > -1 && ind2 > -1)
          {
            String this_go_id = line.substring(ind1+3,ind2).trim();

            line = line.substring(3).trim();
            String aspect = null;
            String term   = null;
            String evidence = null;
            ind1 = -1; 
            ind2 = -1;

            if((ind1 = line.indexOf(" F:"))> -1)
            {
              aspect = "function";
              ind2 = line.indexOf(";",ind1+4);
            }
            else if((ind1 = line.indexOf(" P:"))> -1)
            {
              aspect = "process";
              ind2 = line.indexOf(";",ind1+4);
            }
            else if((ind1 = line.indexOf(" C:"))> -1)
            {
              aspect = "component";
              ind2 = line.indexOf(";",ind1+4);
            }

            if(ind1 > -1 && ind2 > -1)
              term = line.substring(ind1+3,ind2);

            if(ind2 > -1)
              evidence = line.substring(ind2+1).trim();

            StringBuffer goBuff = new StringBuffer();

            goBuff.append("/GO_"+aspect);
            goBuff.append("=\"GO:"+this_go_id+"; ");
            goBuff.append(evidence+"; ");
            goBuff.append(hit.getDB()+":"+hit.getAcc()+";\"<br>");

            goBuff.append("/GO=\"aspect="+aspect+"; ");
            goBuff.append("GOid=GO:"+this_go_id+"; ");
            goBuff.append("term="+term+"; ");
            goBuff.append("evidence="+evidence+"; ");
            goBuff.append("db_xref= ;");
            goBuff.append("with="+hit.getDB()+":"+hit.getAcc()+";\"");

            if(go_id.equals(this_go_id))
              go_ann = goBuff.toString();

            hit.setGoAssociation(this_go_id,goBuff.toString());
            found = true;
//          break;
          }
        }
      }
      catch(IOException ioe) { ioe.printStackTrace(); }
    }
    ann.insert(go_ann,false);

    return go_ann;
  }


  /**
  *
  * @param hit 		HitInfo for a single hit.
  * @param ann		annotation display.
  * @param resultFormat	format of the results (blastp/fasta).
  * @param ortholog	true if ortholog selected.
  *
  */
  private void setAnnotation(HitInfo hit, Annotation ann,
                             String resultFormat, boolean ortholog)
  {
    getzCall(hit,ortholog);

// gene name for orthologs
    StringBuffer orthoText = new StringBuffer();
    if(ortholog)
    {
      //String geneName = hit.getGeneName();

      if(hit.getGeneName() != null &&
         !hit.getGeneName().equals(""))
        orthoText.append("<br>\n/gene=\""+hit.getGeneName()+"\"");

      if(hit.getEC_number() != null)
        orthoText.append("<br>\n/EC_number=\""+hit.getEC_number()+"\"");

      String product = hit.getDescription();

      if(product != null && !product.equals(""))
      {
        if(product.endsWith(".") ||
           product.endsWith(";"))
          product = product.substring(0,product.length()-1);
    
        int ind;
        if((product.startsWith("RecName:") || 
            product.startsWith("SubName:") ||
            product.startsWith("AltName:")) &&
            (ind = product.indexOf("="))>0)
          product = product.substring(ind+1);
        
        orthoText.append("\n<br>\n/product=\""+product.toLowerCase()+"\"");
      }
    }

//  System.out.println("ID "+hit.getID()+"\nOS "+hit.getOrganism()+
//                   "\nDE "+hit.getDescription()+"\nGN "+hit.getGeneName());
    
    StringBuffer buff = new StringBuffer();
    
    if(hit.getDB() != null)
      buff.append(" with="+hit.getDB()+":"+hit.getAcc());
    else
      buff.append(" with=UniProt:"+hit.getAcc());

    if(hit.getEMBL() != null &&
       !hit.getEMBL().equals(""))
      buff.append(" (EMBL:"+hit.getEMBL()+")");
    buff.append(";");

    if(hit.getOrganism() != null &&
       !hit.getOrganism().equals(""))
      buff.append(" "+hit.getOrganism()+";");
    if(hit.getGeneName() != null &&
       !hit.getGeneName().equals(""))
      buff.append(" "+hit.getGeneName()+";");
    if(hit.getDescription() != null)
      buff.append(" "+hit.getDescription()+";");
    if(hit.getLength() != null)
      buff.append(" length="  + hit.getLength()+";");
    if(hit.getIdentity() != null)
      buff.append(" id " + hit.getIdentity()+";");
    if(hit.getUngapped() != null)
      buff.append(" ungapped id " + hit.getUngapped()+";");
    if(hit.getEValue() != null)
      buff.append(" E()="     + hit.getEValue()+";");
    if(hit.getOverlap() != null)
      buff.append(" "+hit.getOverlap()+";");
    if(hit.getQueryRange() != null)
      buff.append(" query "   + hit.getQueryRange()+";");
    if(hit.getSubjectRange() != null)
      buff.append(" subject " + hit.getSubjectRange());
    buff.append("\"");

    ann.insert("\n/similarity=\""+resultFormat+";"+
               buff.toString()+orthoText.toString(), ortholog);
  }


  /**
  *
  * @param hitInfoCollection	Collection of HitInfo.
  *
  */
  private Hashtable getIDHash(Vector hitInfoCollection)
  {
    Hashtable goHash = new Hashtable();
 
    Enumeration hitEnum = hitInfoCollection.elements();
    while(hitEnum.hasMoreElements())
    {
      HitInfo hit = (HitInfo)hitEnum.nextElement();
      Vector gov = hit.getGO();
      if(gov != null)
      {
        Enumeration gov_enum = gov.elements();
        while(gov_enum.hasMoreElements())
        {
          String id = ((String)gov_enum.nextElement()).trim();
          goHash.put(id,"");
        }
      }
    }

    return goHash;
  }
  
  /**
  *
  * @param filename	name of the file containing the GO id and 
  *                     description.
  * @param goHash	Hashtable of GO id's and description.
  *
  */
  protected void getGoHash(String filename, Hashtable goHash)
  {
    try
    {
      String line = null;
      BufferedReader buffRead = new BufferedReader(new FileReader(filename));
      while((line = buffRead.readLine()) != null)
      {
        StringTokenizer tok = new StringTokenizer(line,"\t");
        String ID = tok.nextToken().substring(3);
        String desc  = tok.nextToken();
        
        if(tok.hasMoreTokens())
          desc = desc.concat("; "+tok.nextToken());
        
        if(goHash.containsKey(ID))
          goHash.put(ID,desc);
      }
    }
    catch(FileNotFoundException fnf) { }
    catch(IOException ioe) { ioe.printStackTrace(); }
  }
  
}

