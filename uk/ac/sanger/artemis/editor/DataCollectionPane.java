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

import javax.swing.*;
import java.util.StringTokenizer;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.Dimension;
import java.awt.Color;
import java.awt.Insets;
import java.net.URL;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.util.Vector;
import java.util.Enumeration;
import java.util.Hashtable;

public class DataCollectionPane extends JScrollPane
{

  public DataCollectionPane(Vector hitInfoCollection, String dataFile,
                            final FastaTextPane fastaTextPane,
                            final Annotation ann, final JDesktopPane desktop)
  {
    super();

    Font font = new Font("Monospaced",Font.PLAIN,12);

    Box bdown = Box.createVerticalBox();
    ScrollPanel scrollPanel = new ScrollPanel();
    scrollPanel.add(bdown);

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

    getGoHash("/nfs/disk222/yeastpub/analysis/pathogen/GO/go.flat",goHash);

    final Vector orthoCheckBox = new Vector();
    hitEnum = hitInfoCollection.elements();

    while(hitEnum.hasMoreElements())
    {
      Box bacross = Box.createHorizontalBox();
      final HitInfo hit = (HitInfo)hitEnum.nextElement();

// ortholog / paralog
      final JCheckBox orthoBox = new JCheckBox("ORTH");
      final JCheckBox paraBox  = new JCheckBox("PARA");

      orthoBox.setMargin(new Insets(2,1,1,1));
      paraBox.setMargin(new Insets(2,1,1,1));
  
      orthoBox.setActionCommand(hit.getID());
      orthoCheckBox.add(orthoBox);

      bacross.add(orthoBox);
      orthoBox.setFont(font);
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
              ann.delete(hit.getID(),false);
            }

            setAnnotation(hit,ann,fastaTextPane.getFormat(),true);
          }
          else
            ann.delete(hit.getID(),true);           
        }
      });

      bacross.add(paraBox);
      paraBox.setFont(font);
      paraBox.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          if(paraBox.isSelected())
          {
            if(orthoBox.isSelected())
            {
              orthoBox.setSelected(false);
              ann.delete(hit.getID(),true);
            }
           
            setAnnotation(hit,ann,fastaTextPane.getFormat(),false);
          }
          else
            ann.delete(hit.getID(),false);
        }
      });


// heading
      JLabel hiLabel = new JLabel(hit.getHeader());
      hiLabel.setFont(font);
      hiLabel.setForeground(Color.RED);

// align button
      final JButton selectButt = new JButton("ALIGN");
      selectButt.setMargin(new Insets(1,1,1,1));
      selectButt.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          fastaTextPane.show(hit);
        }
      });
      selectButt.setFont(font);

// retrieve srs entry
      JButton srsButt = new JButton("->SRS");
      srsButt.setMargin(new Insets(1,1,1,1));

      srsButt.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          String srscmd = "srs.sanger.ac.uk/srsbin/cgi-bin/wgetz?-e+";
          if(hit.getID() != null)
          {
            String search = hit.getID();
            srscmd = srscmd.concat("[{swall}-ID:"+search+"*]");
            if(hit.getAcc() != null)
              srscmd = srscmd.concat("|[{swall}-AccNumber:"+hit.getAcc()+"*]");  

            if(BigPane.srsBrowser.isSelected())
              BrowserControl.displayURL(srscmd);

            try
            {
              URL url = new URL("http://"+srscmd);
              if(BigPane.srsTabPane.isSelected())
              {
                if(BigPane.srsFrame == null)
                {
                  BigPane.setUpSRSFrame((2*desktop.getHeight())/3,desktop);
//                int hgt = (2*desktop.getHeight())/3;
//                BigPane.srsFrame = new JInternalFrame("SRS",
//                                                   true, //resizable
//                                                   true, //closable
//                                                   true, //maximizable
//                                                   true);//iconifiable
//                JTabbedPane jtab = new JTabbedPane();
//                BigPane.srsFrame.getContentPane().add(jtab);
//                BigPane.srsFrame.setLocation(0,0);
//                BigPane.srsFrame.setSize(800,hgt);
  
//                JMenuBar menuBar = new JMenuBar();
//                menuBar.add(new CommonMenu(BigPane.srsFrame));
//                BigPane.srsFrame.setJMenuBar(menuBar);

//                desktop.add(BigPane.srsFrame);
                }
                Annotation edPane = new Annotation(url);
                JScrollPane jsp = new JScrollPane(edPane);
                JTabbedPane jtab = (JTabbedPane)BigPane.srsFrame.getContentPane().getComponent(0);
                jtab.insertTab(search,null,jsp,null,0);
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
      });
      srsButt.setFont(font);

      bacross.add(hiLabel);
      bacross.add(selectButt);
      bacross.add(srsButt);
     
      bacross.add(Box.createHorizontalGlue());
      bdown.add(bacross);

      Vector gov = hit.getGO();
      if(gov != null)
      {
        Enumeration gov_enum = gov.elements();
        while(gov_enum.hasMoreElements())
        {
          final String go_id = ((String)gov_enum.nextElement()).trim();
          bacross = Box.createHorizontalBox();
          JButton goButton = new JButton("GO:"+go_id);
          goButton.setFont(font);
          goButton.addActionListener(new ActionListener()
          {
            public void actionPerformed(ActionEvent e)
            {
              String go_cmd = "http://www.godatabase.org/cgi-bin/amigo/go.cgi?query=GO%3A"+go_id;
              BrowserControl.displayURL(go_cmd);
            }
          });
          bacross.add(Box.createHorizontalStrut(10));
          bacross.add(goButton);
          goButton.setMargin(new Insets(0,1,0,1));

          JLabel goLabel = new JLabel((String)goHash.get(go_id));
          goLabel.setFont(font);
          bacross.add(goLabel);
          bdown.add(bacross);
          bacross.add(Box.createHorizontalGlue());
        }
      }
    }

    setViewportView(scrollPanel);
    setPreferredSize(new Dimension(500,300));
  }


  protected static void getzCall(HitInfo hit, boolean ortholog)
  {
    String env[] = { "PATH=/usr/local/pubseq/bin/" };
    if(hit.getOrganism() == null ||
       hit.getDescription()== null)
    {

      String cmd[]   = { "getz", "-f", "org description",
                         "[libs={swall}-acc:"+
                          hit.getAcc()+"]|[libs={swall}-id:"+hit.getID()+"]" };
      ExternalApplication app = new ExternalApplication(cmd,
                                                  env,null);
      String res = app.getProcessStdout();

      int ind1 = res.indexOf("OS ");
      int ind2 = res.indexOf("\n");

      if(ind1 > -1)
        hit.setOrganism(res.substring(ind1+3,ind2).trim());
   
      ind1 = res.indexOf("DE ");
      if(ind1 > -1)
        hit.setDescription(res.substring(ind1+3).trim());
    }

    if(hit.getEMBL() == null)
    {
      String cmd2[]   = { "getz", "-f", "id",
                 "[libs={swall}-id:"+hit.getID()+"]>EMBL" };
      ExternalApplication app = new ExternalApplication(cmd2,env,null);
      String res = app.getProcessStdout();
  
      int ind1 = res.indexOf("ID ");
      if(ind1 > -1)
      {
        StringTokenizer tok = new StringTokenizer(res);
        tok.nextToken();
        hit.setEMBL(tok.nextToken());
      }
    }

    if(ortholog)
    {
      String geneName = hit.getGeneName();
      if(geneName == null)
      {
        String cmd3[]   = { "getz", "-f", "gen",
                       "[libs={swall}-acc:"+
                        hit.getAcc()+"]|[libs={swall}-id:"+hit.getID()+"]" };
        ExternalApplication app = new ExternalApplication(cmd3,
                                                env,null);
        geneName = app.getProcessStdout();
      }

      int ind1 = geneName.indexOf("GN ");
      if(ind1 > -1)
      {
        geneName = geneName.substring(ind1+3).trim();
        if(geneName.startsWith("Name="))
        {
          geneName = geneName.substring(5);
          ind1 = geneName.indexOf(";");
          if(ind1 > -1)
            geneName = geneName.substring(0,ind1);
        }

        if(!geneName.toLowerCase().startsWith("orderedlocusnames="))
          hit.setGeneName(geneName);
      }
    }
  }

  private void setAnnotation(HitInfo hit, Annotation ann,
                             String resultFormat, boolean ortholog)
  {
    getzCall(hit,ortholog);

// gene name for orthologs
    String orthoText = "";
    if(ortholog)
    {
      String geneName = hit.getGeneName();
      if(hit.getGeneName() != null)
        orthoText = "<br>\n/gene=\""+hit.getGeneName()+"\"";
    }

    StringBuffer buff = new StringBuffer();
    
    if(hit.getDB() != null)
      buff.append(hit.getDB()+":"+hit.getID());
    else
      buff.append(" UNIPROT:"+hit.getID());

    if(hit.getEMBL() != null)
      buff.append(" (EMBL:"+hit.getEMBL()+")");
    buff.append(";");

    if(hit.getOrganism() != null)
      buff.append(" "+hit.getOrganism()+";");
    if(hit.getLength() != null)
      buff.append(" length="  + hit.getLength()+";");
    if(hit.getUngapped() != null)
      buff.append("ungapped id=" + hit.getUngapped()+";");
    if(hit.getEValue() != null)
      buff.append(" E()="     + hit.getEValue()+";");
    if(hit.getOverlap() != null)
      buff.append(" "+hit.getOverlap()+";");
    if(hit.getQueryRange() != null)
      buff.append("query "   + hit.getQueryRange()+";");
    if(hit.getSubjectRange() != null)
      buff.append("subject " + hit.getSubjectRange());
    buff.append("\"");

    ann.insert("\n/similarity=\""+resultFormat+";"+
               buff.toString()+orthoText, ortholog);

//  ann.insert("\n/similarity=\""+resultFormat+"; SWALL:"+hit.getID()+
//                   " (EMBL:"+hit.getEMBL()+"); "+
//                   hit.getOrganism()+"; "+
//                   hit.getDescription()+"; "+
//                   "length="  + hit.getLength()+"; "+
//                   "id="      + hit.getIdentity()+"; "+
//                   "ungapped id=" + hit.getUngapped()+"; "+
//                   "E()="     + hit.getEValue()+"; "+
//                   hit.getOverlap()+"; "+
//                   "query "   + hit.getQueryRange()+"; "+
//                   "subject " + hit.getSubjectRange()+"\""+
//                   orthoText, ortholog);
  }

  public void getGoHash(String filename, Hashtable goHash)
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
    catch(IOException ioe) { ioe.printStackTrace(); }
  }


}

