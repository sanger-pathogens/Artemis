/* RunBlastAtNCBI
 *
 * created: 2009
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2009  Genome Research Limited
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
 **/
package uk.ac.sanger.artemis.components;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;

import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

import uk.ac.sanger.artemis.editor.BrowserControl;

/**
 * Used to POST data to the NCBI URLAPI (qBLAST) web services.
 */
class RunBlastAtNCBI extends Thread
{
  private String data;

  public RunBlastAtNCBI(final String data)
  {
    this.data = data;
  }
  
  public void run()
  {
    try
    {
      // Send data
      URL url = new URL("http://blast.ncbi.nlm.nih.gov/Blast.cgi");
      URLConnection conn = url.openConnection();
      conn.setDoOutput(true);
      OutputStreamWriter wr = new OutputStreamWriter(conn.getOutputStream());
      wr.write(data);
      wr.flush();

      // Get the response
      BufferedReader rd = new BufferedReader(new InputStreamReader(conn
          .getInputStream()));
      String urlResults = url + "?RID=";
      String line;
      String eta = "10";
      while ((line = rd.readLine()) != null)
      {
        int index;
        if ((index = line.indexOf("RID =")) > -1)
        {
          line = line.substring(index + 5).trim();
          urlResults = urlResults.concat(line);
          if (line.equals(""))
          {
            JOptionPane.showMessageDialog(null, line, "Blast search error",
                JOptionPane.INFORMATION_MESSAGE);
            return;
          }
        }
        else if ((index = line.indexOf("RTOE =")) > -1)
          eta = line.substring(index + 6).trim();
      }
      wr.close();
      rd.close();

      int waitTime = Integer.parseInt(eta);
      ProgressBarFrame progress = new ProgressBarFrame(waitTime, "BLAST");
      Thread.sleep(waitTime * 1000);
      BrowserControl.displayURL(urlResults + "&CMD=Get&OLD_BLAST=false");

      progress.dispose();
    }
    catch (Exception e)
    {
      e.printStackTrace();
    }
  }
  
  /**
   * Utility for appending key - value pairs to the data used to POST
   * @param data
   * @param key
   * @param value
   * @return
   * @throws UnsupportedEncodingException
   */
  private static String addData(String data, String key, String value) 
          throws UnsupportedEncodingException
  {
    data += "&" +URLEncoder.encode(key, "UTF-8") + "=" + URLEncoder.encode(value, "UTF-8");
    return data;
  }
  
  /**
   * Allow the user to define the blast options.
   * @param programName
   * @param residues
   * @return
   */
  public static String setData(String programName, String residues)
  {
    GridBagLayout layout = new GridBagLayout();
    GridBagConstraints c = new GridBagConstraints();
    JPanel dataPanel = new JPanel(layout);
    int rows = 0;
    
    c.gridx = 0;
    c.gridy = rows;
    c.anchor = GridBagConstraints.EAST;
    JTextField textField = new JTextField("nr",14);
    dataPanel.add(new JLabel("NCBI Database"), c);
    c.gridx = 1;
    dataPanel.add(textField, c);
    
    c.gridx = 0;
    c.gridy = ++rows;
    JTextField hitField = new JTextField("500",14);
    dataPanel.add(new JLabel("Number of hits to keep"), c);
    c.gridx = 1;
    dataPanel.add(hitField, c);
    
    // Gap open and gap close costs
    // Space separated float values
    // ''5 2'' for nuc-nuc, ''11 1'' for proteins, non-affine for megablast
    c.gridx = 0;
    c.gridy = ++rows;
    JTextField gapOpenField = new JTextField((programName.equals("blastn")) ? "5": "11",14);
    dataPanel.add(new JLabel("Gap open cost"), c);
    c.gridx = 1;
    dataPanel.add(gapOpenField, c);
    
    c.gridx = 0;
    c.gridy = ++rows;
    JTextField gapCloseField = new JTextField((programName.equals("blastn")) ? "2": "1",14);
    dataPanel.add(new JLabel("Gap close cost"), c);
    c.gridx = 1;
    dataPanel.add(gapCloseField, c);

    c.gridx = 0;
    c.gridy = ++rows;
    JTextField expectField = new JTextField("10.0",14);
    dataPanel.add(new JLabel("Expect value"), c);
    c.gridx = 1;
    dataPanel.add(expectField, c);
    
    c.gridx = 0;
    c.gridy = ++rows;
    JComboBox filterField = new JComboBox(new String[]{ 
        "None", "Low Complexity", "Human Repeats", "Masked" });
    dataPanel.add(new JLabel("Filter"), c);
    c.gridx = 1;
    c.anchor = GridBagConstraints.WEST;
    dataPanel.add(filterField, c);
    
    
    c.gridx = 0;
    c.gridy = ++rows;
    
    String serviceOptions[];
    if(programName.equals("blastn"))
      serviceOptions = new String[]{"plain", "megablast"};
    else
      serviceOptions = new String[]{"plain"}; // psi blast support?
    
    JComboBox serviceField = new JComboBox(serviceOptions);
    c.anchor = GridBagConstraints.EAST;
    dataPanel.add(new JLabel("Blast service"), c);
    c.gridx = 1;
    c.anchor = GridBagConstraints.WEST;
    dataPanel.add(serviceField, c);
    
    int status = JOptionPane.showConfirmDialog(null, dataPanel, 
        "Options for "+programName, JOptionPane.OK_CANCEL_OPTION);
    
    if(status != JOptionPane.OK_OPTION)
      return null;

    try
    {
      String data = URLEncoder.encode("CMD", "UTF-8") + "=" + 
             URLEncoder.encode("Put", "UTF-8");
      data = addData(data, "QUERY", residues);
      data = addData(data, "DATABASE", textField.getText());
      data = addData(data, "HITLIST_SIZE", hitField.getText());
      if(getFilterOption(filterField) != null)
        data = addData(data, "FILTER", getFilterOption(filterField));
      data = addData(data, "EXPECT", expectField.getText());
      data = addData(data, "FORMAT_TYPE", "HTML");
      data = addData(data, "PROGRAM", programName);
      data = addData(data, "CLIENT", "web");
      data = addData(data, "SERVICE", (String) serviceField.getSelectedItem());
      
      if(((String) serviceField.getSelectedItem()).equals("megablast"))
        data = addData(data, "MEGABLAST", "yes");
      else
        data = addData(data, "GAPCOSTS", gapOpenField.getText()+" "+gapCloseField.getText());
      
      return data;
    }
    catch (UnsupportedEncodingException e)
    {
      e.printStackTrace();
    }
    return null;
  }
  
  /**
   * Blast filter options.
   * @param filterField
   * @return
   */
  private static String getFilterOption(JComboBox filterField)
  {
    String sel = (String) filterField.getSelectedItem();
    if(sel.equals("None"))
      return null;
    else if(sel.equals("Low Complexity"))
      return "L";
    else if(sel.equals("Human Repeats"))
      return "R";
    else
      return "m";
  }
  
  public static void main(String args[])
  {
    String data =
      //BlastAtNCBI.setData("blastp", "ATHIEDLHNITSNQLYETYRTEKLSTSQLLLDSTVXTIDKNLSQHDQVLREDRLR");
      RunBlastAtNCBI.setData("blastn", "aggctgttttccacagatttcacagtattggttcaaatggtcaaaaattgttttaaccagt");
    RunBlastAtNCBI blast= new RunBlastAtNCBI(data);
    blast.start();
  }
}