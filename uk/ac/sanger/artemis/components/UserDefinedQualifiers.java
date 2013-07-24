/* UserDefinedQualifier
 * This file is part of Artemis
 *
 * Copyright(C) 2013  Genome Research Limited
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
package uk.ac.sanger.artemis.components;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.Enumeration;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;

import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.components.genebuilder.JExtendedComboBox;
import uk.ac.sanger.artemis.editor.MultiLineToolTipUI;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.FileDocument;
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.util.URLDocument;

class UserDefinedQualifier extends JPanel
{
  private static final long serialVersionUID = 1L;
  private static org.apache.log4j.Logger logger4j = 
      org.apache.log4j.Logger.getLogger(UserDefinedQualifier.class);
  private QualifierVector qualifiers = new QualifierVector();
  private Box qualifierBox = Box.createHorizontalBox();
  private JButton addButton = new JButton("ADD");
  private ActionListener addListener;
  private QualifierTextArea qualifier_text_area;
  
  UserDefinedQualifier()
  {
    MultiLineToolTipUI.initialize();
    read();
    createComponentToAddCvTerm();
  }

  public void setQualifierTextArea(QualifierTextArea qualifier_text_area)
  {
    this.qualifier_text_area = qualifier_text_area;
  }
  
  public JFrame getFrame()
  {
    return (JFrame) SwingUtilities.getWindowAncestor(this);
  }

  /**
   * Add components to the panel for adding CvTerm's to
   * the annotation.
   */
  private void createComponentToAddCvTerm()
  {
    setLayout(new GridBagLayout());
    final GridBagConstraints c = new GridBagConstraints();

    final JCheckBox ignoreCase = new JCheckBox("Ignore case",true);
    final JComboBox nameCombo = new JComboBox(getQualiferNames());
    final Dimension d  = new Dimension(500, nameCombo.getPreferredSize().height);
    final Dimension d2 = new Dimension(300, nameCombo.getPreferredSize().height);
    
    int row = 0;
    
    c.gridx = 0;
    c.gridy = row;
    c.anchor = GridBagConstraints.EAST;
    add(new JLabel("Name: "), c);
    
    c.gridx = 1;
    c.weightx = 0.5d;
    c.anchor = GridBagConstraints.WEST;
    add(nameCombo, c);
    c.weightx = 0.d;
    
    c.gridx = 2;
    c.weightx = 0.5d;
    add(Box.createHorizontalStrut(150), c);
    
    // keyword
    final JTextField keyWord = new JTextField(45);
    keyWord.setSelectionStart(0);
    keyWord.setSelectionEnd(keyWord.getText().length());
    keyWord.setSelectedTextColor(Color.blue);
    keyWord.setMinimumSize(d2);
    keyWord.addActionListener(new ActionListener(){
      // carry out search when enter key is pressed
      public void actionPerformed(ActionEvent event)
      {
        searchQualifiers(keyWord.getText(), nameCombo, ignoreCase.isSelected(), d);
      }
    });
    c.gridy = ++row;
    c.gridx = 0;
    c.anchor = GridBagConstraints.EAST;
    add(new JLabel("Keywords: "),c);
    c.gridx = 1;
    c.anchor = GridBagConstraints.WEST;
    add(keyWord,c);

    c.gridy = ++row;
    add(ignoreCase,c);
     
    // search button
    c.gridx = 0;
    c.gridy = ++row;
    c.anchor = GridBagConstraints.WEST;
    final JButton search = new JButton("SEARCH");
    search.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent arg0)
      {
        searchQualifiers(keyWord.getText(), nameCombo, ignoreCase.isSelected(), d);
      }
    });
    add(search,c);
    
    c.gridx = 0;
    c.gridy = ++row;
    c.gridwidth = 3;
    qualifierBox.add(Box.createVerticalStrut(25));
    add(qualifierBox, c);

    c.gridx = 0;
    c.gridy = ++row;
    c.gridwidth = 1;
    add(addButton, c);
    
    c.gridx = 1;
    final JButton closeButton = new JButton("CLOSE");
    add(closeButton, c);
    closeButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent arg0)
      {
        final JFrame topFrame = 
            (JFrame) SwingUtilities.getWindowAncestor(UserDefinedQualifier.this);
        topFrame.dispose();
      }
    });
    
    c.gridy = ++row;
    add(Box.createVerticalStrut(15), c);
  }


  /**
   * Search the qualifiers for the selected qualifier name and value.
   */
  private void searchQualifiers(final String keyWord, 
                                final JComboBox nameCombo, 
                                final boolean ignoreCase,
                                final Dimension d)
  {
    final String qName = (String) nameCombo.getSelectedItem();
    final StringVector values = qualifiers.getQualifierByName(qName).getValues();
    
    if(keyWord != null && !keyWord.trim().equals(""))
    {
      final String keyWordLC = keyWord.toLowerCase();
      final StringVector tmp = values.copy();
      for(String val: tmp)
      {
        if(ignoreCase)
        {
          if(val.toLowerCase().indexOf(keyWordLC) == -1)
            values.remove(val);
        }
        else if(val.indexOf(keyWord) == -1)
          values.remove(val);
      }
    }
    
    final JExtendedComboBox valuesList = new JExtendedComboBox(values, true);
    valuesList.getEditor().getEditorComponent().addMouseListener(
        new ComboMouseListener(valuesList));
    
    valuesList.setPreferredSize(d);
    valuesList.setMaximumSize(d);
    
    if(valuesList.getSelectedItem() != null &&
       valuesList.getSelectedItem().equals("") &&
       valuesList.getItemCount() > 2)
      valuesList.setSelectedIndex(1);
    
    qualifierBox.removeAll();
    qualifierBox.add(valuesList);
    
    if(addListener != null)
      addButton.removeActionListener(addListener);
    addListener = new ActionListener(){
      public void actionPerformed(ActionEvent arg0)
      {
        qualifier_text_area.append("/"+qName+"="+valuesList.getSelectedItem()+"\n");
      }
    }; 
    addButton.addActionListener(addListener);

    revalidate();
  }
  
  /**
   * Return a list of the qualifier names.
   * @return
   */
  private Vector<String> getQualiferNames()
  {
    final Vector<String> names = new Vector<String>();
    for(Qualifier q: qualifiers)
      names.add(q.getName());
    return names;
  }
  
  private void read()
  {
    InputStream ins = 
        this.getClass().getClassLoader().getResourceAsStream("etc/artemis.qualifiers");
    try
    {
      if(ins != null)
      {
        loadQualifiers(ins);
        ins.close();
      }
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
    catch (NullPointerException e2) {}
    
    final String [] qualFiles = 
    {
      "artemis.qualifiers",
      System.getProperty("user.home") + File.separator + ".artemis.qualifiers"
    };
    
    for(String fileName: qualFiles)
    {
      final Document doc = new FileDocument(new File(fileName));
      if(doc.readable()) 
      {
        try
        {
          ins = doc.getInputStream();
          loadQualifiers(ins);
          ins.close();
        }
        catch (IOException e)
        {
          e.printStackTrace();
        }
        logger4j.debug("Reading qualifiers from: "+fileName);
      }
    }
  }
  
  /**
   * Load the qualifiers from an input stream.
   * @param ins
   */
  private void loadQualifiers(final InputStream ins)
  {
    final BufferedReader br = new BufferedReader(new InputStreamReader(ins));
    try
    {
      String line;
      while((line = br.readLine()) != null) 
      {
        int idx = line.indexOf("name:");
        if(idx > -1)
        {
          String val = line.substring(idx+5).trim();
          while((line = br.readLine()) != null &&  
                (idx = line.indexOf("namespace:")) > -1)
          {
            String name = line.substring(idx+10).trim();
            logger4j.debug(name + "=" +val);
            qualifiers.addQualifierValues(new Qualifier(name, val));
          }
        }
        else if((idx = line.indexOf("import:")) > -1)
        {
          // import OBO file from a URL
          final URL url = new URL(line.substring(idx+7).trim());
          final Document doc = new URLDocument(url);
          loadQualifiers(doc.getInputStream());
        }
      }
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
  }
  
  class ComboMouseListener extends MouseAdapter 
  {
    private JComboBox cb;
    ComboMouseListener(JComboBox cb)
    {
      this.cb = cb;
    }
    
    public void mouseEntered(MouseEvent me) 
    {
      cb.setToolTipText(getWrappedStr((String) cb.getSelectedItem()));
    }
    
    private String getWrappedStr(String s)
    {
      final StringBuilder buff = new StringBuilder();
      final int lineLen = 60;
      for(int k=0; k<s.length(); k+=lineLen)
      {
        int end = k + lineLen;
        if(end > s.length())
          end = s.length();
        buff.append ( s.substring(k,end) ).append("\n");
      }
      return buff.toString();
    }
  }
  
  public static void main(String args[])
  {
    final javax.swing.plaf.FontUIResource font_ui_resource =
        Options.getOptions().getFontUIResource();

    final Enumeration<Object> keys = UIManager.getDefaults().keys();
    while(keys.hasMoreElements()) 
    {
      Object key = keys.nextElement();
      Object value = UIManager.get(key);
      if(value instanceof javax.swing.plaf.FontUIResource) 
        UIManager.put(key, font_ui_resource);
    }

    final JFrame f = new JFrame("User Defined Qualifiers");
    f.getContentPane().add(new UserDefinedQualifier());
    f.pack();
    f.setVisible(true);
  }
}