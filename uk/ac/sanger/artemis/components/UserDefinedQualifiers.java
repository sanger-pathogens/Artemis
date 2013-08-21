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
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Enumeration;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.UIManager;
import javax.swing.event.MenuEvent;
import javax.swing.event.MenuListener;
import javax.swing.filechooser.FileFilter;

import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.components.genebuilder.JExtendedComboBox;
import uk.ac.sanger.artemis.components.genebuilder.cv.DatePanel;
import uk.ac.sanger.artemis.editor.MultiLineToolTipUI;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.FileDocument;
import uk.ac.sanger.artemis.util.LinePushBackReader;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.util.URLDocument;

/**
 * Load and make available qualifiers from a list in a file 
 * of tag=value pairs and from OBO files.
 */
class UserDefinedQualifiers extends JFrame
{
  private static final long serialVersionUID = 1L;
  private static org.apache.log4j.Logger logger4j = 
      org.apache.log4j.Logger.getLogger(UserDefinedQualifiers.class);
  
  /* collection of loaded files */
  private Vector<QualifierList> obos = new Vector<QualifierList>();
  
  /* component containing list of available values */
  private Box qualifierBox = Box.createHorizontalBox();
  
  /* listener for adding qualifiers **/
  private ActionListener addListener;
  
  /* text area in the feature editor to add qualifiers to */
  private QualifierTextArea qualifier_text_area;
  private Selection selection;
  private JPanel mainPanel;
  
  /* line pattern for tag=value */
  private static Pattern TAG_VALUE_PATTERN = Pattern.compile("^/?(\\w+) ?= ?([^\\n\\r]+)$");
  
  UserDefinedQualifiers()
  {
    super("Qualifier List");
    mainPanel = (JPanel) getContentPane();
    MultiLineToolTipUI.initialize();
    read();
    createComponentToAddCvTerm();
    createMenus();
  }
  
  /**
   * Menu creation
   */
  private void createMenus()
  {
    final JMenuBar menuBar = new JMenuBar();
    final JMenu fileMenu = new JMenu("File");
    menuBar.add(fileMenu);
    final JMenuItem importFile = new JMenuItem("Import qualifiers from file ...");
    importFile.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent arg0)
      {
        importFromFile();
      }
    });
    fileMenu.add(importFile);
    
    final JMenuItem importURL = new JMenuItem("Import qualifiers from URL ...");
    importURL.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent arg0)
      {
        importFromUrl();
      }
    });
    fileMenu.add(importURL);
    
    final JMenuItem pasteQualifiers = new JMenuItem("Paste qualifiers ...");
    pasteQualifiers.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent arg0)
      {
        final QualifierTextArea qualifierText = new QualifierTextArea();
        int status = JOptionPane.showConfirmDialog(UserDefinedQualifiers.this, 
            qualifierText, "Paste in qualifier lists in the format: name = value",
            JOptionPane.OK_CANCEL_OPTION);
        
        if(status == JOptionPane.OK_OPTION)
        {
          QualifierList obo = new QualifierList();
          obo.qualifiers = loadText(new ByteArrayInputStream(qualifierText.getText().getBytes()));
          obos.add(obo);
          createComponentToAddCvTerm();
          mainPanel.revalidate();
          
          writeQualifierList(qualifierText.getText());
        }

      }
    });
    fileMenu.add(pasteQualifiers);
    fileMenu.addSeparator();
    
    final JMenu removeFile = new JMenu("Qualifier source(s)");
    removeFile.addMenuListener(new MenuListener(){
      public void menuCanceled(MenuEvent e){}
      public void menuDeselected(MenuEvent e){}

      public void menuSelected(MenuEvent e)
      {
        removeFile.removeAll();
        for(final QualifierList obo: obos)
        {
          final JCheckBoxMenuItem oboMenu =
              new JCheckBoxMenuItem(obo.doc.getName(), obo.isAcive);
          removeFile.add(oboMenu);
          oboMenu.addActionListener(new ActionListener(){
            public void actionPerformed(ActionEvent e)
            {
              obo.isAcive = oboMenu.isSelected();
              createComponentToAddCvTerm();
              mainPanel.revalidate();
              mainPanel.repaint();
            }
          });
        }
      }
    });
    fileMenu.add(removeFile);
    setJMenuBar(menuBar);
  }

  /**
   * Add components to the panel for adding CvTerm's to
   * the annotation.
   */
  private void createComponentToAddCvTerm()
  {
    mainPanel.removeAll();
    mainPanel.setLayout(new GridBagLayout());
    final GridBagConstraints c = new GridBagConstraints();

    final JCheckBox ignoreCase = new JCheckBox("Ignore case",true);
    final Vector<String> names = getQualiferNames();
    final JComboBox nameCombo = new JComboBox(names);   
    final JButton addButton = new JButton("ADD");
    final JCheckBox addToAllSelected = new JCheckBox("add to all selected features", false);
    final JTextField keyWord = new JTextField(45);
    final Dimension d  = new Dimension(500, ignoreCase.getPreferredSize().height);
    final Dimension d2 = new Dimension(300, ignoreCase.getPreferredSize().height);
    
    int row = 0;
    c.gridx = 0;
    c.gridy = row;
    c.anchor = GridBagConstraints.EAST;
    
    if(names.size() < 1)
    {
      c.gridwidth = 2;
      mainPanel.add(new JLabel("Import a list of qualifiers from the File menu."), c);
      c.gridwidth = 1;
      ignoreCase.setEnabled(false);
      nameCombo.setEnabled(false);
      addButton.setEnabled(false);
      addToAllSelected.setEnabled(false);
      keyWord.setEnabled(false);
    }
    else
    {
      mainPanel.add(new JLabel("Name: "), c);
      c.gridx = 1;
      c.anchor = GridBagConstraints.WEST;
      mainPanel.add(nameCombo, c);
    }
    
    c.gridx = 2;
    mainPanel.add(Box.createHorizontalStrut(150), c);
    
    // keyword
    keyWord.setSelectionStart(0);
    keyWord.setSelectionEnd(keyWord.getText().length());
    keyWord.setSelectedTextColor(Color.blue);
    keyWord.setMinimumSize(d2);
    keyWord.addActionListener(new ActionListener(){
      // carry out search when enter key is pressed
      public void actionPerformed(ActionEvent event)
      {
        searchQualifiers(keyWord.getText(), nameCombo, 
            ignoreCase.isSelected(), addButton, addToAllSelected, d);
      }
    });
    c.gridy = ++row;
    c.gridx = 0;
    c.anchor = GridBagConstraints.EAST;
    mainPanel.add(new JLabel("Keywords: "),c);
    c.gridx = 1;
    c.anchor = GridBagConstraints.WEST;
    mainPanel.add(keyWord,c);

    c.gridy = ++row;
    mainPanel.add(ignoreCase,c);
     
    // search button
    c.gridx = 0;
    c.gridy = ++row;
    c.anchor = GridBagConstraints.WEST;
    final JButton search = new JButton("SEARCH");
    search.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent arg0)
      {
        searchQualifiers(keyWord.getText(), nameCombo, 
            ignoreCase.isSelected(), addButton, addToAllSelected, d);
      }
    });
    mainPanel.add(search,c);
    
    c.gridx = 0;
    c.gridy = ++row;
    c.gridwidth = 3;
    qualifierBox.add(Box.createVerticalStrut(25));
    mainPanel.add(qualifierBox, c);

    c.gridx = 0;
    c.gridy = ++row;
    c.gridwidth = 1;
    mainPanel.add(addButton, c);  

    c.gridx = 1;
    final JButton closeButton = new JButton("CLOSE");
    mainPanel.add(closeButton, c);
    closeButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent arg0)
      {
        dispose();
      }
    });
    
    c.gridx = 0;
    c.gridy = ++row;
    c.gridwidth = 2;
    mainPanel.add(addToAllSelected,c);
    
    c.gridy = ++row;
    mainPanel.add(Box.createVerticalStrut(15), c);
  }

  private StringVector searchOboQualifiers(final QualifierList obo, 
                                           final String qName, 
                                           final String keyWord,  
                                           final boolean ignoreCase)
  {
    final Qualifier q = obo.qualifiers.getQualifierByName(qName);
    if(q == null)
      return null;
    final StringVector values = q.getValues();
    
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
    return values;
  }
  
  /**
   * Search the qualifiers for the selected qualifier name and value.
   */
  private void searchQualifiers(final String keyWord, 
                                final JComboBox nameCombo, 
                                final boolean ignoreCase,
                                final JButton addButton,
                                final JCheckBox addToAllSelected,
                                final Dimension d)
  {
    final String qName = (String) nameCombo.getSelectedItem();
    final StringVector results = new StringVector();
    for(QualifierList obo: obos)
    {
      if(!obo.isAcive)
        continue;
      StringVector oboResult = searchOboQualifiers(obo, qName, keyWord, ignoreCase);
      if(oboResult != null)
        results.add(oboResult);
    }
    
    final JExtendedComboBox valuesList = new JExtendedComboBox(results, true);
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
        String qval = (String) valuesList.getSelectedItem();
        if(qName.equals("biological_process"))
          qval = "/GO=aspect=P;"+qval+";date="+DatePanel.getDate();
        else if(qName.equals("molecular_function"))
          qval = "/GO=aspect=F;"+qval+";date="+DatePanel.getDate();
        else if(qName.equals("cellular_component"))
          qval = "/GO=aspect=C;"+qval+";date="+DatePanel.getDate();
        else
          qval = "/"+qName+"="+qval;
        
        if(addToAllSelected.isSelected())
        {
          final FeatureVector features = selection.getAllFeatures();
          if(!Options.isBlackBeltMode() && features.size() > 1)
          {
            int status = JOptionPane.showConfirmDialog(UserDefinedQualifiers.this, 
              "Add "+qName+" to "+features.size()+" selected features?",
              "Add Qualifier To Selected Features", 
              JOptionPane.OK_CANCEL_OPTION);
            if(status == JOptionPane.CANCEL_OPTION)
              return;
          }
          final int idx = qval.indexOf("=");
          final Qualifier q = new Qualifier(qval.substring(1, idx), qval.substring( idx+1 ));
          
          try
          {
            for(int i=0; i<features.size(); i++)
              features.elementAt(i).addQualifierValues(q);
          }
          catch (ReadOnlyException e)
          {
            e.printStackTrace();
          }
          catch (EntryInformationException e)
          {
            e.printStackTrace();
          }
        }
        else
          qualifier_text_area.append(qval+"\n");
      }
    }; 
    addButton.addActionListener(addListener);

    mainPanel.revalidate();
  }
  
  /**
   * Return a list of the qualifier names.
   * @return
   */
  private Vector<String> getQualiferNames()
  {
    final Vector<String> names = new Vector<String>();
    for(QualifierList obo: obos)
    {
      if(!obo.isAcive)
        continue;
      QualifierVector qualifiers = obo.qualifiers;
      for(Qualifier q: qualifiers)
        if(!names.contains(q.getName()))
           names.add(q.getName());
    }
    return names;
  }
  
  /**
   * Read qualifier lists
   */
  private void read()
  {
    URL url = ClassLoader.getSystemResource("etc/artemis.qualifiers");
    if(url != null)
    {
      Document doc = new URLDocument(
        ClassLoader.getSystemResource("etc/artemis.qualifiers"));
      createObo(doc);
    }
    
    final String [] qualFiles = 
    {
      "artemis.qualifiers",
      System.getProperty("user.home") + File.separator + ".artemis.qualifiers"
    };
    
    for(String fileName: qualFiles)
    {
      Document doc = new FileDocument(new File(fileName));
      if(doc.readable()) 
        createObo(doc);
    }
  }
  
  /**
   * Load the qualifiers from an input stream.
   * @param ins
   * @throws IOException 
   */
  private QualifierVector loadQualifiers(final Document doc) throws IOException
  {
    QualifierVector qualifiers = loadObo(doc.getInputStream());
    if(qualifiers.size() == 0)
      qualifiers = loadText(doc.getInputStream());
    return qualifiers;
  }
  
  /**
   * Load in qualifiers from an input stream that are in OBO format.
   * @param ins
   * @return
   */
  private QualifierVector loadObo(final InputStream ins)
  {
    final LinePushBackReader reader = new LinePushBackReader(
        new BufferedReader(new InputStreamReader(ins)));
    final QualifierVector qualifiers = new QualifierVector();
    try
    {
      String line;
      while((line = reader.readLine()) != null) 
      {
        if(line.startsWith("#"))
          continue;

        if(startOfStanza(line))
        {
          String val = null;
          String name = null;
          String id = null;
          while((line = reader.readLine()) != null &&
                !startOfStanza(line) &&
                 line.indexOf("import:") == -1)
          {
            if(line.startsWith("namespace:"))
              name = line.substring(10).trim();
            else if(line.startsWith("name:"))
              val = line.substring(5).trim();
            else if(line.startsWith("id:"))
              id = line.substring(3).trim();
          }
          
          if(line != null)
            reader.pushBack(line);
          
          if(id != null)
          {
            if(  name != null && 
                (name.equals("biological_process") || 
                 name.equals("molecular_function") || 
                 name.equals("cellular_component")))
              val = "term="+val+"; GOid="+id;
            else
              val = "term="+val+"; id="+id;
          }
          if(name != null)
            qualifiers.addQualifierValues(new Qualifier(name, val));
        }
        else if(line.startsWith("import:"))
        {
          // import OBO file from a URL
          final URL url = new URL(line.substring(7).trim());
          createObo(new URLDocument(url));
        }
      }
      reader.close();
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
    return qualifiers;
  }
  
  /**
   * Each unit of the OBO file starts with one of the three supported stanza 
   * types: [Term], [Typedef], and [Instance]. Parsers/serializers will round-trip 
   * (successfully load and save) unrecognized stanzas. 
   * @param line
   * @return
   */
  private boolean startOfStanza(final String line)
  {
    return line.startsWith("[");
  }
  
  /**
   * Load in qualifiers from an input stream that are in the format:
   * name=value
   * @param ins
   * @return
   */
  private QualifierVector loadText(final InputStream ins)
  {
    final BufferedReader reader = new BufferedReader(new InputStreamReader(ins));
    final QualifierVector qualifiers = new QualifierVector();
    try
    {
      String line;
      while((line = reader.readLine()) != null) 
      {
        Matcher m = TAG_VALUE_PATTERN.matcher(line);
        if(m.matches())
          qualifiers.addQualifierValues(new Qualifier(m.group(1), m.group(2)));
      }
      reader.close();
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
    return qualifiers;
  }
  
  /**
   * Import OBO
   */
  private void importFromFile()
  {
    final StickyFileChooser fc = new StickyFileChooser();
    //final JCheckBox autoImport = new JCheckBox("Automatically import", true);
    //fc.setAccessory(autoImport);
    
    OboFileFilter oboFilter = new OboFileFilter(new String[]{"obo", "OBO"}, "OBO files");
    fc.addChoosableFileFilter(oboFilter);
    fc.setFileFilter(oboFilter);
    int status = fc.showOpenDialog(UserDefinedQualifiers.this);
    if(status != StickyFileChooser.APPROVE_OPTION)
      return;
    
    final FileDocument doc = new FileDocument(fc.getSelectedFile());
    if(!doc.readable())
    {
      JOptionPane.showMessageDialog(UserDefinedQualifiers.this, 
          "Cannot read file "+doc.getName(), 
          "Problem reading file", JOptionPane.WARNING_MESSAGE);
      return;
    }
    
    createObo(doc);
    
    createComponentToAddCvTerm();
    mainPanel.revalidate();
    mainPanel.repaint();
    //if(autoImport.isSelected())
    writeQualifierList(fc.getSelectedFile().getAbsolutePath());
  }
  
  private void importFromUrl()
  {
    String s = (String)JOptionPane.showInputDialog(
        UserDefinedQualifiers.this,
        "URL:", "Import Qualifiers",
        JOptionPane.PLAIN_MESSAGE,
        null, null,
        "http://www.geneontology.org/GO_slims/goslim_generic.obo");
    if(s == null || s.trim().equals("") || s.trim().equals("http://"))
      return;

    try
    {
      String urlStr = s.trim();
      createObo(new URLDocument(new URL(urlStr)));
      createComponentToAddCvTerm();
      mainPanel.revalidate();
      mainPanel.repaint();
      
      writeQualifierList(urlStr);
    }
    catch (MalformedURLException e)
    {
      e.printStackTrace();
    }
  }
  
  private void createObo(Document doc)
  {
    try
    {
      logger4j.debug("Reading qualifiers from: "+doc.getName());
      final QualifierList obo = new QualifierList();
      obo.doc = doc;
      obo.qualifiers = loadQualifiers(doc);
      obos.add(obo);
    }
    catch (IOException e)
    {
      logger4j.error(e.getMessage());
      JOptionPane.showMessageDialog(UserDefinedQualifiers.this, 
          "Problem Loading:\n"+doc.getName()+"\n"+
          e.getMessage(), "Problem Loading", 
          JOptionPane.ERROR_MESSAGE);
    }
  }


  /**
   * Write or append to artemis qualifiers list
   * @param loc
   * @param isOBO
   */
  private static void writeQualifierList(final String loc)
  {
    int status = JOptionPane.showConfirmDialog(
        null, 
        "Automatically import:\n"+loc+"\nbetween session?", 
        "Import Option", JOptionPane.YES_NO_OPTION);
    if(status != JOptionPane.YES_OPTION)
      return;
    
    final boolean isOBO;
    if(loc.toLowerCase().endsWith(".obo"))
      isOBO = true;
    else
    {
      File f = new File(loc);
      if(f.exists())
        isOBO = true;
      else
        isOBO = false;
    }
    
    final String qualPath = System.getProperty("user.home") + File.separator + ".artemis.qualifiers";
    try
    {
      final FileDocument doc = new FileDocument(new File(qualPath));
      if(doc.readable())
      {
        final BufferedReader reader = 
            new BufferedReader(new InputStreamReader(doc.getInputStream()));
        try
        {
          String line;
          while((line = reader.readLine()) != null)
          {
            if(line.startsWith("import: file://"+loc))
            {
              System.out.println("Import line already exists :\n"+line);
              return;
            }
            else if(line.startsWith("import: "+loc))
            {
              System.out.println("Import line already exists :\n"+line);
              return;
            }
          }
        }
        finally
        {
          reader.close();
        }
      }

      final File qualFile = new File(qualPath);
      final boolean exists = qualFile.exists();
      final BufferedWriter bufferedwriter = new BufferedWriter(
          new FileWriter(qualFile, true));
      
      if(!exists)
      {
        bufferedwriter.append("# User defined list of qualifiers in the format of:\n");
        bufferedwriter.append("# name=value\n");
        bufferedwriter.append("# The qualifiers are separated by lines.\n");
        bufferedwriter.append("# OBO files can be used by defining import statements:\n");
        bufferedwriter.append("#import: http://www.geneontology.org/GO_slims/goslim_generic.obo\n\n");
      }
      
      if(isOBO)
      {
        if(loc.startsWith("http://"))
          bufferedwriter.append("\nimport: "+loc+"\n");
        else
          bufferedwriter.append("\nimport: file://"+loc+"\n");
      }
      else
        bufferedwriter.append(loc+"\n");

      bufferedwriter.close();
    }
    catch (FileNotFoundException filenotfoundexception)
    {
      System.err.println(qualPath+" read error");
    }
    catch (IOException e)
    {
      System.err.println(qualPath+" i/o error");
    }
  }
  
  private class OboFileFilter extends FileFilter
  {
    private String[] suffix;
    private String descr;
    OboFileFilter(final String[] suffix, final String descr)
    {
      this.suffix = suffix;
      this.descr = descr;
    }
    
    public boolean accept(final File file)
    {
      if(file.isDirectory())
        return true;

      for(String s: suffix)
      {
        if(file.getName().endsWith("." + s) )
          return true;
      }
      return false;
    }

    public String getDescription()
    {
      return descr;
    }
  }
  
  protected void setQualifierTextArea(final QualifierTextArea qualifier_text_area)
  {
    this.qualifier_text_area = qualifier_text_area;
  }
  
  protected void setSelection(final Selection selection)
  {
    this.selection = selection;
  }
  
  private class ComboMouseListener extends MouseAdapter 
  {
    private JComboBox cb;
    ComboMouseListener(JComboBox cb)
    {
      this.cb = cb;
    }
    
    public void mouseEntered(MouseEvent me) 
    {
      if(cb.getSelectedItem() != null)
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


  private class QualifierList
  {
    Document doc;
    QualifierVector qualifiers;
    boolean isAcive = true;
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

    final JFrame f = new UserDefinedQualifiers();
    f.pack();
    f.setVisible(true);
  }
}