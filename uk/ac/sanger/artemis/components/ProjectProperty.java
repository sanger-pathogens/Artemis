/* ProjectProperty
 * This file is part of Artemis
 *
 * Copyright(C) 2012  Genome Research Limited
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

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.net.URISyntaxException;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Properties;
import java.util.Vector;
import java.util.Map.Entry;
import java.util.Set;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.JToolBar;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.border.Border;
import javax.swing.border.TitledBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.components.database.DatabaseEntrySource;
import uk.ac.sanger.artemis.components.database.DatabaseJPanel;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.FileDocument;
import uk.ac.sanger.artemis.util.DatabaseLocationParser;

/**
 * Project file management system using a properties file.
 * 
 * Example of the syntax for defining a project in the property file:
 * project.Pknowlsei.sequence = Pknowlsei:Pk_strainH_chr01
 * project.Pknowlsei.chado =  genedb-db.sanger.ac.uk:5432/snapshot?genedb_ro
 * project.Pknowlsei.title = Pknowlsei
 */
public class ProjectProperty extends JFrame
{
  private static final long serialVersionUID = 1L;
  private static HashMap<String, HashMap<String, String>> centralProjects;
  private static HashMap<String, HashMap<String, String>> userProjects;
  private Splash splash;
  private DatabaseEntrySource entry_source;
  
  private final static int REFERENCE = 1;
  private final static int ANNOTATION = 2;
  private final static int NEXT_GEN_DATA = 3;
  private final static int CHADO = 4;
  private final static int USERPLOT = 5;
  private final static int LOGUSERPLOT = 6;
  private final static int VCF = 7;
  private static org.apache.log4j.Logger logger4j = 
      org.apache.log4j.Logger.getLogger(ProjectProperty.class);
  
  private final static String[] TYPES = 
    { "title", "sequence", "annotation", "bam", "vcf", "userplot", "log_userplot", "chado" };
  
  public ProjectProperty()
  {
    this(null);
    
    final javax.swing.plaf.FontUIResource font_ui_resource =
        Options.getOptions().getFontUIResource();
    java.util.Enumeration<Object> keys = UIManager.getDefaults().keys();
    while(keys.hasMoreElements()) 
    {
      Object key = keys.nextElement();
      Object value = UIManager.get(key);
      if(value instanceof javax.swing.plaf.FontUIResource) 
        UIManager.put(key, font_ui_resource);
    }
  }
  
  public ProjectProperty(Splash splash)
  {
    super("Project File Manager");
    this.splash = splash;
    InputStream ins = 
        this.getClass().getClassLoader().getResourceAsStream("etc/project.properties");
    
    try
    {
      logger4j.debug("Reading properties from: "+
          this.getClass().getClassLoader().getResource("etc/project.properties").toURI());
    }
    catch (URISyntaxException e1){}
    catch (NullPointerException e2) {}

    final Properties projectProps = new Properties();
    try
    {
      projectProps.load(ins);
      ins.close();
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
    catch (NullPointerException e2) {}
    
    centralProjects = getProjectMap(projectProps);
    projectProps.clear();
    
    final String [] propFiles = 
    {
      "project.properties",
      System.getProperty("user.home") + File.separator + ".artemis.project.properties"
    };
    
    for(int i=0; i<propFiles.length; i++)
    {
      final Document doc =
          new FileDocument(new File(propFiles[i]));
      if(doc.readable()) 
      {
        try
        {
          ins = doc.getInputStream();
          projectProps.load(ins);
        }
        catch (IOException e)
        {
          e.printStackTrace();
        }
        
        logger4j.debug("Reading properties from: "+propFiles[i]);
      }
    }

    userProjects = getProjectMap(projectProps);
    createProjectViewer((JPanel) getContentPane());
    pack();
    setVisible(true);
  }
  
  private void createProjectViewer(JPanel panel)
  {
    final DefaultListModel model = new DefaultListModel();
    
    final JList projectList = new JList(model);
    final JScrollPane jspList = new JScrollPane(projectList);
    Object[] items = centralProjects.keySet().toArray();
    Arrays.sort(items);
    for (int i=0; i<items.length; i++)
      model.add(i, items[i]);
    
    items = userProjects.keySet().toArray();
    Arrays.sort(items);
    for (int i=0; i<items.length; i++)
      model.add(i, items[i]);

    final Box yBox = Box.createVerticalBox();
    final JScrollPane jspProp = new JScrollPane(yBox);
    final LaunchActionListener listener = new LaunchActionListener();
    projectList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    projectList.setVisibleRowCount(-1);
    panel.add(jspList, BorderLayout.WEST);
    panel.add(jspProp, BorderLayout.CENTER);

    // Add / remove project buttons
    
    final JToolBar toolBar = new JToolBar();
    panel.add(toolBar, BorderLayout.PAGE_START);
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    panel.setPreferredSize(new Dimension((int)(screen.width/2.5),screen.height/3));
    panel.setBackground(Color.WHITE);
    final JButton addProjectButton = new JButton("+");
    addProjectButton.setOpaque(false);
    Font font = addProjectButton.getFont().deriveFont(Font.BOLD).deriveFont(14.f);
    addProjectButton.setFont(font);
    addProjectButton.setToolTipText("ADD PROJECT");
    addProjectButton.setForeground(new Color(35, 149, 35));
    addProjectButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent arg0)
      {
        addProject(projectList);
      }
    });
    toolBar.add(addProjectButton);
    
    final JButton removeProjectButton = new JButton("-");
    removeProjectButton.setOpaque(false);
    removeProjectButton.setFont(font);
    removeProjectButton.setToolTipText("REMOVE PROJECT");
    removeProjectButton.setForeground(new Color(149, 35, 35));
    removeProjectButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent arg0)
      {
        removeProject(projectList, yBox, listener);
      }
    });
    toolBar.add(removeProjectButton);
    
    final JButton saveProperties = new JButton("SAVE");
    saveProperties.setFont(font.deriveFont(Font.PLAIN));
    saveProperties.setToolTipText("SAVE PROJECT PROPERTIES");
    saveProperties.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent arg0)
      {
        writeProperties();
      }
    });
    toolBar.add(Box.createHorizontalGlue());
    toolBar.add(saveProperties);
    
    final JButton openArt = new JButton("OPEN");
    openArt.addActionListener(listener);
    final JButton closeButton = new JButton("CLOSE");
    closeButton .addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent arg0)
      {
        dispose();
      }
    });
    
    Box xBox = Box.createHorizontalBox();
    xBox.add(openArt);
    xBox.add(closeButton);
    xBox.add(Box.createHorizontalGlue());
    
    panel.add(xBox, BorderLayout.SOUTH);
    
    final GridBagConstraints c = new GridBagConstraints();
    c.anchor = GridBagConstraints.NORTHWEST;
    c.fill = GridBagConstraints.VERTICAL;
    c.ipadx = 10;
    c.ipady = 10;
    
    
    projectList.addListSelectionListener(new ListSelectionListener()
    {
      public void valueChanged(ListSelectionEvent e)
      {
        if(e.getValueIsAdjusting() == false &&
           projectList.getSelectedIndex() > -1) 
          refreshProperties(projectList, yBox, listener);
      }
    });
  }
  
  private void addProject(final JList projectList)
  {
    DefaultListModel model = (DefaultListModel) projectList.getModel();
    
    String projName = 
        JOptionPane.showInputDialog(ProjectProperty.this,
            "Project Name", "New Project", JOptionPane.QUESTION_MESSAGE);
    if(projName == null)
      return;

    if(model.contains(projName))
    {
      JOptionPane.showMessageDialog(ProjectProperty.this, 
          projName+" is already a project. Please provide a unique project name.", 
          "Project Name", JOptionPane.WARNING_MESSAGE);
      return;
    }
    
    final HashMap<String, String> hMap = new HashMap<String, String>();
    hMap.put("sequence", "");
    userProjects.put(projName, hMap);
    model.add(model.getSize(), projName);
    projectList.repaint();
    projectList.setSelectedIndex(model.getSize()-1);
  }
  
  private void removeProject(final JList projectList, final Box yBox, final LaunchActionListener listener)
  {
    if(projectList.getSelectedValue() == null)
    {
      JOptionPane.showMessageDialog(ProjectProperty.this, 
          "Select a project from the list to be removed.", 
          "Remove", JOptionPane.INFORMATION_MESSAGE);
      return;
    }
    DefaultListModel model = (DefaultListModel) projectList.getModel();
    int status = JOptionPane.showConfirmDialog(
        ProjectProperty.this, "Remove "+projectList.getSelectedValue()+"?",
        "Remove Project", JOptionPane.YES_NO_OPTION);
    if(status != JOptionPane.YES_OPTION)
      return;
    userProjects.remove(projectList.getSelectedValue());
    model.remove(projectList.getSelectedIndex());
    projectList.repaint();
    yBox.removeAll();
    yBox.repaint();
    listener.setSettings(null);
  }
  
  /**
   * Refresh components in the properties panel.
   * @param projectList
   * @param yBox
   * @param listener
   */
  private void refreshProperties(final JList projectList, 
                                 final Box yBox, 
                                 final LaunchActionListener listener)
  {
    yBox.removeAll();

    final HashMap<String, String> projProps;
    if(centralProjects.containsKey(projectList.getSelectedValue()))
      projProps = centralProjects.get(projectList.getSelectedValue());
    else
      projProps = userProjects.get(projectList.getSelectedValue());

    final HashMap<Integer, Vector<JTextField>> settings = new HashMap<Integer, Vector<JTextField>>();
    // order the keys
    Object keys[] = projProps.keySet().toArray();
    Arrays.sort(keys, new TypeComparator());

    for(final Object key: keys)
    {
      final String keyStr = (String) key;
      final Vector<JTextField> vText = new Vector<JTextField>();
      
      Border lineBorder = BorderFactory.createLineBorder(Color.DARK_GRAY);
      TitledBorder title = BorderFactory.createTitledBorder(
          lineBorder, keyStr);
      //title.setTitlePosition(TitledBorder.LEFT);

      final JPanel propPanel = new JPanel(new GridBagLayout());
      GridBagConstraints c = new GridBagConstraints();
      c.gridy = 0;
      
      c.fill = GridBagConstraints.BOTH;
      propPanel.setBorder(title);
      propPanel.setBackground(Color.WHITE);
      
      //
      final Vector<JCheckBox> checkBoxes = new Vector<JCheckBox>();
      final JButton toggle = new JButton("Toggle");
      toggle.setToolTipText("toggle "+keyStr+" selection");
      toggle.addActionListener(new ActionListener(){
        public void actionPerformed(ActionEvent arg0)
        {
           for(JCheckBox cb: checkBoxes)
             cb.setSelected(!cb.isSelected());
        }
      });
               
      final Vector<String> anns = splitLine(projProps.get(keyStr).trim());
      for (int i=0; i<anns.size(); i++)
      {
        c.gridx = 0;
        addProperyToPanel(projectList, propPanel, vText, c, i, anns.get(i), projProps, keyStr, yBox, listener, checkBoxes);
      }
      
      if (!keyStr.equals("title") && !keyStr.equals("chado"))
      {
        Box xBox = Box.createHorizontalBox();
        final JButton selectButton = new JButton(
            keyStr.startsWith("seq") ? "Select " : "Add file");
        selectButton.addActionListener(new ActionListener()
        {
          public void actionPerformed(ActionEvent e)
          {
            StickyFileChooser fileChooser = new StickyFileChooser();
            int status = fileChooser.showOpenDialog(ProjectProperty.this);

            if(status == StickyFileChooser.APPROVE_OPTION)
            {
              if(keyStr.startsWith("seq"))
                vText.get(0).setText(fileChooser.getSelectedFile().getAbsolutePath());
              else
              {
                projProps.put(keyStr, projProps.get(keyStr)+" "+
                            fileChooser.getSelectedFile().getAbsolutePath());
                refreshProperties(projectList, yBox, listener);
              }
            }
          }
        });
        c.gridy = c.gridy+1;
        
        c.gridx = 0;
        xBox.add(selectButton);
        xBox.add(Box.createHorizontalGlue());
        if(checkBoxes.size() > 1) // add toggle option
          xBox.add(toggle, c);
        
        c.gridwidth = 2;
        propPanel.add(xBox, c);
        c.gridwidth = 1;
      }
      
      yBox.add(propPanel);
      
      if(keyStr.startsWith("seq"))
        settings.put(ProjectProperty.REFERENCE, vText);
      else if(keyStr.equals("annotation"))
        settings.put(ProjectProperty.ANNOTATION, vText);
      else if(keyStr.equals("bam"))
        settings.put(ProjectProperty.NEXT_GEN_DATA, vText);
      else if(keyStr.equals("vcf") || keyStr.equals("bcf"))
        settings.put(ProjectProperty.VCF, vText);
      else if(keyStr.equals("chado"))
        settings.put(ProjectProperty.CHADO, vText);
      else if(keyStr.equals("userplot"))
        settings.put(ProjectProperty.USERPLOT, vText);
      else if(keyStr.equals("log_userplot"))
        settings.put(ProjectProperty.LOGUSERPLOT, vText);
    }
    
    // ADD property
    Box xBox = Box.createHorizontalBox();
    final JButton addPropertyButton = new JButton("NEW PROPERTY");
    final JComboBox propertyList = new JComboBox(TYPES);
    xBox.add(addPropertyButton);
    xBox.add(propertyList);
    xBox.add(Box.createHorizontalGlue());
    yBox.add(xBox);
    addPropertyButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent arg0)
      {
        String key = (String) propertyList.getSelectedItem();
        if(!projProps.containsKey(key))
        {
          projProps.put(key, "");
          refreshProperties(projectList, yBox, listener);
        }
      }
    });
    
    //
    yBox.add(Box.createVerticalGlue());
    yBox.revalidate();
    yBox.repaint();
    
    listener.setSettings(settings);
  }

  private Vector<String> splitLine(String line)
  {
    Vector<String> parts = new Vector<String>();
    int index = line.indexOf(" ");
    if(index < 0)
      parts.add(line);
    else
    {
      int startIndex = 0;
      line = line.replaceAll("\\s+", " ");
      while((index = line.indexOf(" ", startIndex)) > -1)
      {
        if(line.charAt(index-1) == '\\')
        {
          startIndex = index+1;
          continue;
        }
        parts.add(line.substring(0, index));
        line = line.substring(index+1);
        startIndex = 0;
      }
      parts.add(line);
    }
    
    return parts;
  }
  
  private void addProperyToPanel(final JList projectList,
                                 final JPanel propPanel,
                                 final Vector<JTextField> vText,
                                 final GridBagConstraints c,
                                 final int index,
                                 final String ann, 
                                 final HashMap<String, String> projProps, 
                                 final String key,
                                 final Box yBox, 
                                 final LaunchActionListener listener,
                                 final Vector<JCheckBox> cbs)
  {
    final JTextField qta = new JTextField(67);
    vText.add(qta);
    
    if(key.equals("title"))
      qta.setText(removeSpaceEscape(ann));
    else
      qta.setText(ann);
    qta.setBorder(BorderFactory.createLineBorder(Color.LIGHT_GRAY));
    qta.getDocument().addDocumentListener(new DocumentListener()
    {
      private void update()
      {
        final String anns[];
        if(key.equals("title")) // only takes one value
          anns = new String[]{ escapeSpace(projProps.get(key).trim()) };
        else
          anns = projProps.get(key).trim().split("\\s+");
        String value = "";
        for(int i=0;i<anns.length;i++)
        {
          if(i == index)
            value += " "+qta.getText();
          else
            value += " "+anns[i];
        }
        
        if(index == anns.length)
          value += " "+qta.getText();
        
        projProps.put(key, value);
      }
      
      public void changedUpdate(DocumentEvent e)
      {
        update();
      }
      public void insertUpdate(DocumentEvent e)
      {
        update();
      }
      public void removeUpdate(DocumentEvent e)
      {
        update();
      }
    });
    
    // REMOVE PROPERTY
    Box xButtons = Box.createHorizontalBox();
    final JButton removeProperty = new JButton("X");
    removeProperty.setOpaque(false);
    Font font = removeProperty.getFont().deriveFont(Font.BOLD);
    removeProperty.setFont(font);
    removeProperty.setToolTipText("REMOVE");
    removeProperty.setForeground(new Color(139, 35, 35));
    removeProperty.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        int status = JOptionPane.showConfirmDialog(
            ProjectProperty.this, "Remove "+key+"?",
            "Remove", JOptionPane.YES_NO_OPTION);
        if(status != JOptionPane.YES_OPTION)
          return;
        final String anns[] = projProps.get(key).trim().split("\\s+");
        String value = "";
        for(int i=0;i<anns.length;i++)
        {
          if(i != index)
            value += " "+anns[i];
        }
        if(value.equals(""))
          projProps.remove(key);
        else
          projProps.put(key, value.trim());
        refreshProperties(projectList, yBox, listener);
      }
    });
    xButtons.add(removeProperty);
    //
    if(!key.equals("title") && !key.startsWith("seq") && !key.equals("chado"))
    {
      final JCheckBox useProperty = new JCheckBox("",true);
      useProperty.addItemListener(new ItemListener(){
        public void itemStateChanged(ItemEvent arg0)
        {
          qta.setEnabled(useProperty.isSelected());
        }
      });

      xButtons.add(useProperty);
      cbs.add(useProperty);
    }
    xButtons.add(Box.createHorizontalGlue());
    
    c.gridy = c.gridy+1;
    propPanel.add(qta, c);
    c.gridx = c.gridx+1;
    propPanel.add(xButtons, c);
  }

  /**
   * Create a project hash of the properties.
   * @param projectProps
   * @return
   */
  private HashMap<String, HashMap<String, String>> getProjectMap(final Properties projectProps)
  {
    final HashMap<String, HashMap<String, String>> projects =
        new HashMap<String, HashMap<String, String>>();
    
    for (Entry<Object, Object> propItem : projectProps.entrySet())
    {
      String key = (String) propItem.getKey();
      String value = (String) propItem.getValue();

      if (key.startsWith("project."))
      {
        key = key.substring(8);
        int ind = key.indexOf(".");
        if (ind > -1)
        {
          String projName = key.substring(0, ind);
          key = key.substring(ind + 1);
          final HashMap<String, String> thisProj;
          if (projects.containsKey(projName))
            thisProj = projects.get(projName);
          else
            thisProj = new HashMap<String, String>();
          thisProj.put(key, value);

          projects.put(projName, thisProj);
        }
      }
    }
    return projects;
  }
  
  /**
  * Write or re-write properties and insert/update the user.dir property
  * @param jemProp      properties file
  * @param uHome        user working directory
  */
  protected static void writeProperties()
  {
    if(userProjects == null)
      return;
    final String prop = System.getProperty("user.home") + File.separator + ".artemis.project.properties";
    File propFile = new File(prop);
    try
    {
      if(userProjects.size() > 0)
      {
        propFile.renameTo(new File(propFile.getAbsolutePath()+".bak"));
        final BufferedWriter bufferedwriter = new BufferedWriter(new FileWriter(propFile));
        for (String project: userProjects.keySet())
        {
          bufferedwriter.write("#");
          bufferedwriter.newLine();
          
          HashMap<String, String> projProps = userProjects.get(project);
          for(final String key: projProps.keySet())
          {
            final String val;
            if(key.equals("title"))
              val = escapeSpace(projProps.get(key).trim());
            else
              val = projProps.get(key).trim().replaceAll("\\s{2,}", " ");

            bufferedwriter.write("project."+project+"."+key+"="+val );
            bufferedwriter.newLine();
          }
          
          // unfortunately Properties.store() adds a timestamp as a comment
          /*myProps.clear();
          HashMap<String, String> projProps = userProjects.get(project);
          for(final String key: projProps.keySet())
            myProps.setProperty("project."+project+"."+key, 
                projProps.get(key).trim().replaceAll("\\s+", " "));
          myProps.store(bufferedwriter, null);*/
        }

        bufferedwriter.close();
      }
      else
        propFile.delete();
    }
    catch (FileNotFoundException filenotfoundexception)
    {
      System.err.println(prop+" read error");
    }
    catch (IOException e)
    {
      System.err.println(prop+" i/o error");
    }
  }
  
  /**
   * Open a database entry
   * @param splash
   * @param featureName gene or sequence ID
   */
  private void openDatabase(final Splash splash, final String featureName)
  {
    String loc = System.getProperty("chado");
    DatabaseLocationParser dlp = new DatabaseLocationParser(loc);
    
    
//    int idx;
//    if((idx = loc.indexOf("?")) > -1 && loc.indexOf("?user=") == -1)
//      loc = loc.substring(0, idx+1) + "user=" + loc.substring(idx+1);

    if( entry_source == null || 
       !entry_source.getLocation().endsWith(loc) )
    {
      entry_source = new DatabaseEntrySource();
      boolean promptUser = true;
      if(System.getProperty("read_only") != null)
      {
        promptUser = false;
        entry_source.setReadOnly(true);
      }
      if(!entry_source.setLocation(promptUser))
        return;
    }
    DatabaseJPanel.getEntryEditFromDatabase(
        entry_source, splash, ProjectProperty.this, featureName);
  }
  
  
  /**
   * Escape the spaces with a double backslash (i.e. '\\ ').
   * @param s
   * @return
   */
  private static String escapeSpace(String s)
  {
    s = removeSpaceEscape(s).replace(" ", "\\\\ ");
    return s;
  }
  
  private static String removeSpaceEscape(String s)
  {
    return s.replace("\\ ", " ");
  }
  
  class LaunchActionListener implements ActionListener
  {
    private HashMap<Integer, Vector<JTextField>> settings;
    
    private void setSettings(HashMap<Integer, Vector<JTextField>> settings)
    {
      this.settings = settings;
    }
    
    private String[] getArgs()
    {
      try
      {
        System.getProperties().remove("bam");
        System.getProperties().remove("chado");
        System.getProperties().remove("userplot");
        System.getProperties().remove("loguserplot");
      }
      catch(Exception e){ e.printStackTrace(); }
      
      boolean seenSequence  = false;
      final Set<Integer> keys = settings.keySet();
      final Vector<String> vargs = new Vector<String>();
      final Vector<String> vann = new Vector<String>();
      
      for(Integer key: keys)
      {
        final Vector<JTextField> vText = settings.get(key);
        switch(key)
        {
          case ProjectProperty.REFERENCE:
            String ref = vText.get(0).getText().trim();
            if(!ref.equals(""))
              seenSequence = true;
            vargs.add( ref );
            break;
          case ProjectProperty.ANNOTATION:
            for(JTextField ann: vText)
              if(ann.isEnabled())
                vann.add( ann.getText().trim() );
            break;
          case ProjectProperty.NEXT_GEN_DATA :
            setBam(vText);
            break;
          case ProjectProperty.VCF :
            setBam(vText);
            break;
          case ProjectProperty.USERPLOT:
            String userplot = "";
            for(JTextField ann: vText)
            {
              if(ann.isEnabled())
                userplot += ","+ann.getText().trim();
            }
            if(!userplot.equals(""))
              System.setProperty("userplot", userplot.replaceFirst(",", "")); 
            break;
          case ProjectProperty.LOGUSERPLOT:
            String loguserplot = "";
            for(JTextField ann: vText)
            {
              if(ann.isEnabled())
                loguserplot += ","+ann.getText().trim();
            }
            if(!loguserplot.equals(""))
              System.setProperty("loguserplot", loguserplot.replaceFirst(",", "")); 
            break;
          case ProjectProperty.CHADO:
            seenSequence = true;
            System.setProperty("chado", vText.get(0).getText().trim());
            break;
        }
      }
      
      if(!seenSequence)
        JOptionPane.showMessageDialog(ProjectProperty.this, 
            "No sequence file entered for this project.", 
            "Sequence Entry Missing", JOptionPane.WARNING_MESSAGE);
      
      String[] args = new String[vargs.size()+(vann.size()*2)];
      for(int i=0; i<vargs.size(); i++)
        args[i] = vargs.get(i);
      for(int i=0; i<vann.size(); i++)
      {
        args[vargs.size()+(i*2)] = "+";
        args[vargs.size()+(i*2)+1] = vann.get(i);
      }
      return args;
    }
    
    private void setBam(final Vector<JTextField> vText)
    {
      String bam = "";
      for(JTextField ann: vText)
        if(ann.isEnabled())
          bam += ","+ann.getText().trim();
      if(!bam.equals(""))
      {
        if(System.getProperty("bam") != null)
          bam += ","+System.getProperty("bam");
        System.setProperty("bam", bam.replaceFirst(",", ""));
      }
    }
    
    public void actionPerformed(ActionEvent arg0)
    {
      if(settings == null)
      {
        JOptionPane.showMessageDialog(ProjectProperty.this, 
            "Select a project.", "No Project", JOptionPane.INFORMATION_MESSAGE);
        return;
      }
      SwingUtilities.invokeLater(new Runnable() 
      {
        public void run() 
        {
          final String[] args = getArgs();
          if(System.getProperty("chado") != null &&
             args != null &&
             args.length == 1 &&
             args[0].indexOf(":") == -1)
          {
            openDatabase(splash, args[0]);
            return;
          }
          
          setCursor(new Cursor(Cursor.WAIT_CURSOR));
          try
          {
            final ArtemisMain main_window;
            if (splash == null)
            {
              main_window = new ArtemisMain(args);
              main_window.setVisible(true);
            }
            else
              main_window = (ArtemisMain) splash;
            main_window.readArgsAndOptions(args, ProjectProperty.this);
          } 
          finally
          {
            setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
          }
        }
      });
    }
  }
  
  class TypeComparator implements Comparator<Object> 
  {
    public int compare(Object o1, Object o2)
    {
      String s1 = (String)o1;
      String s2 = (String)o2;
      if( s1.equals("title") )
        return -1;
      else if(s2.equals("title"))
        return 1;
      
      return s1.compareTo(s2);
    }
  }
  
  public static void main(String args[])
  {
    new ProjectProperty();
  }
}