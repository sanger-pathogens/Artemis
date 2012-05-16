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
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.net.URISyntaxException;

import java.util.Arrays;
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
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.border.Border;
import javax.swing.border.TitledBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import uk.ac.sanger.artemis.components.genebuilder.TextAreaDocumentListener;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.FileDocument;

/**
 * Project file management system using a properties file.
 * 
 * Example of the syntax for defining a project in the property file:
 * project.Pknowlsei.ref = Pknowlsei:Pk_strainH_chr01
 * project.Pknowlsei.chado =  genedb-db.sanger.ac.uk:5432/snapshot?genedb_ro
 * project.Pknowlsei.title = Pknowlsei
 */
public class ProjectProperty extends JFrame
{
  private static final long serialVersionUID = 1L;
  private static HashMap<String, HashMap<String, String>> centralProjects;
  private static HashMap<String, HashMap<String, String>> userProjects;
  private Splash splash;
  
  private final static int REFERENCE = 1;
  private final static int ANNOTATION = 2;
  private final static int NEXT_GEN_DATA = 3;
  private final static int CHADO = 4;
  private static org.apache.log4j.Logger logger4j = 
      org.apache.log4j.Logger.getLogger(ProjectProperty.class);
  
  private final static String[] TYPES = 
    { "title", "ref", "annotation", "bam", "vcf", "chado" };
  
  public ProjectProperty()
  {
    this(null);
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
    centralProjects = getProjectMap(projectProps);
    projectProps.clear();
    
    final String [] propFiles = 
    {
      "project.properties",
      System.getProperty("user.home") + File.separator + ".project.properties"
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
    final LaunchActionListener listener = new LaunchActionListener();
    projectList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    projectList.setVisibleRowCount(-1);
    panel.add(jspList, BorderLayout.WEST);
    panel.add(yBox, BorderLayout.CENTER);
    
    // menus
    final JMenuBar menuBar = new JMenuBar();
    setJMenuBar(menuBar);
    final JMenu fileMenu = new JMenu("File");
    menuBar.add(fileMenu);
    
    final JMenuItem exitMenu = new JMenuItem(splash == null ? "Exit" : "Close");
    fileMenu.add(exitMenu);
    exitMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent arg0)
      {
        if (splash == null)
        {
          int status = JOptionPane.showConfirmDialog(ProjectProperty.this,
              "Exit?", "Exit", JOptionPane.YES_NO_OPTION);
          if (status == JOptionPane.YES_OPTION)
            Splash.exitApp();
        }
        else
          dispose();
      }
    });


    final JMenu editMenu = new JMenu("Edit");
    menuBar.add(editMenu);
    final JMenuItem addProject = new JMenuItem("Add Project");
    addProject.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent arg0)
      {
        String projName = 
            JOptionPane.showInputDialog(ProjectProperty.this,
                "Project Name", "New Project", JOptionPane.QUESTION_MESSAGE);
        if(projName == null)
          return;
        userProjects.put(projName, new HashMap<String, String>());
        model.add(model.getSize(), projName);
        projectList.repaint();
        projectList.setSelectedIndex(model.getSize()-1);
      }
    });
    editMenu.add(addProject);
    final JMenuItem removeProject = new JMenuItem("Remove Project");
    removeProject.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent arg0)
      {
        if(projectList.getSelectedValue() == null)
          return;
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
      }
    });
    editMenu.add(removeProject);
    
    
    //
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    panel.setPreferredSize(new Dimension(screen.width/2,400));
    panel.setBackground(Color.WHITE);
    
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

    final HashMap<Integer, QualifierTextArea> settings = new HashMap<Integer, QualifierTextArea>();
    for(final String key: projProps.keySet())
    {
      Box xBox = Box.createHorizontalBox();
      final QualifierTextArea qta = new QualifierTextArea();
      Border loweredbevel = BorderFactory.createLoweredBevelBorder();
      TitledBorder title = BorderFactory.createTitledBorder(
          loweredbevel, key);
      title.setTitlePosition(TitledBorder.ABOVE_TOP);

      qta.setBorder(title);

      qta.setText(projProps.get(key));
      qta.getDocument().addDocumentListener(
          new TextAreaDocumentListener(qta) {
            public void updateSize(DocumentEvent e)
            {
              projProps.put(key, qta.getText());
              setQualifierTextAreaSize();
            }
          });

      xBox.add(qta);
      
      // REMOVE PROPERTY
      Box yButtons = Box.createVerticalBox();
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
          projProps.remove(key);
          refreshProperties(projectList, yBox, listener);
        }
      });
      yButtons.add(removeProperty);
      //
      if(!key.equals("title") && !key.equals("ref"))
      {
        final JCheckBox useProperty = new JCheckBox("",true);
        useProperty.addActionListener(new ActionListener(){
          public void actionPerformed(ActionEvent arg0)
          {
            qta.setEnabled(useProperty.isSelected());
          }
        });
        yButtons.add(useProperty);
      }
      
      if (!key.equals("title"))
      {
        final JButton selectButton = new JButton("Add file");
        selectButton.addActionListener(new ActionListener()
        {
          public void actionPerformed(ActionEvent e)
          {
            StickyFileChooser fileChooser = new StickyFileChooser();
            int status = fileChooser.showOpenDialog(ProjectProperty.this);

            if(status == StickyFileChooser.APPROVE_OPTION)
            {
              qta.append(System.getProperty("line.separator")+
                  fileChooser.getSelectedFile().getAbsolutePath());
            }
          }
        });
        yButtons.add(selectButton);
      }
      
      xBox.add(yButtons);
      yBox.add(xBox);
      
      if(key.equals("ref"))
        settings.put(ProjectProperty.REFERENCE, qta);
      else if(key.equals("annotation"))
        settings.put(ProjectProperty.ANNOTATION, qta);
      else if(key.equals("bam") || key.equals("vcf") || key.equals("bcf"))
        settings.put(ProjectProperty.NEXT_GEN_DATA, qta);
      else if(key.equals("chado"))
        settings.put(ProjectProperty.CHADO, qta);
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
    final String prop = System.getProperty("user.home") + File.separator + ".project.properties";
    File propFile = new File(prop);
    try
    {
      if(userProjects.size() > 0)
      {
        propFile.delete();
        final BufferedWriter bufferedwriter = new BufferedWriter(new FileWriter(propFile));
        for (String project: userProjects.keySet())
        {
          bufferedwriter.write("#");
          bufferedwriter.newLine();
          
          HashMap<String, String> projProps = userProjects.get(project);
          for(final String key: projProps.keySet())
          {
            bufferedwriter.write("project."+project+"."+key+"="+projProps.get(key));
            bufferedwriter.newLine();
          }
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
  
  class LaunchActionListener implements ActionListener
  {
    private HashMap<Integer, QualifierTextArea> settings;
    
    private void setSettings(HashMap<Integer, QualifierTextArea> settings)
    {
      this.settings = settings;
    }
    
    private String[] getArgs()
    {
      try
      {
        System.getProperties().remove("bam");
        System.getProperties().remove("chado");
      }
      catch(Exception e){ e.printStackTrace(); }
      
      final Set<Integer> keys = settings.keySet();
      final Vector<String> vargs = new Vector<String>();
      final Vector<String> vann = new Vector<String>();
      for(Integer key: keys)
      {
        final QualifierTextArea qta = settings.get(key);
        if(!qta.isEnabled())
          continue;
        switch(key)
        {
          case ProjectProperty.REFERENCE:
            vargs.add( qta.getText().trim() );
            break;
          case ProjectProperty.ANNOTATION:
            vann.add( qta.getText().trim() );
            break;
          case ProjectProperty.NEXT_GEN_DATA:
            System.setProperty("bam", qta.getText().trim().replaceAll("\\s+", ","));
            break;
          case ProjectProperty.CHADO:
            System.setProperty("chado", qta.getText().trim());
            break;
        }
      }
      
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
    
    public void actionPerformed(ActionEvent arg0)
    {
      SwingUtilities.invokeLater(new Runnable() 
      {
        public void run() 
        {
          setCursor(new Cursor(Cursor.WAIT_CURSOR));
          try
          {
            String[] args = getArgs();
            final ArtemisMain main_window;
            if (splash == null)
            {
              main_window = new ArtemisMain(args);
              main_window.setVisible(true);
            }
            else
              main_window = (ArtemisMain) splash;
            main_window.readArgsAndOptions(args);
          } 
          finally
          {
            setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
          }
        }
      });
    }
  }
  
  public static void main(String args[])
  {
    new ProjectProperty();
  }
}