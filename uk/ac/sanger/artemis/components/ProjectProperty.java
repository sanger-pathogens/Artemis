/* ProjectProperty.java
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
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.io.InputStream;

import java.util.HashMap;
import java.util.Properties;
import java.util.Vector;
import java.util.Map.Entry;
import java.util.Set;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.border.Border;
import javax.swing.border.TitledBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import uk.ac.sanger.artemis.components.genebuilder.TextAreaDocumentListener;

class ProjectProperty extends JFrame
{
  private static final long serialVersionUID = 1L;
  private HashMap<String, HashMap<String, String>> projects;
  
  private final static int REFERENCE = 1;
  private final static int ANNOTATION = 2;
  private final static int NEXT_GEN_DATA = 3;
  private final static int CHADO = 4;
  
  ProjectProperty()
  {
    final InputStream ins = 
        this.getClass().getClassLoader().getResourceAsStream("etc/project.properties");
    final Properties projectProps = new Properties();
    try
    {
      projectProps.load(ins);
      ins.close();
      projects = getProjectMap(projectProps);
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
    createProjectViewer((JPanel) getContentPane());
    pack();
    setVisible(true);
  }
  
  private void createProjectViewer(JPanel panel)
  {
    panel.setPreferredSize(new Dimension(650,400));
    
    final JButton openArt = new JButton("OPEN");
    final LaunchActionListener listener = new LaunchActionListener();
    openArt.addActionListener(listener);
    panel.add(openArt, BorderLayout.SOUTH);
    
    final JList projectList = new JList(projects.keySet().toArray());
    projectList.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
    panel.add(projectList, BorderLayout.WEST);
    
    final GridBagConstraints c = new GridBagConstraints();
    c.anchor = GridBagConstraints.NORTHWEST;
    c.fill = GridBagConstraints.VERTICAL;
    c.ipadx = 10;
    c.ipady = 10;
    
    final Box props = Box.createVerticalBox();
    projectList.addListSelectionListener(new ListSelectionListener()
    {
      public void valueChanged(ListSelectionEvent e)
      {
        if(e.getValueIsAdjusting() == false) 
        {
          if(projectList.getSelectedIndex() > -1) 
          {
            props.removeAll();
            
            final HashMap<String, String> projProps = projects.get(projectList.getSelectedValue());
            final HashMap<Integer, QualifierTextArea> settings = new HashMap<Integer, QualifierTextArea>();
            for(String key: projProps.keySet())
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
                  new TextAreaDocumentListener(qta));

              xBox.add(qta);
              props.add(xBox);
              
              if(key.equals("ref"))
                settings.put(ProjectProperty.REFERENCE, qta);
              else if(key.equals("annotation"))
                settings.put(ProjectProperty.ANNOTATION, qta);
              else if(key.equals("bam") || key.equals("vcf") || key.equals("bcf"))
                settings.put(ProjectProperty.NEXT_GEN_DATA, qta);
              else if(key.equals("chado"))
                settings.put(ProjectProperty.CHADO, qta);
             
            }
            props.add(Box.createVerticalGlue());
            props.revalidate();
            props.repaint();
            
            listener.setSettings(settings);
          }
        }
      }
    });
    
    panel.add(props, BorderLayout.CENTER);
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
  
  class LaunchActionListener implements ActionListener
  {
    private String args[] = new String[]{};
    private void setSettings(HashMap<Integer, QualifierTextArea> settings)
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
        switch(key)
        {
          case ProjectProperty.REFERENCE:
            vargs.add( settings.get(key).getText() );
            break;
          case ProjectProperty.ANNOTATION:
            vann.add( settings.get(key).getText() );
            break;
          case ProjectProperty.NEXT_GEN_DATA:
            System.setProperty("bam", settings.get(key).getText().replaceAll("\\s+", ","));
            break;
          case ProjectProperty.CHADO:
            System.setProperty("chado", settings.get(key).getText().trim());
            break;
        }
      }
      
      args = new String[vargs.size()+(vann.size()*2)];
      for(int i=0; i<vargs.size(); i++)
        args[i] = vargs.get(i);
      for(int i=0; i<vann.size(); i++)
      {
        args[vargs.size()+(i*2)] = "+";
        args[vargs.size()+(i*2)+1] = vann.get(i);
      }
    }
    
    public void actionPerformed(ActionEvent arg0)
    {
      SwingUtilities.invokeLater(new Runnable() 
      {
        public void run() 
        {
          final ArtemisMain main_window = new ArtemisMain(args);
          main_window.setVisible(true);
          main_window.readArgsAndOptions(args);
        }
      });
    }
  }
  
  public static void main(String args[])
  {
    new ProjectProperty();
  }
}