/* OrthologPanel.java
 * This file is part of Artemis
 *
 * Copyright (C) 2007  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/ortholog/OrthologPanel.java,v 1.3 2007-04-17 15:19:41 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder.ortholog;

import java.awt.Component;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.io.DocumentEntry;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.components.genebuilder.JExtendedComboBox;

public class OrthologPanel extends JPanel
                      implements FeatureChangeListener
{
  private static final long serialVersionUID = 1L;
  private QualifierVector ortho_para_logQualifiers;
  private static Vector databases;
  
  public OrthologPanel(final Feature feature)
  {
    super(new FlowLayout(FlowLayout.LEFT));
    updateFromFeature(feature);
  }
  
  /**
   * Return true if this is a Ortholog qualifier
   * @param qualifier
   * @return
   */
  public boolean isOrthologTag(final Qualifier qualifier)
  {
    if(qualifier.getName().equals("ortholog"))
      return true;
    return false;
  }
  
  private Component createOrthoParaLogQualifiersComponent(final Feature feature)
  {
    Qualifier orthoQualifier = ortho_para_logQualifiers.getQualifierByName("ortholog");
    Qualifier paraQualifier  = ortho_para_logQualifiers.getQualifierByName("paralog");
    
    if(databases == null)
    {
      DocumentEntry entry = (DocumentEntry)feature.getEmblFeature().getEntry();
      DatabaseDocument doc = (DatabaseDocument)entry.getDocument();
      databases = (Vector)doc.getDatabaseNames();
    }
    
    //
    // ortholog
    Box xBox = Box.createVerticalBox();
    JButton addOrthoButton = new JButton("ADD ORTHOLOG");
    addOrthoButton.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      { 
        JExtendedComboBox dbs = new JExtendedComboBox(databases);
        
        JTextField accession = new JTextField(15);
        
        Box yBox = Box.createHorizontalBox();
        yBox.add(dbs);
        yBox.add(accession);

        int select = JOptionPane.showConfirmDialog(null, 
              yBox, "Add Ortholog",
              JOptionPane.OK_CANCEL_OPTION);
        if(select == JOptionPane.CANCEL_OPTION)
          return;
        
        add("ortholog", ((String)dbs.getSelectedItem())+":"+
                        accession.getText().trim(), feature);
      }
    });
    xBox.add(addOrthoButton);
    
    if(orthoQualifier != null)
    {
      StringVector orthologs = orthoQualifier.getValues();
      for(int i=0; i<orthologs.size(); i++)
      {
        JTextField ortholog = new JTextField( (String)orthologs.get(i) );
        xBox.add(ortholog);
      }
    }
    
    //
    // paralog
    JButton addParaButton = new JButton("ADD PARALOG");
    addParaButton.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      { 
        JExtendedComboBox dbs = new JExtendedComboBox(databases);
        JTextField accession = new JTextField(15);
        
        Box yBox = Box.createHorizontalBox();
        yBox.add(dbs);
        yBox.add(accession);

        int select = JOptionPane.showConfirmDialog(null, 
              yBox, "Add Paralog",
              JOptionPane.OK_CANCEL_OPTION);
        if(select == JOptionPane.CANCEL_OPTION)
          return;
        
        add("paralog", ((String)dbs.getSelectedItem())+":"+
                        accession.getText().trim(), feature);
      }
    });
    xBox.add(addParaButton);
    
    if(paraQualifier != null)
    {
      StringVector paralogs = paraQualifier.getValues();
      for(int i=0; i<paralogs.size(); i++)
      {
        JTextField paralog = new JTextField( (String)paralogs.get(i) );
        xBox.add(paralog);
      }
    }

    return xBox;
  }
  
  /**
   * Update ortho/paralogs for a feature
   * @param feature
   */
  public void updateFromFeature(final Feature feature)
  {
    removeAll();
    if(ortho_para_logQualifiers != null)
      feature.removeFeatureChangeListener(this);
    ortho_para_logQualifiers = feature.getQualifiers().copy();
    
    ortho_para_logQualifiers = new QualifierVector();
    final QualifierVector qualifiers = feature.getQualifiers();  
    for(int i = 0 ; i < qualifiers.size(); ++i) 
    {
      Qualifier qualifier = (Qualifier)qualifiers.elementAt(i);
      if(isOrthologTag(qualifier))
        ortho_para_logQualifiers.addElement(qualifier.copy());
    }
   
    feature.addFeatureChangeListener(this);  
    add(createOrthoParaLogQualifiersComponent(feature));
    repaint();
    revalidate();
  }

  /**
   * Add ortholog/paralog
   * @param name  ortholog or paralog
   * @param value
   */
  public void add(final String name, final String value, final Feature feature)
  {
    final int index;
    
    Qualifier qualifier = ortho_para_logQualifiers.getQualifierByName(name);
    
    if(qualifier == null)
    {
      qualifier = new Qualifier(name);
      index = -1;
    }
    else
     index = ortho_para_logQualifiers.indexOf(qualifier);
       
    qualifier.addValue(value);

    if(index > -1)
    {
      ortho_para_logQualifiers.remove(index);
      ortho_para_logQualifiers.add(index, qualifier);
    }
    else
      ortho_para_logQualifiers.add(qualifier);
    
    removeAll();
    add(createOrthoParaLogQualifiersComponent(feature));
    repaint();
    revalidate();
  }
  
  /**
   * Get the latest (edited) controlled vocab qualifiers
   * @return
   */
  public QualifierVector getOrthologQualifiers()
  {
    return ortho_para_logQualifiers;
  }
  
  public void featureChanged(FeatureChangeEvent event)
  {
    updateFromFeature(event.getFeature());
  }
}