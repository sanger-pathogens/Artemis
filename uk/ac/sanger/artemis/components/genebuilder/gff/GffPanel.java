/* GffPanel.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/gff/GffPanel.java,v 1.8 2007-10-17 14:18:32 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder.gff;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;

import javax.swing.Box;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.StringVector;

public class GffPanel extends JPanel
                      implements FeatureChangeListener
{
  private static final long serialVersionUID = 1L;
  private QualifierVector gffQualifiers;
  private JTextField uniquenameTextField;
  private JTextField timeTextField;
  
  public GffPanel(final Feature feature)
  {
    super(new FlowLayout(FlowLayout.LEFT));
    updateFromFeature(feature);
  }
  
  /**
   * Return true if this is a CV qualifier
   * @param qualifier
   * @return
   */
  public boolean isGffTag(final Qualifier qualifier)
  {
    if(qualifier.getName().equals("ID") ||
       qualifier.getName().equals("feature_id") ||
       qualifier.getName().equals("Parent") ||
       qualifier.getName().equals("Derives_from") ||
       qualifier.getName().equals("feature_relationship_rank") ||
       qualifier.getName().equals("timelastmodified"))
      return true;
    return false;
  }
  
  private Component createGffQualifiersComponent(final Feature feature)
  {
    int nrows = 0;
    
    Qualifier idQualifier = gffQualifiers.getQualifierByName("ID");
    if(idQualifier != null)
      nrows++;
    
    Qualifier parentQualifier = gffQualifiers.getQualifierByName("Parent");
    if(parentQualifier != null)
      nrows += parentQualifier.getValues().size();
     
    Qualifier derivesFromQualifier = gffQualifiers.getQualifierByName("Derives_from");
    if(derivesFromQualifier != null)
      nrows += derivesFromQualifier.getValues().size();
    
    Qualifier timeQualifier = gffQualifiers.getQualifierByName("timelastmodified");
    if(timeQualifier != null)
      nrows += timeQualifier.getValues().size();
    
    Box gffBox = Box.createVerticalBox();
    gffBox.add(Box.createVerticalStrut(10));
    
    GridBagLayout grid = new GridBagLayout();
    GridBagConstraints c = new GridBagConstraints();
    c.ipady = 5;
    JPanel gridPanel = new JPanel(grid);
    gridPanel.setBackground(Color.WHITE);
    
    Dimension cellDimension = null;
    
    int maxLabelWidth = new JLabel("timelastmodified").getPreferredSize().width;
    int featIdWidth   = 0;
    
    nrows = 0;
    if(idQualifier != null)
    {
      final String uniquename = (String)idQualifier.getValues().get(0);
      JLabel idField = new JLabel("ID");
      
      uniquenameTextField = new JTextField(uniquename);
      cellDimension = new Dimension(uniquenameTextField.getPreferredSize().width+10,
                                    idField.getPreferredSize().height+10);
      
      if(feature.getKey().getKeyString().indexOf("exon") > -1)
        uniquenameTextField.setEditable(false);
      uniquenameTextField.setMaximumSize(cellDimension);
      
      c.gridx = 0;
      c.gridy = 0;
      c.ipadx = 5;
      c.anchor = GridBagConstraints.EAST;
      gridPanel.add(idField, c);
      c.gridx = 1;
      c.gridy = 0;
      c.ipadx = 0;
      c.anchor = GridBagConstraints.WEST;
      gridPanel.add(uniquenameTextField, c);
      
      Qualifier featIdQualifier = gffQualifiers.getQualifierByName("feature_id");
      if(featIdQualifier != null)
      {
        JLabel featId = new JLabel("(feature_id="+(String)featIdQualifier.getValues().get(0)+")");
        featIdWidth = featId.getPreferredSize().width;
        c.gridx = 2;
        c.gridy = 0;
        c.ipadx = 0;
        c.anchor = GridBagConstraints.WEST;
        gridPanel.add(featId, c);
      }
      nrows++;
    }
    
    
    if(parentQualifier != null)
    {
      StringVector parents = parentQualifier.getValues();
      JLabel parentField = new JLabel("Parent");
      for(int i=0; i<parents.size(); i++)
      {
        String parent = (String)parents.get(i);
        JTextField parentTextField = new JTextField(parent);
        
        if(cellDimension == null ||
           cellDimension.width < parentTextField.getPreferredSize().width+10)
          cellDimension = new Dimension(parentTextField.getPreferredSize().width+10,
                                        parentField.getPreferredSize().height+10);
        parentTextField.setMaximumSize(cellDimension);
        parentTextField.setEditable(false);
        
        c.gridx = 0;
        c.gridy = nrows;
        c.ipadx = 5;
        c.anchor = GridBagConstraints.EAST;
        gridPanel.add(parentField, c); 
        c.gridx = 1;
        c.gridy = nrows;
        c.ipadx = 0;
        c.anchor = GridBagConstraints.WEST;
        gridPanel.add(parentTextField, c); 
        nrows++;
      }
    }
      
    
    if(derivesFromQualifier != null)
    {
      StringVector derivesFroms = derivesFromQualifier.getValues();
      JLabel derivesFromsField = new JLabel("Derives_from");
      for(int i=0; i<derivesFroms.size(); i++)
      {
        String derivesFrom = (String)derivesFroms.get(i);
        JTextField derivesFromTextField = new JTextField(derivesFrom);
        
        if(cellDimension == null ||
           cellDimension.width < derivesFromTextField.getPreferredSize().width+10)
          cellDimension = new Dimension(derivesFromTextField.getPreferredSize().width+10,
                                        derivesFromsField.getPreferredSize().height+10);
        derivesFromTextField.setMaximumSize(cellDimension);
        derivesFromTextField.setEditable(false);
        
        c.gridx = 0;
        c.gridy = nrows;
        c.ipadx = 5;
        c.anchor = GridBagConstraints.EAST;
        gridPanel.add(derivesFromsField, c); 
        c.gridx = 1;
        c.gridy = nrows;
        c.ipadx = 0;
        c.anchor = GridBagConstraints.WEST;
        gridPanel.add(derivesFromTextField, c); 
        nrows++;
      }
    }
    
    
    if(timeQualifier != null)
    {
      String time = (String)timeQualifier.getValues().get(0);
      
      JLabel timeField = new JLabel("timelastmodified");
      timeField.setPreferredSize(new Dimension(maxLabelWidth,
                       timeField.getPreferredSize().height));
      
      timeTextField = new JTextField(time);
      if(cellDimension == null ||
         cellDimension.width < timeTextField.getPreferredSize().width+10)
         cellDimension = new Dimension(timeTextField.getPreferredSize().width+10,
                                       timeField.getPreferredSize().height+10);
      timeTextField.setMaximumSize(cellDimension);
      timeTextField.setEditable(false);
      
      c.gridx = 0;
      c.gridy = nrows;
      c.ipadx = 5;
      c.anchor = GridBagConstraints.EAST;
      gridPanel.add(timeField, c);
      c.gridx = 1;
      c.gridy = nrows;
      c.ipadx = 0;
      c.anchor = GridBagConstraints.WEST;
      gridPanel.add(timeTextField, c);
      nrows++;
    }  
    

    if(cellDimension != null)
      gridPanel.setPreferredSize(
          new Dimension(maxLabelWidth+cellDimension.width+featIdWidth,
                        (cellDimension.height+5)*nrows));
    gffBox.add(gridPanel);
    return gffBox;
  }
  
  public void updateFromFeature(final Feature feature)
  {
    removeAll();
    if(gffQualifiers != null)
      feature.removeFeatureChangeListener(this);
    gffQualifiers = feature.getQualifiers().copy();
    
    gffQualifiers = new QualifierVector();
    final QualifierVector qualifiers = feature.getQualifiers();  
    for(int i = 0 ; i < qualifiers.size(); ++i) 
    {
      Qualifier qualifier = (Qualifier)qualifiers.elementAt(i);
      if(isGffTag(qualifier))
        gffQualifiers.addElement(qualifier.copy());
    }
   
    feature.addFeatureChangeListener(this);  
    add(createGffQualifiersComponent(feature));
    repaint();
    revalidate();
  }

  /**
   * Get the latest (edited) controlled vocab qualifiers
   * @return
   */
  public QualifierVector getGffQualifiers()
  {
    // check editable components for changes
    
    Qualifier idQualifier = gffQualifiers.getQualifierByName("ID");
    if(!((String)(idQualifier.getValues().get(0))).equals(uniquenameTextField.getText()))
    {
      if(!uniquenameTextField.getText().equals(""))
      {
        gffQualifiers.remove(idQualifier);
        idQualifier = new Qualifier("ID", uniquenameTextField.getText());
        gffQualifiers.addElement(idQualifier);
      }
    }
    
    return gffQualifiers;
  }
  
  public void featureChanged(FeatureChangeEvent event)
  {
    updateFromFeature(event.getFeature());
  }
}