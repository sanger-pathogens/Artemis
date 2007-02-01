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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/gff/GffPanel.java,v 1.2 2007-02-01 11:44:54 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder.gff;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;

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
  private JTextField idTextField;
  
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
       qualifier.getName().equals("Parent"))
      return true;
    return false;
  }
  
  private Component createGffQualifiersComponent()
  {
    int nrows = 0;
    
    Qualifier idQualifier = gffQualifiers.getQualifierByName("ID");
    if(idQualifier != null)
      nrows++;
    
    Qualifier parentQualifier = gffQualifiers.getQualifierByName("Parent");
    if(parentQualifier != null)
      nrows += parentQualifier.getValues().size();
      
    Box gffBox = Box.createVerticalBox();
    gffBox.add(Box.createVerticalStrut(10));
    
    GridBagLayout grid = new GridBagLayout();
    GridBagConstraints c = new GridBagConstraints();
    JPanel gridPanel = new JPanel(grid);
    
    Dimension cellDimension = null;
    JLabel parentField = new JLabel("Parent");
    int maxLabelWidth = parentField.getPreferredSize().width;;
    if(idQualifier != null)
    {
      String featureId = (String)idQualifier.getValues().get(0);
      
      JLabel idField = new JLabel("ID");
      idField.setPreferredSize(new Dimension(maxLabelWidth,
                       idField.getPreferredSize().height));
      
      idTextField = new JTextField(featureId);
      cellDimension = new Dimension(idTextField.getPreferredSize().width+10,
                                    idField.getPreferredSize().height+10);
      idTextField.setMaximumSize(cellDimension);
      
      c.gridx = 0;
      c.gridy = 0;
      c.ipadx = 5;
      c.anchor = GridBagConstraints.EAST;
      gridPanel.add(idField, c);
      c.gridx = 1;
      c.gridy = 0;
      c.ipadx = 0;
      c.anchor = GridBagConstraints.WEST;
      gridPanel.add(idTextField, c);
    }
    
    
    if(parentQualifier != null)
    {
      StringVector parents = parentQualifier.getValues();
      for(int i=0; i<parents.size(); i++)
      {
        String parent = (String)parents.get(i);
        JTextField parentTextField = new JTextField(parent);
        
        if(cellDimension == null ||
           cellDimension.width < parentTextField.getPreferredSize().width+10)
          cellDimension = new Dimension(parentTextField.getPreferredSize().width+10,
                                        parentField.getPreferredSize().height+10);
        parentTextField.setMaximumSize(cellDimension);
        
        c.gridx = 0;
        c.gridy = i+1;
        c.ipadx = 5;
        c.anchor = GridBagConstraints.EAST;
        gridPanel.add(parentField, c); 
        c.gridx = 1;
        c.gridy = i+1;
        c.ipadx = 0;
        c.anchor = GridBagConstraints.WEST;
        gridPanel.add(parentTextField, c); 
      }
    }

    if(cellDimension != null)
      gridPanel.setPreferredSize(new Dimension(maxLabelWidth+cellDimension.width,
                                               cellDimension.height*nrows));
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
    add(createGffQualifiersComponent());
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
    if(!((String)(idQualifier.getValues().get(0))).equals(idTextField.getText()))
    {
      if(!idTextField.getText().equals(""))
      {
        gffQualifiers.remove(idQualifier);
        idQualifier = new Qualifier("ID", idTextField.getText());
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