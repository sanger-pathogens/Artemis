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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/ortholog/OrthologPanel.java,v 1.1 2007-03-15 11:39:07 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder.ortholog;

import java.awt.Component;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;

import javax.swing.Box;
import javax.swing.JLabel;
import javax.swing.JPanel;

import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;

public class OrthologPanel extends JPanel
                      implements FeatureChangeListener
{
  private static final long serialVersionUID = 1L;
  private QualifierVector orthologQualifiers;
  
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
  
  private Component createOrthologQualifiersComponent()
  {
    int nrows = 0;
    
    Qualifier idQualifier = orthologQualifiers.getQualifierByName("ID");
    
    Box gffBox = Box.createVerticalBox();
    gffBox.add(Box.createVerticalStrut(10));
    
    GridBagLayout grid = new GridBagLayout();
    GridBagConstraints c = new GridBagConstraints();
    c.ipady = 5;
    JPanel gridPanel = new JPanel(grid);
    

    if(idQualifier != null)
    {
      final String uniquename = (String)idQualifier.getValues().get(0);
      JLabel idField = new JLabel("ID");
      
      
      c.gridx = 0;
      c.gridy = 0;
      c.ipadx = 5;
      c.anchor = GridBagConstraints.EAST;
      gridPanel.add(idField, c);
      c.gridx = 1;
      c.gridy = 0;
      c.ipadx = 0;
      c.anchor = GridBagConstraints.WEST;
      //gridPanel.add(uniquenameTextField, c);
      

    }
    
    

    gffBox.add(gridPanel);
    return gffBox;
  }
  
  public void updateFromFeature(final Feature feature)
  {
    removeAll();
    if(orthologQualifiers != null)
      feature.removeFeatureChangeListener(this);
    orthologQualifiers = feature.getQualifiers().copy();
    
    orthologQualifiers = new QualifierVector();
    final QualifierVector qualifiers = feature.getQualifiers();  
    for(int i = 0 ; i < qualifiers.size(); ++i) 
    {
      Qualifier qualifier = (Qualifier)qualifiers.elementAt(i);
      if(isOrthologTag(qualifier))
        orthologQualifiers.addElement(qualifier.copy());
    }
   
    feature.addFeatureChangeListener(this);  
    add(createOrthologQualifiersComponent());
    repaint();
    revalidate();
  }

  /**
   * Get the latest (edited) controlled vocab qualifiers
   * @return
   */
  public QualifierVector getOrthologQualifiers()
  {
    // check editable components for changes
    
    // 
    
    return orthologQualifiers;
  }
  
  public void featureChanged(FeatureChangeEvent event)
  {
    updateFromFeature(event.getFeature());
  }
}