/* CVPanel.java
 *
 * created: 2006
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2006  Genome Research Limited
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
 **/

package uk.ac.sanger.artemis.components;

import java.awt.Component;
import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.JComboBox;
import javax.swing.JOptionPane;
import javax.swing.JTextField;
import javax.swing.JPanel;
import javax.swing.JFrame;
import javax.swing.JButton;

import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.chado.ChadoTransactionManager;

/**
 * Panel for display controlled vocabulary terms for Chado
 */
public class CVPanel extends JPanel
             implements FeatureChangeListener
{
  private Feature feature;
  
  public CVPanel(final Feature feature, final JFrame frame)
  {
    super(new BorderLayout());
    this.feature = feature;
    
    add(createCVQualifiersComponent(),
        BorderLayout.CENTER);
    
    feature.addFeatureChangeListener(this);  
  }
  
  /**
   * Return the feature
   * @return
   */
  private Feature getFeature()
  {
    return feature; 
  }
  
  /**
   * Create CV qualifier components
   * @return
   */
  private Component createCVQualifiersComponent()
  {
    Vector cv_tags = new Vector();
    cv_tags.add("GO");
    cv_tags.add("controlled_curation");
    cv_tags.add("product");

    Dimension dimension  = new Dimension(100, 30);
    Dimension dimension2 = new Dimension(150, 30);

    Box cvBox = Box.createVerticalBox();
    
    Box xBox = Box.createHorizontalBox();
    JButton addRemove = new JButton("ADD");
    addRemove.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      { 
        Object[] options = { "CANCEL", "NEXT>"};

        final JComboBox comboCV = new JComboBox(ChadoTransactionManager.cv_tags);
        int select = JOptionPane.showOptionDialog(null, comboCV,
                                    "Add CV",
                                     JOptionPane.YES_NO_CANCEL_OPTION,
                                     JOptionPane.QUESTION_MESSAGE,
                                     null,
                                     options,
                                     options[1]);
        if(select == 0)
          return;
        
        final String cv = ChadoTransactionManager.cv_tags[comboCV.getSelectedIndex()];
        Splash.logger4j.debug("Selected CV is "+cv);
        
        final String options2[] = { "<PREV", "CANCEL", "NEXT>"};
        if(cv.equals("GO"))
        {
          final String aspect[] = { "F", "P", "C" };
          final JComboBox combo = new JComboBox(aspect);
          
          
          select = JOptionPane.showOptionDialog(null, combo,
                                      "Aspect",
                                       JOptionPane.YES_NO_CANCEL_OPTION,
                                       JOptionPane.QUESTION_MESSAGE,
                                       null,
                                       options2,
                                       options2[2]);
        }
        
        
        
      }
    });
    xBox.add(addRemove);
    xBox.add(Box.createHorizontalGlue());
    cvBox.add(xBox);
    
    final QualifierVector qualifiers = getFeature().getQualifiers();

    for(int qualifier_index = 0; qualifier_index < qualifiers.size();
        ++qualifier_index) 
    {
      final Qualifier this_qualifier = (Qualifier)qualifiers.elementAt(qualifier_index);

      if(cv_tags.contains(this_qualifier.getName()))
      {
        final StringVector qualifier_strings = this_qualifier.getValues();
        
        for(int value_index = 0; value_index < qualifier_strings.size();
            ++value_index)
        {
          xBox = Box.createHorizontalBox();
          final String qualifierString = 
            (String)qualifier_strings.elementAt(value_index);
               
          if(this_qualifier.getName().equals("GO"))
          {
            String goId = getField("GOid", qualifierString);
            
            JTextField goTextField = new JTextField(goId);
            goTextField.setPreferredSize(dimension);
            goTextField.setMaximumSize(dimension);
            
            xBox.add(goTextField);
            
            String term = getField("term", qualifierString);
            JTextField termTextField = new JTextField(term);
            termTextField.setPreferredSize(dimension2);
            termTextField.setMaximumSize(dimension2);
            
            xBox.add(termTextField);
            
            String with = getField("with", qualifierString);
            JTextField withTextField = new JTextField(with);
            withTextField.setPreferredSize(dimension2);
            withTextField.setMaximumSize(dimension2);
            
            xBox.add(withTextField);
            
            xBox.add(Box.createHorizontalGlue());
            
            JButton buttRemove = new JButton("REMOVE");
            buttRemove.addActionListener(new ActionListener()
            {
              public void actionPerformed(ActionEvent e)
              { 
              }
            });
            cvBox.add(xBox);
          }
        }
      }
    }
    
    cvBox.add(Box.createVerticalGlue());
    return cvBox;
  }
  
  
  private String getField(final String fieldName, final String qualifierString)
  {
    String field = null;
    
    int ind1 = qualifierString.toLowerCase().indexOf(fieldName.toLowerCase());
    int ind2 = qualifierString.indexOf(";", ind1);
    
    int len = fieldName.length();
    
    if(ind2 > ind1 && ind1 > -1)
      field = qualifierString.substring(ind1+len+1,ind2);
    else if(ind1 > -1)
      field = qualifierString.substring(ind1+len+1);
    
    return field;
  }
  
  /**
   *  Implementation of the FeatureChangeListener interface.  We need to
   *  listen to feature change events from the Features in this object so that
   *  we can update the display.
   *  @param event The change event.
   **/
  public void featureChanged(FeatureChangeEvent event) 
  {
    // re-read the information from the feature
    switch(event.getType()) 
    {
      case FeatureChangeEvent.QUALIFIER_CHANGED:
        /*if(qualifier_text_area.getText().equals(orig_qualifier_text)) 
          updateFromFeature();
        else
        {
          final String message =
            "warning: the qualifiers have changed outside the editor - " +
            "view now?";

          final YesNoDialog yes_no_dialog =
            new YesNoDialog(frame, message);

          if(yes_no_dialog.getResult()) 
            new FeatureViewer(getFeature());
        }*/
        break;
      default:
        break;
    }
  }
}