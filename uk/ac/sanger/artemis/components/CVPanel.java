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

import java.awt.Color;
import java.awt.Component;
import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.JPanel;
import javax.swing.JFrame;
import javax.swing.JButton;
import javax.swing.ListCellRenderer;
import javax.swing.ScrollPaneConstants;
import javax.swing.plaf.basic.BasicComboBoxUI;
import javax.swing.plaf.basic.BasicComboPopup;
import javax.swing.plaf.basic.ComboPopup;

import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.chado.ChadoTransactionManager;

import org.gmod.schema.cv.CvTerm;

/**
 * Panel for display controlled vocabulary terms for Chado
 */
public class CVPanel extends JPanel
             implements FeatureChangeListener
{
  private Feature feature;
  private String[][] evidenceCodes = 
  { 
     {"IC",  "Inferred by Curator"},
     {"IDA", "Inferred from Direct Assay "},
     {"IEA", "Inferred from Electronic Annotation"},
     {"IEP", "Inferred from Expression Pattern"},
     {"IGI", "Inferred from Genetic Interaction"},
     {"IMP", "Inferred from Mutant Phenotype"},
     {"IPI", "Inferred from Physical Interaction"},
     {"ISS", "Inferred from Sequence or Structural Similarity"},
     {"NAS", "Non-traceable Author Statement"},
     {"ND",  "No biological Data available"},
     {"RCA", "inferred from Reviewed Computational Analysis"},
     {"TAS", "Traceable Author Statement"},
     {"NR",  "Not Recorded"}
  };
  
  public CVPanel(final Feature feature)
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
        addCvTerm();
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
            String goId = getField("GOid=", qualifierString);
            
            JTextField goTextField = new JTextField(goId);
            goTextField.setPreferredSize(dimension);
            goTextField.setMaximumSize(dimension);
            
            xBox.add(goTextField);
            
            String term = getField("term=", qualifierString);
            JTextField termTextField = new JTextField(term);
            termTextField.setPreferredSize(dimension2);
            termTextField.setMaximumSize(dimension2);
            
            xBox.add(termTextField);
            
            String with = getField("with=", qualifierString);
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
  
  private void addCvTerm()
  {
    Box xBox = Box.createHorizontalBox();
    final JComboBox comboCV = 
      new ComboBoxPopup(ChadoTransactionManager.cv_tags);
    Dimension d = comboCV.getPreferredSize();
    d = new Dimension(80,(int)d.getHeight());
    comboCV.setPreferredSize(d);
    xBox.add(comboCV);
    
    int step = 0;

    String cv_name = null;
    Vector terms = null;
    while(step < 4)
    {  
      if(step == 0)
      {
        cv_name = prompCv(xBox, comboCV);
        if(cv_name == null)
          return;
        step = 1;
      }
      
      if(step == 1)
      {
        cv_name = promptCvName(xBox, comboCV);
        if(cv_name == null)
          return;
        else if(cv_name.equals(""))
        {
          step = 0;
          continue;
        }
        step = 2;
      }
          
      if(step == 2)
      {
        terms = promptKeyWord(xBox, cv_name);
        if(terms == null)
          return;
        else if(terms.size() < 1)
        {
          if(cv_name.equals("molecular_function") ||
             cv_name.equals("biological_process") ||
             cv_name.equals("cellular_component"))
            step = 1;
          else
            step = 0;
          continue;
        }
        step = 3;
      }
      
      if(step == 3)
      {
        CvTerm cvterm = promptCvTerms(xBox, terms);
        if(cvterm == null)
          return;
        else if(cvterm.getName() == null)
        {
          step = 2;
          continue;
        }
        step = 4;
      }
    }
  }
  
  private String prompCv(final Box xBox, final JComboBox comboCV)
  {
    final Object[] options  = { "CANCEL", "NEXT>"};

    int select = JOptionPane.showOptionDialog(null, 
          xBox, "Add CV",
          JOptionPane.YES_NO_CANCEL_OPTION, 
          JOptionPane.QUESTION_MESSAGE, null,
          options, options[1]);
    if(select == 0)
      return null;
    return ChadoTransactionManager.cv_tags[comboCV.getSelectedIndex()];
  }
  
  private String promptCvName(final Box xBox, final JComboBox comboCV)
  {
    final String options[] = { "<PREV", "CANCEL", "NEXT>"};
     
    String cv_name = ChadoTransactionManager.cv_tags[comboCV
          .getSelectedIndex()];
    Splash.logger4j.debug("Selected CV is " + cv_name);

    if(cv_name.equals("GO"))
    {
      final String aspect[] = { "F", "P", "C" };
      final JComboBox combo = new JComboBox(aspect);
      
      xBox.removeAll();
      xBox.add(comboCV);
      xBox.add(combo);

      int select = JOptionPane.showOptionDialog(null, 
            xBox, "Aspect",
            JOptionPane.YES_NO_CANCEL_OPTION, 
            JOptionPane.QUESTION_MESSAGE,
            null, options, options[2]);
      if(select == 1)
        return null;
      else if(select == 0)
      {
        xBox.remove(combo);
        return "";
      }

      if(((String) combo.getSelectedItem()).equals("F"))
          cv_name = "molecular_function";
      else if(((String) combo.getSelectedItem()).equals("P"))
        cv_name = "biological_process";
      else if(((String) combo.getSelectedItem()).equals("C"))
        cv_name = "cellular_component";
    }
    else if(cv_name.equals("product"))
      cv_name = DatabaseDocument.PRODUCTS_TAG_CVNAME;
    else if(cv_name.equals("controlled_curation"))
      cv_name = DatabaseDocument.CONTROLLED_CURATION_TAG_CVNAME;
    
    return cv_name;
  }
  
  private Vector promptKeyWord(final Box xBox, final String cv_name)
  {
    final String options[] = { "<PREV", "CANCEL", "NEXT>"};
    final JTextField tfield = new JTextField(12);
    tfield.setSelectionStart(0);
    tfield.setSelectionEnd(tfield.getText().length());
    tfield.setSelectedTextColor(Color.blue);
    xBox.add(tfield);
    
    int select = JOptionPane.showOptionDialog(null, xBox,
        "keyword term selection",
         JOptionPane.YES_NO_CANCEL_OPTION,
         JOptionPane.QUESTION_MESSAGE,
         null,
         options,
         options[2]);
    
    if(select == 1)
      return null;
    else if(select == 0)
    {
      xBox.remove(tfield);
      return new Vector(0);
    }
    
    xBox.remove(tfield);
    
    return DatabaseDocument.getCvterms(tfield.getText().trim(), cv_name);
  }
  
  private CvTerm promptCvTerms(final Box xBox, final Vector terms)
  {
    final String options[] = { "<PREV", "CANCEL", "NEXT>"};
    JComboBox term_list = new ComboBoxPopup(terms);
    
    Dimension d = term_list.getPreferredSize();
    d = new Dimension(280,(int)d.getHeight());
    term_list.setPreferredSize(d);
    
    term_list.setRenderer(new CVCellRenderer());
    xBox.add(term_list);
    
    int select = JOptionPane.showOptionDialog(null, xBox,
        "CV term selection",
         JOptionPane.YES_NO_CANCEL_OPTION,
         JOptionPane.QUESTION_MESSAGE,
         null,
         options,
         options[2]);
    
    if(select == 1)
      return null;
    else if(select == 0)
    {
      xBox.remove(term_list);
      return new CvTerm();
    }
    
    return (CvTerm)term_list.getSelectedItem();
  }
  
  /**
   * Strip out the value of a field of interest from a qualifier string
   * 
   * @param fieldName
   * @param qualifierString
   * @return
   */
  private String getField(final String fieldName, final String qualifierString)
  {
    String field = null;
    
    int ind1 = qualifierString.toLowerCase().indexOf(fieldName.toLowerCase());
    int ind2 = qualifierString.indexOf(";", ind1);
    
    int len = fieldName.length();

    if(ind2 > ind1 && ind1 > -1)
      field = qualifierString.substring(ind1+len,ind2);
    else if(ind1 > -1)
      field = qualifierString.substring(ind1+len);
    
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
  
  
  class CVCellRenderer extends JLabel implements ListCellRenderer 
  {
    public CVCellRenderer() 
    {
      setOpaque(true);
    }

    public Component getListCellRendererComponent(JList list,
                                                  Object value,
                                                  int index,
                                                  boolean isSelected,
                                                  boolean cellHasFocus) 
    {
      setText(((CvTerm)value).getName());
      return this;
    }
  }
  
  
  class ComboBoxPopup extends JComboBox
  {
    
    /**
    *  Creates a JComboBox that contains the elements in 
    *  the specified array.
    *  @param items       array of objects to insert into the combo box
    */
    public ComboBoxPopup(Object[] items)
    {
      super(items);
      setUI(new myComboUI());
      setBorder(BorderFactory.createEtchedBorder());
    }

    /**
    *  Creates a JComboBox that contains the elements in
    *  the vextor
    *  @param items       vector of objects to insert into the combo box
    */
    public ComboBoxPopup(Vector items)
    {
      super(items);
      setUI(new myComboUI());
      setBorder(BorderFactory.createEtchedBorder());
    }
            
    /**
    * Create scroller
    */                             
    public class myComboUI extends BasicComboBoxUI
    {
      protected ComboPopup createPopup()
      {
        BasicComboPopup popup = new BasicComboPopup(comboBox)
        {
          protected JScrollPane createScroller() 
          {
            return new JScrollPane( list, ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED,
                                 ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED );
          }
        };
        return popup;
      }
    }
  }

  
}