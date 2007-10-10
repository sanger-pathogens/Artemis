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

package uk.ac.sanger.artemis.components.genebuilder.cv;

import java.awt.Color;
import java.awt.Component;
import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Collections;
import java.util.Comparator;
import java.util.Vector;


import javax.swing.Box;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JTextField;
import javax.swing.JPanel;
import javax.swing.JButton;
import javax.swing.ListCellRenderer;

import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.chado.ChadoTransactionManager;
//import uk.ac.sanger.artemis.components.Splash;
import uk.ac.sanger.artemis.components.genebuilder.GeneEditorPanel;
import uk.ac.sanger.artemis.components.genebuilder.JExtendedComboBox;

import org.gmod.schema.cv.CvTerm;

/**
 * Panel for display controlled vocabulary terms for Chado
 */
public class CVPanel extends JPanel
             implements FeatureChangeListener
{

  /** */
  private static final long serialVersionUID = 1L;

  private QualifierVector cvQualifiers;
  private Vector editableComponents;
  public static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(CVPanel.class);
  
  private JExtendedComboBox evidenceList;
  
  private JButton hide_show_CC;
  private JButton hide_show_GO;
  private Feature feature;
  
  public CVPanel(final Feature feature)
  {
    super(new BorderLayout());

    updateFromFeature(feature);
  }
  
  /**
   * Return true if this is a CV qualifier
   * @param qualifier
   * @return
   */
  public boolean isCvTag(final Qualifier qualifier)
  {
    return isCvTag(qualifier.getName());
  }
  
  /**
   * Return true if this is a CV qualifier
   * @param qualifierName
   * @return
   */
  public boolean isCvTag(String qualifierName)
  {
    if(qualifierName.startsWith("/"))
      qualifierName = qualifierName.substring(1);
    
    if(qualifierName.equals("product") ||
       qualifierName.equals("controlled_curation") ||
       qualifierName.equals("GO") ||
       qualifierName.equals("class"))
      return true;
    return false;
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
    cv_tags.add("class");
    
    editableComponents = new Vector();
    
    final Dimension dimension  = new Dimension(100, 
        (new JTextField()).getPreferredSize().height);
    final Dimension dimension4 = new Dimension(dimension.width*4, dimension.height);

    Box cvBox = Box.createVerticalBox();
    
    Box xBox = Box.createHorizontalBox();
    JButton addRemove = new JButton("ADD");
    addRemove.setOpaque(false);
    addRemove.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      { 
        addCvTerm();
      }
    });
    xBox.add(addRemove);
    JButton lookUp = new JButton("LOOK UP");
    lookUp.setOpaque(false);
    lookUp.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      { 
        DatabaseDocumentEntry entry =
           (DatabaseDocumentEntry)feature.getEmblFeature().getEntry();
        DatabaseDocument doc = (DatabaseDocument)entry.getDocument();
        
        doc.showCvTermLookUp();
      }
    });
    xBox.add(lookUp);
    xBox.add(Box.createHorizontalGlue());
    cvBox.add(xBox);

    Dimension go_dimension = null;
    int n = 0;
    for(int qualifier_index = 0; qualifier_index < cvQualifiers.size();
        ++qualifier_index) 
    {
      Qualifier this_qualifier = (Qualifier)cvQualifiers.elementAt(qualifier_index);

      if(this_qualifier.getName().equals("GO"))
      {
        final Box yBox = Box.createVerticalBox();
        if(hide_show_GO == null)
          hide_show_GO = new JButton("-");
        
        addHideShowButton(yBox, hide_show_GO);
        
        Box xHeadings = Box.createHorizontalBox();
        cvBox.add(xHeadings);
        
        n++;
        final StringVector qualifier_strings = this_qualifier.getValues();
        
        for(int value_index = 0; value_index < qualifier_strings.size();
            ++value_index)
        {
          final int v_index = value_index;
          
          final String qualifierString = 
            (String)qualifier_strings.elementAt(value_index);
                      
          GoBox go_box = new GoBox(this_qualifier, qualifierString, value_index, 
                                   go_dimension, dimension);
          go_dimension = go_box.getGoDimension();
          editableComponents.add(go_box);
          
          xBox = go_box.getBox();
          xBox.add(Box.createHorizontalGlue());
          xBox.add(getRemoveButton(this_qualifier, v_index));
          
          yBox.add(xBox);
        }
        
        // add column headings
        JLabel lab = new JLabel("GO terms");
        lab.setPreferredSize(go_dimension);
        lab.setFont(lab.getFont().deriveFont(Font.BOLD));
        xHeadings.add(lab);
        final JLabel withLabel = new JLabel("WITH/FROM");
        withLabel.setPreferredSize(dimension);
        xHeadings.add(withLabel);
        
        final JLabel dbxrefLabel = new JLabel("Dbxref");
        dbxrefLabel.setPreferredSize(dimension);
        xHeadings.add(dbxrefLabel);
        
        final JLabel evidenceLabel = new JLabel("Evidence");
        evidenceLabel.setPreferredSize(GoBox.getEvidenceListDimension());
        xHeadings.add(evidenceLabel);
        
        final JLabel qualLabel = new JLabel("Qualifier");
        qualLabel.setPreferredSize(dimension);
        xHeadings.add(qualLabel);
        
        final JLabel dateLabel = new JLabel("Date");
        xHeadings.add(dateLabel);
        
        xHeadings.add(Box.createHorizontalGlue());
        xHeadings.add(hide_show_GO);
              
        // add go rows
        cvBox.add(yBox);
        if(hide_show_GO.getText().equals("+"))
          yBox.setVisible(false);
      }
    }
    
    
    if(n > 0)
      GeneEditorPanel.addLightSeparator(cvBox);
    
    
    n = 0;
    for(int qualifier_index = 0; qualifier_index < cvQualifiers.size();
        ++qualifier_index) 
    {
      final Qualifier this_qualifier = (Qualifier)cvQualifiers.elementAt(qualifier_index);
      if(this_qualifier.getName().equals("controlled_curation"))
      {
        final StringVector qualifier_strings = this_qualifier.getValues();
        
        final Box yBox = Box.createVerticalBox();
        if(hide_show_CC == null)
          hide_show_CC = new JButton("-");
        
        addHideShowButton(yBox, hide_show_CC);
        
        if(n == 0)
        {
          final Box xLabel = Box.createHorizontalBox();
          JLabel lab = new JLabel("Controlled Curation");
          lab.setFont(lab.getFont().deriveFont(Font.BOLD));
          xLabel.add(lab);
          xLabel.add(Box.createHorizontalGlue());
          xLabel.add(hide_show_CC);
          cvBox.add(xLabel);
        
        
          final Box xHeadings = Box.createHorizontalBox();
          yBox.add(xHeadings);
//        add column headings
          final JLabel termLabel = new JLabel("Term");
          
          if(go_dimension != null)
            termLabel.setPreferredSize(
              new Dimension(go_dimension.width+dimension.width,
                            dimension.height));
          else
            termLabel.setPreferredSize(
                new Dimension(dimension.width,
                              dimension.height));
          xHeadings.add(termLabel);
          
          final JLabel dbxrefLabel = new JLabel("Dbxref");
          dbxrefLabel.setPreferredSize(dimension);
          xHeadings.add(dbxrefLabel);
          
          final JLabel evidenceLabel = new JLabel("Evidence");
          evidenceLabel.setPreferredSize(GoBox.getEvidenceListDimension());
          xHeadings.add(evidenceLabel);
          
          final JLabel qualLabel = new JLabel("Qualifier");
          qualLabel.setPreferredSize(dimension);
          xHeadings.add(qualLabel);
          
          final JLabel dateLabel = new JLabel("Date");
          xHeadings.add(dateLabel);
          
          xHeadings.add(Box.createHorizontalGlue());
        }
        
        for(int value_index = 0; value_index < qualifier_strings.size();
            ++value_index)
        {
          n++;
          xBox = Box.createHorizontalBox();
          final String qualifierString = 
             (String)qualifier_strings.elementAt(value_index);
          
          final ControlledCurationBox ccBox = new ControlledCurationBox(
                  this_qualifier,
                  qualifierString, value_index, 
                  dimension, go_dimension);
          editableComponents.add(ccBox);
          
          xBox = ccBox.getBox();
          xBox.add(Box.createHorizontalGlue());
          xBox.add(getRemoveButton(this_qualifier, value_index));         
          yBox.add(xBox);
        }
  
        // add CC rows
        cvBox.add(yBox); 
        if(hide_show_CC.getText().equals("+"))
          yBox.setVisible(false);
      }
    }
    
    if(n > 0)
      GeneEditorPanel.addLightSeparator(cvBox);
    
    //
    // RILEY
    //
    
    if(go_dimension == null)
      go_dimension = new JLabel("product ").getPreferredSize();
    n = 0;
    for(int qualifier_index = 0; qualifier_index < cvQualifiers.size(); ++qualifier_index)
    {
      final Qualifier this_qualifier = (Qualifier) cvQualifiers
          .elementAt(qualifier_index);
      if(this_qualifier.getName().equals("class"))
      {
        final StringVector qualifier_strings = this_qualifier.getValues();
        n++;
        
        for(int value_index = 0; value_index < qualifier_strings.size(); ++value_index)
        {
          final int v_index = value_index;
          
          xBox = Box.createHorizontalBox();
          final String qualifierString = (String) qualifier_strings.elementAt(value_index);

          JLabel label = new JLabel("class");
          if(go_dimension != null)
            label.setPreferredSize(go_dimension);
          xBox.add(label);

          int index = qualifierString.indexOf("::");
          String classStr = qualifierString.substring(0, index);
          JTextField termTextField = new JTextField(classStr);
          termTextField.setOpaque(false);
          termTextField.setEditable(false);
          
          termTextField
              .setToolTipText(DatabaseDocument.getCvTermByCvTermId(
                  Integer.parseInt(qualifierString.substring(index + 2)), feature.getEmblFeature())
                  .getName());
          termTextField.setPreferredSize(dimension);
          termTextField.setMaximumSize(dimension);
          termTextField.setCaretPosition(0);

          xBox.add(termTextField);
          xBox.add(Box.createHorizontalGlue());
          xBox.add(getRemoveButton(this_qualifier, v_index));
          cvBox.add(xBox);

          // Splash.logger4j.debug(this_qualifier.getName());
        }
      }
    }
    
    
    if(n > 0)
      GeneEditorPanel.addLightSeparator(cvBox);
    
    for(int qualifier_index = 0; qualifier_index < cvQualifiers.size();
        ++qualifier_index) 
    {
      final Qualifier this_qualifier = (Qualifier)cvQualifiers.elementAt(qualifier_index);
      if(cv_tags.contains(this_qualifier.getName()))
      {
        final StringVector qualifier_strings = this_qualifier.getValues();

        for(int value_index = 0; value_index < qualifier_strings.size(); ++value_index)
        {
          final int v_index = value_index;

          xBox = Box.createHorizontalBox();
          final String qualifierString = (String) qualifier_strings
              .elementAt(value_index);
          if(this_qualifier.getName().equals("product"))
          {
            JLabel label = new JLabel("product");
            if(go_dimension != null)
              label.setPreferredSize(go_dimension);
            xBox.add(label);
            
            JTextField termTextField = new JTextField(qualifierString);
            termTextField.setOpaque(false);
            termTextField.setEditable(false);
            termTextField.setToolTipText("term column");
            termTextField.setPreferredSize(dimension4);
            termTextField.setMaximumSize(dimension4);
            termTextField.setCaretPosition(0);
            
            xBox.add(termTextField);
            xBox.add(Box.createHorizontalGlue());
            xBox.add(getRemoveButton(this_qualifier, v_index));
            cvBox.add(xBox);
          }
          //Splash.logger4j.debug(this_qualifier.getName());
        }
      }
    }
    
    cvBox.add(Box.createVerticalGlue());
    validate();
    return cvBox;
  }

  /**
   * Add hide/show button to CV section
   * @param box
   */
  private void addHideShowButton(final Box box, final JButton hide_show)
  {
    hide_show.setOpaque(false);
    
    // remove any old listeners
    ActionListener l[] = hide_show.getActionListeners();
    if(l != null)
      for(int i=0;i<l.length;i++)
        hide_show.removeActionListener(l[i]);
    
    hide_show.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        if(hide_show.getText().equals("-"))
        {
          hide_show.setText("+");
          box.setVisible(false);
        }
        else
        {
          hide_show.setText("-");
          box.setVisible(true);
        }
      }
    });
  }
  

  
  private JButton getRemoveButton(final Qualifier this_qualifier, 
                                  final int v_index)
  {
    JButton buttRemove = new JButton("X");
    buttRemove.setOpaque(false);
    Font font = buttRemove.getFont().deriveFont(Font.BOLD);
    buttRemove.setFont(font);
    buttRemove.setToolTipText("REMOVE");
    buttRemove.setForeground(new Color(139,35,35));

    buttRemove.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      { 
        removeCvTerm(this_qualifier.getName(), v_index);
      }
    }); 
    return buttRemove;
  }
  
  public void updateFromFeature(final Feature feature)
  {
    this.feature = feature;
    if(cvQualifiers != null)
      feature.removeFeatureChangeListener(this);
    //cvQualifiers = feature.getQualifiers().copy();
    
    cvQualifiers = new QualifierVector();
    final QualifierVector qualifiers = feature.getQualifiers();  
    for(int i = 0 ; i < qualifiers.size(); ++i) 
    {
      Qualifier this_qualifier = (Qualifier)qualifiers.elementAt(i);
      
      if(this_qualifier.getName().equals("GO"))
      {
        final StringVector qualifier_strings = this_qualifier.getValues();
        // sort by aspect (molecular_function, biological_process, cellular_component)
        Collections.sort(qualifier_strings, new StringVectorComparator());
        this_qualifier = new Qualifier("GO", qualifier_strings);
      }
      
      if(isCvTag(this_qualifier))
        cvQualifiers.addElement(this_qualifier.copy());
    }
   
    feature.addFeatureChangeListener(this);  
    
    removeAll();
    add(createCVQualifiersComponent(),
        BorderLayout.CENTER);
  }
  
  public void updateFromQualifiers(final QualifierVector qualifiers)
  {
    cvQualifiers = new QualifierVector();
   
    for(int i = 0 ; i < qualifiers.size(); ++i) 
    {
      Qualifier this_qualifier = (Qualifier)qualifiers.elementAt(i);
      
      if(this_qualifier.getName().equals("GO"))
      {
        final StringVector qualifier_strings = this_qualifier.getValues();
        // sort by aspect (molecular_function, biological_process, cellular_component)
        Collections.sort(qualifier_strings, new StringVectorComparator());
        this_qualifier = new Qualifier("GO", qualifier_strings);
      }
      
      if(isCvTag(this_qualifier))
        cvQualifiers.addElement(this_qualifier.copy());
    }
    
    removeAll();
    add(createCVQualifiersComponent(),
        BorderLayout.CENTER);
  }
  
  /**
   * Add a CV term to the feature
   */
  private void addCvTerm() 
  {
    final Box xBox = Box.createHorizontalBox();
    final JExtendedComboBox comboCV = 
      new JExtendedComboBox(ChadoTransactionManager.cv_tags);
 
    final java.util.List cvNames = 
      DatabaseDocument.getCvControledCurationNames();
    comboCV.addItem(JExtendedComboBox.SEPARATOR);
    for(int i=0; i<cvNames.size(); i++)
      comboCV.addItem(cvNames.get(i));
    
    JExtendedComboBox term_list = null;
    
    Dimension d = comboCV.getPreferredSize();
    d = new Dimension(80,(int)d.getHeight());
    comboCV.setPreferredSize(d);
    xBox.add(comboCV);
    
    int step = 0;

    CvTerm cvterm  = null;
    String cv_name = null;
    String cv_type = null;
    Vector terms   = null;
    while(step < 5)
    {  
      if(step == 0)
      {
        cv_type = prompCv(xBox, comboCV);
        if(cv_type == null)
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
        final Object obj = promptKeyWord(xBox, cv_name);
        
        if(obj == null)                   // CANCEL
          return;
        else if(!(obj instanceof Vector)) // PREV
        {
          step = 1;
          continue;
        }

        terms = (Vector)obj;
        if(terms.size() < 1)
        {
          JOptionPane.showMessageDialog(this, 
              "No terms found for in the \""+cv_name+"\" controlled vocabulary", 
              "Terms Not Found", JOptionPane.INFORMATION_MESSAGE);
          logger4j.warn("No CV terms found for "+cv_name+
                        " (these may not have been loaded into chado)");
          return;
        }
        
        step = 3;
      }
      
      if(step == 3)
      {
        cvterm = promptCvTerms(xBox, terms, term_list);
        if(cvterm == null)                // CANCEL
          return;
        else if(cvterm.getName() == null) // PREV
        {
          step = 2;
          continue;
        }
        
        if(cv_name.equals("molecular_function") ||
           cv_name.equals("biological_process") ||
           cv_name.equals("cellular_component"))
          step = 4;
        else
          step = 5;
      }
      
      if(step == 4)
      {
        int select = promptEvidence(xBox);
        if(select == 0)       // PREV
        {
          step = 3;
          continue;
        }
        else if(select == 1)  // CANCEL
          return;
        
        step = 5;
      }

    }

    if(!cv_type.equals("GO") &&
       !cv_type.equals("controlled_curation") && 
       !cv_type.equals("product") &&
       !cv_type.equals("class"))
      cv_type = "controlled_curation";
    
    Qualifier cv_qualifier = cvQualifiers.getQualifierByName(cv_type);
    
    final int index;
    if(cv_qualifier == null)
    {
      cv_qualifier = new Qualifier(cv_type);
      index = -1;
    }
    else
     index = cvQualifiers.indexOf(cv_qualifier);
       
    if(cv_type.equals("GO"))
      cv_qualifier.addValue("GOid=GO:"+cvterm.getDbXRef().getAccession()+";"+
          "aspect="+cv_name+";"+
          "term="+cvterm.getName()+";"+
          "evidence="+ GoBox.evidenceCodes[2][ evidenceList.getSelectedIndex() ]);
    else if(cv_type.equals("controlled_curation"))
      cv_qualifier.addValue("term="+cvterm.getName());
    else if(cv_type.equals("product"))
      cv_qualifier.addValue(cvterm.getName());
    else if(cv_type.equals("class"))
      cv_qualifier.addValue(cvterm.getDbXRef().getAccession()+
                            "::"+cvterm.getCvTermId());
    
    if(index > -1)
    {
      cvQualifiers.remove(index);
      cvQualifiers.add(index, cv_qualifier);
    }
    else
      cvQualifiers.add(cv_qualifier);
    
    removeAll();
    add(createCVQualifiersComponent(),
        BorderLayout.CENTER);
    revalidate();
    repaint();
  }
  
  private void removeCvTerm(final String qualifier_name, 
                            final int value_index)
  {
    Qualifier qual = cvQualifiers.getQualifierByName(qualifier_name);
    StringVector values = qual.getValues();
    values.remove(value_index);
    
    cvQualifiers.removeQualifierByName(qualifier_name);
    
    if(values.size() > 0)
    {
      qual = new Qualifier(qualifier_name, values);
      cvQualifiers.addQualifierValues(qual);
    }
    
    removeAll();

    add(createCVQualifiersComponent(),
        BorderLayout.CENTER);
    validate();
  }
  
  /**
   * Get the latest (edited) controlled vocab qualifiers
   * @return
   */
  public QualifierVector getCvQualifiers()
  {
    // check editable components for changes
    
    if(editableComponents != null)
    {
      logger4j.debug("CV checking......");
      for(int i=0; i<editableComponents.size(); i++)
      {
        AbstractCvBox cvBox = (AbstractCvBox)editableComponents.get(i);
        if(cvBox.isQualifierChanged())
        {
          logger4j.debug("CV QUALIFIER CHANGED");
          cvBox.updateQualifier(cvQualifiers);
        }
      }
    }
    
    return cvQualifiers;
  }
  
  /**
   * Prompt for the CV type
   * @param xBox
   * @param comboCV
   * @return
   */
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
    return (String)comboCV.getSelectedItem(); // ChadoTransactionManager.cv_tags[comboCV.getSelectedIndex()];
  }
  
  /**
   * Prompt for the name of the CV
   * @param xBox
   * @param comboCV
   * @return
   */
  private String promptCvName(final Box xBox, final JComboBox comboCV)
  {
    final String options[] = { "<PREV", "CANCEL", "NEXT>"};
     
    String cv_name = (String)comboCV.getSelectedItem();
       //ChadoTransactionManager.cv_tags[comboCV.getSelectedIndex()];
    logger4j.debug("Selected CV is " + cv_name);

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
    else if(cv_name.equals("class"))
      cv_name = DatabaseDocument.RILEY_TAG_CVNAME;
    
    return cv_name;
  }
  
  /**
   * Keyword search of CV terms
   * @param xBox
   * @param cv_name
   * @return
   */
  private Object promptKeyWord(final Box xBox, final String cv_name)
  {
    final String options[] = { "<PREV", "CANCEL", "NEXT>"};
    final JTextField tfield = new JTextField(10);
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
      return new Object();
    }
    
    xBox.remove(tfield);
    
    logger4j.debug("CvTerm cache lookup: "+tfield.getText().trim()+" from "+cv_name);
    return DatabaseDocument.getCvterms(tfield.getText().trim(), cv_name);
  }
  
  /**
   * Prompt for CV term
   * @param xBox
   * @param terms
   * @return
   */
  private CvTerm promptCvTerms(final Box xBox, final Vector terms,
                               JExtendedComboBox term_list)
  {
    final String options[] = { "<PREV", "CANCEL", "NEXT>"};   
    
    if(term_list == null)
    {
      Collections.sort(terms, new CvTermsComparator());
      term_list = new JExtendedComboBox(terms);
    
      Dimension d = term_list.getPreferredSize();
      d = new Dimension(280,(int)d.getHeight());
      term_list.setPreferredSize(d);
    
      term_list.setRenderer(new CVCellRenderer());
    }
    else
      xBox.remove(term_list);
    
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
  
  private int promptEvidence(final Box xBox)
  {
    final String options[] = { "<PREV", "CANCEL", "NEXT>"};
    evidenceList = new JExtendedComboBox(GoBox.evidenceCodes[1]);
    evidenceList.setSelectedItem("NR \t:: not recorded");
    xBox.add(evidenceList);
    
    int select = JOptionPane.showOptionDialog(null, xBox,
        "CV term selection",
         JOptionPane.YES_NO_CANCEL_OPTION,
         JOptionPane.QUESTION_MESSAGE,
         null,
         options,
         options[2]);
    
    if(select == 0)
      xBox.remove(evidenceList);
    
    return select;
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
  
  /**
   * Cell renderer for CvTerm's
   */
  class CVCellRenderer extends JLabel implements ListCellRenderer 
  {
    /** */
    private static final long serialVersionUID = 1L;

    public CVCellRenderer() 
    {
      super();
    }

    public Component getListCellRendererComponent(JList list,
                                                  Object value,
                                                  int index,
                                                  boolean isSelected,
                                                  boolean cellHasFocus) 
    {
      setText(((CvTerm)value).getName());
      //setToolTipText(((CvTerm)value).getCv().getName());
      return this;
    }
  }
  
  class StringVectorComparator implements Comparator
  {
    public StringVectorComparator()
    {
    }

    public final int compare(Object a, Object b)
    {
      int result;

      String strA = getField("aspect", (String)a);
      String strB = getField("aspect", (String)b);
      result = strA.compareTo(strB);

      return result;
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
      String field = "";
      
      int ind1 = qualifierString.toLowerCase().indexOf(fieldName.toLowerCase());
      int ind2 = qualifierString.indexOf(";", ind1);
      
      int len = fieldName.length();

      if(ind2 > ind1 && ind1 > -1)
        field = qualifierString.substring(ind1+len,ind2);
      else if(ind1 > -1)
        field = qualifierString.substring(ind1+len);
      
      return field;
    }
  }
  
}
