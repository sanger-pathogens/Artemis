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
import java.awt.FontMetrics;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Collections;
import java.util.Comparator;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.JPanel;
import javax.swing.JButton;

import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.chado.ChadoTransactionManager;
import uk.ac.sanger.artemis.components.genebuilder.GeneEditorPanel;

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
  private Vector<AbstractCvBox> editableComponents;
  public static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(CVPanel.class);
  
  private JButton hide_show_CC;
  private JButton hide_show_GO;
  private Feature feature;
  
  private DatabaseDocument doc;
  //used to test if cv panel has contents
  private boolean empty = true;
  
  public CVPanel(final Feature feature)
  {
    super(new BorderLayout());
    updateFromFeature(feature);
    doc = (DatabaseDocument)
      ((GFFStreamFeature)feature.getEmblFeature()).getDocumentEntry().getDocument();
  }
  
  /**
   * Return true if this is a CV qualifier
   * @param qualifier
   * @return
   */
  public static boolean isCvTag(final Qualifier qualifier)
  {
    return isCvTag(qualifier.getName());
  }
  
  /**
   * Return true if this is a CV qualifier
   * @param qualifierName
   * @return
   */
  public static boolean isCvTag(String qualifierName)
  {
    if(qualifierName.startsWith("/"))
      qualifierName = qualifierName.substring(1);
    
    return ChadoTransactionManager.isCvTag(qualifierName);
  }
  
  /**
   * Create CV qualifier components
   * @return
   */
  private Component createCVQualifiersComponent()
  {
    empty = true;
    Vector<String> cv_tags = new Vector<String>();
    cv_tags.add("GO");
    cv_tags.add("controlled_curation");
    cv_tags.add("product");
    cv_tags.add("class");
    
    editableComponents = new Vector<AbstractCvBox>();
    
    final Dimension dimension  = new Dimension(100, 
        (new JTextField()).getPreferredSize().height);

    Box cvBox = Box.createVerticalBox();
    
    Box xBox = Box.createHorizontalBox();
    JButton addRemove = new JButton("ADD");
    addRemove.setOpaque(false);
    addRemove.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {        
        CvTermSelector cvTermSelector = new CvTermSelector();
        if(cvTermSelector.showOptions(null, doc))
        {
          addCvTermQualifier(cvTermSelector.getCvTerm(),
              cvTermSelector.getCvName(),
              cvTermSelector.getCv(),
              cvTermSelector.getEvidenceCode());
        }
        
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
    
    JButton addHistory = new JButton("ADD HISTORY");
    addHistory.setOpaque(false);
    addHistory.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      { 
        addHistoryQualifier();
      }
    });
    if(ChadoTransactionManager.HISTORY_CV != null)
      xBox.add(addHistory);
    
    xBox.add(Box.createHorizontalGlue());
    cvBox.add(xBox);

    Dimension go_dimension = null;
    Box goBox = null;
    Box goHeadings = null;
    for(Qualifier this_qualifier: cvQualifiers) 
    {
      if(this_qualifier.getName().equals("GO"))
      {
        empty = false;

        if(hide_show_GO == null)
          hide_show_GO = new JButton("-");
        
        if(goBox == null)
          goBox = Box.createVerticalBox();
        
        addHideShowButton(goBox, hide_show_GO);
        
        if(goHeadings == null)
          goHeadings = Box.createHorizontalBox();

        final StringVector qualifier_strings = this_qualifier.getValues();
        for(int i = 0; i < qualifier_strings.size(); ++i)
        {
          GoBox go_box = new GoBox(this_qualifier, qualifier_strings.elementAt(i), i, 
                                   go_dimension, dimension);
          go_dimension = go_box.getGoDimension();
          editableComponents.add(go_box);
          
          xBox = go_box.getBox();
          xBox.add(Box.createHorizontalGlue());
          xBox.add(getRemoveButton(this_qualifier, i));
          
          goBox.add(xBox);
        }
        
        // add column headings
        final JLabel lab = new JLabel("GO terms");
        lab.setPreferredSize(go_dimension);
        lab.setFont(lab.getFont().deriveFont(Font.BOLD));
        goHeadings.add(lab);
        final JLabel withLabel = new JLabel("WITH/FROM");
        withLabel.setPreferredSize(dimension);
        goHeadings.add(withLabel);
        
        final JLabel dbxrefLabel = new JLabel("Dbxref");
        dbxrefLabel.setPreferredSize(dimension);
        goHeadings.add(dbxrefLabel);
        
        final JLabel evidenceLabel = new JLabel("Evidence");
        evidenceLabel.setPreferredSize(GoBox.getEvidenceListDimension());
        goHeadings.add(evidenceLabel);
        
        final JLabel qualLabel = new JLabel("Qualifier");
        qualLabel.setPreferredSize(dimension);
        goHeadings.add(qualLabel);
        
        final JLabel sourceLabel = new JLabel("Source");
        sourceLabel.setPreferredSize(GoBox.getSourceListDimension());
        goHeadings.add(sourceLabel);
        
        goHeadings.add(new JLabel("Date"));
        
        goHeadings.add(Box.createHorizontalGlue());
        goHeadings.add(hide_show_GO);
              
        if(hide_show_GO.getText().equals("+"))
          goBox.setVisible(false);
      }
    }
    
    int n = 0;
    for(Qualifier this_qualifier: cvQualifiers) 
    {
      if(this_qualifier.getName().equals("product"))
      {
        final StringVector qualifier_strings = this_qualifier.getValues();
        for(int i = 0; i < qualifier_strings.size(); ++i)
        {
          empty = false;
          final ProductBox productBox = new ProductBox(
                    this_qualifier, qualifier_strings.elementAt(i), 
                    i, dimension, go_dimension);
          editableComponents.add(productBox);

          if(qualifier_strings.size() == 1)
            productBox.getRecommended().setEnabled(false);
          xBox = productBox.getBox();
          xBox.add(Box.createHorizontalGlue());
          xBox.add(getRemoveButton(this_qualifier, i));         
          n++;
          cvBox.add(productBox.getHeadingsBox()); 
          cvBox.add(xBox); 
        }
      }
    }
    
    if(n > 0)
      GeneEditorPanel.addLightSeparator(cvBox);
    
    // history field
    n = 0;
    for(Qualifier this_qualifier: cvQualifiers) 
    {
      if(this_qualifier.getName().equals("history"))
      {
        final StringVector qualifier_strings = this_qualifier.getValues();

        for(int i = 0; i < qualifier_strings.size(); ++i)
        {
          empty = false;
           
          final HistoryBox historyBox = new HistoryBox(
                    this_qualifier,
                    qualifier_strings.elementAt(i), i, 
                    dimension, go_dimension);
          editableComponents.add(historyBox);
            
          xBox = historyBox.getBox();
          xBox.add(Box.createHorizontalGlue());
          xBox.add(getRemoveButton(this_qualifier, i));
          
          if(n == 0)
          {
            final Box xLabel = Box.createHorizontalBox();
            JLabel lab = new JLabel("History");
            lab.setFont(lab.getFont().deriveFont(Font.BOLD));
            xLabel.add(lab);
            xLabel.add(Box.createHorizontalGlue());
            cvBox.add(xLabel); 
          }
          n++;
          cvBox.add(xBox); 
        }
      }
    }
    
    if(n > 0)
      GeneEditorPanel.addLightSeparator(cvBox);
    
    //
    if(goBox != null)
    {
      cvBox.add(goHeadings);
      cvBox.add(goBox);
      GeneEditorPanel.addLightSeparator(cvBox);
    }


    n = 0;
    for(Qualifier this_qualifier: cvQualifiers) 
    {
      if(this_qualifier.getName().equals("controlled_curation"))
      {
        empty = false;

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
                new Dimension(getWidthOfGoField()+dimension.width,
                              dimension.height));
          xHeadings.add(termLabel);


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
        }

        final StringVector qualifier_strings = this_qualifier.getValues();
        for(int i = 0; i < qualifier_strings.size(); ++i)
        {
          n++;

          final ControlledCurationBox ccBox = new ControlledCurationBox(
                  this_qualifier,
                  qualifier_strings.elementAt(i), i, 
                  dimension, go_dimension);
          editableComponents.add(ccBox);

          xBox = ccBox.getBox();
          xBox.add(Box.createHorizontalGlue());
          xBox.add(getRemoveButton(this_qualifier, i));         
          yBox.add(xBox);
          yBox.add(Box.createVerticalStrut(2));
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
    for(Qualifier this_qualifier: cvQualifiers)
    {
      if(this_qualifier.getName().equals("class"))
      {
        empty = false;
        final StringVector qualifier_strings = this_qualifier.getValues();
        n++;
        
        for(int i = 0; i < qualifier_strings.size(); ++i)
        {
          xBox = Box.createHorizontalBox();
          final String qualifierStr = qualifier_strings.elementAt(i);

          JLabel label = new JLabel("class");
          if(go_dimension != null)
            label.setPreferredSize(go_dimension);
          xBox.add(label);

          int index = qualifierStr.indexOf("::");
          String classStr = qualifierStr.substring(0, index);
          JLabel termTextField = new JLabel(classStr);
          termTextField.setOpaque(false);

          termTextField
              .setToolTipText(DatabaseDocument.getCvTermByCvTermId(
                  Integer.parseInt(qualifierStr.substring(index + 2)), feature.getEmblFeature())
                  .getDefinition());
          termTextField.setPreferredSize(dimension);
          termTextField.setMaximumSize(dimension);

          xBox.add(termTextField);
          xBox.add(Box.createHorizontalGlue());
          xBox.add(getRemoveButton(this_qualifier, i));
          cvBox.add(xBox);
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


  private static int getWidthOfGoField()
  {
    JTextField textField = new JTextField();
    FontMetrics fm  = textField.getFontMetrics(textField.getFont());
    return fm.stringWidth("GO:0001234 [F] ");  
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
   * Add a history qualifier
   */
  private void addHistoryQualifier()
  {
    cvQualifiers = getCvQualifiers();
    Qualifier cv_qualifier = cvQualifiers.getQualifierByName("history");
    
    final int index;
    if(cv_qualifier == null)
    {
      cv_qualifier = new Qualifier("history");
      index = -1;
    }
    else
     index = cvQualifiers.indexOf(cv_qualifier);
    
    cv_qualifier.addValue(
        "term="+HistoryBox.getDefaultTerm().getName()+";"+
        "curatorName="+doc.getUserName()+";"+
        "date="+ DatePanel.getDate());
    
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

  /**
   * Add CvTerm qualifier
   * @param cvterm
   * @param go_name
   * @param cv_type
   */
  private void addCvTermQualifier(CvTerm cvterm, 
                                  final String go_name,
                                  String cv_type,
                                  final String evidenceCode)
  {
    if(!cv_type.equals("GO") &&
       !cv_type.equals("controlled_curation") && 
       !cv_type.equals("product") &&
       !cv_type.equals("class"))
      cv_type = "controlled_curation";

    cvQualifiers = getCvQualifiers();
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
          "aspect="+go_name+";"+
          "term="+cvterm.getName()+";"+
          "evidence="+ evidenceCode+";"+
          "date="+ DatePanel.getDate());
    else if(cv_type.equals("controlled_curation"))
      cv_qualifier.addValue("term="+cvterm.getName()+";cv="+cvterm.getCv().getName());
    else if(cv_type.equals("product"))
      cv_qualifier.addValue("term="+cvterm.getName());
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
    cvQualifiers = getCvQualifiers();
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

  public static String getDescription()
  {
    final StringBuffer buff = new StringBuffer(); 
    buff.append("Section for controlled vocabulary qualifiers:\n");
    buff.append("GO, controlled_curation, product, Riley class");
    
    return buff.toString();
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
  
  
  class StringVectorComparator implements Comparator<String>
  {
    public final int compare(String a, String b)
    {
      int result;

      String strA = getField("aspect", a);
      String strB = getField("aspect", b);
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

  public boolean isEmpty()
  {
    return empty;
  }

  public void setEmpty(boolean empty)
  {
    this.empty = empty;
  }

}
