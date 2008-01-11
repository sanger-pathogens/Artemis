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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/gff/GffPanel.java,v 1.13 2008-01-11 14:17:13 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder.gff;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.SwingConstants;

import org.gmod.schema.cv.CvTerm;


import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.chado.ChadoTransactionManager;
import uk.ac.sanger.artemis.components.genebuilder.GeneEditorPanel;
import uk.ac.sanger.artemis.components.genebuilder.JExtendedComboBox;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.StringVector;

public class GffPanel extends JPanel
                      implements FeatureChangeListener
{
  private static final long serialVersionUID = 1L;
  private QualifierVector gffQualifiers;
  private JTextField uniquenameTextField;
  private JTextField timeTextField;
  private Feature feature;
  private boolean empty = true;
  
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
  public boolean isPropertiesTag(final Qualifier qualifier, final Feature feature)
  {
    if(qualifier.getName().equals("ID") ||
       qualifier.getName().equals("feature_id") ||
       qualifier.getName().equals("Parent") ||
       qualifier.getName().equals("Derives_from") ||
       qualifier.getName().equals("feature_relationship_rank") ||
       qualifier.getName().equals("timelastmodified") ||
       ChadoTransactionManager.isSynonymTag(qualifier.getName(), 
           (GFFStreamFeature)feature.getEmblFeature()))
      return true;
    return false;
  }
  
  private Component createGffQualifiersComponent()
  {
    empty = true;
    int nrows = 0;
    
    Qualifier idQualifier          = gffQualifiers.getQualifierByName("ID");
    Qualifier parentQualifier      = gffQualifiers.getQualifierByName("Parent");
    Qualifier derivesFromQualifier = gffQualifiers.getQualifierByName("Derives_from");
    Qualifier timeQualifier        = gffQualifiers.getQualifierByName("timelastmodified");
    
    Box gffBox = Box.createVerticalBox();
    gffBox.add(Box.createVerticalStrut(5));
    
    GridBagLayout grid = new GridBagLayout();
    GridBagConstraints c = new GridBagConstraints();
    
    int maxLabelWidth = new JLabel("temporary_systematic_id ").getPreferredSize().width;

    c.ipady = 3;
    JPanel gridPanel = new JPanel(grid);
    gridPanel.setBackground(Color.WHITE);
    
    //
    // create components for synonym qualifiers
    int maxSynonymWidth = 0;
    for(int i=0; i<gffQualifiers.size(); i++)
    {
      final Qualifier qualifier = (Qualifier)gffQualifiers.get(i);
      if( ChadoTransactionManager.isSynonymTag(qualifier.getName(), 
          (GFFStreamFeature)feature.getEmblFeature()) &&
          isSystematicId(qualifier.getName()))
      {
        addSynonymComponent(qualifier, c, gridPanel, nrows,
            maxLabelWidth, maxSynonymWidth);  
        nrows++;
      }
    }
    
    boolean firstFound = true;
    for(int i=0; i<gffQualifiers.size(); i++)
    {
      final Qualifier qualifier = (Qualifier)gffQualifiers.get(i);
      if( ChadoTransactionManager.isSynonymTag(qualifier.getName(), 
          (GFFStreamFeature)feature.getEmblFeature()) &&
          !isSystematicId(qualifier.getName()))
      {
        if(firstFound && nrows > 0)
        {
          addJSeparator(gridPanel, nrows, maxSynonymWidth);
          nrows++;
        }
        addSynonymComponent(qualifier, c, gridPanel, nrows, 
            maxLabelWidth, maxSynonymWidth);  
        nrows++;
        firstFound = false;
      }
    }
    
    Dimension cellDimension = null;
    nrows = 0;
    if(idQualifier != null)
    {
      final String uniquename = (String)idQualifier.getValues().get(0);
      JLabel idField = new JLabel("Internal ID");
      
      uniquenameTextField = new JTextField(uniquename);
      cellDimension = new Dimension(uniquenameTextField.getPreferredSize().width+10,
                                    idField.getPreferredSize().height+10);
      
      if(feature.getKey().getKeyString().indexOf("exon") > -1)
        uniquenameTextField.setEditable(false);
      uniquenameTextField.setMaximumSize(cellDimension);
      
      c.gridx = 3;
      c.gridy = 0;
      c.ipadx = 5;
      c.fill = GridBagConstraints.NONE;
      c.anchor = GridBagConstraints.NORTHEAST;
      gridPanel.add(idField, c);
      c.gridx = 4;
      c.gridy = 0;
      c.ipadx = 0;
      c.fill = GridBagConstraints.HORIZONTAL;
      c.anchor = GridBagConstraints.NORTHWEST;
      gridPanel.add(uniquenameTextField, c);
      
      Qualifier featIdQualifier = gffQualifiers.getQualifierByName("feature_id");
      if(featIdQualifier != null)
      {
        
        idField.setToolTipText("feature_id="+(String)featIdQualifier.getValues().get(0));
        uniquenameTextField.setToolTipText("feature_id="+(String)featIdQualifier.getValues().get(0));
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
        
        c.gridx = 3;
        c.gridy = nrows;
        c.ipadx = 5;
        c.fill = GridBagConstraints.NONE;
        c.anchor = GridBagConstraints.NORTHEAST;
        gridPanel.add(parentField, c); 
        c.gridx = 4;
        c.gridy = nrows;
        c.ipadx = 0;
        c.fill = GridBagConstraints.HORIZONTAL;
        c.anchor = GridBagConstraints.NORTHWEST;
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
        
        c.gridx = 3;
        c.gridy = nrows;
        c.ipadx = 5;
        c.fill = GridBagConstraints.NONE;
        c.anchor = GridBagConstraints.NORTHEAST;
        gridPanel.add(derivesFromsField, c); 
        c.gridx = 4;
        c.gridy = nrows;
        c.ipadx = 0;
        c.fill = GridBagConstraints.HORIZONTAL;
        c.anchor = GridBagConstraints.NORTHWEST;
        gridPanel.add(derivesFromTextField, c); 
        nrows++;
      }
    }
    
    
    if(timeQualifier != null)
    {
      String time = (String)timeQualifier.getValues().get(0);
      
      JLabel timeField = new JLabel("timelastmodified");
      //timeField.setPreferredSize(new Dimension(timeField.getPreferredSize().width+10,
      //                 timeField.getPreferredSize().height));
      
      timeTextField = new JTextField(time);
      if(cellDimension == null ||
         cellDimension.width < timeTextField.getPreferredSize().width+10)
         cellDimension = new Dimension(timeTextField.getPreferredSize().width+10,
                                       timeField.getPreferredSize().height+10);
      timeTextField.setMaximumSize(cellDimension);
      timeTextField.setEditable(false);
      
      c.gridx = 3;
      c.gridy = nrows;
      c.ipadx = 5;
      c.fill = GridBagConstraints.NONE;
      c.anchor = GridBagConstraints.NORTHEAST;
      gridPanel.add(timeField, c);
      c.gridx = 4;
      c.gridy = nrows;
      c.ipadx = 0;
      c.fill = GridBagConstraints.HORIZONTAL;
      c.anchor = GridBagConstraints.NORTHWEST;
      gridPanel.add(timeTextField, c);
      nrows++;
    }  
    
    gffBox.add(gridPanel);
    
    
    // add/remove buttons
    Box xBox = Box.createHorizontalBox();

    final JButton addSynonym = new JButton("ADD");
    addSynonym.setOpaque(false);
    addSynonym.setToolTipText("Add id or synonym");
    addSynonym.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      { 
        addSynonym();
      }
    });
    xBox.add(addSynonym);
    xBox.add(Box.createHorizontalStrut(10));
    
    final JButton removeSynonym = new JButton("REMOVE");
    removeSynonym.setOpaque(false);
    removeSynonym.setToolTipText("Remove id or synonym");
    removeSynonym.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      { 
        removeSynonym();
      }
    });
    xBox.add(removeSynonym);
    
    xBox.add(Box.createHorizontalGlue());
    gffBox.add(xBox);
    
    return gffBox;
  }
  
  public void updateFromFeature(final Feature feature)
  {
    this.feature = feature;
    removeAll();
    if(gffQualifiers != null)
      feature.removeFeatureChangeListener(this);
    gffQualifiers = feature.getQualifiers().copy();
    
    gffQualifiers = new QualifierVector();
    final QualifierVector qualifiers = feature.getQualifiers();  
    for(int i = 0 ; i < qualifiers.size(); ++i) 
    {
      Qualifier qualifier = (Qualifier)qualifiers.elementAt(i);
      if(isPropertiesTag(qualifier, feature))
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
  public QualifierVector getGffQualifiers(final Feature feature)
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
  
  private boolean isSystematicId(final String synonymType)
  {
    if(synonymType.indexOf("systematic_id") > -1)
      return true;
    return false;
  }
  
  private void removeSynonym()
  {
    Vector synonymTypes = new Vector();
    for(int i=0; i<gffQualifiers.size(); i++)
    {
      Qualifier qualifier = (Qualifier) gffQualifiers.get(i);
      if(ChadoTransactionManager.isSynonymTag(qualifier.getName(), 
         (GFFStreamFeature)feature.getEmblFeature()))
        synonymTypes.add(qualifier.getName());
    }
    
    final JExtendedComboBox list = new JExtendedComboBox(synonymTypes);
    final String options[] = { "CANCEL", "NEXT>"};   
    
    int select = JOptionPane.showOptionDialog(null, list,
        "Select synonym type",
         JOptionPane.YES_NO_CANCEL_OPTION,
         JOptionPane.QUESTION_MESSAGE,
         null,
         options,
         options[1]);
    
    if(select == 0)
      return;
    
    Box xBox = Box.createHorizontalBox();
    final String synonymName = (String) list.getSelectedItem();
    final JLabel name = new JLabel( synonymName );
    xBox.add(name);
    
    Qualifier qualifier = gffQualifiers.getQualifierByName(synonymName);
    StringVector values = qualifier.getValues();
    final JExtendedComboBox valueList = new JExtendedComboBox(values);
    xBox.add(valueList);
    
    select = JOptionPane.showConfirmDialog(null, xBox, 
        "Input name", JOptionPane.OK_CANCEL_OPTION);
    
    if(select == JOptionPane.CANCEL_OPTION)
      return;
    
    if(values.size()==1)
      gffQualifiers.removeQualifierByName(synonymName);
    else
    {
      int index = gffQualifiers.indexOfQualifierWithName(synonymName);
      values.remove(valueList.getSelectedIndex());
      gffQualifiers.remove(index); 
      gffQualifiers.add(index, new Qualifier(synonymName, values));
    }
    
    removeAll();
    add(createGffQualifiersComponent());
    revalidate();
  }
  
  private void addSynonym()
  {
    final Vector synonyms = DatabaseDocument.getCvterms("", "genedb_synonym_type");
    final JExtendedComboBox list = new JExtendedComboBox(synonyms);
    final String options[] = { "CANCEL", "NEXT>"};   
    
    int select = JOptionPane.showOptionDialog(null, list,
        "Select synonym type",
         JOptionPane.YES_NO_CANCEL_OPTION,
         JOptionPane.QUESTION_MESSAGE,
         null,
         options,
         options[1]);
    
    if(select == 0)
      return;
    
    Box xBox = Box.createHorizontalBox();
    final String synonymName = ((CvTerm)list.getSelectedItem()).getName();
    final JLabel name = new JLabel( synonymName );
    xBox.add(name);
    
    final JTextField newSynonym = new JTextField(15);
    xBox.add(newSynonym);
    
    final JCheckBox current = new JCheckBox("make current", true);
    
    if(isSystematicId(synonymName))
      xBox.add(current);
    
    select = JOptionPane.showConfirmDialog(null, xBox, 
        "Input name", JOptionPane.OK_CANCEL_OPTION);
    
    if(select == JOptionPane.CANCEL_OPTION || newSynonym.getText().equals(""))
      return;
    
    Qualifier synonymQualifier = gffQualifiers.getQualifierByName(synonymName);
    
    String newSynonymValue = newSynonym.getText();
    if(isSystematicId(synonymName) && !current.isSelected())
      newSynonymValue = newSynonymValue + ";current=false";
    
    if(synonymQualifier == null)
    {
      synonymQualifier = new Qualifier(synonymName, newSynonymValue);
      gffQualifiers.add(synonymQualifier);
    }
    else
      synonymQualifier.addValue(newSynonymValue);
    
    final StringVector values = synonymQualifier.getValues();
    final StringVector newValues = new StringVector();
    for(int i=0; i<values.size(); i++)
    {
      String thisValue = (String) values.get(i);
      StringVector str = StringVector.getStrings(thisValue, ";");
      String synonymValue = (String) str.get(0);
      
      if(isSystematicId(synonymName) && !synonymValue.equals(newSynonymValue))
        if(current.isSelected() && !thisValue.endsWith(";current=false"))
          thisValue = thisValue + ";current=false";
      
      newValues.add(thisValue);
    }

    int index = gffQualifiers.indexOfQualifierWithName(synonymName);
    if(index == -1)
      gffQualifiers.setQualifier(new Qualifier(synonymName,newValues));
    else
    {
      gffQualifiers.remove(index);
      gffQualifiers.add(index, new Qualifier(synonymName,newValues));
    }
    removeAll();
    add(createGffQualifiersComponent());
    revalidate();
  }
  
  private void addJSeparator(final JPanel gridPanel, 
                             final int nrows, 
                             final int maxLabelWidth)
  {
    GridBagConstraints c = new GridBagConstraints();
    c.gridx = 0;
    c.gridy = nrows;
    c.ipady = 6;
    c.gridwidth = 2;
    c.anchor = GridBagConstraints.CENTER;
    c.fill   = GridBagConstraints.HORIZONTAL;
    JSeparator separator = GeneEditorPanel.getSeparator(gridPanel, false);
    separator.setPreferredSize(new Dimension(maxLabelWidth, separator
        .getPreferredSize().height));
    gridPanel.add(separator, c);
  }
  
  /**
   * Add the synonym components
   * @param qualifier
   * @param c
   * @param gridPanel
   * @param nrows
   * @param maxLabelWidth
   * @param maxSynonymWidth
   * @return
   */
  private void addSynonymComponent(final Qualifier qualifier, 
                                   final GridBagConstraints c, 
                                   final JPanel gridPanel, 
                                   final int nrows, 
                                   final int maxLabelWidth,
                                   int maxSynonymWidth)
  {
    empty = false;
    int current = 0;
    final StringVector values = qualifier.getValues();
    final Vector featureSynonym = new Vector();
    for(int j=0; j<values.size(); j++)
    {
      String value = (String) values.get(j);
      StringVector strings = StringVector.getStrings(value, ";");
      
      if(strings.size()>1)
        featureSynonym.add((String) strings.get(0));
      else
        featureSynonym.add(0, value);
    }
    
    final JLabel sysidField = new JLabel(qualifier.getName()+" ");
    sysidField.setPreferredSize(new Dimension(maxLabelWidth, 
                     sysidField.getPreferredSize().height));
    sysidField.setHorizontalAlignment(SwingConstants.RIGHT);
    
    c.gridx = 0;
    c.gridy = nrows;
    c.ipadx = 5;
    c.fill = GridBagConstraints.NONE;
    c.anchor = GridBagConstraints.NORTHEAST;
    gridPanel.add(sysidField, c);

    c.gridx = 1;
    c.gridy = nrows;
    c.ipadx = 5;
    c.fill = GridBagConstraints.HORIZONTAL;
    c.anchor = GridBagConstraints.NORTHWEST;
    
    if(isSystematicId(qualifier.getName()))
    {
      final JExtendedComboBox comboBox = new JExtendedComboBox(
                                               featureSynonym);
      comboBox.setHighLightCurrent(true);
      comboBox.setSelectedIndex(current);
      comboBox.setCurrent(current);
      comboBox.setOpaque(false);
      comboBox.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          comboBoxAction(comboBox, qualifier);
        }
      });
      
      final String tt = "current "+qualifier.getName()+" is shown";
      sysidField.setToolTipText(tt);
      comboBox.setToolTipText(tt);
      
      gridPanel.add(comboBox, c);
    }
    else
    { 
      final DefaultListModel listModel = new DefaultListModel();
      for(int k = featureSynonym.size()-1; k > -1; k--)
        listModel.addElement(featureSynonym.get(k));
      
      final JList name_list = new JList(listModel);
      
      if(featureSynonym.size() > 5)
      {
        final JScrollPane jsp = new JScrollPane(name_list);
        gridPanel.add(jsp, c);
      }
      else
        gridPanel.add(name_list, c);
    }
    
    c.gridx = 2;
    c.gridy = nrows;
    
    gridPanel.add(Box.createHorizontalStrut(25), c);
  }
  
  private void comboBoxAction(final JExtendedComboBox comboBox,
                              final Qualifier qualifier)
  {
    if(comboBox.getSelectedIndex() == comboBox.getCurrent())
      return;
    
    int select = JOptionPane.showConfirmDialog(null, 
        "Change the current value of "+
        qualifier.getName()+" to "+comboBox.getSelectedItem()+"?", 
        qualifier.getName()+" change", 
        JOptionPane.OK_CANCEL_OPTION);
    if(select == JOptionPane.CANCEL_OPTION)
      comboBox.setSelectedIndex(comboBox.getCurrent());
    else
    {
      StringVector values = qualifier.getValues();
      StringVector newValues = new StringVector();
      String selectedValue = (String) comboBox.getSelectedItem();
      
      for(int i =0; i<values.size(); i++)
      {
        String value = (String) values.get(i);
        StringVector strings = StringVector.getStrings(value, ";");
        if(selectedValue.equals(strings.get(0)))
          value = selectedValue;
        else
          value = strings.get(0)+";current=false";
        
        newValues.add(value);
      }
      int index = gffQualifiers.indexOfQualifierWithName(qualifier.getName());
      gffQualifiers.remove(index);
      gffQualifiers.add(index, new Qualifier(qualifier.getName(), newValues));
      comboBox.setCurrent(comboBox.getSelectedIndex());
    }  
  }
  
  public void featureChanged(FeatureChangeEvent event)
  {
    updateFromFeature(event.getFeature());
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