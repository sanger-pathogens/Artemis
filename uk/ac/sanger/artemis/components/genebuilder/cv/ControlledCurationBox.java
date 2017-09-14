/* ControlledCurationBox.java
 *
 * created: 2007
 *
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
 **/

package uk.ac.sanger.artemis.components.genebuilder.cv;

import java.awt.Dimension;
import java.awt.FontMetrics;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.JTextField;

import org.gmod.schema.cv.CvTerm;

import uk.ac.sanger.artemis.components.genebuilder.JExtendedComboBox;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.StringVector;

class ControlledCurationBox extends AbstractCvBox
{
  private Box xBox;
  private int value_index;
  private WrapTextArea termTextField;
  private JTextField withTextField;
  private JTextField dbxrefTextField;
  private JExtendedComboBox evidenceList;
  private JTextField qualfTextField;
  private DatePanel dateField;
  private String origQualifierString;
  private Qualifier origQualifier;
  
  public ControlledCurationBox(final Qualifier qualifier,
                               final String qualifierString,
                               final int value_index,
                               final Dimension dimension,
                               final Dimension go_dimension)
  {
    this.origQualifier = qualifier;
    this.origQualifierString = qualifierString;
    this.value_index  = value_index;
    this.xBox = Box.createHorizontalBox();
    
    //JLabel cclabel = new JLabel("controlled_curation");
    //if(go_dimension != null)
    //  cclabel.setPreferredSize(go_dimension);
      
    //xBox.add(cclabel);
    
    final String term = getField("term=", qualifierString);
    final String cvName = getField("cv=", qualifierString);

    CvTerm cvTerm = DatabaseDocument.getCvTermByCvPartAndCvTerm(term,cvName);
    if(cvTerm == null)
      cvTerm = DatabaseDocument.getCvTermByCvPartAndCvTerm(term,"CC");
    termTextField = new WrapTextArea(term, go_dimension, dimension.width);
    termTextField.setOpaque(false);
    termTextField.setEditable(false);
    termTextField.setToolTipText(cvTerm.getCv().getName());
    
    final Dimension d;
    if(go_dimension != null)
    {
      d = new Dimension(go_dimension.width+dimension.width,
          termTextField.getPreferredSize().height);
    }
    else
    {
      FontMetrics fm  = termTextField.getFontMetrics(termTextField.getFont());
      int width= fm.stringWidth("GO:0001234 [F] ");
      d = new Dimension(width+dimension.width,
          termTextField.getPreferredSize().height);
    }
    termTextField.setPreferredSize(d);
    termTextField.setMaximumSize(d);
    
    termTextField.setCaretPosition(0);
    xBox.add(termTextField);
    
    //
    String with = getField("with=", qualifierString);
    withTextField = new JTextField(with);
    withTextField.setToolTipText("with/from column");
    withTextField.setPreferredSize(dimension);
    withTextField.setMaximumSize(dimension);
    withTextField.setActionCommand("with=");
    xBox.add(withTextField);
    
    String dbxref = getField("db_xref=", qualifierString);
    dbxrefTextField = new JTextField(dbxref);
    dbxrefTextField.setToolTipText("dbxref column");
    dbxrefTextField.setPreferredSize(dimension);
    dbxrefTextField.setMaximumSize(dimension);
    dbxrefTextField.setActionCommand("db_xref=");
    xBox.add(dbxrefTextField);
 
    // feature_cvterm_prop's
    String evidence = getField("evidence=", qualifierString);
    
    evidenceList = new JExtendedComboBox(GoBox.evidenceCodes[1]);
    evidenceList.setOpaque(false);
    evidenceList.setToolTipText("evidence column");
    evidenceList.setSelectedIndex( GoBox.getEvidenceIndex(evidence) );
  
    Dimension d2 = evidenceList.getPreferredSize();
    d2 = new Dimension(90,(int)d2.getHeight());
    evidenceList.setPreferredSize(d2);
    evidenceList.setMaximumSize(d2);
    evidenceList.setActionCommand("evidence=");
    xBox.add(evidenceList);
    
    String qual = getField("qualifier=", qualifierString);
    qualfTextField = new JTextField(qual);      
    qualfTextField.setToolTipText("qualifier column");
    qualfTextField.setPreferredSize(dimension);
    qualfTextField.setMaximumSize(dimension);
    qualfTextField.setActionCommand("qualifier=");
    xBox.add(qualfTextField);
    

    dateField = new DatePanel(getField("date=", qualifierString),
                              dimension.height);

    xBox.add(dateField);
  }


  
  protected boolean isQualifierChanged()
  {
    String old = getField("with=", origQualifierString);
    if(!old.equals(withTextField.getText().trim()))
      return true;
    
    old = getField("db_xref=", origQualifierString);
    if(!old.equals(dbxrefTextField.getText().trim()))
      return true;
    
    old = getField("evidence=", origQualifierString);
    
    if(!(old.equals("") && evidenceList.getSelectedIndex() == -1) )
      if(!old.equals(GoBox.evidenceCodes[2][ evidenceList.getSelectedIndex() ]) &&
         !old.equals(GoBox.evidenceCodes[0][ evidenceList.getSelectedIndex() ]))
        return true;
    
    old = getField("qualifier=", origQualifierString);
    if(!old.equals(qualfTextField.getText()))
      return true;
    
    old = getField("date=", origQualifierString);
    if(!old.equals(dateField.getText()))
      return true;
    
    return false;
  }

  protected void updateQualifier(QualifierVector qv)
  {
    int index = qv.indexOfQualifierWithName(origQualifier.getName());
    Qualifier oldQualifier = qv.getQualifierByName(origQualifier.getName());
    
    final String term = getField("term=", origQualifierString);
    
    StringVector oldValues = oldQualifier.getValues();
    Vector<Integer> values_index = new Vector<Integer>();
    for(int i=0; i<oldValues.size(); i++)
    {
      String oldValue = (String)oldValues.get(i);
      String newTerm = getField("term=", oldValue);
      if(newTerm.equals(term))
        values_index.add(new Integer(i));
    }
  
    if(values_index.size() > 0)
    { 
      String oldValue = (String) oldValues.get(value_index);
      String oldTermId  = getField("term=", oldValue);
      
      if(!term.equals(oldTermId))
      {
        if(values_index.size() == 1)
        value_index = ((Integer)values_index.get(0)).intValue();
        else
        {
          final String with = getField("with=", origQualifierString);
          final String evidence = getField("evidence=", origQualifierString);
          final String dbxref = getField("dbxref=", origQualifierString);
          for(int i=0; i<values_index.size(); i++)
          {
            int ind = ((Integer)values_index.get(i)).intValue();
            value_index = ind;
            String value = (String) oldValues.get(ind);

            if(!with.equals(""))
            {
              String thisWith = getField("with=", value);
              if(thisWith.equals(with))
                break;
            }
          
            if(!dbxref.equals(""))
            {
              String thisDbxref = getField("dbxref=", value);
              if(thisDbxref.equals(dbxref))
                break;
            }
          
            String thisEvidence = getField("evidence=", value);
            if(thisEvidence.equals(evidence))
              break;
          }
        }
      }
    }
    else
      value_index = -99;

    if(value_index > -1)
      oldValues.remove(value_index);
    
    String updatedQualifierString = updateQualifierString();
    oldValues.add(value_index, updatedQualifierString);
    
    origQualifier = new Qualifier(origQualifier.getName(), oldValues);
    qv.remove(index);
    qv.add(index, origQualifier);
  }
  
  private String updateQualifierString()
  {
    String newQualifierString = origQualifierString;
    
    String old = getField("with=", origQualifierString);
    if(!old.equals(withTextField.getText().trim()))
    {
      newQualifierString = changeField("with=", withTextField.getText().trim(), 
                                       newQualifierString);
    }
    
    old = getField("date=", origQualifierString);
    if(!old.equals(dateField.getText()))
    {
      newQualifierString = changeField("date=", dateField.getText().trim(), 
                                       newQualifierString);
    }
    
    old = getField("db_xref=", origQualifierString);
    if(!old.equals(dbxrefTextField.getText().trim()))
    {    
      newQualifierString = changeField("db_xref=", dbxrefTextField.getText().trim(), 
                                       newQualifierString);
    }
    
    old = getField("evidence=", origQualifierString);
    if(evidenceList.getSelectedIndex() > -1 &&
       !old.equals(GoBox.evidenceCodes[2][ evidenceList.getSelectedIndex() ]))
    {
      newQualifierString = changeField("evidence=", GoBox.evidenceCodes[2][ evidenceList.getSelectedIndex() ], 
                                       newQualifierString);
    }
    
    old = getField("qualifier=", origQualifierString);
    if(!old.equals(qualfTextField.getText()))
    {
      newQualifierString = changeField("qualifier=", qualfTextField.getText().trim(), 
                                       newQualifierString);
    }
    
    return newQualifierString;
  }
  
  protected Box getBox()
  {
    return xBox;
  }
}
