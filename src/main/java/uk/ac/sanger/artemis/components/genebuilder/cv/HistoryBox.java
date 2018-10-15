/* HistoryBox
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

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JTextField;

import org.gmod.schema.cv.CvTerm;

import uk.ac.sanger.artemis.chado.ChadoTransactionManager;
import uk.ac.sanger.artemis.components.QualifierTextArea;
import uk.ac.sanger.artemis.components.genebuilder.JExtendedComboBox;
import uk.ac.sanger.artemis.components.genebuilder.TextAreaDocumentListener;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.StringVector;

public class HistoryBox extends AbstractCvBox
{
  private Box xBox;
  private int value_index;
  private JExtendedComboBox termCombo;
  private JTextField curatorNameField;
  private QualifierTextArea qualfTextField;
  private DatePanel dateField;
  private String origQualifierString;
  private Qualifier origQualifier;
  
  /**
   * Display history qualifiers
   * @param qualifier
   * @param qualifierString
   * @param value_index
   * @param dimension
   * @param go_dimension
   */
  public HistoryBox(final Qualifier qualifier,
                    final String qualifierString,
                    final int value_index,
                    final Dimension dimension,
                    final Dimension go_dimension)
  {
    this.origQualifier = qualifier;
    this.origQualifierString = qualifierString;
    this.value_index  = value_index;
    this.xBox = Box.createHorizontalBox();
    
    Box yBox = Box.createVerticalBox();

    Box lineBox = Box.createHorizontalBox();
    final String term = getField("term=", qualifierString);
    
    termCombo = 
      new JExtendedComboBox(getCvTermStrings());
    
    final Dimension d;
    if(go_dimension != null)
    {
      d = new Dimension(go_dimension.width+dimension.width,
          termCombo.getPreferredSize().height);
    }
    else
    {
      FontMetrics fm  = termCombo.getFontMetrics(termCombo.getFont());
      int width= fm.stringWidth("GO:0001234 [F] ");
      d = new Dimension(width+dimension.width,
          termCombo.getPreferredSize().height);
    }
    termCombo.setPreferredSize(d);
    termCombo.setMaximumSize(d);
    
    termCombo.setSelectedItem(term);
    lineBox.add(termCombo);

    // feature_cvterm_prop's
    Dimension dimension2 = new Dimension(d.width*2, d.height);
    String qual = getField("curatorName=", qualifierString);
    curatorNameField = new JTextField(qual);
    curatorNameField.setToolTipText("qualifier column");
    curatorNameField.setPreferredSize(dimension2);
    curatorNameField.setMaximumSize(dimension2);
    curatorNameField.setActionCommand("curatorName=");
    curatorNameField.setCaretPosition(0);
    lineBox.add(curatorNameField);

    dateField = new DatePanel(getField("date=", qualifierString),
                              dimension.height);

    lineBox.add(dateField);
    lineBox.add(Box.createHorizontalGlue());
    yBox.add(lineBox);
    
    lineBox = Box.createHorizontalBox();
    lineBox.add(Box.createHorizontalStrut(5));
    Dimension dimension4 = new Dimension(
        termCombo.getPreferredSize().width+
        curatorNameField.getPreferredSize().width+
        dateField.getPreferredSize().width-5, dimension.height*20);
    qual = getFieldIgnoreSeparator("qualifier", qualifierString);
    qualfTextField = new QualifierTextArea();
    qualfTextField.setUseHyperlinks(false);
    qualfTextField.setBorder(BorderFactory.createLineBorder(Color.gray));

    qualfTextField.setText(qual);      
    qualfTextField.setToolTipText("qualifier column");
    qualfTextField.setPreferredSize(dimension4);
    qualfTextField.setMaximumSize(dimension4);
    qualfTextField.getDocument().addDocumentListener(
        new TextAreaDocumentListener(qualfTextField));
    qualfTextField.setCaretPosition(0);
    lineBox.add(qualfTextField);
    lineBox.add(Box.createHorizontalGlue());
    yBox.add(lineBox);
    
    xBox.add(yBox);
  }

  protected boolean isQualifierChanged()
  {
    String old = getField("term=", origQualifierString);
    if(!old.equals(termCombo.getSelectedItem()))
      return true;
    
    old = getFieldIgnoreSeparator("qualifier", origQualifierString);
    if(!old.equals(qualfTextField.getText()))
      return true;
    
    old = getField("date=", origQualifierString);
    if(!old.equals(dateField.getText()))
      return true;
    
    old = getField("curatorName=", origQualifierString);
    if(!old.equals(curatorNameField.getText()))
      return true;
    
    return false;
  }
  
  public static boolean contains(StringVector oldQualValues, String newQualString)
  {
	for(int i=0; i<oldQualValues.size(); i++)
	{
	  String oldQualString = (String) oldQualValues.get(i);
      String oldStr = getField("term=", oldQualString);
      String newStr = getField("term=", newQualString);
      if(!oldStr.equals(newStr))
        continue;
    
      oldStr = getFieldIgnoreSeparator("qualifier", oldQualString);
      newStr = getFieldIgnoreSeparator("qualifier", newQualString);
      if(!oldStr.equals(newStr))
        continue;
    
      oldStr = getField("date=", oldQualString);
      newStr = getField("date=", newQualString);
      if(!oldStr.equals(newStr))
        continue;
    
      oldStr = getField("curatorName=", oldQualString);
      newStr = getField("curatorName=", newQualString);
      if(oldStr.equals(newStr))
        return true;
	}
	
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
    String old = getField("term=", origQualifierString);
    if(!old.equals(dateField.getText()))
    {
      newQualifierString = changeField("term=", (String) termCombo.getSelectedItem(), 
                                       newQualifierString);
    }
    
    old = getField("date=", origQualifierString);
    if(!old.equals(dateField.getText()))
    {
      newQualifierString = changeField("date=", dateField.getText().trim(), 
                                       newQualifierString);
    }
    
    old = getField("curatorName=", origQualifierString);
    if(!old.equals(curatorNameField.getText().trim()))
    {
      newQualifierString = changeField("curatorName=", curatorNameField.getText().trim(), 
                                       newQualifierString);
    }
    
    old = getFieldIgnoreSeparator("qualifier", origQualifierString);
    if(!old.equals(qualfTextField.getText()))
    {
      newQualifierString = newQualifierString.replaceAll(old, "");
      newQualifierString = newQualifierString.replaceAll("qualifier=[^;]+", "");
      newQualifierString = newQualifierString.replaceAll("qualifier=", "");
      newQualifierString = newQualifierString.replaceAll("[;]{2,}", ";");
      newQualifierString += ";qualifier="+qualfTextField.getText().trim();
    }

    return newQualifierString;
  }
  
  protected Box getBox()
  {
    return xBox;
  }
  
  protected static Vector<CvTerm> getCvTerms()
  {
    return DatabaseDocument.getCvterms("", ChadoTransactionManager.HISTORY_CV, false);
  }
  
  public static Vector<String> getCvTermStrings()
  {
    Vector<CvTerm> cvTerms = getCvTerms();
    Vector<String> cvTermStr = new Vector<String>();
    for(int i=0; i<cvTerms.size(); i++)
    {
      CvTerm cvTerm = cvTerms.get(i);
      cvTermStr.add(cvTerm.getName());
    }
    return cvTermStr;
  }
  
  public static CvTerm getDefaultTerm()
  {
    CvTerm cvterm =
      (CvTerm) DatabaseDocument.getCvterms("", ChadoTransactionManager.HISTORY_CV, false).get(0);
    return cvterm;
  }
}
