/* GoBox.java
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
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JTextField;

import org.gmod.schema.cv.CvTerm;

import uk.ac.sanger.artemis.components.Splash;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.StringVector;

class GoBox extends CvBoxA
{
  protected static String[][] evidenceCodes = 
  { 
     {"IC", "IDA", "IEA", "IEP", "IGC", "IGI", "IMP", "IPI", "ISS", 
      "NAS", "ND", "RCA", "TAS", "NR" }, 
     {"IC \t:: inferred by curator",
      "IDA\t:: inferred from direct assay",
      "IEA\t:: inferred from electronic annotation",
      "IEP\t:: inferred from expression pattern",
      "IGC\t:: inferred from genomic context",
      "IGI\t:: inferred from genetic interaction",
      "IMP\t:: inferred from mutant phenotype",
      "IPI\t:: inferred from physical interaction",
      "ISS\t:: inferred from sequence or structural similarity",
      "NAS\t:: non-traceable author statement",
      "ND \t:: no biological data available",
      "RCA\t:: inferred from reviewed computational analysis",
      "TAS\t:: traceable author statement",
      "NR \t:: not recorded"},
      {"inferred by curator",
        "inferred from direct assay",
        "inferred from electronic annotation",
        "inferred from expression pattern",
        "inferred from genomic context",
        "inferred from genetic interaction",
        "inferred from mutant phenotype",
        "inferred from physical interaction",
        "inferred from sequence or structural similarity",
        "non-traceable author statement",
        "no biological data available",
        "inferred from reviewed computational analysis",
        "traceable author statement",
        "not recorded"}
  };
  
  private Dimension go_dimension;
  private Box xBox;
  private int value_index;
  private JTextField withTextField;
  private JTextField dbxrefTextField;
  private JExtendedComboBox evidenceList;
  private JTextField qualfTextField;
  private JFormattedTextField dateField;
  private String origQualifierString;
  private Qualifier origQualifier;
  
  protected GoBox(final Qualifier qualifier,
                  final String qualifierString,
                  final int value_index,
                  Dimension go_dimension, final Dimension dimension)
  {
    this.origQualifier = qualifier;
    this.origQualifierString = qualifierString;
    this.go_dimension = go_dimension;
    this.value_index  = value_index;
    this.xBox = Box.createHorizontalBox();
    
    Vector editable = new Vector(5);
    
    String goId = getField("GOid=", qualifierString);
    final String term = getField("term=", qualifierString);
    CvTerm cvTerm = DatabaseDocument.getCvTermByCvTermName(term);

    if(cvTerm.getCv().getName().equals("molecular_function"))
      goId = goId+"  aspect=F ";
    else if(cvTerm.getCv().getName().equals("biological_process"))
      goId = goId+"  aspect=P ";
    else if(cvTerm.getCv().getName().equals("cellular_component"))
      goId = goId+"  aspect=C ";
    
    JLabel goTermField = new JLabel(goId);
    if(go_dimension == null)
      this.go_dimension = new Dimension(goTermField.getPreferredSize().width+10,
                                        goTermField.getPreferredSize().height);
    goTermField.setPreferredSize(this.go_dimension);
    goTermField.setToolTipText(term);
    xBox.add(goTermField);
    
    // the WITH column is associated with one or more FeatureCvTermDbXRef
    String with = getField("with=", qualifierString);
    withTextField = new JTextField(with);
    withTextField.setToolTipText("with/from column");
    withTextField.setPreferredSize(dimension);
    withTextField.setMaximumSize(dimension);
    withTextField.setActionCommand("with=");
    editable.add(withTextField);
    xBox.add(withTextField);
 

    // N.B. for /GO the db_xref is a Pub (for primary pubs) 
    //      or FeatureCvTermPub (for others) in /GO
    String dbxref = getField("db_xref=", qualifierString);
    dbxrefTextField = new JTextField(dbxref);
    dbxrefTextField.setToolTipText("dbxref column");
    dbxrefTextField.setPreferredSize(dimension);
    dbxrefTextField.setMaximumSize(dimension);
    dbxrefTextField.setActionCommand("db_xref=");
    editable.add(dbxrefTextField);
    xBox.add(dbxrefTextField);
 
    // feature_cvterm_prop's
    String evidence = getField("evidence=", qualifierString);
    
    evidenceList = new JExtendedComboBox(evidenceCodes[1]);
    evidenceList.setToolTipText("evidence column");
    evidenceList.setSelectedIndex( getEvidenceIndex(evidence) );
  
    Dimension d = evidenceList.getPreferredSize();
    d = new Dimension(80,(int)d.getHeight());
    evidenceList.setPreferredSize(d);
    evidenceList.setMaximumSize(d);
    evidenceList.setActionCommand("evidence=");
    editable.add(evidenceList);
    xBox.add(evidenceList);
    
    String qual = getField("qualifier=", qualifierString);
    qualfTextField = new JTextField(qual);      
    qualfTextField.setToolTipText("qualifier column");
    qualfTextField.setPreferredSize(dimension);
    qualfTextField.setMaximumSize(dimension);
    qualfTextField.setActionCommand("qualifier=");
    editable.add(qualfTextField);
    xBox.add(qualfTextField);
    
    String date = getField("date=", qualifierString);
    Date this_date = getDate(date);
    dateField = new JFormattedTextField(new SimpleDateFormat("yyyyMMdd"))
    {
      protected int getColumnWidth()
      {
        return dateField.getFontMetrics(getFont()).charWidth('0');
      }
    };
    dateField.setInputVerifier(new DateVerifier());
    dateField.setValue(this_date);     
    dateField.setToolTipText("date column");
    dateField.setColumns(8);
    dateField.setMaximumSize(dimension);
    dateField.setActionCommand("date=");
    editable.add(dateField);
    xBox.add(dateField);
  }
  
  private int getEvidenceIndex(String evidence)
  {
    for(int i=0; i<evidenceCodes[2].length; i++)
    {
      if(evidenceCodes[2][i].equals(evidence))
        return i;
    }
    return -1;
  }
  
  protected Dimension getGoDimension()
  {
    return go_dimension;
  }
  
  protected Box getBox()
  {
    return xBox;
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
    if(!old.equals(evidenceCodes[2][ evidenceList.getSelectedIndex() ]))
      return true;
    
    old = getField("qualifier=", origQualifierString);
    if(!old.equals(qualfTextField.getText()))
      return true;
    
    old = getField("date=", origQualifierString);
    if(!old.equals(dateField.getText()))
      return true;
    
    return false;
  }
  
  protected int getValueIndex()
  {
    return value_index;  
  }
  
  protected void updateQualifier(final QualifierVector qv)
  {
    StringVector values = origQualifier.getValues();
    values.remove(value_index);
    String updatedQualifierString = updateQualifierString();
    
    Splash.logger4j.debug(origQualifierString);
    Splash.logger4j.debug(updatedQualifierString);
    values.add(value_index, updatedQualifierString);
    
    int index = qv.indexOfQualifierWithName(origQualifier.getName());
    
    origQualifier = new Qualifier(origQualifier.getName(), values);
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
    
    old = getField("db_xref=", origQualifierString);
    if(!old.equals(dbxrefTextField.getText().trim()))
    {    
      newQualifierString = changeField("db_xref=", dbxrefTextField.getText().trim(), 
                                       newQualifierString);
    }
    
    old = getField("evidence=", origQualifierString);
    if(!old.equals(evidenceCodes[2][ evidenceList.getSelectedIndex() ]))
    {
      newQualifierString = changeField("evidence=", evidenceCodes[2][ evidenceList.getSelectedIndex() ], 
                                       newQualifierString);
    }
    
    old = getField("qualifier=", origQualifierString);
    if(!old.equals(qualfTextField.getText()))
    {
      newQualifierString = changeField("qualifier=", qualfTextField.getText().trim(), 
                                       newQualifierString);
    }
    
    old = getField("date=", origQualifierString);
    if(!old.equals(dateField.getText()))
    {
      newQualifierString = changeField("date=", dateField.getText().trim(), 
                                       newQualifierString);
    }
    
    return newQualifierString;
  }
  
}