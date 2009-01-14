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

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.JLabel;
import javax.swing.JTextField;

import org.gmod.schema.cv.CvTerm;

import uk.ac.sanger.artemis.components.genebuilder.JExtendedComboBox;
import uk.ac.sanger.artemis.components.Splash;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.StringVector;

public class GoBox extends AbstractCvBox
{
  protected static String[][] evidenceCodes = 
  { 
     {"EXP", "IC", "IDA", "IEA", "IEP", "IGC", "IGI", 
      "IMP", "IPI", "ISA", "ISM", "ISO", "ISS", 
      "NAS", "ND", "RCA", "TAS", "NR" }, 
     {"EXP\t:: inferred from experiment",
      "IC \t:: inferred by curator",
      "IDA\t:: inferred from direct assay",
      "IEA\t:: inferred from electronic annotation",
      "IEP\t:: inferred from expression pattern",
      "IGC\t:: inferred from genomic context",
      "IGI\t:: inferred from genetic interaction",
      "IMP\t:: inferred from mutant phenotype",
      "IPI\t:: inferred from physical interaction",
      "ISA\t:: inferred from sequence alignment",
      "ISM\t:: inferred from sequence model",
      "ISO\t:: inferred from sequence ontology",
      "ISS\t:: inferred from sequence or structural similarity",
      "NAS\t:: non-traceable author statement",
      "ND \t:: no biological data available",
      "RCA\t:: inferred from reviewed computational analysis",
      "TAS\t:: traceable author statement",
      "NR \t:: not recorded"},
      { "inferred from experiment",
        "inferred by curator",
        "inferred from direct assay",
        "inferred from electronic annotation",
        "inferred from expression pattern",
        "inferred from genomic context",
        "inferred from genetic interaction",
        "inferred from mutant phenotype",
        "inferred from physical interaction",
        "inferred from sequence alignment",
        "inferred from sequence model",
        "inferred from sequence ontology",
        "inferred from sequence or structural similarity",
        "non-traceable author statement",
        "no biological data available",
        "inferred from reviewed computational analysis",
        "traceable author statement",
        "not recorded"}
  };
  
  private Dimension go_dimension;
  private static Dimension evidenceListDimension;
  
  private Box xBox;
  private int value_index;
  private JTextField withTextField;
  private JTextField dbxrefTextField;
  private JExtendedComboBox evidenceList;
  private JTextField qualfTextField;
  private DatePanel dateField;
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

    JLabel goTermField = new JLabel(goId);
    JLabel goAspect = null;
    
    Font font = goTermField.getFont().deriveFont(Font.BOLD);
    
    if(cvTerm.getCv().getName().indexOf("molecular_function")>-1)
    {
      goAspect = new JLabel(" [F] ");
      goAspect.setForeground(Color.RED);
      goAspect.setFont(font);
    }
    else if(cvTerm.getCv().getName().indexOf("biological_process")>-1)
    {
      goAspect = new JLabel(" [P] ");
      goAspect.setForeground(Color.GREEN);
      goAspect.setFont(font);
    }
    else if(cvTerm.getCv().getName().indexOf("cellular_component")>-1)
    {
      goAspect = new JLabel(" [C] ");
      goAspect.setForeground(Color.BLUE);
      goAspect.setFont(font);
    }
    else
    {
      goAspect = new JLabel(" [?] ");
      goAspect.setForeground(Color.BLACK);
      goAspect.setFont(font);
    }
    
    if(go_dimension == null)
      this.go_dimension = new Dimension(goTermField.getPreferredSize().width+
                                        goAspect.getPreferredSize().width,
                                        goTermField.getPreferredSize().height);
    goTermField.setToolTipText(term);
    xBox.add(goTermField);
    xBox.add(goAspect);
    
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
    evidenceList.setOpaque(false);
    evidenceList.setToolTipText("evidence column");
    evidenceList.setSelectedIndex( getEvidenceIndex(evidence) );
  
    evidenceListDimension = evidenceList.getPreferredSize();
    evidenceListDimension = new Dimension(90,(int)evidenceListDimension.getHeight());
    evidenceList.setPreferredSize(evidenceListDimension);
    evidenceList.setMaximumSize(evidenceListDimension);
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
    
    dateField = new DatePanel( getField("date=", qualifierString), 
                                        dimension.height); 
    
    editable.add(dateField);
    xBox.add(dateField.getDateSpinner());
  }
  
  protected static int getEvidenceIndex(String evidence)
  {
    for(int i=0; i<evidenceCodes[2].length; i++)
    {
      if(evidenceCodes[2][i].equalsIgnoreCase(evidence))
        return i;
    }
    
    // this is mainly to catch RCA
    // - reviewed computational analysis (inferred from missing)
    /*for(int i=0; i<evidenceCodes[2].length; i++)
    {
      if(evidenceCodes[2][i].indexOf(evidence) > -1)
        return i;
    }*/
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
    
    if(evidenceList.getSelectedIndex() > -1 &&
       !old.equalsIgnoreCase(evidenceCodes[2][ evidenceList.getSelectedIndex() ]))
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
    int index = qv.indexOfQualifierWithName(origQualifier.getName());
    Qualifier newQualifier = qv.getQualifierByName(origQualifier.getName());
    
    final String goId = getField("GOid=", origQualifierString);
    
    StringVector values = newQualifier.getValues();
    
    int value_index = -10;
    
    for(int i=0; i<values.size(); i++)
    {
      String newGoId = getField("GOid=", (String)values.get(i));
      if(newGoId.equals(goId))
      {
        value_index = i;
        break;
      }
    }
    
    if(value_index > -1)
      values.remove(value_index);
    
    String updatedQualifierString = updateQualifierString();
    
    Splash.logger4j.debug(origQualifierString);
    Splash.logger4j.debug(updatedQualifierString);
    values.add(value_index, updatedQualifierString);
    
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

  public static Dimension getEvidenceListDimension()
  {
    if(evidenceListDimension == null)
    {
      JExtendedComboBox evidenceList = new JExtendedComboBox(evidenceCodes[1]);
      evidenceListDimension = evidenceList.getPreferredSize();
      evidenceListDimension = new Dimension(80,(int)evidenceListDimension.getHeight());
    }
    
    return evidenceListDimension;
  }
  
  /**
   * Given the string:
   * aspect=F;GOid=GO:0003674;term=molecular_function;evidence=No biological Data available
   * return the string:
   * aspect=F;GOid=GO:0003674;term=molecular_function;evidence=ND
   * @param goText
   * @return
   */
  public static String getEvidenceCodeGoTextFromText(String goText)
  {
    final String oldEvidence = getField("evidence=", goText); 
    String newEvidence = oldEvidence;
    
    for(int i=0; i<evidenceCodes[2].length; i++)
    {
      if(evidenceCodes[2][i].equalsIgnoreCase(oldEvidence.toLowerCase()))
      {
        newEvidence = evidenceCodes[0][i];
        break;
      }
    }
    
    if(!oldEvidence.equals(newEvidence))
      goText = goText.replaceAll(oldEvidence, newEvidence);
    return goText;
  }
  
}