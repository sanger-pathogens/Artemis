/* ProductBox.java
 *
 * created: 2008
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2008  Genome Research Limited
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
import java.awt.Font;

import javax.swing.Box;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JTextField;

import uk.ac.sanger.artemis.components.Splash;
import uk.ac.sanger.artemis.components.genebuilder.JExtendedComboBox;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.StringVector;

/**
 * Product display box. 
 * 
 * Where there are multiple products associated with a feature assume 
 * the product is the recommended product unless rank=1 which means it 
 * is an alternative product.
 */
class ProductBox extends AbstractCvBox
{
  private Box xBox;
  private int value_index;
  private WrapTextArea termTextField;
  private JTextField withTextField;
  private JTextField dbxrefTextField;
  private JExtendedComboBox evidenceList;
  private JCheckBox recommended = new JCheckBox();
  private String origQualifierStr;
  private Qualifier origQualifier;
  private Box xHeadings = Box.createHorizontalBox();
  
  public ProductBox(final Qualifier qualifier,
                    final String qualifierStr,
                    final int value_index,
                    final Dimension dimension,
                    final Dimension go_dimension)
  {
    this.origQualifier = qualifier;
    this.origQualifierStr = qualifierStr;
    this.value_index  = value_index;
    this.xBox = Box.createHorizontalBox();
    
    String term = getField("term=", qualifierStr);

    // this may not be stored as a CV (product_cv=no?)
    if(term.equals(""))
      term = qualifierStr;
   
    termTextField = new WrapTextArea(term, go_dimension, dimension.width*2);
    
    xBox.add(termTextField);
    
    // the WITH column is associated with one or more FeatureCvTermDbXRef
    String with = getField("with=", qualifierStr);
    withTextField = new JTextField(with);
    withTextField.setToolTipText("with/from column");
    withTextField.setPreferredSize(dimension);
    withTextField.setMaximumSize(dimension);
    withTextField.setActionCommand("with=");
    
    xBox.add(withTextField);

    // N.B. for /GO the db_xref is a Pub (for primary pubs) 
    //      or FeatureCvTermPub (for others) in /GO
    String dbxref = getField("db_xref=", qualifierStr);
    dbxrefTextField = new JTextField(dbxref);
    dbxrefTextField.setToolTipText("dbxref column");
    dbxrefTextField.setPreferredSize(dimension);
    dbxrefTextField.setMaximumSize(dimension);
    dbxrefTextField.setActionCommand("db_xref=");
    
    xBox.add(dbxrefTextField); 
    
    // feature_cvterm_prop's
    String evidence = getField("evidence=", qualifierStr);
    
    evidenceList = new JExtendedComboBox(GoBox.evidenceCodes[1]);
    evidenceList.setOpaque(false);
    evidenceList.setToolTipText("evidence column");
    evidenceList.setSelectedIndex( GoBox.getEvidenceIndex(evidence) );
  
    Dimension de = evidenceList.getPreferredSize();
    de = new Dimension(90,(int)de.getHeight());
    evidenceList.setPreferredSize(de);
    evidenceList.setMaximumSize(de);
    evidenceList.setActionCommand("evidence=");
    xBox.add(evidenceList);
    
    // check box for recommended product
    final String rank = getField("rank=", qualifierStr);
    recommended.setSelected( (rank.length() == 0 || rank.equals("0")) ? true : false );
    evidenceList.setActionCommand("rank=");
    xBox.add(recommended);
    
    if (value_index == 0)
    {
      JLabel lab = new JLabel("Product");
      lab.setFont(lab.getFont().deriveFont(Font.BOLD));
      xHeadings.add(lab);

      xHeadings.add(Box.createRigidArea(new Dimension(termTextField
          .getPreferredSize().width - lab.getPreferredSize().width, 0)));
      JLabel withLabel = new JLabel("WITH/FROM");
      withLabel.setPreferredSize(dimension);
      withLabel.setMaximumSize(dimension);
      xHeadings.add(withLabel);

      JLabel dbxrefLabel = new JLabel("Dbxref");
      dbxrefLabel.setPreferredSize(dimension);
      dbxrefLabel.setMaximumSize(dimension);
      xHeadings.add(dbxrefLabel);

      xHeadings.add(Box.createHorizontalStrut(de.width));
      JLabel recommendedLabel = new JLabel("Recommended");
      recommendedLabel.setPreferredSize(dimension);
      recommendedLabel.setMaximumSize(dimension);
      xHeadings.add(recommendedLabel);

      xHeadings.add(Box.createHorizontalGlue());
    }
  }

  /**
   * Check if this qualifier has been changed
   */
  protected boolean isQualifierChanged()
  {
    String old = getField("with=", origQualifierStr);
    if(!old.equals(withTextField.getText().trim()))
      return true;
    
    old = getField("db_xref=", origQualifierStr);
    if(!old.equals(dbxrefTextField.getText().trim()))
      return true;
    
    old = getField("evidence=", origQualifierStr);
    if(evidenceList.getSelectedIndex() > -1 &&
       !old.equalsIgnoreCase(GoBox.evidenceCodes[2][ evidenceList.getSelectedIndex() ]))
      return true;
    
    // test rank change for recommended/alternative product
    if(recommended.isEnabled())
    {
      old = getField("rank=", origQualifierStr);
      if( !recommended.isSelected() && (old.equals("0") || old.equals("")) )
        return true;
      else if(recommended.isSelected() && !old.equals("0"))
        return true;
    }
    return false;
  }

  protected void updateQualifier(QualifierVector qv)
  {
    Qualifier oldQualifier = qv.getQualifierByName(origQualifier.getName());
    StringVector values = oldQualifier.getValues();

    values.remove(value_index);
    String updatedQualifierString = updateQualifierString();
    
    Splash.logger4j.debug(origQualifierStr);
    Splash.logger4j.debug(updatedQualifierString);
    values.add(value_index, updatedQualifierString);
    
    int index = qv.indexOfQualifierWithName(origQualifier.getName());
    
    origQualifier = new Qualifier(origQualifier.getName(), values);
    qv.remove(index);
    qv.add(index, origQualifier);
  }
  
  private String updateQualifierString()
  {
    String newQualifierString = origQualifierStr;
    
    String old = getField("with=", origQualifierStr);
    if(!old.equals(withTextField.getText().trim()))
    {
      newQualifierString = changeField("with=", withTextField.getText().trim(), 
                                       newQualifierString);
    }
    
    old = getField("db_xref=", origQualifierStr);
    if(!old.equals(dbxrefTextField.getText().trim()))
    {    
      newQualifierString = changeField("db_xref=", dbxrefTextField.getText().trim(), 
                                       newQualifierString);
    }
    
    old = getField("evidence=", origQualifierStr);
    if(evidenceList.getSelectedIndex() > -1 &&
       !old.equals(GoBox.evidenceCodes[2][ evidenceList.getSelectedIndex() ]))
    {
      newQualifierString = changeField("evidence=", 
                   GoBox.evidenceCodes[2][ evidenceList.getSelectedIndex() ], 
                                       newQualifierString);
    }
    
    // test rank change
    old = getField("rank=", origQualifierStr);
    if( !recommended.isSelected() )
    {
      if(old.equals("0") || old.equals(""))
        newQualifierString = changeField("rank=", "1", newQualifierString);
    }
    else
    {
      if(recommended.isEnabled() && !old.equals("0") && !old.equals(""))
        newQualifierString = changeField("rank=", "0", newQualifierString);
    }

    return newQualifierString;
  }
  
  protected Box getBox()
  {
    return xBox;
  }
  
  protected Box getHeadingsBox()
  {
    return xHeadings;
  }

  /**
   * @return the recommended
   */
  protected JCheckBox getRecommended()
  {
    return recommended;
  }
}
