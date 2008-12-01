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
import java.awt.FontMetrics;

import javax.swing.Box;
import javax.swing.JLabel;
import javax.swing.JTextArea;
import javax.swing.JTextField;

import uk.ac.sanger.artemis.components.Splash;
import uk.ac.sanger.artemis.components.genebuilder.JExtendedComboBox;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.StringVector;

class ProductBox extends AbstractCvBox
{
  private Box xBox;
  private int value_index;
  private ProductTextArea termTextField;
  private JTextField withTextField;
  private JTextField dbxrefTextField;
  private JExtendedComboBox evidenceList;
  private String origQualifierString;
  private Qualifier origQualifier;
  
  public ProductBox(final Qualifier qualifier,
                    final String qualifierString,
                    final int value_index,
                    final Dimension dimension,
                    final Dimension go_dimension)
  {
    this.origQualifier = qualifier;
    this.origQualifierString = qualifierString;
    this.value_index  = value_index;
    this.xBox = Box.createHorizontalBox();
    
    JLabel label = new JLabel("product");
    
    final String term = getField("term=", qualifierString);
    termTextField = new ProductTextArea(term, go_dimension, dimension);
    
    label.setPreferredSize(
        new Dimension(termTextField.getLabelWidth(), label.getPreferredSize().height));
    xBox.add(label);
    
    xBox.add(termTextField);
    
    // the WITH column is associated with one or more FeatureCvTermDbXRef
    String with = getField("with=", qualifierString);
    withTextField = new JTextField(with);
    withTextField.setToolTipText("with/from column");
    withTextField.setPreferredSize(dimension);
    withTextField.setMaximumSize(dimension);
    withTextField.setActionCommand("with=");
    
    xBox.add(withTextField);

    // N.B. for /GO the db_xref is a Pub (for primary pubs) 
    //      or FeatureCvTermPub (for others) in /GO
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
  
    Dimension de = evidenceList.getPreferredSize();
    de = new Dimension(80,(int)de.getHeight());
    evidenceList.setPreferredSize(de);
    evidenceList.setMaximumSize(de);
    evidenceList.setActionCommand("evidence=");
    xBox.add(evidenceList);
  }

  /**
   * Check if this qualifier has been changed
   */
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
       !old.equalsIgnoreCase(GoBox.evidenceCodes[2][ evidenceList.getSelectedIndex() ]))
      return true;
    
    return false;
  }

  protected void updateQualifier(QualifierVector qv)
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
    if(!old.equals(GoBox.evidenceCodes[2][ evidenceList.getSelectedIndex() ]))
    {
      newQualifierString = changeField("evidence=", 
                   GoBox.evidenceCodes[2][ evidenceList.getSelectedIndex() ], 
                                       newQualifierString);
    }

    return newQualifierString;
  }
  
  protected Box getBox()
  {
    return xBox;
  }
  
  class ProductTextArea extends JTextArea
  {
    private static final long serialVersionUID = 1L;
    private int labelWidth;
    
    public ProductTextArea(final String text, 
                           final Dimension go_dimension,
                           final Dimension dimension)
    {
      super(text);
      setOpaque(false);
      setEditable(false);
      setLineWrap(true);
      setWrapStyleWord(true);
      FontMetrics fm  = getFontMetrics(getFont());
      int stringWidth = fm.stringWidth(text);
      
      if(go_dimension != null)
        labelWidth = go_dimension.width;
      else
        labelWidth= fm.stringWidth("GO:0001234 [F] ");
      
      int width = labelWidth+(dimension.width*2);

      int rows = Math.round((stringWidth/width)+.5f);
      
      int indexSpace = text.indexOf(' ');
      if(indexSpace > -1 && fm.stringWidth(text.substring(0, indexSpace)) > width)
        rows++;
      
      int rowHeight = getTextRowHeight();
      Dimension d = new Dimension(width, (int) (rowHeight*rows) );
      setPreferredSize(d);
      setMaximumSize(d);
    }
    
    protected int getTextRowHeight()
    {
      return super.getRowHeight();
    }
    
    protected int getLabelWidth()
    {
      return labelWidth;
    }
  }
}
