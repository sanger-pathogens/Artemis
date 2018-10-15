/* ReferencesPanel.java
 * This file is part of Artemis
 *
 * Copyright (C) 2009  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/gff/PropertiesPanel.java,v 1.11 2009-08-17 12:50:42 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder;

import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;

import javax.swing.BorderFactory;
import javax.swing.JLabel;
import javax.swing.JPanel;

import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.components.QualifierTextArea;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.StringVector;

public class ReferencesPanel extends JPanel
{
  private static final long serialVersionUID = 1L;
  private QualifierTextArea literatureTextArea;
  private QualifierTextArea dbxrefTextArea;
  
  public ReferencesPanel(final Feature feature)
  {
    super(new FlowLayout(FlowLayout.LEFT));
    updateFromFeature(feature);
  }

  /**
   * Return true if this is a literature or dbxref qualifier
   * @param qualifier
   * @return
   */
  public static boolean isReferenceTag(final Qualifier qualifier)
  {
    if(qualifier.getName().equals("literature") ||
       qualifier.getName().equalsIgnoreCase("Dbxref"))
      return true;
    return false;
  }
  
  public void updateFromFeature(Feature feature)
  {
    super.removeAll();
    GridBagConstraints c = new GridBagConstraints();
    JPanel gridPanel = new JPanel(new GridBagLayout());
    gridPanel.setBackground(Color.white);
    
    //
    // literature & dbxref
    literatureTextArea = new QualifierTextArea();
    literatureTextArea.setBorder(BorderFactory.createLineBorder(Color.gray));
    dbxrefTextArea = new QualifierTextArea();
    dbxrefTextArea.setBorder(BorderFactory.createLineBorder(Color.gray));
    
    literatureTextArea.getDocument().addDocumentListener(
        new TextAreaDocumentListener(literatureTextArea));
    
    dbxrefTextArea.getDocument().addDocumentListener(
        new TextAreaDocumentListener(dbxrefTextArea));
    
    final QualifierVector qualifiers = feature.getQualifiers().copy();  
    final StringBuffer litBuffer = new StringBuffer();
    final StringBuffer dbxrefBuffer = new StringBuffer();

    for(int i = 0 ; i < qualifiers.size(); ++i) 
    {
      Qualifier this_qualifier = (Qualifier)qualifiers.elementAt(i);
      if(this_qualifier.getName().equals("literature"))
        appendToBuffer(this_qualifier.getValues(), litBuffer);
      else if(this_qualifier.getName().equalsIgnoreCase("Dbxref"))
        appendToBuffer(this_qualifier.getValues(), dbxrefBuffer);
    }
    
    c.gridx = 0;
    c.gridy = 0;
    c.ipadx = 5;
    c.ipady = 5;
    c.anchor = GridBagConstraints.NORTHEAST;
    c.fill = GridBagConstraints.NONE;
    JLabel litLab = new JLabel("Literature");
    litLab.setToolTipText("Comma separated list, e.g. PMID:12345, PMID:56789...");
    gridPanel.add(litLab,c);
    c.gridx = 1;
    c.anchor = GridBagConstraints.NORTHWEST;
    gridPanel.add(literatureTextArea,c);
    
    c.gridx = 0;
    c.gridy = 1;
    c.anchor = GridBagConstraints.NORTHEAST;
    JLabel dbLab = new JLabel("Dbxref");
    dbLab.setToolTipText("Comma separated list, e.g. UniProtKB:Q9NFB6, ...");
    gridPanel.add(dbLab,c);
    c.gridx = 1;
    c.anchor = GridBagConstraints.NORTHWEST;
    gridPanel.add(dbxrefTextArea,c);

    add(gridPanel);

    literatureTextArea.setText(litBuffer.toString()+"\n");
    dbxrefTextArea.setText(dbxrefBuffer.toString()+"\n");
  }
  
  private void appendToBuffer(StringVector qualifierStr, StringBuffer buff)
  {
    for (int j = 0; j < qualifierStr.size(); ++j)
    {
      String str = (String) qualifierStr.elementAt(j);
      if(j!=0)
        buff.append(", "+str);
      else
        buff.append(str);
    }
  }
  
  /**
   * Get the latest (edited) literature/dbxref qualifiers
   * @return
   */
  public QualifierVector getQualifiers()
  {
    QualifierVector referenceQualifier = null;
    String literatureTxt = literatureTextArea.getText().trim();
    
    if(!literatureTxt.equals(""))
    {
      referenceQualifier = new QualifierVector();
      String[] lits = getValues(literatureTxt);
      StringVector litValues = new StringVector(lits);
      Qualifier literature = new Qualifier("literature",litValues); 
      referenceQualifier.setQualifier(literature);
    }
    
    String dbxrefTxt = dbxrefTextArea.getText().trim();
    if(!dbxrefTxt.equals(""))
    {
      if(referenceQualifier == null)
        referenceQualifier = new QualifierVector();
      String[] dbxrefs = getValues(dbxrefTxt);
      StringVector dbxrefsValues = new StringVector(dbxrefs);
      Qualifier dbxrefsQualifier = new Qualifier("Dbxref",dbxrefsValues); 
      referenceQualifier.setQualifier(dbxrefsQualifier);
    }
    
    return referenceQualifier;
  }
  
  protected boolean isEmpty()
  {
    QualifierVector qualifiers = getQualifiers();
    if(qualifiers == null)
      return true;
    return false;
  }
  
  private String[] getValues(String txt)
  {
    String delim = "[\t\n\f\r,;]";
    String[] lits = txt.split(delim);
    for(int i=0; i<lits.length; i++)
      lits[i] = lits[i].trim();
    return lits;
  }
  
}