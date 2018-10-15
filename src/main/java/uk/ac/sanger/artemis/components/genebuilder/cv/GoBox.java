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
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JTextField;

import org.gmod.schema.cv.CvTerm;

import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.components.genebuilder.JExtendedComboBox;
import uk.ac.sanger.artemis.components.SwingWorker;
import uk.ac.sanger.artemis.editor.BrowserControl;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.StringVector;

public class GoBox extends AbstractCvBox
{
  public static String[][] evidenceCodes = 
  { 
     {"EXP", "IC", "IDA", "IEA", "IEP", "IGC", 
       "IBA", "IBD", "IKR", "IRD", "IGI", 
      "IMP", "IPI", "ISA", "ISM", "ISO", "ISS", 
      "NAS", "ND", "RCA", "TAS", "NR" }, 
     {"EXP\t:: Inferred from Experiment",
      "IC \t:: Inferred by Curator",
      "IDA\t:: Inferred from Direct Assay",
      "IEA\t:: Inferred from Electronic Annotation",
      "IEP\t:: Inferred from Expression Pattern",
      "IGC\t:: Inferred from Genomic Context",
      "IBA\t:: Inferred from Biological aspect of Ancestor",
      "IBD\t:: Inferred from Biological aspect of Descendent",
      "IKR\t:: Inferred from Key Residues",
      "IRD\t:: Inferred from Rapid Divergence",
      "IGI\t:: Inferred from Genetic Interaction",
      "IMP\t:: Inferred from Mutant Phenotype",
      "IPI\t:: Inferred from Physical Interaction",
      "ISA\t:: Inferred from Sequence Alignment",
      "ISM\t:: Inferred from Sequence Model",
      "ISO\t:: Inferred from Sequence Orthology",
      "ISS\t:: Inferred from Sequence or Structural Similarity",
      "NAS\t:: Non-traceable Author Statement",
      "ND \t:: No biological Data available",
      "RCA\t:: inferred from Reviewed Computational Analysis",
      "TAS\t:: Traceable Author Statement",
      "NR \t:: Not Recorded"},
      { "Inferred from Experiment",
        "Inferred by Curator",
        "Inferred from Direct Assay",
        "Inferred from Electronic Annotation",
        "Inferred from Expression Pattern",
        "Inferred from Genomic Context",
        "Inferred from Biological aspect of Ancestor",
        "Inferred from Biological aspect of Descendent",
        "Inferred from Key Residues",
        "Inferred from Rapid Divergence",
        "Inferred from Genetic Interaction",
        "Inferred from Mutant Phenotype",
        "Inferred from Physical Interaction",
        "Inferred from Sequence Alignment",
        "Inferred from Sequence Model",
        "Inferred from Sequence Orthology",
        "Inferred from Sequence or Structural Similarity",
        "Non-traceable Author Statement",
        "No biological Data available",
        "inferred from Reviewed Computational Analysis",
        "Traceable Author Statement",
        "Not Recorded"}
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
  private static Cursor cbusy = new Cursor(Cursor.WAIT_CURSOR);
  private static Cursor cdone = new Cursor(Cursor.DEFAULT_CURSOR);
  private static Cursor chand = new Cursor(Cursor.HAND_CURSOR);
  
  public static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(GoBox.class);
  
  private static String AMIGOURL;
  
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
    
    String goId = getField("GOid=", qualifierString);
    final String term = getField("term=", qualifierString);
    CvTerm cvTerm = getGOCvTerm(term);
    
    final JLabel goTermField = new JLabel(goId);
    addGoLabelLiteners(goTermField);
    
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
    
    evidenceList = new JExtendedComboBox(evidenceCodes[1]);
    evidenceList.setOpaque(false);
    evidenceList.setToolTipText("evidence column");
    evidenceList.setSelectedIndex( getEvidenceIndex(evidence) );
    evidenceList.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        if(((String)evidenceList.getSelectedItem()).startsWith("NR \t::"))
          JOptionPane.showMessageDialog(null, 
              "This evicence code is obsolete:\n"+
              evidenceList.getSelectedItem(), 
              "Obsolete Evidence Code", JOptionPane.WARNING_MESSAGE);
      } 
    });
    evidenceListDimension = evidenceList.getPreferredSize();
    evidenceListDimension = new Dimension(90,(int)evidenceListDimension.getHeight());
    evidenceList.setPreferredSize(evidenceListDimension);
    evidenceList.setMaximumSize(evidenceListDimension);
    evidenceList.setActionCommand("evidence=");
    xBox.add(evidenceList);
    
    String qual = getField("qualifier=", qualifierString);
    qualfTextField = new JTextField(qual);      
    qualfTextField.setToolTipText("qualifier column");
    qualfTextField.setPreferredSize(dimension);
    qualfTextField.setMaximumSize(dimension);
    qualfTextField.setActionCommand("qualifier=");
    xBox.add(qualfTextField);
    
    dateField = new DatePanel( getField("date=", qualifierString), 
                                        dimension.height); 
    xBox.add(dateField);
  }
  
  public static CvTerm getGOCvTerm(String term)
  {
    CvTerm cvTerm = DatabaseDocument.getCvTermByCvTermName(term);
    
    if(cvTerm.getCv().getName().indexOf("molecular_function") < 0 &&
       cvTerm.getCv().getName().indexOf("biological_process") < 0 &&
       cvTerm.getCv().getName().indexOf("cellular_component") < 0)
    {
      CvTerm thisCvTerm = DatabaseDocument.getCvTermByCvAndCvTerm(term,
                                                 "molecular_function");
      
      if(thisCvTerm == null)
        thisCvTerm = DatabaseDocument.getCvTermByCvAndCvTerm(term,
                                            "biological_process");
      
      if(thisCvTerm == null)
        thisCvTerm = DatabaseDocument.getCvTermByCvAndCvTerm(term,
                                            "cellular_component");
      if(thisCvTerm != null)
        cvTerm = thisCvTerm;
    }
    return cvTerm;
  }
  
  /**
   * Add GO listeners for opening Amigo.
   * @param goTermField
   */
  private void addGoLabelLiteners(final JLabel goTermField)
  {
    setAmigoUrl();
    goTermField.addMouseListener(new MouseAdapter()
    {
      public void mouseClicked(final MouseEvent e)
      {
        SwingWorker browserLaunch = new SwingWorker()
        {
          public Object construct()
          {
            if(e.getClickCount() == 1)
            {
              goTermField.setCursor(cbusy);
              BrowserControl.displayURL((AMIGOURL+goTermField.getText()).replaceFirst("GO:GO:", "GO:"));
              goTermField.setCursor(cdone);
            }
            return null;
          }
        };
        browserLaunch.start();
      }
    });
    
    goTermField.addMouseMotionListener(new MouseMotionAdapter()
    {
      public void mouseMoved(MouseEvent e)
      {
        goTermField.setCursor(chand);
      }
    });
  }
  
  /**
   * Set the Amigo URL for hyperlinking GO terms.
   */
  private void setAmigoUrl()
  {
    if(AMIGOURL == null)
    {
      StringVector dbsLinks = Options.getOptions().getOptionValues("hyperlinks");
      for(int i=0; i<dbsLinks.size(); i+=2)
      {
        if(dbsLinks.get(i).equals("GO"))
        {
          AMIGOURL = dbsLinks.get(i+1);
          return;
        }
      }
      AMIGOURL = "http://amigo.geneontology.org/amigo/term/GO:";
    }
  }
  
  protected static int getEvidenceIndex(String evidence)
  {
    for(int i=0; i<evidenceCodes[2].length; i++)
    {
      // look for full text or abbreviation of the code
      if(evidenceCodes[2][i].equalsIgnoreCase(evidence) ||
         evidenceCodes[0][i].equalsIgnoreCase(evidence))
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
  
  /**
   * Update the qualifier from the GO form.
   */
  protected void updateQualifier(final QualifierVector qv)
  {
    int index = qv.indexOfQualifierWithName(origQualifier.getName());
    Qualifier oldQualifier = qv.getQualifierByName(origQualifier.getName());
    
    final String goId = getField("GOid=", origQualifierString);
    
    StringVector oldValues = oldQualifier.getValues();
    Vector<Integer> values_index = new Vector<Integer>();
    for(int i=0; i<oldValues.size(); i++)
    {
      String oldValue = oldValues.get(i);
      String newGoId = getField("GOid=", oldValue);
      if(newGoId.equals(goId))
        values_index.add(new Integer(i));
    }
  
    if(values_index.size() > 0)
    { 
      String oldValue = oldValues.get(value_index);
      String oldGoId  = getField("GOid=", oldValue);
      
      if(!goId.equals(oldGoId))
      {
        if(values_index.size() == 1)
          value_index = values_index.get(0).intValue();
        else
        {
          final String with = getField("with=", origQualifierString);
          final String evidence = getField("evidence=", origQualifierString);
          final String dbxref = getField("dbxref=", origQualifierString);
          for(int i=0; i<values_index.size(); i++)
          {
            int ind = values_index.get(i).intValue();
            value_index = ind;
            String value = oldValues.get(ind);

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
    logger4j.debug(origQualifierString);
    logger4j.debug(updatedQualifierString);
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
      newQualifierString = changeField("with=", withTextField.getText().trim(), 
                                       newQualifierString);
    
    old = getField("db_xref=", origQualifierString);
    if(!old.equals(dbxrefTextField.getText().trim()))
      newQualifierString = changeField("db_xref=", dbxrefTextField.getText().trim(), 
                                       newQualifierString);
    
    old = getField("evidence=", origQualifierString);
    if(!old.equals(evidenceCodes[2][ evidenceList.getSelectedIndex() ]))
      newQualifierString = changeField("evidence=", evidenceCodes[2][ evidenceList.getSelectedIndex() ], 
                                       newQualifierString);
    
    old = getField("qualifier=", origQualifierString);
    if(!old.equals(qualfTextField.getText()))
      newQualifierString = changeField("qualifier=", qualfTextField.getText().trim(), 
                                       newQualifierString);
    
    old = getField("date=", origQualifierString);
    if(!old.equals(dateField.getText()))
      newQualifierString = changeField("date=", dateField.getText().trim(), 
                                       newQualifierString);
    return newQualifierString;
  }

  protected static Dimension getEvidenceListDimension()
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
