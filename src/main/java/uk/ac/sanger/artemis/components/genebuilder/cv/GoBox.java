/* GoBox.java
 *
 * created: 2007
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2019  Genome Research Limited
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

/**
 * 
 * Widget for adding/editing/deleting gene GO terms.
 * 
 */
public class GoBox extends AbstractCvBox
{
  public static String[][] evidenceCodes = 
  { 
     {"EXP", "HDA", "HEP", "HGI", "HMP", "HTP", "IC", "IDA",  
      "IEA", "IEP", "IGC", "IBA", "IBD", "IKR", "IRD", "IGI", 
      "IMP", "IPI", "ISA", "ISM", "ISO", "ISS", 
      "NAS", "ND", "RCA", "TAS", "NR" }, 
     {"EXP\t:: Inferred from Experiment",
      "HDA\t:: Inferred from High Throughput Direct Assay",
      "HEP\t:: Inferred from High Throughput Expression Pattern",
      "HGI\t:: Inferred from High Throughput Genetic Interaction",
      "HMP\t:: Inferred from High Throughput Mutant Phenotype",
      "HTP\t:: Inferred from High Throughput Experiment",
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
    	"Inferred from High Throughput Direct Assay",
    	"Inferred from High Throughput Expression Pattern",
    	"Inferred from High Throughput Genetic Interaction",
    	"Inferred from High Throughput Mutant Phenotype",
    	"Inferred from High Throughput Experiment",
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
  
  /** 
   * Default option for the source combo (index 0). 
   * Unfortunately, we have to make this a space else it won't display
   */
  public static final String SOURCE_DEFAULT_OPTION = " ";
  
  /** 
   * List of source databases obtained from application properties, 
   * e.g. GeneDB.<br/>
   * This list is not generated from the Chado db table as it would need significant 
   * filtering to be usable.
   */
  public static final StringVector sources = generateSourcesList();
  
  private Dimension go_dimension;
  private static Dimension evidenceListDimension;
  private static Dimension sourceListDimension;
  
  private Box xBox;
  
  /** The index of this GO Box row within the "table". */
  private int value_index;
  
  private JTextField withTextField;
  private JTextField dbxrefTextField;
  private JExtendedComboBox evidenceList;
  private JTextField qualfTextField;
  private JExtendedComboBox sourceList;
  private DatePanel dateField;
  private String origQualifierString;
  private Qualifier origQualifier;
  private static Cursor cbusy = new Cursor(Cursor.WAIT_CURSOR);
  private static Cursor cdone = new Cursor(Cursor.DEFAULT_CURSOR);
  private static Cursor chand = new Cursor(Cursor.HAND_CURSOR);
  
  public static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(GoBox.class);
  
  private static String AMIGOURL;
  
  /**
   * Constructor
   * @param qualifier Qualifier
   * @param qualifierString Qualifier
   * @param value_index int
   * @param go_dimension Dimension
   * @param dimension Dimension
   */
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
    CvTerm cvTerm = getGOCvTermFromChado(term);
    
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
    
    // ========================================
    // Source combo - RT ticket 621414
    // ========================================
    
    // Get source database name from the qualifier string
    String source = getField(ASSIGNEDBY_QUALIFIER, qualifierString);
    
    sourceList = new JExtendedComboBox(sources);
    sourceList.setOpaque(false);
    sourceList.setToolTipText("source column");
    
    if ("".equals(source) || getSourceIndex(source) == -1)
    {
      // Default option
      sourceList.setSelectedIndex(0);
    }
    else
    {
      sourceList.setSelectedIndex( getSourceIndex(source) );
    }
    
    sourceListDimension = sourceList.getPreferredSize();
    sourceListDimension = new Dimension(90,(int)sourceListDimension.getHeight());
    sourceList.setPreferredSize(sourceListDimension);
    sourceList.setMaximumSize(sourceListDimension);
    sourceList.setActionCommand(ASSIGNEDBY_QUALIFIER);
    xBox.add(sourceList);
    
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
   * A Non-static version of getGOCvTerm.
   * @param term String
   * @return CvTerm
   */
  public CvTerm getGOCvTermFromChado(String term)
  {
	  return getGOCvTerm(term);
  }
  
  /**
   * Add GO listeners for opening Amigo.
   * @param goTermField JLabel
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
  
  /**
   * Determine the integer index of a given item, in the 
   * evidence list.
   * @param source String
   * @return int
   */
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
  
  /**
   * Determine the integer index of a given item, in the 
   * sources list.
   * @param source String
   * @return int
   */
  protected static int getSourceIndex(final String source) 
  {
	final int numSources = sources.size();
	for(int i = 0; i < numSources; i++)
	{
	  // look for full text in the list of sources
	  if(sources.elementAt(i).equalsIgnoreCase(source))
	    return i;
	}
	
	return -1;
  }

  /**
   * Determine if the user has changed any of the GO fields.
   * @return boolean
   */
  protected boolean isQualifierChanged()
  {
    final String origQual = getOrigQualifierString();
    
    String old = getField("with=", origQual);
    if(!old.equals(getWithTextFieldValue()))
      return true;
    
    old = getField("db_xref=", origQual);
    if(!old.equals(getDbxrefTextFieldValue()))
      return true;
    
    old = getField("evidence=", origQual);
    if(getEvidenceListSelectedIndex() > -1 &&
       !old.equalsIgnoreCase(evidenceCodes[2][ getEvidenceListSelectedIndex() ]))
      return true;
    
    old = getField("qualifier=", origQual);
    if(!old.equals(getQualifierTextFieldValue()))
      return true;
    
    old = getField(ASSIGNEDBY_QUALIFIER, origQual);
    if(getSourceListSelectedIndex() > -1 &&
    	!old.equalsIgnoreCase(getSourceListSelectedValue()) )
      return true;
    
    old = getField("date=", origQual);
    if(!old.equals(getDateFieldValue()))
      return true;
    
    return false;
  }
  
  /**
   * Update the qualifier from the GO form.
   * @param  qv QualifierVector - Vector of controlled vocab qualifiers
   */
  protected void updateQualifier(final QualifierVector qv)
  {
	int index = qv.indexOfQualifierWithName(origQualifier.getName());
	
	// Get the qualifier containing all the cached GO terms strings
    Qualifier oldQualifier = qv.getQualifierByName(origQualifier.getName());
    
    // Get the GO id
    final String origQualifierString = getOrigQualifierString();
    final String goId = getField("GOid=", origQualifierString);
    
    // Get the cached GO term strings from the GO qualifier
    StringVector oldValues = oldQualifier.getValues();
    Vector<Integer> values_index = new Vector<Integer>();
    
    // Find only the cached GO term rows that match our goId
    // and add their indexes to the values_index list.
    for(int i=0; i<oldValues.size(); i++)
    {
      String oldValue = oldValues.get(i);
      String newGoId = getField("GOid=", oldValue);
      if(newGoId.equals(goId))
      {
    	// We have a match, add to our list
        values_index.add(Integer.valueOf(i));
      }
    }
  
    if(values_index.size() > 0)
    { 
      // We are here, because there were some cached GO rows
      // that matched the GO id of this term.
    	
      // Get the GO id at this term's current index from the cached terms list
      String oldValue = oldValues.get(value_index);
      String oldGoId  = getField("GOid=", oldValue);
      
      // Really not sure what the point of this code block is....
      // Potentially remove in the future.
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

    // Remove this GO term from the cached list of GO table rows
    if(value_index > -1)
      oldValues.remove(value_index);
    
    // Get the latest field values and add in
    
    String updatedQualifierString = updateQualifierString();
    logger4j.debug(origQualifierString);
    logger4j.debug(updatedQualifierString);
    oldValues.add(value_index, updatedQualifierString);
    
    origQualifier = new Qualifier(origQualifier.getName(), oldValues);
    qv.remove(index);
    qv.add(index, origQualifier);
  }
  
  /**
   * Update the qualifier string with any new values
   * from the GUI fields.
   * @return String - new qualifier string
   */
  protected String updateQualifierString()
  {
    String newQualifierString = getOrigQualifierString();
    final String origQual = newQualifierString;
    
    String old = getField("with=", origQual);
    if(!old.equals(getWithTextFieldValue()))
      newQualifierString = changeField("with=", getWithTextFieldValue(), 
                                       newQualifierString);
    
    old = getField("db_xref=", origQual);
    if(!old.equals(getDbxrefTextFieldValue()))
      newQualifierString = changeField("db_xref=", getDbxrefTextFieldValue(), 
                                       newQualifierString);
    
    old = getField("evidence=", origQual);
    if(!old.equals(evidenceCodes[2][ getEvidenceListSelectedIndex() ]))
      newQualifierString = changeField("evidence=", evidenceCodes[2][ getEvidenceListSelectedIndex() ], 
                                       newQualifierString);
    
    old = getField("qualifier=", origQual);
    if(!old.equals(getQualifierTextFieldValue()))
      newQualifierString = changeField("qualifier=", getQualifierTextFieldValue(), 
                                       newQualifierString);
    
    old = getField(ASSIGNEDBY_QUALIFIER, origQual);
    String newSource = getSourceListSelectedValue();
    if( !old.equals(newSource) )
        newQualifierString = changeField(ASSIGNEDBY_QUALIFIER, newSource, 
                                       newQualifierString);
    
    old = getField("date=", origQual);
    if(!old.equals(getDateFieldValue()))
      newQualifierString = changeField("date=", getDateFieldValue(), 
                                       newQualifierString);
    return newQualifierString;
  }

  /**
   * Get the dimension associated with the evidence list combo.
   * @return Dimension
   */
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
   * Get the dimension associated with the source list combo.
   * @return Dimension
   */
  protected static Dimension getSourceListDimension()
  {
    if(sourceListDimension == null)
    {
      JExtendedComboBox sourceList = new JExtendedComboBox(sources);
      sourceListDimension = sourceList.getPreferredSize();
      sourceListDimension = new Dimension(80,(int)sourceListDimension.getHeight());
    }
    return sourceListDimension;
  }
  
  /**
   * Given the string:
   * aspect=F;GOid=GO:0003674;term=molecular_function;evidence=No biological Data available
   * return the string:
   * aspect=F;GOid=GO:0003674;term=molecular_function;evidence=ND
   * @param goText
   * @return String
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
  
  /**
   * Initialise the static list of sources.
   * A "source" corresponds to the "assigned by" field
   * in a GAF file.
   * 
   * @return StringVector of sources
   */
  public static StringVector generateSourcesList()
  {
	  // Get list of sources from system properties...
	  StringVector sourcesList = Options.getOptions().getOptionValues("GO_evidence_sources"); 
	  if (sourcesList == null)
	  {
	    	// Just create an empty list
		  sourcesList = new StringVector();
	  }
	  
	  // Add blank default element to the start of the sources list
	  sourcesList.insertElementAt(SOURCE_DEFAULT_OPTION, 0);
	  
	  return sourcesList;
  }
  
  // ========================================================
  // Basic getters and setters
  // ========================================================
  
  /**
   * Get the GO box dimension object
   * @return Dimension
   */
  protected Dimension getGoDimension()
  {
    return go_dimension;
  }
  
  /** 
   * Get the GoBox window widget.
   * @return Box
   */
  protected Box getBox()
  {
    return xBox;
  }
  
  /**
   * Get the value index
   * @return int
   */
  protected int getValueIndex()
  {
    return value_index;  
  }
  
  /**
   * Getter method for the origQualifierString property.
   * @return String
   */
  public String getOrigQualifierString()
  {
	  return origQualifierString;
  }
  
  /**
   * Return the trimmed contents of the "With" text field.
   * @return String
   */
  public String getWithTextFieldValue()
  {
	  return withTextField.getText().trim();
  }
  
  /**
   * Return the trimmed contents of the "Dbxref" text field.
   * @return String
   */
  public String getDbxrefTextFieldValue()
  {
	  return dbxrefTextField.getText().trim();
  }
  
  /**
   * Return the integer index of the currently selected 
   * Evidence combo item.
   * @return int
   */
  public int getEvidenceListSelectedIndex()
  {
	  return evidenceList.getSelectedIndex();
  }
  
  /**
   * Return the trimmed contents of the "Qualifier" text field.
   * @return String
   */
  public String getQualifierTextFieldValue()
  {
	  return qualfTextField.getText().trim();
  }
  
  /**
   * Return the integer index of the currently selected 
   * Source combo item.
   * @return int
   */
  public int getSourceListSelectedIndex()
  {
	  return sourceList.getSelectedIndex();
  }
  
  /**
   * Return the value of the currently selected 
   * Source combo item.
   * Obviously, if you call this method with nothing selected
   * you'll get an ArrayIndexOutOfBoundsException.
   * @return String
   */
  public String getSourceListSelectedValue()
  {
	  return sources.elementAt(getSourceListSelectedIndex()).trim();
  }
  
  /**
   * Return the trimmed contents of the "Date" text field.
   * @return String
   */
  public String getDateFieldValue()
  {
	  return dateField.getText().trim();
  }
}
