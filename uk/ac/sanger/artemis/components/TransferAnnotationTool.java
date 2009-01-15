/* TransferAnnotationTool.java
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
 */

package uk.ac.sanger.artemis.components;

import java.awt.Cursor;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextArea;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeaturePredicate;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.chado.ChadoTransactionManager;
import uk.ac.sanger.artemis.components.genebuilder.GeneEdit;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.PartialSequence;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;

class TransferAnnotationTool extends JFrame
{
  private static final long serialVersionUID = 1L;
  private static String[] NON_TRANSFERABLE_QUALIFIERS =
  {
    "ID",
    "feature_id",
    "Derives_from",
    "feature_relationship_rank",
    "Parent",
    "isObsolete",
    "timelastmodified",
    "orthologous_to",
    "paralogous_to",
    "fasta_file",
    "blastp_file",
    "blastn_file",
    "systematic_id",
    "previous_systematic_id"
  };
  
  public TransferAnnotationTool(final Feature feature, 
  		                          final EntryGroup entryGroup)
  {
    super("Transfer Annotation Tool :: "
        + feature.getIDString());

    JPanel pane = (JPanel) getContentPane();
    pane.setLayout(new GridBagLayout());
    GridBagConstraints c = new GridBagConstraints();
    int nrows = 0;

    c.gridx = 0;
    c.gridy = 0;
    JLabel geneLabel = new JLabel(feature.getIDString() + " Qualifiers");
    geneLabel.setFont(geneLabel.getFont().deriveFont(Font.BOLD));
    pane.add(geneLabel, c);

    c.gridx = 1;
    JLabel label = new JLabel("Gene List");
    label.setFont(label.getFont().deriveFont(Font.BOLD));
    pane.add(label, c);

    c.gridx = 0;
    c.gridy = ++nrows;
    c.anchor = GridBagConstraints.WEST;

    final Vector qualifierCheckBoxes = new Vector();
    final QualifierVector qualifiers = feature.getQualifiers();
    for(int i = 0; i < qualifiers.size(); i++)
    {
      Qualifier qualifier = ((Qualifier) qualifiers.get(i));

      if(isNonTransferable(qualifier.getName()))
          continue;

      JCheckBox checkBox = new JCheckBox(qualifier.getName(), true);
      pane.add(checkBox, c);
      qualifierCheckBoxes.add(checkBox);
      c.gridy = ++nrows;
    }

    
    
    c.gridx = 1;
    c.gridy = 1;
    c.gridheight = nrows;
    c.fill = GridBagConstraints.BOTH;
    final JTextArea geneNameTextArea = new JTextArea("gene1");
    geneNameTextArea.setEditable(true);
    pane.add(geneNameTextArea, c);

    c.gridy = ++nrows;
    c.gridheight = 1;
    c.fill = GridBagConstraints.NONE;
    c.gridx = 0;
    final JButton toggle = new JButton("Toggle Selection");
    toggle.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        for(int i = 0; i < qualifierCheckBoxes.size(); i++)
        {
          JCheckBox cb = (JCheckBox) qualifierCheckBoxes.get(i);
          cb.setSelected(!cb.isSelected());
        }
      }
    });
    pane.add(toggle, c);

    final JCheckBox sameKeyCheckBox = new JCheckBox("Add to feature of same key", true);

    final JButton transfer = new JButton(">>TRANSFER");
    transfer.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        transferAnnotation(qualifierCheckBoxes, geneNameTextArea, feature, 
        		entryGroup, sameKeyCheckBox.isSelected());
      }
    });
    c.gridx = 1;
    pane.add(transfer, c);
    
    c.gridy = ++nrows;
    c.gridx = 1;
    pane.add(sameKeyCheckBox, c);

    final JButton close = new JButton("CLOSE");
    close.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        dispose();
      }
    });
    c.gridy = ++nrows;
    pane.add(close, c);
    pack();
    setVisible(true);
  }

  /**
   * Returns true if this qualifier is non-transferable
   * @param qualifierName
   * @return
   */
  private boolean isNonTransferable(String qualifierName)
  {
    for(int i=0; i<NON_TRANSFERABLE_QUALIFIERS.length; i++)
    {
      if(NON_TRANSFERABLE_QUALIFIERS[i].equals(qualifierName))
        return true;
    }
    return false;
  }
  
  /**
   * Transfer selected qualifiers to the list of features defined
   * by the names in the JTextArea 
   * @param qualifierCheckBoxes - list of qualifier check boxes
   * @param geneNameTextArea - text with a list of feature names to transfer to
   * @param feature - feature to copy from 
   * @param entryGroup
   */
  private void transferAnnotation(final Vector qualifierCheckBoxes, 
  		                            final JTextArea geneNameTextArea,
  		                            final Feature orginatingFeature,
  		                            final EntryGroup entryGroup,
  		                            final boolean sameKey)
  {
  	setCursor(new Cursor(Cursor.WAIT_CURSOR));
    // transfer selected annotation to genes
  	final QualifierVector qualifiers = orginatingFeature.getQualifiers();
  	final QualifierVector qualifiersToTransfer = new QualifierVector();
  	for(int i = 0; i < qualifierCheckBoxes.size(); i++)
    {
  		JCheckBox cb = (JCheckBox) qualifierCheckBoxes.get(i);
  		if(cb.isSelected())
  		{
  			qualifiersToTransfer.addElement(
  					qualifiers.getQualifierByName(cb.getText()).copy());
  		}
    }
  	
  	String geneNames[] = geneNameTextArea.getText().split("\\s");

  	final String key = orginatingFeature.getKey().getKeyString();
  	final FeatureVector features = entryGroup.getAllFeatures();

  	// transfer selected annotation
  	entryGroup.getActionController().startAction();
  	geneNames = transfer(features, qualifiersToTransfer, key, 
  			                 sameKey, GeneUtils.isDatabaseEntry(entryGroup), geneNames);
  	entryGroup.getActionController().endAction();
  	
  	//
  	// Commit changes to genes not in Artemis but in the database
  	//
    DatabaseDocumentEntry db_entry = 
    	(DatabaseDocumentEntry) orginatingFeature.getEntry().getEMBLEntry();
    DatabaseDocument doc = (DatabaseDocument) db_entry.getDocument();
    Vector genesNotFound = null;
    
  	for(int i=0; i<geneNames.length; i++)
  	{
  		DatabaseDocumentEntry newDbEntry = 
  				GeneEdit.makeGeneEntry(null, geneNames[i], doc, null);
  		
  		if(newDbEntry == null)
  		{
  			if(genesNotFound == null)
  				genesNotFound = new Vector();
  			genesNotFound.add(geneNames[i]);
  		  continue;	
  		}
  		
  		char[] c = new char[1];
      PartialSequence ps = new PartialSequence(c, 100, 0, null, null);
      newDbEntry.setPartialSequence(ps);
      Entry entry = null;
      try
      {
        entry = new Entry(newDbEntry);
      }
      catch(Exception e) { e.printStackTrace(); }
      
      SimpleEntryGroup entry_group = new SimpleEntryGroup();
      entry_group.addElement(entry);
      
      ChadoTransactionManager ctm = new ChadoTransactionManager();
      entry_group.addFeatureChangeListener(ctm);
      entry_group.addEntryChangeListener(ctm);
      ctm.setEntryGroup(entry_group);
      
  		transfer(entry.getAllFeatures(), qualifiersToTransfer, key, 
          sameKey, true, geneNames);
  		ChadoTransactionManager.commit(
  				(DatabaseDocument)newDbEntry.getDocument(), false, ctm);
  		
  		entry_group.removeFeatureChangeListener(ctm);
  		entry_group.removeEntryChangeListener(ctm);
      //if(newDbEntry != null)
      //  GeneEdit.showGeneEditor(null, geneNames[i], newDbEntry);
  	}
  	
  	setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
  	
  	if(genesNotFound != null)
  		JOptionPane.showMessageDialog(this, 
  				"Gene(s) Not Found:\n"+genesNotFound.toString(), 
  				"Gene(s) Not Found", JOptionPane.WARNING_MESSAGE);
  }
  
  /**
   * 
   * @param features
   * @param qualifiersToTransfer
   * @param key
   * @param sameKey
   * @param isDatabaseEntry
   * @param geneNames
   * @return
   */
  private String[] transfer(final FeatureVector features,
  		                      final QualifierVector qualifiersToTransfer,
  		                      final String key, 
  		                      final boolean sameKey,
  		                      final boolean isDatabaseEntry,
  		                      String[] geneNames)
	{
		final TransferFeaturePredicate predicate = new TransferFeaturePredicate(
				key, sameKey, isDatabaseEntry, geneNames);

		for (int i = 0; i < features.size(); i++)
		{
			Feature thisFeature = features.elementAt(i);
			if (predicate.testPredicate(thisFeature))
			{
				for (int j = 0; j < qualifiersToTransfer.size(); j++)
				{
					Qualifier qualifier = (Qualifier) qualifiersToTransfer.elementAt(j);
					try
					{
						thisFeature.addQualifierValues(qualifier);
					} catch (Exception e1)
					{
						e1.printStackTrace();
					}
				}
				geneNames = removeArrayElement(geneNames, predicate.getGeneName());
			}
		}
		return geneNames;
	}
  
  /**
   * Remove a string from an array of strings. If the string appears multiple 
   * times in the array this method will delete all occurrences.
   * @param strArr
   * @param str
   * @return
   */
  private String[] removeArrayElement(final String strArr[], final String str)
  {
  	String[] newarray = new String[strArr.length - 1];
  	int count = 0;
  	for(int i=0;i<strArr.length; i++)
  	{
  		if(strArr[i].equals(str))
  			continue;
  		
  	  // not found str return original array
  		if(count>=newarray.length) 
  			return strArr;
  		newarray[count] = strArr[i];
  		count++;
  	}
  	
  	if(count < newarray.length)
  	{
  		String[] tmparray = new String[count];
  		System.arraycopy(newarray, 0, tmparray, 0, count);
  		newarray = tmparray;
  	}
  	
    return newarray;
  }
   
  /**
   * Test if the feature is nominated to have annotation transferred
   * to it.
   */
  class TransferFeaturePredicate implements FeaturePredicate
  {
  	private String geneName;
  	private String key;
  	private boolean sameKey;
  	private boolean isDatabaseEntry;
  	private String[] geneNames;
  	
  	public TransferFeaturePredicate(final String key, 
  			                            final boolean sameKey,
  			                            final boolean isDatabaseEntry,
  			                            final String[] geneNames)
  	{
  		this.key             = key;
  		this.sameKey         = sameKey;
  		this.isDatabaseEntry = isDatabaseEntry;
  		this.geneNames       = geneNames;
  	}
  	
		public boolean testPredicate(Feature targetFeature)
		{
			String targetKey = targetFeature.getKey().getKeyString();
			if(!sameKey || !targetKey.equals(key))
				return false;
			
     	String chadoGeneName = null;
     	if(isDatabaseEntry)
     	{	
     		GFFStreamFeature gffFeature = ((GFFStreamFeature)targetFeature.getEmblFeature());
     		if(gffFeature.getChadoGene() != null)
     		  chadoGeneName = gffFeature.getChadoGene().getGeneUniqueName();
     	}
      	
			String thisFeatureSystematicName = targetFeature.getSystematicName();
			
			for(int i=0;i<geneNames.length;i++)
			{
				if( geneNames[i].equals(thisFeatureSystematicName) ||
						(chadoGeneName != null && geneNames[i].equals(chadoGeneName)) )
				{
					geneName = geneNames[i];
					return true;
				}
			}
			return false;
		}
		
		public String getGeneName()
		{
			return geneName;
		}
  }
}
