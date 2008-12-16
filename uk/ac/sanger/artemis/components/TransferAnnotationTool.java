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
import javax.swing.JPanel;
import javax.swing.JTextArea;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeaturePredicate;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;

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
    "blastp_file"
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

    final JButton transfer = new JButton(">>TRANSFER");
    transfer.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        transferAnnotation(qualifierCheckBoxes, geneNameTextArea, feature, entryGroup);
      }
    });
    c.gridx = 1;
    pane.add(transfer, c);

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
  		                            final Feature feature,
  		                            final EntryGroup entryGroup)
  {
    // transfer selected annotation to genes
  	final QualifierVector qualifiers = feature.getQualifiers();
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
  	
  	final String geneNames[] = geneNameTextArea.getText().split("\\s");
  	final String key = feature.getKey().getKeyString();
  	final FeatureVector features = entryGroup.getAllFeatures();

  	final FeaturePredicate predicate = new FeaturePredicate()
  	{
			public boolean testPredicate(Feature feature)
			{
	     	String chadoGeneName = null;
	     	if(GeneUtils.isDatabaseEntry(entryGroup))
	     	{	
	     		GFFStreamFeature gffFeature = ((GFFStreamFeature)feature.getEmblFeature());
	     		if(gffFeature.getChadoGene() != null)
	     		  chadoGeneName = gffFeature.getChadoGene().getGeneUniqueName();
	     	}
	      	
				String thisFeatureSystematicName = feature.getSystematicName();
				for(int i=0;i<geneNames.length;i++)
				{
					if(feature.getKey().getKeyString().equals(key))
					{
						if( geneNames[i].equals(thisFeatureSystematicName) ||
								(chadoGeneName != null && geneNames[i].equals(chadoGeneName)) )
						{
							return true;
						}
					}
				}
				return false;
			}
  	};
  	
  	// transfer selected annotation
  	entryGroup.getActionController().startAction();
  	for(int i=0; i<features.size(); i++)
  	{
  	  Feature thisFeature = features.elementAt(i);
  	  if(predicate.testPredicate(thisFeature))
  	  {
  	  	for(int j=0; j<qualifiersToTransfer.size(); j++)
  	  	{
  	  		Qualifier qualifier = (Qualifier) qualifiersToTransfer.elementAt(j);
  	      try
					{
  	      	thisFeature.addQualifierValues(qualifier);
					} 
  	      catch (Exception e1)
					{
						e1.printStackTrace();
					} 
  	  	}
  	  }
  	}
  	entryGroup.getActionController().endAction();
  }
}
