/* TransferAnnotationTool.java
 *
 * created: 2009
 *
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
 */

package uk.ac.sanger.artemis.components;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.SwingConstants;
import javax.swing.border.EtchedBorder;
import javax.swing.border.TitledBorder;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeaturePredicate;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.chado.ChadoTransactionManager;
import uk.ac.sanger.artemis.components.filetree.LocalAndRemoteFileManager;
import uk.ac.sanger.artemis.components.genebuilder.GeneEdit;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.components.genebuilder.ortholog.MatchPanel;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.PartialSequence;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.StringVector;

public class TransferAnnotationTool extends JFrame
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
    "Start_range",
    "End_range",
    "timelastmodified",
    "cytoplasm_location",
    "cytoplasmic_polypeptide_region",
    "membrane_structure",
    "non_cytoplasm_location",
    "non_cytoplasmic_polypeptide_region",
    "orthologous_to",
    "paralogous_to",
    "pepstats_file",
    "PlasmoAP_score",
    "polypeptide_domain",
    "fasta_file",
    "blastp_file",
    "blastn_file",
    "systematic_id",
    "transmembrane",
    "transmembrane_polypeptide_region",
    "previous_systematic_id"
  };
  private MatchPanel matchPanel;
  
  private static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(TransferAnnotationTool.class);
  
  protected static Color STEEL_BLUE = new Color(25, 25, 112);
  
  /**
   * Tool for transferring annotation from one feature to other feature(s)
   * @param feature
   * @param entryGroup
   * @param geneNames
   */
  public TransferAnnotationTool(final Feature feature, 
  		                        final EntryGroup entryGroup,
  		                        final MatchPanel matchPanel)
  {
    super("Transfer Annotation Tool :: "+ feature.getIDString());
    this.matchPanel = matchPanel;
    
    List<String> geneNames = null;
    if(matchPanel != null)
      geneNames = matchPanel.getGeneNameList();

    JPanel panel = new JPanel(new FlowLayout(FlowLayout.LEFT));
    JPanel pane = new JPanel(new GridBagLayout());
    JScrollPane jsp = new JScrollPane(panel);
    panel.setBackground(Color.white);
    pane.setBackground(Color.white);
    panel.add(pane);
    
    JPanel framePanel = (JPanel)getContentPane();
    framePanel.add(jsp, BorderLayout.CENTER);
    framePanel.setPreferredSize(new Dimension(600,600));
    
    final Vector<JCheckBox> geneNameCheckBoxes = new Vector<JCheckBox>();
    final Vector<QualifierPanel> qualifierPanels = new Vector<QualifierPanel>();

    addMainPanel(feature, pane, qualifierPanels, 
                 geneNameCheckBoxes, geneNames);
    addBottomButtons(qualifierPanels, geneNameCheckBoxes, 
                     framePanel, entryGroup);
    pack();
    setVisible(true);
  }
  
  /**
   * Construct the panel for setting up the gene list and the
   * list of qualifiers to transfer.
   * @param feature
   * @param pane
   * @param qualifierCheckBoxes
   * @param geneNameCheckBoxes
   * @param geneNames
   */
  private void addMainPanel(final Feature feature, 
                            final JPanel pane, 
                            final Vector<QualifierPanel> qualifierPanels, 
                            final Vector<JCheckBox> geneNameCheckBoxes,
                            final List<String> geneNames)
  {
    GridBagConstraints c = new GridBagConstraints();
    int nrows = 0;

    c.anchor = GridBagConstraints.NORTHWEST;
    c.gridx = 2;
    c.gridy = 0;
    c.ipadx = 50;
    
    JLabel geneLabel = new JLabel("Qualifier(s)");
    geneLabel.setFont(geneLabel.getFont().deriveFont(Font.BOLD));
    pane.add(geneLabel, c);
    
    c.gridy = 0;
    c.gridx = 0;
    JLabel label = new JLabel("Gene List");
    label.setFont(label.getFont().deriveFont(Font.BOLD));
    pane.add(label, c);

    nrows+=3;
    c.gridx = 2;
    c.gridy = nrows;
    c.anchor = GridBagConstraints.WEST;
    
    addQualifierPanel(feature, qualifierPanels, c, nrows, pane);
    nrows+=2;
    
    if(feature.getEmblFeature() instanceof GFFStreamFeature)
    {
      GFFStreamFeature gffFeature = 
        ((GFFStreamFeature) feature.getEmblFeature());
      if (gffFeature.getChadoGene() != null)
      {
        String id = GeneUtils.getUniqueName(gffFeature);
        ChadoCanonicalGene chadoGene = gffFeature.getChadoGene();
        Feature gene = (Feature) chadoGene.getGene().getUserData();
        
        if(!id.equals( GeneUtils.getUniqueName(((GFFStreamFeature)chadoGene.getGene())) ))
         addQualifierPanel(gene, qualifierPanels, c, nrows, pane);
        
        nrows+=2;
       
        String transcriptName = 
          chadoGene.getTranscriptFromName(GeneUtils.getUniqueName(gffFeature));
        
        if(transcriptName != null)
        {
          GFFStreamFeature transcript = 
                 (GFFStreamFeature) chadoGene.getFeatureFromId(transcriptName);
          addQualifierPanel((Feature)transcript.getUserData(), 
              qualifierPanels, c, nrows, pane);
          nrows+=2;
          
          Set<uk.ac.sanger.artemis.io.Feature> children = chadoGene.getChildren(transcript);
          Iterator<uk.ac.sanger.artemis.io.Feature> it = children.iterator();
          
          while(it.hasNext())
          {
            GFFStreamFeature kid = (GFFStreamFeature)it.next();
            if(id.equals( GeneUtils.getUniqueName(((GFFStreamFeature)kid)) ))
              continue;
            addQualifierPanel((Feature)kid.getUserData(), qualifierPanels,
                              c, nrows, pane);
            nrows+=2;
          }
        }
      }
    }
    
    c.gridx = 0;
    c.gridy = 3;
    c.gridheight = nrows;
    c.fill = GridBagConstraints.BOTH;

    final Box geneNameBox = Box.createVerticalBox();
    pane.add(geneNameBox, c);
    
    if(geneNames != null)
    {
      for(int i = 0; i < geneNames.size(); i++)
      {
        JCheckBox cb = new JCheckBox((String) geneNames.get(i),true);
        geneNameBox.add(cb);
        geneNameCheckBoxes.add(cb);
      }
    }

    c.gridy = 1;
    c.gridheight = 1;
    c.fill = GridBagConstraints.NONE;
    c.gridx = 2;
    final JButton toggle = new JButton("Toggle");
    toggle.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        for(int i=0; i<qualifierPanels.size(); i++)
        {
          QualifierPanel qP = qualifierPanels.get(i);
          Enumeration<JCheckBox> enumQualifiers = qP.getQualifierCheckBoxes().keys();
          while(enumQualifiers.hasMoreElements())
          {
            JCheckBox cb = enumQualifiers.nextElement();
            cb.setSelected(!cb.isSelected());
          }
        }
      }
    });
    pane.add(toggle, c);
      
    Box xBox = Box.createHorizontalBox();
    final JButton toggleGeneList = new JButton("Toggle");
    toggleGeneList.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        for(int i = 0; i < geneNameCheckBoxes.size(); i++)
        {
          JCheckBox cb = geneNameCheckBoxes.get(i);
          cb.setSelected(!cb.isSelected());
        }
        geneNameBox.repaint();
      }
    });
    xBox.add(toggleGeneList); 
    
    final JButton addGenes = new JButton("Add");
    addGenes.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        JTextArea geneNameTextArea = new JTextArea();
        geneNameTextArea.setEditable(true);
        JScrollPane jsp = new JScrollPane(geneNameTextArea);
        
        int res = JOptionPane.showConfirmDialog(TransferAnnotationTool.this,
                 jsp, "Paste Feature Names to Add", 
                 JOptionPane.OK_CANCEL_OPTION);
        if(res == JOptionPane.CANCEL_OPTION)
          return;
        
        String geneNames[] = geneNameTextArea.getText().split("\\s");
        for(int i=0;i<geneNames.length; i++)
        {
          if(geneNames[i] == null || geneNames[i].equals(""))
            continue;
           JCheckBox cb = new JCheckBox(geneNames[i],true);
           geneNameBox.add(cb);
           geneNameCheckBoxes.add(cb);
        }
        pane.revalidate();
      }
    });
    xBox.add(addGenes);
    c.gridx = 0;
    pane.add(xBox, c);  
    
    final List<String> clusterList = (matchPanel == null ? null : matchPanel.getGeneNameList(true));
    if(clusterList != null && !geneNames.contains(clusterList.get(0)))
    {
      final JButton importCluster = new JButton("Import Cluster Names");
      importCluster.addActionListener(new ActionListener()
      {
          public void actionPerformed(ActionEvent e)
          {
            for(String n: clusterList)
            {
              if(n == null || n.equals(""))
                continue;
               JCheckBox cb = new JCheckBox(n,true);
               geneNameBox.add(cb);
               geneNameCheckBoxes.add(cb);
            }
            importCluster.setEnabled(false);
            pane.revalidate();
          }
        });
      c.gridy = 2;
      pane.add(importCluster, c);  
    }
  }
  
  /**
   * Add a panel to display a given features qualifiers.
   * @param f
   * @param qualifierPanels
   * @param c
   * @param nrows
   * @param pane
   */
  private void addQualifierPanel(Feature f,
                                 Vector<QualifierPanel> qualifierPanels,
                                 GridBagConstraints c,
                                 int nrows,
                                 JPanel pane)
  {
    QualifierPanel qPanel = new QualifierPanel(f,f.getKey().getKeyString());
    if(qPanel.nrows == 0)
      return;

    c.fill = GridBagConstraints.HORIZONTAL;
    c.anchor = GridBagConstraints.WEST;
    c.weightx = 100;
    qualifierPanels.add(qPanel);
    c.gridy = ++nrows;

    JLabel l = new JLabel(f.getIDString());
    l.setFont(l.getFont().deriveFont(Font.BOLD));
    l.setForeground(STEEL_BLUE);
    pane.add(l, c);

    c.gridy = ++nrows;
    pane.add(qPanel, c);
    c.weightx = 0.d;
  }


  /**
   * Add panel for the transfer and close button.
   * @param qualifierCheckBoxes
   * @param geneNameCheckBoxes
   * @param framePanel
   * @param feature
   * @param entryGroup
   */
  private void addBottomButtons(final Vector<QualifierPanel> qualifierPanels,
                                final Vector<JCheckBox> geneNameCheckBoxes,
                                final JPanel framePanel, 
                                final EntryGroup entryGroup)
  {
    final JCheckBox sameKeyCheckBox = new JCheckBox("Add to feature of same key", true);
    
    final JCheckBox overwriteCheckBox = new JCheckBox("Overwrite", false);
    overwriteCheckBox.setToolTipText("overwrite rather than append values");

    final JCheckBox cvCheckBox = new JCheckBox("Set evidence as ISO and link to source in WITH/FROM", false);
    cvCheckBox.setToolTipText("for GO and Product qualifiers set the evidence as ISO (Inferred from\n"+
                              "Sequence Orthology) and add a link to the source in the WITH/FROM field");

    Box buttonBox = Box.createHorizontalBox();
    final JButton transfer = new JButton(">>TRANSFER");
    transfer.setToolTipText("transfer annotation");
    transfer.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        if(overwriteCheckBox.isSelected())
        {
          int res = JOptionPane.showConfirmDialog(TransferAnnotationTool.this, 
              "Overwrite selected annotation?", "Overwrite", JOptionPane.OK_CANCEL_OPTION);
          if(res == JOptionPane.CANCEL_OPTION)
            return;
        }
        
        final boolean autoHistorySetting = LocalAndRemoteFileManager.isAutomaticHistory();
        StringBuffer summary = new StringBuffer();
        try
        {
          LocalAndRemoteFileManager.setAutomaticHistory(false);
          setCursor(new Cursor(Cursor.WAIT_CURSOR));
          StringBuffer buff = new StringBuffer();
        
          for(int i = 0; i < qualifierPanels.size(); i++)
          {
            QualifierPanel qP = qualifierPanels.get(i);
            int res = transferAnnotation(qP.getQualifierCheckBoxes(), 
              geneNameCheckBoxes, qP.getFeature(), entryGroup, 
              sameKeyCheckBox.isSelected(),
              overwriteCheckBox.isSelected(),
              cvCheckBox.isSelected(),
              buff, summary);
            if(res == -1)
              break;
          }
        
          if(buff.length() > 0)
            logger4j.debug("TRANSFERRED ANNOTATION SUMMARY:\n"+buff.toString());
          setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
          
          if(summary.length()>0)
          {
            final JTextArea list = new JTextArea(summary.toString());
            final JScrollPane jsp = new JScrollPane(list);
            jsp.setPreferredSize(new Dimension(300,200));
            JOptionPane.showMessageDialog(
                TransferAnnotationTool.this, jsp, 
                "Summary of Genes Changed",
                JOptionPane.INFORMATION_MESSAGE);
          }
        }
        finally
        {
          LocalAndRemoteFileManager.setAutomaticHistory(autoHistorySetting);
        }
      }
    });
    Box yBox = Box.createVerticalBox();
    yBox.add(transfer);
    yBox.add(sameKeyCheckBox);
    yBox.add(overwriteCheckBox);
    yBox.add(cvCheckBox);
    buttonBox.add(yBox);
    
    final JButton close = new JButton("CLOSE");
    close.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        dispose();
      }
    });
    yBox = Box.createVerticalBox();
    yBox.add(close);
    yBox.add(Box.createVerticalGlue());
    buttonBox.add(yBox);
    buttonBox.add(Box.createHorizontalGlue());
    framePanel.add(buttonBox, BorderLayout.SOUTH);
  }
  
  
  /**
   * Returns true if this qualifier is non-transferable
   * @param qualifierName
   * @return
   */
  protected static boolean isNonTransferable(String qualifierName)
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
   * by the selected names.
   * @param qualifierCheckBoxes - list of qualifier check boxes
   * @param geneNameTextArea - text with a list of feature names to transfer to
   * @param feature - feature to copy from 
   * @param entryGroup
   * @param sameKey
   * @param overwrite
   */
  protected static int transferAnnotation(
             final Hashtable<JCheckBox, Vector<JCheckBox>> qualifierCheckBoxes, 
  		     final Vector<JCheckBox> geneNameCheckBoxes,
  		     final Feature orginatingFeature,
  		     final EntryGroup entryGroup,
  		     final boolean sameKey,
  		     final boolean overwrite,
  		     final boolean setEvidenceAndWithFrom,
  		     final StringBuffer buff,
  		     final StringBuffer genesUpdated)
  {
    // transfer selected annotation to genes
    final QualifierVector qualifiers = orginatingFeature.getQualifiers();
    final QualifierVector qualifiersToTransfer = new QualifierVector();

    Enumeration<JCheckBox> enumQualifiers = qualifierCheckBoxes.keys();
    while(enumQualifiers.hasMoreElements())
    {
      JCheckBox cb = enumQualifiers.nextElement();
      if (cb.isSelected())
      {
        Vector<JCheckBox> qualifierValuesCheckBox = qualifierCheckBoxes.get(cb);
        final StringVector values = qualifiers.getQualifierByName(cb.getText()).getValues();
        StringVector valuesToTransfer = new StringVector(values);
        
        logger4j.debug("TRANSFER "+cb.getText());
        for(int i=0; i<qualifierValuesCheckBox.size(); i++)
        {
          JCheckBox valuesCb = qualifierValuesCheckBox.get(i);
          if(!valuesCb.isSelected())
          {
            valuesToTransfer.remove(valuesCb.getText());
            logger4j.debug("NOT TRANSFERING "+valuesCb.getText());
          }
        }

        if(valuesToTransfer.size() < 1)
          continue;

        valuesToTransfer = new StringVector( getTransferValues(
            setEvidenceAndWithFrom, orginatingFeature, cb.getText(), valuesToTransfer) );

        qualifiersToTransfer.addElement(new Qualifier(cb.getText(), valuesToTransfer));
      }
    }

    int count = 0;
  	for(int i =0; i<geneNameCheckBoxes.size(); i++)
  	{
  	  if( geneNameCheckBoxes.get(i).isSelected() )
  	    count++;
  	}
  	
  	if(count < 1)
  	{
  	  JOptionPane.showMessageDialog(null, "No genes selected.", 
        "Warning", JOptionPane.WARNING_MESSAGE);
  	  return -1;
  	} 
  	
  	String geneNames[] = new String[count];
  	count = 0;
  	for(int i =0; i<geneNameCheckBoxes.size(); i++)
    {
  	  JCheckBox cb = geneNameCheckBoxes.get(i);
      if( cb.isSelected() )
      {
        geneNames[count] = cb.getText();
        logger4j.debug("TRANSFER ANNOTATION TO "+geneNames[count]);
        count++;
      }
    }

  	final String key = orginatingFeature.getKey().getKeyString();
  	final FeatureVector features = entryGroup.getAllFeatures();

  	// transfer selected annotation
  	entryGroup.getActionController().startAction();
  	geneNames = transfer(features, qualifiersToTransfer, key, sameKey, overwrite,
  			             GeneUtils.isDatabaseEntry(entryGroup), geneNames, genesUpdated);
  	entryGroup.getActionController().endAction();
  	
  	//
  	// Commit changes to genes not in Artemis but in the database
  	//
    Vector<String> genesNotFound = null;
    if (geneNames != null &&
        orginatingFeature.getEntry().getEMBLEntry() instanceof DatabaseDocumentEntry)
    {
      DatabaseDocumentEntry db_entry =
        (DatabaseDocumentEntry) orginatingFeature.getEntry().getEMBLEntry();
      DatabaseDocument doc = (DatabaseDocument) db_entry.getDocument();

      for (int i = 0; i < geneNames.length; i++)
      {
        DatabaseDocumentEntry newDbEntry = GeneEdit.makeGeneEntry(null,
            geneNames[i], doc, null);

        if (newDbEntry == null)
        {
          if (genesNotFound == null)
            genesNotFound = new Vector<String>();
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
        catch (Exception e)
        {
          e.printStackTrace();
        }

        SimpleEntryGroup entry_group = new SimpleEntryGroup();
        entry_group.addElement(entry);

        ChadoTransactionManager ctm = new ChadoTransactionManager();
        entry_group.addFeatureChangeListener(ctm);
        entry_group.addEntryChangeListener(ctm);
        ctm.setEntryGroup(entry_group);

        transfer(entry.getAllFeatures(), qualifiersToTransfer, key, sameKey,
            overwrite, true, geneNames, genesUpdated);
        
        for(int j=0; j<ctm.getTransactionCount(); j++)
          buff.append(ctm.getTransactionAt(j).getLogComment()+"\n");
        ChadoTransactionManager.commit((DatabaseDocument) newDbEntry
            .getDocument(), false, ctm);

        entry_group.removeFeatureChangeListener(ctm);
        entry_group.removeEntryChangeListener(ctm);
        // if(newDbEntry != null)
        // GeneEdit.showGeneEditor(null, geneNames[i], newDbEntry);
      }
    }

  	if(genesNotFound != null)
  		JOptionPane.showMessageDialog(null, 
  				"Gene(s) Not Found:\n"+genesNotFound.toString(), 
  				"Gene(s) Not Found", JOptionPane.WARNING_MESSAGE);
  	return 0;
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
  private static String[] transfer(final FeatureVector features,
                            final QualifierVector qualifiersToTransfer, 
                            final String key,
                            final boolean sameKey,
                            final boolean overwrite,
                            final boolean isDatabaseEntry, 
                            String[] geneNames,
                            final StringBuffer genesUpdated)
  {
    final TransferFeaturePredicate predicate = new TransferFeaturePredicate(
        key, sameKey, isDatabaseEntry, geneNames);

    for (int i = 0; i < features.size(); i++)
    {
      Feature thisFeature = features.elementAt(i);
      if (predicate.testPredicate(thisFeature))
      {
        StringBuffer qualifierBuffer = new StringBuffer();
        for (int j = 0; j < qualifiersToTransfer.size(); j++)
        {
          Qualifier newQualifier = qualifiersToTransfer.elementAt(j);
          String qualifierName = newQualifier.getName();
          try
          {
            if(overwrite)
            {
              thisFeature.setQualifier(newQualifier);
              qualifierBuffer.append("  "+qualifierName+" (overwritten)\n"+
                  parseStringVector(newQualifier.getValues()));
            }
            else
            {
              final StringVector oldValues;
              if (thisFeature.getQualifierByName(newQualifier.getName()) == null)
                oldValues = null;
              else
                oldValues = thisFeature.getQualifierByName(
                    newQualifier.getName()).getValues();

              final Qualifier newQualifierTmp = getQualifierWithoutDuplicateValues(
                  newQualifier, oldValues);
              if (newQualifierTmp == null)
                continue;
              thisFeature.addQualifierValues(newQualifierTmp);
              qualifierBuffer.append("  "+qualifierName+" (added)\n"+
                  parseStringVector(newQualifier.getValues()));
            }
          }
          catch (Exception e1)
          {
            e1.printStackTrace();
          }
        }

        geneNames = removeArrayElement(geneNames, predicate.getGeneName());
        if(qualifierBuffer.length() > 0)
          genesUpdated.append(thisFeature.getSystematicName()+
                              " ("+key+")\n"+qualifierBuffer);
      }
    }
    return geneNames;
  }
  
  /**
   * Get a StringBuffer representation of the values in a StringVector
   * @param v
   * @return
   */
  private static StringBuffer parseStringVector(final StringVector v)
  {
    StringBuffer buff = new StringBuffer();
    for(int i=0; i<v.size(); i++)
      buff.append("    "+v.elementAt(i)+"\n");
    return buff;
  }
  
  /**
   * Return a qualifier copy of the qualifier provided that does not contain
   * any of the values in the StringVector.
   * @param newQualifier
   * @param oldValues
   * @return
   * @throws InvalidRelationException
   */
  protected static Qualifier getQualifierWithoutDuplicateValues( 
      final Qualifier qualifier,
      final StringVector values) throws InvalidRelationException
  {
    final Qualifier newQualifier;
    if (values == null || values.size() < 1)
      newQualifier = qualifier;
    else
    {
      StringVector newValues =  qualifier.getValues();
      StringVector valuesToAdd = new StringVector();

      for (int k = 0; k < newValues.size(); k++)
      {
        if(!values.contains(newValues.get(k)))
        {
          if(qualifier.getName().equals("history"))
          {
        	if(!uk.ac.sanger.artemis.components.genebuilder.cv.HistoryBox.contains(values, (String)newValues.get(k)))
        	  valuesToAdd.add(newValues.get(k));
          }
          else
            valuesToAdd.add(newValues.get(k));
        }
      }

      if(valuesToAdd.size() == 0)
        return null;
      newQualifier = new Qualifier(qualifier.getName(), valuesToAdd);
    }
    return newQualifier;
  }
  
  /**
   * Remove a string from an array of strings. If the string appears multiple 
   * times in the array this method will delete all occurrences.
   * @param strArr
   * @param str
   * @return
   */
  private static String[] removeArrayElement(final String strArr[], final String str)
  {
    if(strArr.length == 1)
    {
      if(strArr[0].equals(str))
        return null;
      return strArr;
    }
    String[] newarray = new String[strArr.length - 1];
    int count = 0;
    for (int i = 0; i < strArr.length; i++)
    {
      if (strArr[i].equals(str))
        continue;

      // not found str return original array
      if (count >= newarray.length)
        return strArr;
      newarray[count] = strArr[i];
      count++;
    }

    if (count < newarray.length)
    {
      String[] tmparray = new String[count];
      System.arraycopy(newarray, 0, tmparray, 0, count);
      newarray = tmparray;
    }
    return newarray;
  }
  
  /**
   * Optionally transfer GO fields with evidence code ISO and link back to
   * the original source in the WITH/FROM column.
   * @param setEvidenceAndWithFrom
   * @param feature
   * @param qName
   * @param values
   * @return
   */
  private static StringVector getTransferValues(final boolean setEvidenceAndWithFrom,
                                                final Feature feature, 
                                                final String qName, 
                                                final StringVector values)
  {
    if(!setEvidenceAndWithFrom)
      return values;

    if(qName.equals("GO") || qName.equals("product"))
    {
      final StringVector tvalues = new StringVector();
      final String gene = getGeneName(feature);
      for (int i = 0; i < values.size(); i++)
      {
        String val =  changeField("evidence=", 
              "Inferred from Sequence Orthology",
              null, values.get(i));

        if(gene != null)
          val =  changeField("with=", 
              "GeneDB:"+gene, "|", val);
        tvalues.add(val);
      }
      return tvalues;
    }
    return values;
  }
  
  private static String getGeneName(Feature feature)
  {
    try
    {
       return
        ((GFFStreamFeature)feature.getEmblFeature()).getChadoGene().getGeneUniqueName();
    }
    catch(Exception e){}
    return null;
  }
  
  /**
   * Replace or add the value of a field in a qualifier string
   * @param fieldName
   * @param newFieldStr
   * @param separator
   * @param qualStr
   * @return
   */
  private static String changeField(final String fieldName, 
                             final String newFieldStr,
                             final String separator,
                             String qualStr)
  {
    int idx1 = qualStr.toLowerCase().indexOf(fieldName.toLowerCase());
    int idx2 = qualStr.indexOf(";", idx1);
    int len  = fieldName.length();
    if(idx2 > idx1 && idx1 > -1)
    {
      if(separator != null)
        qualStr = qualStr.substring(0, idx2) + separator + newFieldStr +
                  qualStr.substring(idx2);
      else
        qualStr = qualStr.substring(0, idx1+len) + newFieldStr +
                  qualStr.substring(idx2);
    }
    else if(idx1 > -1)
    {
      if(separator != null)
        qualStr = qualStr + separator + newFieldStr;
      else
        qualStr = qualStr.substring(0, idx1+len) + newFieldStr;
    }
    else if(!newFieldStr.equals(""))
        qualStr = qualStr + ";" + 
                  fieldName + newFieldStr;
    return qualStr;
  }
}

class QualifierPanel extends JPanel
{
  private static final long serialVersionUID = 1L;
  private Hashtable<JCheckBox, Vector<JCheckBox>> qualifierCheckBoxes = new Hashtable<JCheckBox, Vector<JCheckBox>>();
  private Feature feature;
  protected int nrows = 0;
  
  public QualifierPanel(Feature feature, String title)
  {
    super(new GridBagLayout());
    
    this.feature = feature;
    
    TitledBorder titleBorder = BorderFactory.createTitledBorder(
        BorderFactory.createEtchedBorder(EtchedBorder.LOWERED), title);
    titleBorder.setTitleJustification(TitledBorder.LEFT);
    titleBorder.setTitleColor(TransferAnnotationTool.STEEL_BLUE);
    setBorder(titleBorder);
    
    GridBagConstraints c = new GridBagConstraints();
    c.anchor = GridBagConstraints.WEST;
    c.ipadx = 0;
    final QualifierVector qualifiers = feature.getQualifiers();

    for(int i = 0; i < qualifiers.size(); i++)
    {
      nrows = 
        addQualifierComponents(qualifiers.get(i), 
                              qualifierCheckBoxes, c, nrows);
    }
    
    setMinimumSize(new Dimension(
        titleBorder.getMinimumSize(this).width, 
                   getMinimumSize().height));
  }
  
  /**
   * Add a qualifier to the list of transferable annotation
   * @param qualifier
   * @param qualifierCheckBoxes
   * @param c
   * @param nrows
   * @return
   */
  private int addQualifierComponents(
      final Qualifier qualifier, 
      final Hashtable<JCheckBox, Vector<JCheckBox>> qualifierCheckBoxes,
      final GridBagConstraints c,
      int nrows)
  {
    if(TransferAnnotationTool.isNonTransferable(qualifier.getName()))
      return nrows;
    
    final JCheckBox qualifierNameCheckBox = new JCheckBox(qualifier.getName(), false);
    c.gridx = 1;
    c.fill = GridBagConstraints.NONE;
    c.weightx = 0.d;
    Box qualifierValueBox = Box.createVerticalBox();
    
    JButton detailsShowHide = new JButton("+");
    final Vector<JCheckBox> qualifierValuesCheckBox = setExpanderButton(detailsShowHide,
        qualifier, qualifierValueBox, qualifierNameCheckBox);
    
    qualifierNameCheckBox.addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent e)
      {
        if(qualifierNameCheckBox.isSelected())
        {
          for(int i=0; i<qualifierValuesCheckBox.size(); i++)
          {
            JCheckBox cb = qualifierValuesCheckBox.get(i);
            if(cb.isSelected())
              return;
          }
        }
        
        for(int i=0; i<qualifierValuesCheckBox.size(); i++)
        {
          JCheckBox cb = qualifierValuesCheckBox.get(i);
          cb.setSelected(qualifierNameCheckBox.isSelected());
        }
      }        
    });
    add(detailsShowHide, c);
    
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weightx = 100; 
    c.gridx = 2;

    add(qualifierNameCheckBox, c);
    qualifierCheckBoxes.put(qualifierNameCheckBox, qualifierValuesCheckBox);
    c.gridy = ++nrows;
    add(qualifierValueBox, c);
    c.gridy = ++nrows;
    return nrows;
  }
  
  /**
   * Set up the expander button to display qualifier values.
   * @param butt - expander button
   * @param qualifier - the qualifer that is being displayed
   * @param qualifierValueBox - Box containing the values
   * @param qualifierNameCheckBox - JCheckBox for the given qualifier
   * @param pane
   * @return
   */
  private Vector<JCheckBox> setExpanderButton(final JButton butt,
                                   final Qualifier qualifier, 
                                   final Box qualifierValueBox,
                                   final JCheckBox qualifierNameCheckBox)
  {
    butt.setMargin(new Insets(0, 0, 0, 0));
    butt.setHorizontalAlignment(SwingConstants.RIGHT);
    butt.setHorizontalTextPosition(SwingConstants.RIGHT);
    butt.setBorderPainted(false);
    butt.setFont(butt.getFont().deriveFont(Font.BOLD));
    butt.setForeground(TransferAnnotationTool.STEEL_BLUE);
    
    butt.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        if (butt.getText().equals("+"))
          butt.setText("-");
        else
          butt.setText("+");

        qualifierValueBox.setVisible(butt.getText().equals("-"));
        revalidate();
      }
    });

    // set-up qualifier values list
    qualifierValueBox.setVisible(false);
    final Vector<JCheckBox> qualifierValuesCheckBox = new Vector<JCheckBox>();
    final StringVector values = qualifier.getValues();
    if(values != null)
    {
      for (int i = 0; i < values.size(); i++)
      {
        final JCheckBox cb = new JCheckBox(values.get(i),
          qualifierNameCheckBox.isSelected());
        cb.setFont(cb.getFont().deriveFont(Font.ITALIC));
        qualifierValueBox.add(cb);
        qualifierValuesCheckBox.add(cb);
      }
    }
    return qualifierValuesCheckBox;
  }

  protected Hashtable<JCheckBox, Vector<JCheckBox>> getQualifierCheckBoxes()
  {
    return qualifierCheckBoxes;
  }
  
  protected Feature getFeature()
  {
    return feature;
  }
}

/**
 * Test if the feature is nominated to have annotation transferred
 * to it. For genes in a chado database it looks at the gene name
 * and transcript name.
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
    this.key = key;
    this.sameKey = sameKey;
    this.isDatabaseEntry = isDatabaseEntry;
    this.geneNames = geneNames;
  }

  public boolean testPredicate(Feature targetFeature)
  {
    String targetKey = targetFeature.getKey().getKeyString();
    if (sameKey && !targetKey.equals(key))
      return false;

    Vector<String> chadoNames = null;
    if (isDatabaseEntry)
    {
      GFFStreamFeature gffFeature = 
        ((GFFStreamFeature) targetFeature.getEmblFeature());
      if (gffFeature.getChadoGene() != null)
      {
        chadoNames = new Vector<String>();
        
        ChadoCanonicalGene chadoGene = gffFeature.getChadoGene();
        chadoNames.add(chadoGene.getGeneUniqueName());
        List<uk.ac.sanger.artemis.io.Feature> transcripts = chadoGene.getTranscripts();
        for(int i=0;i<transcripts.size();i++)
        {
          GFFStreamFeature feature = (GFFStreamFeature) transcripts.get(i);
          chadoNames.add(GeneUtils.getUniqueName(feature));
        }
      }
    }

    String thisFeatureSystematicName = targetFeature.getSystematicName();
    for (int i = 0; i < geneNames.length; i++)
    {
      if(geneNames[i].equals(thisFeatureSystematicName) ||
         (chadoNames != null && chadoNames.contains(geneNames[i])))
      {
        geneName = geneNames[i];
        return true;
      }
    }
    return false;
  }

  protected String getGeneName()
  {
    return geneName;
  }
}
