/* OrthologPanel.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/ortholog/MatchPanel.java,v 1.31 2009-05-27 15:00:10 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder.ortholog;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTable;
import javax.swing.JTextField;

import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.sequence.FeatureCvTerm;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.DocumentEntry;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.PartialSequence;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierLazyLoading;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.chado.ChadoTransactionManager;
import uk.ac.sanger.artemis.components.SwingWorker;
import uk.ac.sanger.artemis.components.genebuilder.GeneEdit;
import uk.ac.sanger.artemis.components.genebuilder.GeneEditorPanel;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.components.genebuilder.JExtendedComboBox;

/**
 * For similarity, orthologue, paralogue qualifiers
 */
public class MatchPanel extends JPanel
                      implements FeatureChangeListener
{
  private static final long serialVersionUID = 1L;
  private QualifierVector matchQualifiers;
  private static Vector databases;
  private SimilarityTable similarityTable;
  private OrthoParalogTable orthoparaLogTable;
  private OrthoParalogTable clusterTable;

  private Vector editableComponents;
  private JButton hide_show_ortho;
  private JButton hide_show_cluster;
  private JButton hide_show_sim;
  public static String ORTHOLOG = "orthologous_to";
  public static String PARALOG  = "paralogous_to";
  //public static String CLUSTER  = "cluster";
  public static String SIMILARITY = "similarity";
  
  private static String[] SO_CLUSTER_NAMES = { 
          ORTHOLOG, 
          PARALOG, 
    //    CLUSTER, 
          SIMILARITY };
  private DocumentEntry entry;
  
  // used to test if match panel has contents
  private boolean empty = true;
  private static boolean LOADING = false;
  private Feature feature;
  
  public MatchPanel(final Feature feature, 
                    final DocumentEntry entry)
  {
    super(new BorderLayout());
    this.entry = entry;
    this.feature = feature;
  }
  
  /**
   * Return true if this is a Ortholog qualifier
   * @param qualifier
   * @return
   */
  public static boolean isMatchTag(final Qualifier qualifier)
  {
    return isMatchTag(qualifier.getName());
  }
  
  
  /**
   * Return true if this is a match qualifier
   * @param qualifierName
   * @return
   */
  public static boolean isMatchTag(final String qualifierName)
  {
    for(int i=0; i<SO_CLUSTER_NAMES.length;i++)
      if(qualifierName.equals(SO_CLUSTER_NAMES[i]) || 
         qualifierName.startsWith("/"+SO_CLUSTER_NAMES[i]+"="))
        return true;

    return false;
  }
  
  /**
   * Return true if this is a cluster, ortholog, paralog qualifier
   * @param qualifierName
   * @return
   */
  public static boolean isClusterTag(final String qualifierName)
  {
    for(int i=0; i<SO_CLUSTER_NAMES.length-1;i++)
      if(qualifierName.equals(SO_CLUSTER_NAMES[i]) || 
         qualifierName.startsWith("/"+SO_CLUSTER_NAMES[i]+"="))
        return true;
    return false;
  }
  
  
  /**
   * Create components for ortholog, paralog, cluster and similarity
   * @param feature
   * @return
   */
  private Component createMatchQualifiersComponent(final Feature feature)
  {
    empty = true;
    editableComponents = new Vector();
    final Qualifier orthoQualifier   = matchQualifiers.getQualifierByName(ORTHOLOG);
    final Qualifier paraQualifier    = matchQualifiers.getQualifierByName(PARALOG);
    final Qualifier simQualifier     = matchQualifiers.getQualifierByName(SIMILARITY);
    
    //DocumentEntry entry = (DocumentEntry)feature.getEmblFeature().getEntry();
    final DatabaseDocument doc = (DatabaseDocument)entry.getDocument();
    if(databases == null)
    {
      databases = (Vector) doc.getOrganismNames();
    }
    
    //
    // ortholog / paralog / cluster
    Box matchVerticalBox = Box.createVerticalBox();
    JButton addOrthoButton = new JButton("ADD ORTHOLOG/PARALOG");
    addOrthoButton.setOpaque(false);
    addOrthoButton.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      { 
        addOrthoParalog(feature, doc);
      }
    });
    Box xBox = Box.createHorizontalBox();
    xBox.add(addOrthoButton);
    xBox.add(Box.createHorizontalGlue());
    matchVerticalBox.add(xBox);
    
    
    ///
    if(orthoQualifier != null || paraQualifier != null)
    {
      empty = false;
      if(orthoQualifier != null && orthoQualifier instanceof QualifierLazyLoading)
        ((QualifierLazyLoading)orthoQualifier).setForceLoad(true);
      
      if(paraQualifier != null && paraQualifier instanceof QualifierLazyLoading)
        ((QualifierLazyLoading)paraQualifier).setForceLoad(true);
      
      if(OrthoParalogTable.hasOrthoParlaog(orthoQualifier, paraQualifier, 
                                          (GFFStreamFeature)feature.getEmblFeature()))
      {
        if(hide_show_ortho == null)
          hide_show_ortho = new JButton("-");
        
        orthoparaLogTable = new OrthoParalogTable(doc, orthoQualifier,
            paraQualifier, feature, false);
        
        if(orthoparaLogTable.getTable().getRowCount() > 0)
        {
          addHideShowButton(orthoparaLogTable.getTable(), hide_show_ortho);
          xBox.add(hide_show_ortho);
          editableComponents.add(orthoparaLogTable);

          Box horizontalBox = Box.createHorizontalBox();
          horizontalBox.add(orthoparaLogTable.getTable().getTableHeader());
          horizontalBox.add(Box.createHorizontalGlue());
          matchVerticalBox.add(horizontalBox);
      
      
          horizontalBox = Box.createHorizontalBox();
          horizontalBox.add(orthoparaLogTable.getTable());
          horizontalBox.add(Box.createHorizontalGlue());
          matchVerticalBox.add(horizontalBox);
        }
      }
      
      //
      // clusters
      //
      if(OrthoParalogTable.hasCluster(orthoQualifier, paraQualifier,
          (GFFStreamFeature)feature.getEmblFeature()))
      {
        empty = false;
        if(OrthoParalogTable.hasOrthoParlaog(orthoQualifier, paraQualifier,
            (GFFStreamFeature)feature.getEmblFeature()))
          GeneEditorPanel.addLightSeparator(matchVerticalBox);
       
        if(hide_show_cluster == null)
          hide_show_cluster = new JButton("-");
        
        clusterTable = new OrthoParalogTable(doc, orthoQualifier,
            paraQualifier, feature, true);
        addHideShowButton(clusterTable.getTable(), hide_show_cluster);
        Box horizontalBox = Box.createHorizontalBox();
        horizontalBox.add(Box.createHorizontalGlue());
        horizontalBox.add(hide_show_cluster);
        matchVerticalBox.add(horizontalBox);
        editableComponents.add(clusterTable);

        horizontalBox = Box.createHorizontalBox();
        horizontalBox.add(clusterTable.getTable().getTableHeader());
        horizontalBox.add(Box.createHorizontalGlue());
        matchVerticalBox.add(horizontalBox);

        horizontalBox = Box.createHorizontalBox();
        horizontalBox.add(clusterTable.getTable());
        horizontalBox.add(Box.createHorizontalGlue());
        matchVerticalBox.add(horizontalBox);
      }
    }

    
    //
    // similarity
    GeneEditorPanel.addLightSeparator(matchVerticalBox);
    
    /*
    JButton addSimButton = new JButton("ADD SIMILARITY");
    addSimButton.setOpaque(false);
    addSimButton.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      { 
        JTextField accession = new JTextField(15);
        
        Box yBox = Box.createHorizontalBox();
        yBox.add(accession);

        int select = JOptionPane.showConfirmDialog(null, 
              yBox, "Add Similarity",
              JOptionPane.OK_CANCEL_OPTION);
        if(select == JOptionPane.CANCEL_OPTION)
          return;
        
        add(SIMILARITY, accession.getText().trim(), feature);
      }
    });
    */
    
    
    if(simQualifier != null)
    {
      xBox = Box.createHorizontalBox();
      final JLabel simLabel = new JLabel("Similarity:");
      simLabel.setFont( simLabel.getFont().deriveFont(Font.BOLD ));
      xBox.add(simLabel);
      xBox.add(Box.createHorizontalGlue());
      matchVerticalBox.add(xBox);
      
      empty = false;
      if(simQualifier instanceof QualifierLazyLoading)
        ((QualifierLazyLoading)simQualifier).setForceLoad(true);
      
      similarityTable = new SimilarityTable(simQualifier,doc);
      if(hide_show_sim == null)
        hide_show_sim = new JButton("-");
      addHideShowButton(similarityTable.getTable(), hide_show_sim);
      editableComponents.add(similarityTable);
      
      xBox.add(similarityTable.getInfoLevelButton());
      xBox.add(hide_show_sim);
      matchVerticalBox.add(xBox);
      matchVerticalBox.add(similarityTable.getTable().getTableHeader());
      matchVerticalBox.add(similarityTable.getTable());
    }
    
    return matchVerticalBox;
  }
  
  /**
   * Add an ortholog or paralog to the table
   * @param feature
   */
  private void addOrthoParalog(final Feature feature,
                               final DatabaseDocument doc)
  {
    JExtendedComboBox dbs = new JExtendedComboBox(databases);
    JTextField geneField = new JTextField(15);
    JRadioButton ortho = new JRadioButton(ORTHOLOG, true);
    JRadioButton para  = new JRadioButton(PARALOG, false);
    ButtonGroup group = new ButtonGroup();
    group.add(ortho);
    group.add(para);
    
    Box xBox = Box.createVerticalBox();
    Box yBoxRef = Box.createHorizontalBox();
    yBoxRef.add(dbs);
    yBoxRef.add(geneField);
    yBoxRef.add(Box.createHorizontalGlue());
    xBox.add(yBoxRef);
    
    Box yBoxType = Box.createHorizontalBox();
    yBoxType.add(ortho);
    yBoxType.add(para);
    yBoxType.add(Box.createHorizontalGlue());
    xBox.add(yBoxType);

    boolean found = false;
    JComboBox polypepList = null;
    String uniqueName = null;
    String label = "Add Ortholog/Paralog";
    int select;
    Vector polypeptides = new Vector();
    while(!found)
    {
      select = JOptionPane.showConfirmDialog(null, xBox,
          label, JOptionPane.OK_CANCEL_OPTION);
      if(select == JOptionPane.CANCEL_OPTION)
        return;

      try
      {
        uniqueName = geneField.getText().trim();
        polypeptides = doc.getPolypeptideFeatures(uniqueName);
        
        if(polypeptides == null || polypeptides.size() == 0)
          polypeptides = doc.getPartOfFeatures(uniqueName);

        final Vector polypeptideNames = new Vector();
        
        for(int i=0; i<polypeptides.size(); i++)
        {
          org.gmod.schema.sequence.Feature ppfeature = 
            (org.gmod.schema.sequence.Feature)polypeptides.get(i);
          
          if(!ppfeature.getOrganism().getCommonName().equals(dbs.getSelectedItem()))
            JOptionPane.showMessageDialog(null, "Found in "+
                ppfeature.getOrganism().getCommonName(), "Organism Mismatch", 
                JOptionPane.WARNING_MESSAGE);
          
          polypeptideNames.add(ppfeature.getUniqueName());
        }
        
        polypepList = new JComboBox(polypeptideNames);
        found = true;
      }
      catch(NullPointerException npe)
      {
        found = false;
        label = "Gene : "+uniqueName+"  not found! Try again!";
      }
    }
    
    Box yBoxPeptide = Box.createHorizontalBox();
    yBoxPeptide.add(polypepList);
    yBoxPeptide.add(new JLabel("Add annotation to selected feature"));
    yBoxPeptide.add(Box.createHorizontalGlue());
    xBox.add(yBoxPeptide);
    
    // if there are existing links the optionally transfer annotation link
    final Qualifier orthoQualifier   = matchQualifiers.getQualifierByName(ORTHOLOG);
    final Qualifier paraQualifier    = matchQualifiers.getQualifierByName(PARALOG);
    JCheckBox transferToList = null;
    if(orthoQualifier != null || paraQualifier != null)
    {
      transferToList = 
        new JCheckBox("Add to the existing list of ortho/paralogs", false);
      yBoxPeptide = Box.createHorizontalBox();
      yBoxPeptide.add(transferToList);
      yBoxPeptide.add(Box.createHorizontalGlue());
      xBox.add(yBoxPeptide);
    }
    
    select = JOptionPane.showConfirmDialog(null, 
        xBox, "Add Ortholog/Paralog",
        JOptionPane.OK_CANCEL_OPTION);
    if(select == JOptionPane.CANCEL_OPTION)
      return;
    
    final String type;
    if(ortho.isSelected())
      type = MatchPanel.ORTHOLOG;
    else
      type = MatchPanel.PARALOG;
    
    int rank = 0;
    if(orthoparaLogTable != null)
      rank = orthoparaLogTable.getTable().getRowCount();
    
    // find product
    String product = getProductFromFeatures(polypeptides, 
                                            (String)polypepList.getSelectedItem());

    String qualifierStr = ((String)dbs.getSelectedItem())+":"+
                                uniqueName+" link="+
                                polypepList.getSelectedItem()+" type="+
                                type+"; rank="+rank;
    
    if(product != null)
      qualifierStr = qualifierStr.concat("; product="+product);
    
    if(ortho.isSelected())
      add(ORTHOLOG, qualifierStr, feature);
    else
      add(PARALOG, qualifierStr, feature);
    
    if(transferToList != null && transferToList.isSelected())
    {
      try
      {
        addToExistingList((DatabaseDocumentEntry) feature.getEntry().getEMBLEntry(), 
                       (String)dbs.getSelectedItem()+":"+uniqueName,
                       (String)polypepList.getSelectedItem(),
                       (String)dbs.getSelectedItem());
      }
      catch (Exception e)
      {
        e.printStackTrace();
      }
    }
  }
  
  /**
   * Add a link to any existing ortho/paralog genes.
   * @param db_entry
   * @param orthoOrParaLog   new ortholog/paralog gene (e.g. Pfalciparum:PFA0010c)
   * @param orthoOrParaLogLink this is the name of the feature (usualy polypep)
   * that stores the ortholog/paralog link (e.g. Pfalciparum:PFA0010c:pep)
   * @param organismCommonName
   * @throws ReadOnlyException
   * @throws EntryInformationException
   */
  private void addToExistingList(
                         final DatabaseDocumentEntry db_entry,
                         final String orthoOrParaLog,
                         final String orthoOrParaLogLink,
                         final String organismCommonName) 
          throws ReadOnlyException, EntryInformationException
  {
    DatabaseDocument doc = (DatabaseDocument) db_entry.getDocument();
    String TYPES[] = { ORTHOLOG, PARALOG };
    for(int i=0; i<TYPES.length; i++)
    {
      final Qualifier qualifier = matchQualifiers.getQualifierByName(TYPES[i]);
      if(qualifier != null)
        addToList(doc, qualifier, orthoOrParaLog, 
                  orthoOrParaLogLink, organismCommonName);
    }
  }
  
  /**
   * Add a link to any existing ortho/paralog genes.
   * @param doc
   * @param thisQualifier
   * @param orthoOrParaLog
   * @param orthoOrParaLogLink
   * @throws ReadOnlyException
   * @throws EntryInformationException
   */
  private void addToList(final DatabaseDocument doc, 
                         final Qualifier thisQualifier, 
                         final String orthoOrParaLog,
                         final String orthoOrParaLogLink,
                         final String organismCommonName) 
          throws ReadOnlyException, EntryInformationException
  {
    StringVector values = thisQualifier.getValues();
    
    for (int i = 0; i < values.size(); i++)
    {
      String link = AbstractMatchTable.getField("link", (String) values.get(i));
      link = link.split(" ")[0];
      if(link.equals(orthoOrParaLogLink))
        continue;
      
      DatabaseDocumentEntry newDbEntry = GeneEdit.makeGeneEntry(null,
          link, doc, null);

      if (newDbEntry == null)
        continue;

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
      ctm.setEntryGroup(entry_group);

      
      FeatureVector features = entry_group.getAllFeatures();     
      Feature existingLinkedFeature = features.elementAt(0);
      
      org.gmod.schema.sequence.Feature f =
        ((DatabaseDocument)newDbEntry.getDocument()).getFeatureByUniquename(
              GeneUtils.getUniqueName(existingLinkedFeature.getEmblFeature()));
      
      final String TYPE;
      if(f.getOrganism().getCommonName().equals(organismCommonName))
        TYPE = PARALOG;
      else
        TYPE = ORTHOLOG;
     
      int rank;
      QualifierVector qualifiers = existingLinkedFeature.getQualifiers();
      Qualifier qualifier = qualifiers.getQualifierByName(TYPE);
      
      if(qualifier == null || qualifier.getValues().size() < 1)
        rank = 0;
      else
      {
        ((QualifierLazyLoading)qualifier).setForceLoad(true);
        rank = qualifier.getValues().size();
      }
      String qualifierStr = orthoOrParaLog+" link="+orthoOrParaLogLink+
                            " type="+TYPE+"; rank="+rank;

      int index;
      if(qualifier == null)
      {
        qualifier = new Qualifier(TYPE);
        index = -1;
      }
      else
       index = qualifiers.indexOf(qualifier);
         

      StringVector sv = qualifier.getValues();
      if(sv == null)
        sv = new StringVector();
      sv.add(qualifierStr);
      
      existingLinkedFeature.setQualifier(new Qualifier(TYPE, sv));
      // ADD SQL TO EXISTING ChadoTransactionManager?
      ChadoTransactionManager.commit((DatabaseDocument) newDbEntry
          .getDocument(), false, ctm);

      entry_group.removeFeatureChangeListener(ctm);
    }
  }
  
  /**
   * Return the product for the feature with the given uniquename 
   * from a list of features.
   * @param polypeptides
   * @param selectedPolyPetideName
   * @return
   */
  private String getProductFromFeatures(final Vector features, 
                                        final String uniqueName)
  {
    for(int i=0; i<features.size(); i++)
    {
      org.gmod.schema.sequence.Feature ppfeature = 
        (org.gmod.schema.sequence.Feature)features.get(i);
      if(ppfeature.getUniqueName().equals(uniqueName))
      {
        Collection fc = ppfeature.getFeatureCvTerms();
        Iterator it = fc.iterator();
        while(it.hasNext())
        {
          FeatureCvTerm featureCvTerm = (FeatureCvTerm)it.next();
          CvTerm cvTerm = featureCvTerm.getCvTerm();
          
          if(cvTerm.getCv().getName().equals(
              uk.ac.sanger.artemis.chado.ChadoTransactionManager.PRODUCT_CV))
          {
            return cvTerm.getName();
          }
        }
      }
    }
    return null;
  }
  
  /**
   * Add hide/show button 
   * @param box
   */
  private void addHideShowButton(final JTable table, final JButton hide_show)
  {
    hide_show.setOpaque(false);
    
    // remove any old listeners
    ActionListener l[] = hide_show.getActionListeners();
    if(l != null)
      for(int i=0;i<l.length;i++)
        hide_show.removeActionListener(l[i]);
    
    hide_show.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        if(hide_show.getText().equals("-"))
        {
          hide_show.setText("+");
          table.setVisible(false);
          table.getTableHeader().setVisible(false);
        }
        else
        {
          hide_show.setText("-");
          table.setVisible(true);
          table.getTableHeader().setVisible(true);
        }
      }
    });
  }
  
  /**
   * Update ortho/paralogs for a feature
   * @param feature
   */
  public synchronized void updateFromFeature(final Feature feature)
  {   
    if (matchQualifiers != null)
      feature.removeFeatureChangeListener(MatchPanel.this);
    
    matchQualifiers = new QualifierVector();
    final QualifierVector qualifiers = feature.getQualifiers();
    for (int i = 0; i < qualifiers.size(); ++i)
    {
      final Qualifier qualifier = (Qualifier) qualifiers.elementAt(i);
      if (isMatchTag(qualifier))
        matchQualifiers.addElement(qualifier.copy());
    }

    if(matchQualifiers.size() < 1)
      LOADING = false;
    else
      LOADING = true;
    
    SwingWorker createMatchWorker = new SwingWorker()
    {
      public Object construct()
      {
        removeAll();
        final JLabel loadingData = new JLabel("LOADING DATA...");
        loadingData.setForeground(Color.red);
        add(loadingData, BorderLayout.NORTH);

        feature.addFeatureChangeListener(MatchPanel.this);
        Component matchComponent = createMatchQualifiersComponent(feature);
        LOADING = false;
        
        add(matchComponent);
       
        remove(loadingData);
        repaint();
        revalidate();
        return null;
      }
    };
    createMatchWorker.start();
  }
  
  public void updateFromQualifiers(final QualifierVector qualfiers,
                                   final Feature feature)
  {
    removeAll();
    matchQualifiers = qualfiers;
    add(createMatchQualifiersComponent(feature));
    repaint();
    revalidate();
  }

  /**
   * Add ortholog/paralog
   * @param name  ortholog or paralog
   * @param value
   * @param feature
   */
  private void add(final String name, final String value, final Feature feature)
  {
    final int index;
    
    Qualifier qualifier = matchQualifiers.getQualifierByName(name);
    
    if(qualifier == null)
    {
      qualifier = new Qualifier(name);
      index = -1;
    }
    else
     index = matchQualifiers.indexOf(qualifier);
       
    StringVector sv = qualifier.getValues();
    if(sv == null)
      sv = new StringVector();
    sv.add(value);
    
    qualifier = new Qualifier(name, sv);
    if(index > -1)
    {
      matchQualifiers.remove(index);
      matchQualifiers.add(index, qualifier);
    }
    else
      matchQualifiers.add(qualifier);
    
    removeAll();
    add(createMatchQualifiersComponent(feature));
    repaint();
    revalidate();
  }

  
  /**
   * Get the latest (edited) controlled vocab qualifiers
   * @return
   */
  public QualifierVector getMatchQualifiers()
  {
    if(editableComponents != null)
    {
      for(int i=0; i<editableComponents.size(); i++)
      {
        AbstractMatchTable matchTable = (AbstractMatchTable)editableComponents.get(i);
        //System.out.println("CHECKING MATCHES "+i);
        if(matchTable.isQualifierChanged())
        {
          //System.out.println("UPDATING MATCHES "+i);
          matchTable.updateQualifier(matchQualifiers);
        }
      }
    }
    return matchQualifiers;
  }
  
  public List<String> getGeneNameList()
  {
	return getGeneNameList(false);
  }
  
  /**
   * Get a list of the gene names from the ortholog table or cluster table.
   * @param cluster if true then retrieve names from the cluster table
   * @return
   */
  public List<String> getGeneNameList(boolean cluster)
  {
    final OrthoParalogTable table;
  	if(orthoparaLogTable == null || orthoparaLogTable.getTable().getRowCount()<1 || cluster)
  	{
  	   if(clusterTable == null || clusterTable.getTable().getRowCount()<1)
  	      return null;
  	  table = clusterTable;
  	}
  	else
  	  table = orthoparaLogTable;

  	int columnIndex = table.getColumnIndex(OrthoParalogTable.GENE_COL);
  	List<String> geneNames = new Vector<String>(table.getTable().getRowCount());
  	for(int row=0; row<table.getTable().getRowCount(); row++)
  	{
  		String name[] = 
  			((String)table.getTable().getValueAt(row, columnIndex)).split(":");
  		geneNames.add(name[name.length-1]);
  	}
  	
  	return geneNames;
  }
  
  public static String getDescription()
  {
    return "Ortholog/Paralog/Similarity";  
  }
  
  public void featureChanged(FeatureChangeEvent event)
  {
    updateFromFeature(event.getFeature());
  }

  public boolean isEmpty()
  {
    if(LOADING)
      return false;
    return empty;
  }

  public void setEmpty(boolean empty)
  {
    this.empty = empty;
  }
  
}