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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/ortholog/MatchPanel.java,v 1.23 2008-01-10 14:17:18 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder.ortholog;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTable;
import javax.swing.JTextField;

import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.io.DocumentEntry;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierLazyLoading;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.components.genebuilder.GeneEditorPanel;
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
  
  public MatchPanel(final Feature feature, final DocumentEntry entry)
  {
    super(new BorderLayout());
    this.entry = entry;
    updateFromFeature(feature);
  }
  
  /**
   * Return true if this is a Ortholog qualifier
   * @param qualifier
   * @return
   */
  public boolean isMatchTag(final Qualifier qualifier)
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
   * Return true if this is a similarity qualifier
   * @param qualifierName
   * @return
   */
  public static boolean isSimilarityTag(final String qualifierName)
  {
    if(qualifierName.equals(SIMILARITY) || 
       qualifierName.startsWith("/"+SIMILARITY+"="))
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
    xBox = Box.createHorizontalBox();
    //xBox.add(addSimButton);
    final JLabel simLabel = new JLabel("Similarity:");
    simLabel.setFont( simLabel.getFont().deriveFont(Font.BOLD ));
    xBox.add(simLabel);
    xBox.add(Box.createHorizontalGlue());
    matchVerticalBox.add(xBox);
    
    if(simQualifier != null)
    {
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
    while(!found)
    {
      select = JOptionPane.showConfirmDialog(null, xBox,
          label, JOptionPane.OK_CANCEL_OPTION);
      if(select == JOptionPane.CANCEL_OPTION)
        return;

      try
      {
        uniqueName = geneField.getText().trim();
        final Vector polypep = doc.getPolypeptideNames(uniqueName);
        polypepList = new JComboBox(polypep);
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
    
    
    final String qualifierStr = ((String)dbs.getSelectedItem())+":"+
                                uniqueName+" link="+
                                polypepList.getSelectedItem()+" type="+
                                type+"; rank="+rank;
    if(ortho.isSelected())
      add(ORTHOLOG, qualifierStr, feature);
    else
      add(PARALOG, qualifierStr, feature);
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
  public void updateFromFeature(final Feature feature)
  {
    removeAll();
    if(matchQualifiers != null)
      feature.removeFeatureChangeListener(this);
    //matchQualifiers = feature.getQualifiers().copy();
    
    matchQualifiers = new QualifierVector();
    final QualifierVector qualifiers = feature.getQualifiers();  
    for(int i = 0 ; i < qualifiers.size(); ++i) 
    {
      Qualifier qualifier = (Qualifier)qualifiers.elementAt(i);
      if(isMatchTag(qualifier))
        matchQualifiers.addElement(qualifier.copy());
    }
   
    feature.addFeatureChangeListener(this);  
    add(createMatchQualifiersComponent(feature));
    repaint();
    revalidate();
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
  public void add(final String name, final String value, final Feature feature)
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
  
  public void featureChanged(FeatureChangeEvent event)
  {
    updateFromFeature(event.getFeature());
  }

  public boolean isEmpty()
  {
    return empty;
  }

  public void setEmpty(boolean empty)
  {
    this.empty = empty;
  }
  
}