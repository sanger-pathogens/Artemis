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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/ortholog/MatchPanel.java,v 1.13 2007-08-02 08:54:10 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder.ortholog;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTable;
import javax.swing.JTextField;

import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.io.DocumentEntry;
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
  private OrthoParalogTable orthologTable;
  private OrthoParalogTable paralogTable;
  private Vector editableComponents;
  private JButton hide_show_ortho;
  private JButton hide_show_sim;
  public static String ORTHOLOG = "orthologous_to";
  public static String PARALOG  = "paralogous_to";
  public static String SIMILARITY = "similarity";
  private static String[] SO_CLUSTER_NAMES =  
            { ORTHOLOG, PARALOG, SIMILARITY };
  
  public MatchPanel(final Feature feature)
  {
    super(new BorderLayout());
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
  
  private Component createMatchQualifiersComponent(final Feature feature)
  {
    editableComponents = new Vector();
    final Qualifier orthoQualifier = matchQualifiers.getQualifierByName(ORTHOLOG);
    final Qualifier paraQualifier  = matchQualifiers.getQualifierByName(PARALOG);
    final Qualifier simQualifier   = matchQualifiers.getQualifierByName(SIMILARITY);
    
    if(databases == null)
    {
      DocumentEntry entry = (DocumentEntry)feature.getEmblFeature().getEntry();
      DatabaseDocument doc = (DatabaseDocument)entry.getDocument();
      databases = (Vector) doc.getOrganismNames();
    }
    
    //
    // ortholog
    Box matchVerticalBox = Box.createVerticalBox();
    JButton addOrthoButton = new JButton("ADD ORTHOLOG");
    addOrthoButton.setOpaque(false);
    addOrthoButton.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      { 
        JExtendedComboBox dbs = new JExtendedComboBox(databases);
        
        JTextField accession = new JTextField(15);
        
        Box yBox = Box.createHorizontalBox();
        yBox.add(dbs);
        yBox.add(accession);

        int select = JOptionPane.showConfirmDialog(null, 
              yBox, "Add Ortholog",
              JOptionPane.OK_CANCEL_OPTION);
        if(select == JOptionPane.CANCEL_OPTION)
          return;
        
        add(ORTHOLOG, ((String)dbs.getSelectedItem())+":"+
                        accession.getText().trim()+"; rank="+
                        orthologTable.getTable().getRowCount(), feature);
      }
    });
    Box xBox = Box.createHorizontalBox();
    xBox.add(addOrthoButton);
    xBox.add(Box.createHorizontalGlue());
    matchVerticalBox.add(xBox);
    
    
    ///
    /// temp
    ///
    if(orthoQualifier != null)
    {
      if(orthoQualifier instanceof QualifierLazyLoading)
        ((QualifierLazyLoading)orthoQualifier).setForceLoad(true);
    //}
      if(hide_show_ortho == null)
        hide_show_ortho = new JButton("-");
      
      DocumentEntry entry = (DocumentEntry)feature.getEmblFeature().getEntry();
      DatabaseDocument doc = (DatabaseDocument)entry.getDocument();
      
      orthologTable = new OrthoParalogTable(doc, orthoQualifier, feature);
      addHideShowButton(orthologTable.getTable(), hide_show_ortho);
      xBox.add(hide_show_ortho);
      editableComponents.add(orthologTable);
      matchVerticalBox.add(orthologTable.getTable().getTableHeader());
      matchVerticalBox.add(orthologTable.getTable());
    }

    
    //
    // paralog
    GeneEditorPanel.addLightSeparator(matchVerticalBox);
    JButton addParaButton = new JButton("ADD PARALOG");
    addParaButton.setOpaque(false);
    addParaButton.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      { 
        JExtendedComboBox dbs = new JExtendedComboBox(databases);
        JTextField accession = new JTextField(15);
        
        Box yBox = Box.createHorizontalBox();
        yBox.add(dbs);
        yBox.add(accession);

        int select = JOptionPane.showConfirmDialog(null, 
              yBox, "Add Paralog",
              JOptionPane.OK_CANCEL_OPTION);
        if(select == JOptionPane.CANCEL_OPTION)
          return;
        
        add(PARALOG, ((String)dbs.getSelectedItem())+":"+
                        accession.getText().trim(), feature);
      }
    });
    xBox = Box.createHorizontalBox();
    xBox.add(addParaButton);
    xBox.add(Box.createHorizontalGlue());
    matchVerticalBox.add(xBox);
    
    if(paraQualifier != null)
    {
      if(paraQualifier instanceof QualifierLazyLoading)
        ((QualifierLazyLoading)paraQualifier).setForceLoad(true);

      if(hide_show_ortho == null)
        hide_show_ortho = new JButton("-");
      
      DocumentEntry entry = (DocumentEntry)feature.getEmblFeature().getEntry();
      DatabaseDocument doc = (DatabaseDocument)entry.getDocument();
      
      paralogTable = new OrthoParalogTable(doc, paraQualifier, feature);
      addHideShowButton(paralogTable.getTable(), hide_show_ortho);
      xBox.add(hide_show_ortho);
      editableComponents.add(paralogTable);
      matchVerticalBox.add(paralogTable.getTable().getTableHeader());
      matchVerticalBox.add(paralogTable.getTable());
    }
    
    
    //
    // similarity
    GeneEditorPanel.addLightSeparator(matchVerticalBox);
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
    xBox = Box.createHorizontalBox();
    xBox.add(addSimButton);
    xBox.add(Box.createHorizontalGlue());
    matchVerticalBox.add(xBox);
    
    if(simQualifier != null)
    {
      if(simQualifier instanceof QualifierLazyLoading)
        ((QualifierLazyLoading)simQualifier).setForceLoad(true);
      
      similarityTable = new SimilarityTable(simQualifier);
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
  
}