/* GeneComponentTree.java
 *
 * created: 2006
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2006  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or(at your option) any later version.
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/GeneComponentTree.java,v 1.23 2009-06-12 13:50:35 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder;

import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.DatabaseInferredFeature;
import uk.ac.sanger.artemis.io.Feature;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.FeatureVector;

import javax.swing.*;
import javax.swing.event.*;

import javax.swing.tree.*;
import java.util.List;
import java.util.Enumeration;
import java.util.Vector;

/**
 * Tree to display a gene hierarchy.
 */
public class GeneComponentTree extends JTree
{
  
  /** */
  private static final long serialVersionUID = 1L;
  private ChadoCanonicalGene chado_gene;
  private GeneBuilderFrame gene_builder;
  private Selection selection;
  private GeneTreeSelectionListener selection_listener;
  
  public GeneComponentTree(final ChadoCanonicalGene chado_gene,
                           final GeneBuilderFrame gene_builder,
                           final Selection selection)
  {
    this.chado_gene   = chado_gene;
    this.gene_builder = gene_builder;
    this.selection    = selection;
    
    final Feature gene = (Feature)chado_gene.getGene();
    final String gene_id;
    try
    {
      gene_id = (String)gene.getQualifierByName("ID").getValues().get(0);
      DefaultMutableTreeNode top =
           new DefaultMutableTreeNode(gene_id);

      createNodes(top, chado_gene);
      DefaultTreeModel model = new DefaultTreeModel(top);
      setModel(model);
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
    
    //Listen for when a file is selected
    selection_listener = new GeneTreeSelectionListener();
    addTreeSelectionListener(selection_listener);
    setSelection(selection);
  }
  
  /**
   * Select the tree nodes based on the features selected.
   * @param selection
   */
  protected void setSelection(Selection selection)
  {
    removeTreeSelectionListener(selection_listener);
    FeatureVector features = selection.getAllFeatures();
    TreePath path[] = new TreePath[features.size()]; 
    for(int i=0; i<features.size(); i++)
    {
      uk.ac.sanger.artemis.Feature feature = 
        (uk.ac.sanger.artemis.Feature)features.elementAt(i);

      try
      {
        String id = 
          (String)feature.getQualifierByName("ID").getValues().get(0);
        DefaultMutableTreeNode node = getNodeFromName(id);
        if(node == null)
          return;
        
        path[i] = new TreePath(node.getPath());
      }
      catch(InvalidRelationException e)
      {
        e.printStackTrace();
      }
    }
    setSelectionPaths(path);
    addTreeSelectionListener(selection_listener);
  }
  
  /**
   * Build the node hierarchy for this gene, i.e. transcript, exons, 
   * CDS, polypeptide.
   * @param gene_node
   * @param chado_gene
   * @throws InvalidRelationException
   */
  private void createNodes(final DefaultMutableTreeNode gene_node,
                           final ChadoCanonicalGene chado_gene) 
          throws InvalidRelationException 
  {
    List transcripts = chado_gene.getTranscripts();
    
    Feature transcript;
    Feature exon;
    Object  protein;
    String transcript_id;
    String exon_id;
    String protein_id;
    DefaultMutableTreeNode transcript_node;
    DefaultMutableTreeNode exon_node;
    DefaultMutableTreeNode protein_node;
    
    for(int i=0; i<transcripts.size(); i++)
    {
      transcript = (Feature)transcripts.get(i);
      transcript_id = 
        (String)transcript.getQualifierByName("ID").getValues().get(0);
      transcript_node = new DefaultMutableTreeNode(transcript_id);
      gene_node.add(transcript_node);
      
      // exons
      List exons = chado_gene.getSplicedFeaturesOfTranscript(transcript_id);
      if(exons != null)
      {
        for(int j = 0; j < exons.size(); j++)
        {
          exon = (Feature) exons.get(j);
          
          if(exon instanceof DatabaseInferredFeature)
            continue;
          exon_id = (String) exon.getQualifierByName("ID").getValues().get(0);

          exon_node = new DefaultMutableTreeNode(exon_id);
          transcript_node.add(exon_node);
        }
      }
      
      // utr & other features
      List utrs3  = chado_gene.get3UtrOfTranscript(transcript_id);
      List utrs5  = chado_gene.get5UtrOfTranscript(transcript_id);
      List others = chado_gene.getOtherFeaturesOfTranscript(transcript_id);
      
      List utrs  = new Vector();
      if(utrs3 != null)
        utrs.addAll(utrs3);
      if(utrs5 != null)
        utrs.addAll(utrs5);
      if(others != null)
        utrs.addAll(others);
      
      for(int j=0; j<utrs.size(); j++)
      {
        Feature utr = (Feature)utrs.get(j);
        if(utr.getQualifierByName("ID") == null)
          continue;
        String utr_id = (String)utr.getQualifierByName("ID").getValues().get(0);
        transcript_node.add(new DefaultMutableTreeNode(utr_id));
      }    
      
      // protein
      protein = chado_gene.getProteinOfTranscript(transcript_id);
      if(protein == null)
        continue;
      
      protein_id = 
          (String)((Feature)protein).getQualifierByName("ID").getValues().get(0);

      
      protein_node = new DefaultMutableTreeNode(protein_id);
      
      transcript_node.add(protein_node);
    }
  }
  
  /**
   * Change a node name.
   * @param old_id  the old uniquename of the feature
   * @param new_id  the new uniquename of the feature
   */
  protected void changeNode(String old_id, String new_id)
  {
    DefaultMutableTreeNode change_node = getNodeFromName(old_id);
    if(change_node != null)
    {
      change_node.setUserObject(new_id);
      repaint();
      revalidate();
    }
  }
  
  /**
   * Get the node in the tree with a given name. If it is not part
   * of the tree return null.
   * @param name
   * @return
   */
  private DefaultMutableTreeNode getNodeFromName(final String name)
  {
    DefaultMutableTreeNode root = (DefaultMutableTreeNode)getModel().getRoot();
    
    if(name.equals((String)root.getUserObject()))
      return root;
        
    DefaultMutableTreeNode change_node = searchChildren(root, name);
    if(change_node != null)
      return change_node;
    
    Enumeration root_children = root.children();
    while(root_children.hasMoreElements())
    {
      DefaultMutableTreeNode child = 
           (DefaultMutableTreeNode)root_children.nextElement();
      
      change_node = searchChildren(child, name);
      if(change_node != null)
        return change_node;
    }
    
    return null;
  }
  
  /**
   * Delete a node and all its descendents from the tree.
   * @param id  the uniquename of the feature being removed
   */
  protected void deleteNode(final String id)
  {
    DefaultMutableTreeNode root = (DefaultMutableTreeNode)getModel().getRoot();
    deleteChildNode(id,root);
  }
  
  /**
   * Recursively check children of the parent node against the ID
   * and delete the child that matches the ID.
   * @param id
   * @param parentNode
   */
  private void deleteChildNode(final String id,
                               final DefaultMutableTreeNode parentNode)
  {
    Enumeration root_children = parentNode.children();
    while(root_children.hasMoreElements())
    {
      DefaultMutableTreeNode child = 
           (DefaultMutableTreeNode)root_children.nextElement();
      
      // check children of this node!!
      if(child.getSiblingCount() > 0)
        deleteChildNode(id, child);
      
      if(id.equals((String)child.getUserObject()))
      {
        ((DefaultTreeModel)getModel()).removeNodeFromParent(child);
        return;
      }
    }
  }
  
  /**
   * Add a new node for the selected feature. It uses the Parent
   * qualifier information to work out where it should be added.
   * @param feature
   */
  protected void addNode(final uk.ac.sanger.artemis.Feature feature)
  {
    try
    {
      final String parent;
      
      if(feature.getQualifierByName("Parent") != null)
        parent = (String)feature.getQualifierByName("Parent").getValues().get(0);
      else
        parent = (String)feature.getQualifierByName("Derives_from").getValues().get(0);
      
      String name = 
        (String)feature.getQualifierByName("ID").getValues().get(0);
      
      if(getNodeFromName(name) != null)
        return;
      
      DefaultMutableTreeNode parentNode = getNodeFromName(parent);
      
      if(parentNode !=  null)
      {
        ((DefaultTreeModel)getModel()).insertNodeInto(
            new DefaultMutableTreeNode(name), parentNode, 
            parentNode.getChildCount());
      }
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
  }
  
  private DefaultMutableTreeNode searchChildren(
        final DefaultMutableTreeNode node,
        final String id)
  {
    Enumeration root_children = node.children();
    while(root_children.hasMoreElements())
    {
      DefaultMutableTreeNode child = 
           (DefaultMutableTreeNode)root_children.nextElement();

      if(id.equals((String)child.getUserObject()))
        return child;
    }
    return null;
  }
  
  class GeneTreeSelectionListener implements TreeSelectionListener
  {
    public void valueChanged(TreeSelectionEvent e)
    {
      DefaultMutableTreeNode node = 
        (DefaultMutableTreeNode)GeneComponentTree.this.getLastSelectedPathComponent();
      if(node == null)
        return;

      Feature embl_feature = 
            (Feature)chado_gene.getFeatureFromId((String)node.getUserObject());
      
      final uk.ac.sanger.artemis.Feature feature;
      
      if(embl_feature.getUserData() == null)
        feature = new uk.ac.sanger.artemis.Feature(embl_feature);    
      else
        feature =
          (uk.ac.sanger.artemis.Feature)embl_feature.getUserData();
      boolean isSet = true;
      if(feature.isReadOnly())
        isSet = false;
      
      gene_builder.setActiveFeature(feature, isSet);
      if(selection != null)
        selection.set(feature);
    }
  }
        
}
