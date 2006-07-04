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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/GeneComponentTree.java,v 1.3 2006-07-04 15:57:57 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder;

import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.Feature;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.chado.ChadoFeature;

import javax.swing.*;
import javax.swing.event.*;

import javax.swing.tree.*;
import java.util.List;

/**
 * Tree to display a gene hierarchy.
 */
public class GeneComponentTree extends JTree
{
  
  public GeneComponentTree(final ChadoCanonicalGene chado_gene,
                           final GeneBuilderFrame gene_builder)
  {
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
    addTreeSelectionListener(new TreeSelectionListener()
    {
      public void valueChanged(TreeSelectionEvent e)
      {
        DefaultMutableTreeNode node = 
          (DefaultMutableTreeNode)GeneComponentTree.this.getLastSelectedPathComponent();
        if(node == null)
          return;
       
        Feature feature = (Feature)chado_gene.getFeatureFromId((String)node.getUserObject());
        gene_builder.setActiveFeature((uk.ac.sanger.artemis.Feature)feature.getUserData());
      }
    });
    
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
    List exons;
    Feature transcript;
    Feature exon;
    Feature protein;
    String transcript_id;
    String exon_id;
    String protein_id;
    DefaultMutableTreeNode transcript_node;
    DefaultMutableTreeNode exon_node;
    DefaultMutableTreeNode protein_node;
    
    for(int i=0; i<transcripts.size(); i++)
    {
      transcript = (Feature)transcripts.get(i);
      transcript_id = (String)transcript.getQualifierByName("ID").getValues().get(0);
      transcript_node = new DefaultMutableTreeNode(transcript_id);
      gene_node.add(transcript_node);
      
      exons = chado_gene.getExonsOfTranscript(transcript_id);
      if(exons == null)
        continue;

      for(int j=0; j<exons.size(); j++)
      {
        if(exons.get(j) instanceof Feature)
        {
          exon = (Feature)exons.get(j);
          exon_id = (String)exon.getQualifierByName("ID").getValues().get(0);
        }
        else  // ChadoFeature
          exon_id = ((ChadoFeature)exons.get(j)).getUniquename();
        
        exon_node = new DefaultMutableTreeNode(exon_id);
        transcript_node.add(exon_node);
      }
      
      protein = (Feature)chado_gene.getProteinOfTranscript(transcript_id);
      if(protein == null)
        continue;
      
      protein_id = (String)protein.getQualifierByName("ID").getValues().get(0);
      protein_node = new DefaultMutableTreeNode(protein_id);
      
      transcript_node.add(protein_node);
    }
  }
  
}
