/* ProteinMapPanel.java
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
 **/

package uk.ac.sanger.artemis.components.genebuilder;

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.util.List;

import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.io.Feature;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;

public class BasicProteinMapPanel extends ProteinMapPanel
{
  private static final long serialVersionUID = 1L;
  private BasicGeneBuilderFrame gbFrame;
  
  public BasicProteinMapPanel(final Feature feature, 
                         final ChadoCanonicalGene chado_gene,
                         final Selection selection,
                         final BasicGeneBuilderFrame gbFrame) 
  {
    super(feature, chado_gene, selection);
    this.gbFrame = gbFrame;
    setFont(uk.ac.sanger.artemis.Options.getOptions().getFont());
    setPreferredSize(new Dimension(getSize().width, 100));
  }

  public void paintComponent(Graphics g)
  {
    super.paintComponent(g);
    Graphics2D g2d = (Graphics2D)g;
    toolTipPositions.clear();

    Feature embl_gene = (Feature)chado_gene.getGene();
    uk.ac.sanger.artemis.Feature gene =
      (uk.ac.sanger.artemis.Feature)embl_gene.getUserData();

    int ypos = border;
    int geneStart = embl_gene.getFirstBase();
    int geneEnd   = embl_gene.getLastBase();
    float fraction = (float)(getSize().width - (2*border))/
                     (float)(geneEnd-geneStart);

    uk.ac.sanger.artemis.Feature transcript = gbFrame.getSelectedTranscriptFeature();
    if(transcript == null)
      return;
    String transcriptName = GeneUtils.getUniqueName(
        transcript.getEmblFeature()); 
    Feature protein_embl_feature = 
       (Feature)chado_gene.getProteinOfTranscript(transcriptName);

    if(protein_embl_feature == null)
      return;

    uk.ac.sanger.artemis.Feature protein = 
      (uk.ac.sanger.artemis.Feature)protein_embl_feature.getUserData();

    List<Feature> exons = chado_gene.getSpliceSitesOfTranscript(transcriptName, 
                                          DatabaseDocument.EXONMODEL);
    if(exons == null || exons.size() == 0)
      exons = chado_gene.getSpliceSitesOfTranscript(transcriptName,
                                                 "pseudogenic_exon");

    if(exons == null || exons.size() == 0)
      return;

    Feature exon_embl_feature = exons.get(0);
    uk.ac.sanger.artemis.Feature exon = 
      (uk.ac.sanger.artemis.Feature)exon_embl_feature.getUserData();

    int ppLength = exon.getTranslationBasesLength()/3;
    int emblStart = protein_embl_feature.getFirstBase();
    int emblEnd   = protein_embl_feature.getLastBase();

    // draw protein
    g2d.drawString(protein.getIDString()+
        ( GeneUtils.isObsolete((GFFStreamFeature)embl_gene) ? " (obsolete)" : "") , border, ypos);

    int ppStart = border+(int)((emblStart-geneStart)*fraction);
    int ppEnd   = border+(int)((emblEnd-geneStart)*fraction);

    drawFeature(g2d, ppStart, ppEnd, 
                ypos, protein.getColour(), 1, 
                selection.contains(gene), 2.f, getFontHeight());

    //record position for tooltip
    Rectangle r = new Rectangle(ppStart, ypos, ppEnd-ppStart, 1*getFontHeight());
    toolTipPositions.put(r, protein.getIDString()+" length="+ppLength);

    ypos += border*2;    
    ypos = drawDomain(protein_embl_feature, 
            g2d, ypos, ppStart, ppEnd, ppLength);
      
    QualifierVector qualifiers = protein_embl_feature.getQualifiers();
    if(qualifiers.getQualifierByName(TMHMM[0]) != null)
    {
      g2d.drawString("Transmembrane Domains:", ppStart, ypos);
      drawPrediction(protein_embl_feature, 
                     g2d, ypos, ppStart, ppEnd, ppLength, TMHMM);
      ypos += border * 2;
    }
      
    if(qualifiers.getQualifierByName(SIGNALP[0]) != null)
    {
      g2d.drawString("SignalP:", ppStart, ypos);
      drawPrediction(protein_embl_feature, 
                   g2d, ypos, ppStart, ppEnd, ppLength, 
                   new String[]{ SIGNALP[0] });
      ypos += border * 2;
    }

    Qualifier gpiAnchor;
    if((gpiAnchor=qualifiers.getQualifierByName(GPI_ANCHORED)) != null)
    {
      g2d.drawString("GPI anchor cleavage site:", ppStart, ypos);
      drawGPIArrow(g2d, gpiAnchor, ppStart, ppEnd, ppLength, ypos);
    }
  }
}