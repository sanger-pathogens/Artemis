/* BasicPropertiesPanel
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/gff/PropertiesPanel.java,v 1.11 2009-08-17 12:50:42 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder.gff;

import java.awt.FlowLayout;

import javax.swing.Box;
import javax.swing.JPanel;

import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.components.genebuilder.BasicGeneBuilderFrame;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.QualifierVector;

public class BasicPropertiesPanel extends JPanel
{
  private static final long serialVersionUID = 1L;
  private PropertiesPanel genePropPanel;
  private PropertiesPanel transcriptPropPanel;
  private PropertiesPanel pepPropPanel;
  private Feature gene;
  
  public BasicPropertiesPanel(final ChadoCanonicalGene chadoGene,
		                      BasicGeneBuilderFrame gbFrame)
  {
    super(new FlowLayout(FlowLayout.LEFT));
    
    Box yBox = Box.createVerticalBox();
    gene = (Feature) chadoGene.getGene().getUserData();
    genePropPanel = new PropertiesPanel(gene, false, true, false);
    yBox.add(genePropPanel);
    
    Feature transcript = gbFrame.getSelectedTranscriptFeature();
    transcriptPropPanel = new PropertiesPanel(transcript, false, false, false);
    yBox.add(transcriptPropPanel);
    
    String transcriptName = GeneUtils.getUniqueName(transcript.getEmblFeature());
    uk.ac.sanger.artemis.io.Feature pep =
      chadoGene.getProteinOfTranscript(transcriptName);
    
    if(pep != null)
    {
      Feature protein = (Feature) pep.getUserData();
      pepPropPanel = new PropertiesPanel(protein, false, false, false);
      yBox.add(pepPropPanel);
    }
    
    add(yBox);
  }
  
  public QualifierVector getProteinProperties(Feature feature)
  {
    return pepPropPanel.getGffQualifiers(feature);
  }
  
  public QualifierVector getGeneProperties(Feature feature)
  {
    return genePropPanel.getGffQualifiers(feature);
  }
  
  public QualifierVector getTranscriptProperties(Feature feature)
  {
    return transcriptPropPanel.getGffQualifiers(feature);
  }
  
  public void updateObsoleteSettings()
  {
    if(!genePropPanel.obsoleteChanged)
      return;
    
    genePropPanel.obsoleteChanged = false;
    PropertiesPanel.updateObsoleteSettings((GFFStreamFeature)gene.getEmblFeature());
  }

  public void setObsoleteChanged(boolean obs, FeatureVector features)
  {
    PropertiesPanel.updateObsoleteSettings((GFFStreamFeature)gene.getEmblFeature());
    //obsoleteField.setSelected(obsoleteChanged);
  }
}