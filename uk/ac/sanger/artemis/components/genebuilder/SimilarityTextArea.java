/* SimilarityTextArea
 *
 * created: 2007
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2007  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/SimilarityTextArea.java,v 1.1 2007-02-09 15:17:55 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder;

import javax.swing.JSplitPane;
import javax.swing.JTextArea;

import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.io.GFFStreamFeature;

public class SimilarityTextArea extends JTextArea
             implements FeatureChangeListener
{

  /** */
  private static final long serialVersionUID = 1L;
  private JSplitPane splitPane;
  
  public SimilarityTextArea(final Feature gffFeature,
                            final JSplitPane splitPane)
  {
    super();
    setLineWrap (true);
    this.splitPane = splitPane;
    updateFromFeature(gffFeature);
  }
  
  public void updateFromFeature(final Feature feature)
  {
    feature.removeFeatureChangeListener(this);
    GFFStreamFeature gffFeature = (GFFStreamFeature)feature.getEmblFeature();
    final String text = uk.ac.sanger.artemis.chado.Similarity.getSimilarityString(gffFeature);
    setText(text);
    if(gffFeature.getSimilarityFeatures() == null)
    {
      setRows(0);
      splitPane.setDividerLocation(1.0);
    }
    else
    {
      splitPane.setDividerLocation(0.5);
      setRows((gffFeature.getSimilarityFeatures().size()*2)+1);
    }
    feature.addFeatureChangeListener(this);
  }
  

  public void featureChanged(FeatureChangeEvent event)
  {
    updateFromFeature( event.getFeature() );
  }
  
}