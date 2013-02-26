/* ValidateViewer
 *
 * created: 2013
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2013 Genome Research Limited
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
 */
package uk.ac.sanger.artemis.components;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.EntryGroupChangeEvent;
import uk.ac.sanger.artemis.EntryGroupChangeListener;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.ValidateFeature;


class ValidateViewer extends FileViewer implements EntryGroupChangeListener
{
  private static final long serialVersionUID = 1L;
  private EntryGroup entryGrp;
  
  /**
   * Viewer to display validation results
   * @param entryGrp
   * @param features
   */
  public ValidateViewer(final EntryGroup entryGrp,
                        final FeatureVector features)
  {
    super("Validation Report :: "+features.size()+
        " feature(s)", false, false, true);
    this.entryGrp = entryGrp;

    update(features);
    setVisible(true);

    final JButton fixButton = new JButton("Auto-Fix");
    fixButton.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e) 
      {
        try
        {
          entryGrp.getActionController().startAction();
          for(int i=0; i<features.size(); i++)
            fixBoundary(features.elementAt(i));
        }
        finally
        {
          //update(features);
          entryGrp.getActionController().endAction();
        }
      }
    });
    button_panel.add(fixButton);
 
    entryGrp.addEntryGroupChangeListener(new EntryGroupChangeListener(){
      public void entryGroupChanged(EntryGroupChangeEvent event)
      {
        switch(event.getType())
        {
          case EntryGroupChangeEvent.DONE_GONE:
            entryGrp.removeEntryGroupChangeListener(this);
            dispose();
            break;
        }
      }
    });
  }
  
  private void update(final FeatureVector features)
  {
    super.setText("");
    final ValidateFeature gffTest = new ValidateFeature(entryGrp);
    for(int i=0; i<features.size(); i++)
      gffTest.featureValidate(features.elementAt(i).getEmblFeature(), this);
  }
  
  private void fixBoundary(final Feature feature)
  {
    if( ! (feature.getEmblFeature() instanceof GFFStreamFeature) )
      return;
    
    final GFFStreamFeature gffFeature = (GFFStreamFeature) feature.getEmblFeature();
    final ChadoCanonicalGene gene = gffFeature.getChadoGene();
    if(gene != null && !gene.getGene().isReadOnly() && GeneUtils.isBoundaryOK(gene) > 0)
    {
      //updatedFeatures.add(feature);
      
      GeneUtils.checkGeneBoundary(gene, false);
    }
  }

  public void entryGroupChanged(EntryGroupChangeEvent event)
  {
    entryGrp.removeEntryGroupChangeListener(ValidateViewer.this);
    dispose();
  }
}