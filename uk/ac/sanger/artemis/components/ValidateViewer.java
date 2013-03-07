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
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import javax.swing.JButton;
import javax.swing.JCheckBox;

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
  private JCheckBox showFailedFeatures = new JCheckBox("Show only failed features", true);
  
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

    if( entryGrp == null || GeneUtils.isGFFEntry( entryGrp ) )
    {
      final JButton fixButton = new JButton("Auto-Fix Boundaries");
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
            entryGrp.getActionController().endAction();
            update(features);
          }
        }
      });
      button_panel.add(fixButton);
    }
    
    button_panel.add(showFailedFeatures);
    showFailedFeatures.addItemListener(new ItemListener(){
      public void itemStateChanged(ItemEvent arg0)
      {
        update(features);
      }
    });
    
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
    int nfail = 0;
    for(int i=0; i<features.size(); i++)
    {
      if(!gffTest.featureValidate(features.elementAt(i).getEmblFeature(), 
          this, showFailedFeatures.isSelected()))
        nfail++;
    }
    setTitle("Validation Report :: "+features.size()+
        " feature(s) Pass: "+(features.size()-nfail)+" Failed: "+nfail);
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