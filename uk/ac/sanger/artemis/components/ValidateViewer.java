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

import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.EntryGroupChangeEvent;
import uk.ac.sanger.artemis.EntryGroupChangeListener;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.ValidateFeature;
import uk.ac.sanger.artemis.util.ReadOnlyException;


class ValidateViewer extends FileViewer implements EntryGroupChangeListener
{
  private static final long serialVersionUID = 1L;
  private EntryGroup entryGrp;
  private FeatureVector selectedFeatures;
  private JCheckBox showFailedFeatures = new JCheckBox("Show only failed features", true);
  private boolean inAutoFix = false;
  /**
   * Viewer to display validation results
   * @param entryGrp
   * @param features
   */
  public ValidateViewer(final EntryGroup entryGrp,
                        final FeatureVector selectedFeatures)
  {
    super("Validation Report :: ", false, false, true);
    this.entryGrp = entryGrp;
    this.selectedFeatures = selectedFeatures;
    
    final boolean allFeatures = (selectedFeatures == null);

    update();
    setVisible(true);

    final JButton fixButton = new JButton("Auto-Fix");
    fixButton.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e) 
      {
        final JPanel options = new JPanel( new GridLayout(2,1) );
        final JCheckBox extendToNextStop = new JCheckBox("Fix stop codon", true);
        options.add(extendToNextStop);
        final JCheckBox boundary = new JCheckBox("Gene Model boundaries", true);
        options.add(boundary);
        if( entryGrp != null && !GeneUtils.isGFFEntry( entryGrp ) )
        {
          boundary.setSelected(false);
          boundary.setEnabled(false);
        }
        
        int status = 
            JOptionPane.showConfirmDialog(ValidateViewer.this, options, 
                "Auto-Fix", JOptionPane.OK_CANCEL_OPTION);
        if(status == JOptionPane.CANCEL_OPTION)
          return;
          
        try
        {
          ValidateFeature validate = new ValidateFeature(entryGrp);
          entryGrp.getActionController().startAction();
          inAutoFix = true;
          final FeatureVector features = getFeatures();
          for(int i=0; i<features.size(); i++)
          {
            final Feature f = features.elementAt(i);
            if( extendToNextStop.isSelected() && 
                !validate.hasValidStop(f.getEmblFeature()) )
              f.fixStopCodon();
            
            if(boundary.isSelected())
              fixBoundary(f);
          }
        }
        catch (ReadOnlyException e1)
        {
          JOptionPane.showMessageDialog(ValidateViewer.this, 
              "Read only entry", "Error", JOptionPane.ERROR_MESSAGE);
        }
        finally
        {
          inAutoFix = false;
          entryGrp.getActionController().endAction();
          update();
        }
      }
    });
    button_panel.add(fixButton);

    
    button_panel.add(showFailedFeatures);
    showFailedFeatures.addItemListener(new ItemListener(){
      public void itemStateChanged(ItemEvent arg0)
      {
        update();
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
    
    entryGrp.addFeatureChangeListener(new FeatureChangeListener(){
      public void featureChanged(FeatureChangeEvent event)
      {
        if(inAutoFix)
          return;

        inAutoFix = true;
        ValidateViewer.this.selectedFeatures = null;

        SwingUtilities.invokeLater(new Runnable() 
        {
          public void run() 
          {
            update();
            inAutoFix = false;
          }
        });
      }
    });
  }
  
  private void update()
  {
    final FeatureVector features = getFeatures();
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
  
  private FeatureVector getFeatures()
  {
    if(selectedFeatures == null)
      selectedFeatures = entryGrp.getAllFeatures();
    return selectedFeatures;
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