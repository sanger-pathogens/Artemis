package uk.ac.sanger.artemis.components;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;

import uk.ac.sanger.artemis.EntryChangeEvent;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.EntryGroupChangeEvent;
import uk.ac.sanger.artemis.EntryGroupChangeListener;
import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.ValidateFeature;


class ValidateViewer extends FileViewer implements EntryGroupChangeListener
{
  private static final long serialVersionUID = 1L;
  private EntryGroup entryGrp;
  
  public ValidateViewer(final String label, 
                        final boolean visible, 
                        final boolean showClearButton,
                        final boolean showSaveButton, 
                        final EntryGroup entryGrp,
                        final FeatureVector features)
  {
    super(label, visible, showClearButton, showSaveButton);
    this.entryGrp = entryGrp;

    update(features);
    setVisible(true);
    
    
    final JButton fixButton = new JButton("Auto-Fix");
    fixButton.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e) 
      {
        //entryGrp.getActionController().startAction();
        for(int i=0; i<features.size(); i++)
        {
          if(features.elementAt(i).getEmblFeature() instanceof GFFStreamFeature)
          {
            fixBoundary((GFFStreamFeature)features.elementAt(i).getEmblFeature());
          }
        }
        
        //entryGrp.entryChanged(new EntryChangeEvent());
        //entryGrp.getActionController().endAction();
      }
    });
    button_panel.add(fixButton);

    entryGrp.addFeatureChangeListener(new FeatureChangeListener(){
      public void featureChanged(FeatureChangeEvent event)
      {
        System.out.println("HEREREE");
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
    for(int i=0; i<features.size(); i++)
      gffTest.featureValidate(features.elementAt(i).getEmblFeature(), this);
  }
  
  private void fixBoundary(final GFFStreamFeature gffFeature)
  {
    final ChadoCanonicalGene gene = gffFeature.getChadoGene();
    if(gene != null && !gene.getGene().isReadOnly() && GeneUtils.isBoundaryOK(gene) > 0)
    {
      System.out.println(GeneUtils.getUniqueName(gffFeature));
      
      GeneUtils.checkGeneBoundary(gene);
    }
  }

  public void entryGroupChanged(EntryGroupChangeEvent event)
  {
    entryGrp.removeEntryGroupChangeListener(ValidateViewer.this);
    dispose();
  }
}