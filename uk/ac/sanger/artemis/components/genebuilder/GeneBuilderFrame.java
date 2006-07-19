/* GeneBuilderFrame.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/GeneBuilderFrame.java,v 1.4 2006-07-19 16:05:17 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder;

import javax.swing.*;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.ItemListener;
import java.util.List;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryChangeEvent;
import uk.ac.sanger.artemis.EntryChangeListener;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.GotoEventSource;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.components.FeatureDisplay;
import uk.ac.sanger.artemis.components.FeatureEdit;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.FeatureVector;

public class GeneBuilderFrame extends JFrame
       implements EntryChangeListener, FeatureChangeListener
{
  
  private Feature active_feature; 
  private FeatureEdit feature_editor;
  
  public GeneBuilderFrame(final Feature feature,
                          final EntryGroup entry_group,
                          final Selection selection,
                          final GotoEventSource goto_event_source)
  {
    super("Artemis Gene Builder: " + feature.getIDString() +
          (feature.isReadOnly() ?
          "  -  (read only)" :
          ""));
    
    this.active_feature = feature;


    GFFStreamFeature gff_feature = (GFFStreamFeature)feature.getEmblFeature();
    final ChadoCanonicalGene chado_gene = gff_feature.getChadoGene();
    GeneComponentTree tree = new GeneComponentTree(chado_gene,
                                                   this);
    
    JScrollPane jsp_tree = new JScrollPane(tree);
    jsp_tree.setPreferredSize( new Dimension(150, jsp_tree.getPreferredSize().height) );
    getContentPane().add(jsp_tree, BorderLayout.WEST);
    
    GeneViewerPanel viewer = new GeneViewerPanel(gff_feature.getChadoGene());

    Box xBox = Box.createHorizontalBox();
    xBox.add(buildCheckBoxes(viewer, chado_gene, selection));
    xBox.add(viewer);
    
    JScrollPane jsp_viewer = new JScrollPane(xBox);
    jsp_viewer.setPreferredSize( viewer.getPreferredSize() );
    getContentPane().add(jsp_viewer, BorderLayout.CENTER);

    if(entry_group != null)
    {
      JTabbedPane tabpane = new JTabbedPane();
    
      feature_editor = new FeatureEdit(feature, entry_group,
                                       selection, goto_event_source, this);
      tabpane.addTab("Annotation", feature_editor);
      getContentPane().add(tabpane, BorderLayout.SOUTH);
    }
    
    addWindowListener(new WindowAdapter() 
    {
      public void windowClosing(WindowEvent event) 
      {
        stopListening();
        dispose();
      }
    });
    
    pack();
    setVisible(true);
  }
  
  /**
   * Create a Box component with the checkboxes used to determine if
   * a transcripts exons are visible in artemis.
   * @param viewer
   * @param chado_gene
   * @return
   */
  private Box buildCheckBoxes(final GeneViewerPanel viewer,
                              final ChadoCanonicalGene chado_gene,
                              final Selection selection)
  {
    Box yBox = Box.createVerticalBox();
    yBox.add(Box.createVerticalStrut(viewer.getViewerBorder()*3));
    
    java.util.List transcripts = chado_gene.getTranscripts();
    
    for(int i=0; i<transcripts.size(); i++)
    {
      final GFFStreamFeature transcript = 
        (GFFStreamFeature)((uk.ac.sanger.artemis.io.Feature)transcripts.get(i));
      
      JCheckBox cb = new JCheckBox();
      cb.setSelected(true);
      cb.addItemListener(new ItemListener()
      {
        public void itemStateChanged(ItemEvent e)
        {
          List exons = chado_gene.getExonsOfTranscript(
              (String)transcript.getQualifierByName("ID").getValues().get(0));
          
          boolean visible = true;
          if(e.getStateChange() == ItemEvent.DESELECTED)
            visible = false;        
          
          for(int j=0; j<exons.size(); j++)
          {
            GFFStreamFeature embl_exon = (GFFStreamFeature)exons.get(j);
            embl_exon.setVisible(visible);
          }
          
          FeatureVector fv = selection.getAllFeatures();
          selection.clear();
          selection.set(fv);
        }
      });
      yBox.add(cb);
      yBox.add(Box.createVerticalStrut(
          viewer.getTranscriptSize() - cb.getPreferredSize().height)); 
    }
    yBox.add(Box.createVerticalGlue());
    return yBox;
  }
  
  protected void setActiveFeature(final Feature active_feature)
  {
    if(this.active_feature != null)
      stopListening();
    
    this.active_feature = active_feature;
    feature_editor.setActiveFeature(active_feature);
  }

  /**
   *  Remove this object as a feature and entry change listener.
   **/
  private void stopListening() 
  {
    getEntry().removeEntryChangeListener(this);
    active_feature.removeFeatureChangeListener(this);
  }

  /**
   *  Implementation of the EntryChangeListener interface.  We listen to
   *  EntryChange events so we can notify the user if of this component if the
   *  feature gets deleted.
   **/
  public void entryChanged(EntryChangeEvent event) 
  {
    switch(event.getType())
    {
      case EntryChangeEvent.FEATURE_DELETED:
        if(event.getFeature() == active_feature) 
        {
          stopListening();
          dispose();
        }
        break;
      default:
        // do nothing;
        break;
    }
  }
  
  /**
   *  Implementation of the FeatureChangeListener interface.  We need to
   *  listen to feature change events from the Features in this object so that
   *  we can update the display.
   *  @param event The change event.
   **/
  public void featureChanged(FeatureChangeEvent event) 
  {
    active_feature.resetColour();
    /*
    // re-read the information from the feature
    switch(event.getType()) 
    {
      case FeatureChangeEvent.LOCATION_CHANGED:
        updateLocation();
        break;
      case FeatureChangeEvent.KEY_CHANGED:
        updateKey();
        break;
      case FeatureChangeEvent.QUALIFIER_CHANGED:
        if(qualifier_text_area.getText().equals(orig_qualifier_text)) 
          updateFromFeature();
        else
        {
          final String message =
            "warning: the qualifiers have changed outside the editor - " +
            "view now?";

          final YesNoDialog yes_no_dialog =
            new YesNoDialog(FeatureEdit.this, message);

          if(yes_no_dialog.getResult()) 
            new FeatureViewer(gene_feature);
        }
        break;
      default:
        updateFromFeature();
        break;
    }
    */
  }
  
  
  private Entry getEntry()
  {
    return active_feature.getEntry();
  }
  
  
}