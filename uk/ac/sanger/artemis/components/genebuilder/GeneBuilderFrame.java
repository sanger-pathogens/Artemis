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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/GeneBuilderFrame.java,v 1.3 2006-07-06 15:10:14 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder;

import javax.swing.*;

import java.awt.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryChangeEvent;
import uk.ac.sanger.artemis.EntryChangeListener;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.GotoEventSource;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.components.FeatureEdit;
import uk.ac.sanger.artemis.io.GFFStreamFeature;


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
    GeneComponentTree tree = new GeneComponentTree(gff_feature.getChadoGene(),
                                                   this);
    
    JScrollPane jsp_tree = new JScrollPane(tree);
    jsp_tree.setPreferredSize( new Dimension(150, jsp_tree.getPreferredSize().height) );
    getContentPane().add(jsp_tree, BorderLayout.WEST);
    
    GeneViewerPanel viewer = new GeneViewerPanel(gff_feature.getChadoGene());
    JScrollPane jsp_viewer = new JScrollPane(viewer);
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