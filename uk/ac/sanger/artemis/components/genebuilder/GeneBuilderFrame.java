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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/GeneBuilderFrame.java,v 1.6 2006-07-25 10:50:20 tjc Exp $
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
import uk.ac.sanger.artemis.chado.ChadoFeature;
import uk.ac.sanger.artemis.components.FeatureDisplay;
import uk.ac.sanger.artemis.components.FeatureEdit;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.FeatureVector;

public class GeneBuilderFrame extends JFrame
       implements EntryChangeListener, FeatureChangeListener
{
  
  private Feature active_feature; 
  private FeatureEdit feature_editor;
  private GeneViewerPanel viewer;
  private GeneComponentTree tree;
  
  public GeneBuilderFrame(final Feature feature,
                          final EntryGroup entry_group,
                          Selection selection,
                          final GotoEventSource goto_event_source)
  {
    super("Artemis Gene Builder: " + feature.getIDString() +
          (feature.isReadOnly() ?
          "  -  (read only)" :
          ""));
    
    this.active_feature = feature;


    GFFStreamFeature gff_feature = (GFFStreamFeature)feature.getEmblFeature();
    final ChadoCanonicalGene chado_gene = gff_feature.getChadoGene();
    
   
    try
    {
      addListeners(chado_gene);
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
    
    
    tree = new GeneComponentTree(chado_gene, this);
    
    JScrollPane jsp_tree = new JScrollPane(tree);
    jsp_tree.setPreferredSize( new Dimension(150, jsp_tree.getPreferredSize().height) );
    getContentPane().add(jsp_tree, BorderLayout.WEST);
    
    if(selection == null)
      selection = new Selection(null);
    
    viewer = new GeneViewerPanel(
                gff_feature.getChadoGene(), selection);

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
    
    // add menus
    JMenuBar menuBar = new JMenuBar();
    
    JMenu editMenu = new JMenu("Edit");
    menuBar.add(editMenu);
    viewer.createMenus(editMenu);
    setJMenuBar(menuBar);
    
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
          
          if(exons == null)
            return;
          
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
        QualifierVector qualifiers = event.getFeature().getQualifiers();
        tree.deleteNode( (String)qualifiers.getQualifierByName("ID").getValues().get(0) );
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
   
    // re-read the information from the feature
    switch(event.getType()) 
    {
      case FeatureChangeEvent.LOCATION_CHANGED:
        break;
      case FeatureChangeEvent.KEY_CHANGED:
        break;
      case FeatureChangeEvent.QUALIFIER_CHANGED:
        break;
      case FeatureChangeEvent.SEGMENT_CHANGED:
        Feature feature = event.getFeature();    
        QualifierVector old_qualifiers = event.getOldQualifiers();
        QualifierVector new_qualifiers = feature.getQualifiers();
      
        tree.changeNode( 
            (String)old_qualifiers.getQualifierByName("ID").getValues().get(0),
            (String)new_qualifiers.getQualifierByName("ID").getValues().get(0));
      default:
        break;
    }
    viewer.repaint();
  }
  
  
  private Entry getEntry()
  {
    return active_feature.getEntry();
  }
  
  /**
   * Add feature listeners for each artemis feature.
   * @throws InvalidRelationException
   */
  private void addListeners(final ChadoCanonicalGene chado_gene) 
               throws InvalidRelationException
  {
    // add feature listeners
    uk.ac.sanger.artemis.io.Feature embl_gene = 
      (uk.ac.sanger.artemis.io.Feature)chado_gene.getGene();
    Feature gene = (Feature)embl_gene.getUserData();
    gene.addFeatureChangeListener(this);
    gene.getEntry().addEntryChangeListener(this);
    
    List transcripts = chado_gene.getTranscripts();
    for(int i=0; i<transcripts.size(); i++)
    {
      uk.ac.sanger.artemis.io.Feature transcript =
         (uk.ac.sanger.artemis.io.Feature)transcripts.get(i);
      Feature trans = (Feature)transcript.getUserData();
       
       if(trans == null)
         trans = new Feature(transcript);
       
       trans.addFeatureChangeListener(this);
       trans.getEntry().addEntryChangeListener(this);
       
       List exons = chado_gene.getExonsOfTranscript(
           (String)trans.getQualifierByName("ID").getValues().get(0));
       
       if(exons == null)
         continue;
       
       if(exons.get(0) instanceof ChadoFeature)
         return;
       
       for(int j=0; j<exons.size(); j++)
       {
         uk.ac.sanger.artemis.io.Feature embl_exon = 
            (uk.ac.sanger.artemis.io.Feature)exons.get(j);

         Feature exon = (Feature)embl_exon.getUserData();
         
         if(exon == null)
           exon = new Feature(embl_exon);
         exon.addFeatureChangeListener(this);
         exon.getEntry().addEntryChangeListener(this);
       }
    }
  }
  
}