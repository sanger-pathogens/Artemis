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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/GeneBuilderFrame.java,v 1.24 2007-04-27 13:57:19 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder;

import javax.swing.*;
import javax.swing.border.Border;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.ItemListener;
import java.util.List;
import java.util.Hashtable;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryChangeEvent;
import uk.ac.sanger.artemis.EntryChangeListener;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.SelectionChangeEvent;
import uk.ac.sanger.artemis.SelectionChangeListener;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.GotoEventSource;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.components.FeatureEdit;
import uk.ac.sanger.artemis.components.Utilities;
import uk.ac.sanger.artemis.io.GFFEntryInformation;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.FeatureVector;

public class GeneBuilderFrame extends JFrame
       implements EntryChangeListener, FeatureChangeListener
{
  
  /** */
  private static final long serialVersionUID = 1L;
  private Feature active_feature; 
  private FeatureEdit feature_editor;
  private GeneViewerPanel viewer;
  private GeneComponentTree tree;
  private Box yBox;
  private Hashtable transcriptBoxes = new Hashtable();
  private Component glue = Box.createVerticalGlue();
  private Selection selection;
  private ChadoCanonicalGene chado_gene;
  private JLabel status_line = new JLabel("");
  
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
    this.selection = selection;
    
    if(selection != null)
      selection.addSelectionChangeListener(new SelectionChangeListener()
      {
        public void selectionChanged(SelectionChangeEvent event)
        {
          viewer.repaint();
          tree.setSelection(GeneBuilderFrame.this.selection);
        }
      });
    
    
    GFFStreamFeature gff_feature = (GFFStreamFeature)feature.getEmblFeature();
    chado_gene = gff_feature.getChadoGene();
    
    try
    {
      addListeners(chado_gene);
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
    
    tree = new GeneComponentTree(chado_gene, this, selection);
    
    JScrollPane jsp_tree = new JScrollPane(tree);
    jsp_tree.setPreferredSize( new Dimension(150, jsp_tree.getPreferredSize().height) );
    
    if(selection == null)
      selection = new Selection(null);
    
    viewer = new GeneViewerPanel(
                gff_feature.getChadoGene(), selection, entry_group, this, status_line);

    Box xBox = Box.createHorizontalBox();
    xBox.add(buildCheckBoxes(viewer, chado_gene));
    xBox.add(viewer);
    
    JScrollPane jsp_viewer = new JScrollPane(xBox);
    jsp_viewer.getViewport().setBackground(Color.white);
    jsp_viewer.setPreferredSize( viewer.getPreferredSize() );
    
    ///
    status_line.setFont(Options.getOptions().getFont());
    final FontMetrics fm =
      this.getFontMetrics(status_line.getFont());

    final int font_height = fm.getHeight()+10;

    status_line.setMinimumSize(new Dimension(100, font_height));
    status_line.setPreferredSize(new Dimension(100, font_height));

    Border loweredbevel = BorderFactory.createLoweredBevelBorder();
    Border raisedbevel = BorderFactory.createRaisedBevelBorder();
    Border compound = BorderFactory.createCompoundBorder(raisedbevel,loweredbevel);
    status_line.setBorder(compound);
    jsp_viewer.setColumnHeaderView(status_line);
    ///
    
    JSplitPane top = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, jsp_tree, jsp_viewer);
    JSplitPane all = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
    all.setTopComponent(top);
    
    if(entry_group != null)
      feature_editor = new FeatureEdit(feature, entry_group,
                        selection, goto_event_source, this);
    else
      feature_editor = new FeatureEdit(feature, null,
          null, null, this, new GFFEntryInformation());
    
    JTabbedPane tabpane = new JTabbedPane();
    tabpane.addTab("Annotation", feature_editor);
    all.setBottomComponent(tabpane);
    
    getContentPane().add(all);
    
    addWindowListener(new WindowAdapter() 
    {
      public void windowClosing(WindowEvent event) 
      {
        try
        {
          stopListeningAll();
         
        }
        catch(InvalidRelationException e)
        {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }
        dispose();
      }
    });
    
    // add menus
    JMenuBar menuBar = new JMenuBar();
    
    JMenu fileMenu = new JMenu("File");
    menuBar.add(fileMenu);
    JMenuItem close = new JMenuItem("Close");
    close.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        dispose();
      }
    });
    fileMenu.add(close);
    
    JMenu editMenu = new JMenu("Edit");
    menuBar.add(editMenu);
    viewer.createMenus(editMenu, entry_group);
    setJMenuBar(menuBar);
    
    pack();
    
    Utilities.centreFrame(this);
    setVisible(true);
    all.setDividerLocation(0.35);
  }
  
  /**
   * Create a Box component with the checkboxes used to determine if
   * a transcripts exons are visible in artemis.
   * @param viewer
   * @param chado_gene
   * @return
   */
  private Box buildCheckBoxes(final GeneViewerPanel viewer,
                              final ChadoCanonicalGene chado_gene)
  {
    yBox = Box.createVerticalBox();
    yBox.add(Box.createVerticalStrut(viewer.getViewerBorder()*3));
    
    java.util.List transcripts = chado_gene.getTranscripts();
    
    for(int i=0; i<transcripts.size(); i++)
    {
      final GFFStreamFeature transcript = 
        (GFFStreamFeature)((uk.ac.sanger.artemis.io.Feature)transcripts.get(i));
      final String transcript_name =
        (String)transcript.getQualifierByName("ID").getValues().get(0);
      
      Box transcriptBox =  createTranscriptBox(transcript, transcript_name);
      yBox.add(transcriptBox);
    }
    yBox.add(glue);
    return yBox;
  }
  
  /**
   * Create a Box containing the checkbox for a transcript
   * @param transcript
   * @param transcript_name
   * @return
   */
  private Box createTranscriptBox(final GFFStreamFeature transcript,
                                  final String transcript_name)
  {
    JCheckBox cb = new JCheckBox();
    cb.setOpaque(false);
    cb.setSelected(transcript.isVisible());
    cb.addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent e)
      {
        List exons = chado_gene.getSpliceSitesOfTranscript(transcript_name, "exon");
        
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
    
    Box transcriptBox = Box.createVerticalBox();

    transcriptBox.add(cb);
    transcriptBox.add(Box.createVerticalStrut(
        viewer.getTranscriptSize() - cb.getPreferredSize().height)); 
    transcriptBoxes.put(transcript_name, transcriptBox);
    return transcriptBox;
  }
  
  protected void setActiveFeature(final Feature active_feature,
                                  final boolean isSet)
  {  
    this.active_feature = active_feature;
    feature_editor.setActiveFeature(active_feature, isSet);
  }

  /**
   *  Remove this object as a feature and entry change listener.
   **/
  private void stopListening(final Feature feature) 
  {
    if(getEntry() != null)
      getEntry().removeEntryChangeListener(this);
    feature.removeFeatureChangeListener(this);
  }
  
  private void startListening(final Feature feature) 
  {
    if(getEntry() != null)
      getEntry().addEntryChangeListener(this);
    feature.addFeatureChangeListener(this);
  }

  /**
   *  Implementation of the EntryChangeListener interface.  We listen to
   *  EntryChange events so we can notify the user if of this component if the
   *  feature gets deleted.
   **/
  public void entryChanged(EntryChangeEvent event) 
  {
    Feature feature = event.getFeature();
    stopListening(feature);
    
    QualifierVector qualifiers = feature.getQualifiers();
    String name = (String)qualifiers.getQualifierByName("ID").getValues().get(0);
    switch(event.getType())
    {
      case EntryChangeEvent.FEATURE_DELETED:
        feature.removeFeatureChangeListener(this);
        tree.deleteNode( name );
        
        if(chado_gene.isTranscript(name))
        {
          Box transBox = (Box)transcriptBoxes.get(name);
          yBox.remove(transBox);
          yBox.revalidate();
        }
        break;
      case EntryChangeEvent.FEATURE_ADDED:
        Qualifier parent_qualifier = 
          qualifiers.getQualifierByName("Parent");
        
        if(parent_qualifier == null)
          return;
        tree.addNode(event.getFeature());
        feature.addFeatureChangeListener(this);
        
        String parent = (String)parent_qualifier.getValues().get(0);
        String gene_name = null;
        try
        {
          gene_name = 
            (String)chado_gene.getGene().getQualifierByName("ID").getValues().get(0);
        }
        catch(InvalidRelationException e)
        {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }
        
        if(parent.equals(gene_name) &&
           !transcriptBoxes.containsKey(name))
        {
          Box transcriptBox = createTranscriptBox(
                   (GFFStreamFeature)feature.getEmblFeature(), name);
          yBox.remove(glue);
          yBox.add(transcriptBox);
          yBox.add(glue);
          yBox.revalidate();
        }
      default:
        // do nothing;
        break;
    }
    
    startListening(feature);
  }
  
  /**
   *  Implementation of the FeatureChangeListener interface.  We need to
   *  listen to feature change events from the Features in this object so that
   *  we can update the display.
   *  @param event The change event.
   **/
  public void featureChanged(FeatureChangeEvent event) 
  {
    Feature feature = event.getFeature();
    stopListening(feature);
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
        QualifierVector old_qualifiers = event.getOldQualifiers();
        QualifierVector new_qualifiers = feature.getQualifiers();
        tree.changeNode( 
            (String)old_qualifiers.getQualifierByName("ID").getValues().get(0),
            (String)new_qualifiers.getQualifierByName("ID").getValues().get(0));
      default:
        break;
    }
    viewer.repaint();
    startListening(feature);
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
    
    if(gene.getEntry() != null)
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
       
       if(trans.getEntry() != null)
         trans.getEntry().addEntryChangeListener(this);
       
       List exons = chado_gene.getSpliceSitesOfTranscript(
           (String)trans.getQualifierByName("ID").getValues().get(0), "exon");
       
       if(exons == null || exons.size() < 1)
         continue;
       
       if(exons.get(0) instanceof org.gmod.schema.sequence.Feature)
         return;
       
       for(int j=0; j<exons.size(); j++)
       {
         uk.ac.sanger.artemis.io.Feature embl_exon = 
            (uk.ac.sanger.artemis.io.Feature)exons.get(j);

         Feature exon = (Feature)embl_exon.getUserData();
         
         if(exon == null)
           exon = new Feature(embl_exon);
         exon.addFeatureChangeListener(this);
         
         if(exon.getEntry() != null)
           exon.getEntry().addEntryChangeListener(this);
       }
    }
  }
  
  private void stopListeningAll() throws InvalidRelationException
  {
    uk.ac.sanger.artemis.io.Feature embl_gene = 
      (uk.ac.sanger.artemis.io.Feature)chado_gene.getGene();
    Feature gene = (Feature)embl_gene.getUserData();   
    stopListening(gene);
    
    List transcripts = chado_gene.getTranscripts();
    for(int i=0; i<transcripts.size(); i++)
    {
      uk.ac.sanger.artemis.io.Feature transcript =
         (uk.ac.sanger.artemis.io.Feature)transcripts.get(i);
      Feature trans = (Feature)transcript.getUserData();
       
      stopListening(trans);
      
      List exons = chado_gene.getSpliceSitesOfTranscript(
          (String)trans.getQualifierByName("ID").getValues().get(0), "exon");
      
      if(exons == null)
        continue;
      
      for(int j=0; j<exons.size(); j++)
      {
        uk.ac.sanger.artemis.io.Feature embl_exon = 
           (uk.ac.sanger.artemis.io.Feature)exons.get(j);

        Feature exon = (Feature)embl_exon.getUserData();
        stopListening(exon);      
      }
    }
  }

}
