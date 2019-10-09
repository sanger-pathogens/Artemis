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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/GeneBuilderFrame.java,v 1.47 2009-09-24 15:01:27 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder;


import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.border.Border;

import java.awt.Color;
import java.awt.Component;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Toolkit;
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
import uk.ac.sanger.artemis.chado.ChadoTransactionManager;
import uk.ac.sanger.artemis.components.FeatureEdit;
import uk.ac.sanger.artemis.components.Utilities;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;
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
  private GeneBuilderSelectionChangeListener geneBuilderSelectionChangeListener;
  private JTabbedPane tabpane;
  private ChadoTransactionManager chadoTransactionManager;
  private EntryGroup entry_group;
  private GotoEventSource goto_event_source;
  private Hashtable geneBuilderHash;
  
  public GeneBuilderFrame(Feature feature,
                          final EntryGroup entry_group,
                          final Selection selection,
                          final GotoEventSource goto_event_source)
  {
    this(feature, entry_group, selection, goto_event_source, null);
  }
  
  public GeneBuilderFrame(Feature feature,
      final EntryGroup entry_group,
      final Selection selection,
      final GotoEventSource goto_event_source,
      final ChadoTransactionManager chadoTransactionManager)
  {
    super();
    
    this.entry_group = entry_group;
    this.goto_event_source = goto_event_source;
    
    // set title
    final String title = "Artemis Gene Builder: " + 
          ((Feature)((GFFStreamFeature)feature.getEmblFeature()).getChadoGene().getGene().getUserData()).getIDString() +
           (feature.isReadOnly() ?
          "  -  (read only)" : "");

    setTitle(title);
    
    this.selection = selection;
    if(feature.getKey().getKeyString().equals(DatabaseDocument.EXONMODEL) ||
       (DatabaseDocument.CHADO_INFER_CDS && feature.getKey().getKeyString().equals("CDS")))
    {
      Feature proteinFeature = getProteinFeature(feature, selection);
      if(proteinFeature != null)
        feature = proteinFeature;
    }
    
    this.active_feature = feature;
    this.chadoTransactionManager = chadoTransactionManager;
 
    if(selection != null)
    {
      geneBuilderSelectionChangeListener = new GeneBuilderSelectionChangeListener();
      selection.addSelectionChangeListener(geneBuilderSelectionChangeListener);
    }
    
    final GFFStreamFeature gff_feature = (GFFStreamFeature)feature.getEmblFeature();
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
    final JScrollPane jsp_tree = new JScrollPane(tree);
    jsp_tree.setPreferredSize( new Dimension(150, jsp_tree.getPreferredSize().height) );
    
    viewer = new GeneViewerPanel(
                gff_feature.getChadoGene(), selection, 
                entry_group, this, status_line);

    Box xBox = Box.createHorizontalBox();
    xBox.add(buildCheckBoxes(viewer, chado_gene));
    xBox.add(viewer);
    
    final JTabbedPane tabMapViewer = new JTabbedPane();
    
    final JScrollPane jsp_viewer = new JScrollPane(xBox);
    jsp_viewer.getViewport().setBackground(Color.white);
    jsp_viewer.setPreferredSize( viewer.getPreferredSize() );
    tabMapViewer.addTab("Gene Map", jsp_viewer);
    
    //
    // protein map
    List proteins = ProteinMapPanel.getProteinsWithProteinMapElement(gff_feature);
    if(proteins != null)
    {
      final ProteinMapPanel proteinMap = 
        new ProteinMapPanel((GFFStreamFeature) proteins.get(0), 
                            gff_feature.getChadoGene(), selection);
      
      final JScrollPane jsp_ppviewer = new JScrollPane(proteinMap);
      jsp_ppviewer.getViewport().setBackground(Color.white);
      jsp_ppviewer.setPreferredSize( proteinMap.getPreferredSize() );
      tabMapViewer.addTab("Protein Map", jsp_ppviewer);
    }
    
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
    
    JSplitPane top = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, jsp_tree, tabMapViewer);
    JSplitPane all = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
    all.setTopComponent(top);
    
    feature_editor = new FeatureEdit(feature, entry_group,
                        selection, goto_event_source, this);
   
    feature_editor.getQualifierTextArea().getDocument().addDocumentListener(
        new TextAreaDocumentListener(feature_editor.getQualifierTextArea()));

    tabpane = new JTabbedPane();
    tabpane.addTab("Annotation", feature_editor);
    setTabTitle();
    all.setBottomComponent(tabpane);
    
    getContentPane().add(all);
    
    addWindowListener(new WindowAdapter() 
    {
      public void windowClosing(WindowEvent event) 
      {
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
    
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();

    int height = getPreferredSize().height;
    // Default size adjustment so all components are visible.
    int width  = (int)(getPreferredSize().width*1.18);
    if(height > (screen.height*0.9))
    {    
      setSize(width, (int)(screen.height*0.9));
    }
    else
    {
      setSize(width, height);
    }

    Utilities.centreFrame(this);
    setVisible(true);
    all.setDividerLocation(0.30);
  }
  
  
  private Feature getProteinFeature(final Feature feature, final Selection selection)
  {
    Feature proteinFeature = null;
    try
    {
      Qualifier idQualifier = feature.getQualifierByName("ID");
      if(idQualifier != null)
      {
        final ChadoCanonicalGene chadoGene = ((GFFStreamFeature)feature.getEmblFeature()).getChadoGene();
        final String ID = (String)idQualifier.getValues().get(0);
        final String transcriptName = chadoGene.getTranscriptFromName(ID);
        proteinFeature = (Feature)chadoGene.getProteinOfTranscript(transcriptName).getUserData();
        selection.clear();
        selection.add(proteinFeature);
      }
    }
    catch(Exception e){}
    return proteinFeature;
  }
  
  public void dispose()
  {
    dispose(false);
  }
  
  /**
   * Override to ensure the GeneBuilderFrame removes all listeners
   */
  public void dispose(final boolean reopen)
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
    
    if(chadoTransactionManager != null &&
       chadoTransactionManager.hasTransactions())
    { 
      int select = JOptionPane.showConfirmDialog(this, 
          "Commit changes back to the database?", 
          "Commit", JOptionPane.YES_NO_OPTION);
      
      if(select == JOptionPane.YES_OPTION)
      {
        DatabaseDocument dbDoc =
          (DatabaseDocument)(
            (GFFStreamFeature)active_feature.getEmblFeature()).getDocumentEntry().getDocument();
        ChadoTransactionManager.commit(dbDoc, false,chadoTransactionManager);
      }
    }
    
    if(geneBuilderHash != null)
    {
      geneBuilderHash.remove(chado_gene.getGeneUniqueName());
    }
    
    if(reopen)
    {
      final GeneBuilderFrame gbFrame =
        new GeneBuilderFrame(active_feature, entry_group, 
                             selection, goto_event_source, 
                             chadoTransactionManager);
      if(geneBuilderHash != null)
      {
        gbFrame.addGeneBuilderHash(geneBuilderHash);
        geneBuilderHash.put(chado_gene.getGeneUniqueName(), gbFrame);
      }
    }
    else 
    {
      if(!chado_gene.getGene().isReadOnly() && GeneUtils.isBoundaryOK(chado_gene) > 0)
      {
        int result = JOptionPane.showConfirmDialog(this, 
          "Gene model boundary needs fixing.\nFix this now?", 
          "Gene Boundary", JOptionPane.YES_NO_OPTION);
        if(result == JOptionPane.YES_OPTION)
          GeneUtils.checkGeneBoundary(chado_gene);
      }
      
      if(!GeneUtils.isStrandOK(chado_gene))
      {
        JOptionPane.showMessageDialog(this, 
          "All gene features should be on the same strand.\n"+
          "Check the strand of each of the features (gene,\n"+
          "transcript, CDS, polypeptide).", 
          "Check Strand", JOptionPane.WARNING_MESSAGE);
      }
    }
    
    super.dispose();
  }
  
  /**
   * Add the hash from EditMenu. This records what gene builders are
   * open.
   * @param geneBuilderFocusHash
   */
  public void addGeneBuilderHash(final Hashtable geneBuilderHash)
  {
    this.geneBuilderHash = geneBuilderHash;
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
    final JCheckBox cb = new JCheckBox();

    List exons = chado_gene.getSpliceSitesOfTranscript(transcript_name, 
        DatabaseDocument.EXONMODEL);
    if(exons != null && exons.size() > 0)
      cb.setSelected(((GFFStreamFeature)exons.get(0)).isVisible());
    else
      cb.setSelected(transcript.isVisible());
    
    cb.setOpaque(false);
    cb.addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent e)
      { 
        final List exons = chado_gene.getSpliceSitesOfTranscript(transcript_name, 
            DatabaseDocument.EXONMODEL);
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
    setCursor(new Cursor(Cursor.WAIT_CURSOR));
    this.active_feature = active_feature;
    feature_editor.stopListening();
    feature_editor.setActiveFeature(active_feature, isSet);
    setTabTitle();
    setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
  }
  
  private void setTabTitle()
  {
    try
    {
      tabpane.setTitleAt(0, "Annotation :: "+
          (String)active_feature.getQualifierByName("ID").getValues().get(0));
    }
    catch(InvalidRelationException e){}
  }

  /**
   *  Remove this object as a feature and entry change listener.
   **/
  private void stopListening(final Feature feature) 
  {
    if(getEntry() != null)
      getEntry().removeEntryChangeListener(this);
    feature.removeFeatureChangeListener(this);
    
    if(selection != null && geneBuilderSelectionChangeListener != null)
      selection.removeSelectionChangeListener(geneBuilderSelectionChangeListener);
  }
  
  private void startListening(final Feature feature) 
  {
    if(getEntry() != null)
      getEntry().addEntryChangeListener(this);
    feature.addFeatureChangeListener(this);
    
    if(selection != null)
    {
      geneBuilderSelectionChangeListener = new GeneBuilderSelectionChangeListener();
      selection.addSelectionChangeListener(geneBuilderSelectionChangeListener);
    }
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
        
        Qualifier derives_from_qualifier = 
          qualifiers.getQualifierByName("Derives_from");
        
        if(parent_qualifier == null && derives_from_qualifier == null)
          return;

        tree.addNode(event.getFeature());
        feature.addFeatureChangeListener(this);
        
        // if polypeptide added then we are done
        if(derives_from_qualifier != null)
          return;
           
        final String parent = (String)parent_qualifier.getValues().get(0);

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
  
  protected FeatureEdit getFeatureEdit()
  {
    return feature_editor;
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
           (String)trans.getQualifierByName("ID").getValues().get(0), 
           DatabaseDocument.EXONMODEL);
       
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
          (String)trans.getQualifierByName("ID").getValues().get(0), 
          DatabaseDocument.EXONMODEL);
      
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

  class GeneBuilderSelectionChangeListener implements SelectionChangeListener
  {
    public void selectionChanged(SelectionChangeEvent event)
    {
      viewer.repaint();
      tree.setSelection(GeneBuilderFrame.this.selection);
    } 
  }
}
