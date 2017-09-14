/* BasicGeneBuilderFrame
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2009  Genome Research Limited
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
package uk.ac.sanger.artemis.components.genebuilder;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.FontMetrics;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.Iterator;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.border.Border;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryChangeEvent;
import uk.ac.sanger.artemis.EntryChangeListener;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.chado.ChadoTransactionManager;
import uk.ac.sanger.artemis.components.KeyChoice;
import uk.ac.sanger.artemis.components.MessageDialog;
import uk.ac.sanger.artemis.components.QualifierTextArea;
import uk.ac.sanger.artemis.components.Utilities;
import uk.ac.sanger.artemis.components.genebuilder.cv.CVPanel;

import uk.ac.sanger.artemis.components.genebuilder.gff.BasicPropertiesPanel;
import uk.ac.sanger.artemis.components.genebuilder.gff.PropertiesPanel;
import uk.ac.sanger.artemis.components.genebuilder.ortholog.MatchPanel;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.DocumentEntry;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.LocationParseException;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierInfo;
import uk.ac.sanger.artemis.io.QualifierLazyLoading;
import uk.ac.sanger.artemis.io.QualifierParseException;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.RangeVector;
import uk.ac.sanger.artemis.io.StreamQualifier;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.util.StringVector;


public class BasicGeneBuilderFrame extends JFrame
       implements EntryChangeListener, FeatureChangeListener
{
  private static final long serialVersionUID = 1L;
  private Feature activeFeature;
  private ChadoCanonicalGene chadoGene;
  private ChadoTransactionManager chadoTransactionManager;
  private EntryGroup entry_group;
  private Selection selection;
  private CVPanel cvPanel;
  private MatchPanel matchForm;
  private BasicPropertiesPanel propertiesPanel;
  private ReferencesPanel refPanel;
  private QualifierTextArea qualifier_text_area;
  private JTabbedPane tabPane = new JTabbedPane();
  private JTextField locationText = new JTextField(60);
  private KeyChoice keyChoice;
  private int lastTabSelected = 0;
  private Border empty = new EmptyBorder(0,0,0,0);
  private JScrollPane jsp = new JScrollPane();
  private DatabaseDocument dbDoc;

  public BasicGeneBuilderFrame(Feature feature,
                          final EntryGroup entry_group,
                          final Selection selection)
  {
    this(feature, entry_group, selection, null);
  }
  
  public BasicGeneBuilderFrame(Feature feature,
      final EntryGroup entry_group,
      final Selection selection,
      final ChadoTransactionManager chadoTransactionManager)
  {
    super();
    
    this.activeFeature = feature;
    this.entry_group = entry_group;
    this.chadoTransactionManager = chadoTransactionManager;
    this.selection = selection;
    
    
    dbDoc =
      (DatabaseDocument)(
        (GFFStreamFeature)activeFeature.getEmblFeature()).getDocumentEntry().getDocument();
    
    chadoGene = 
      ((GFFStreamFeature)feature.getEmblFeature()).getChadoGene();
    
    // set title
    final String title = "Artemis Gene Builder: " + 
          ((Feature)chadoGene.getGene().getUserData()).getIDString() +
           (feature.isReadOnly() ?
          "  -  (read only)" : "");

    final List<uk.ac.sanger.artemis.io.Feature> transcripts =
      chadoGene.getTranscripts();
    
    JPanel mainPanel = (JPanel) getContentPane();
    mainPanel.setLayout(new BorderLayout());
    
    JTabbedPane tabMapViewer = new JTabbedPane();
    JPanel northPanel = new JPanel(new BorderLayout());
    JLabel statusLine = createStatusBar();

    final BasicGeneViewerPanel nucViewer = new BasicGeneViewerPanel(
        this, chadoGene, selection, 
        entry_group, this, statusLine);
    
    northPanel.add(nucViewer, BorderLayout.CENTER);
    northPanel.add(statusLine, BorderLayout.NORTH);
    
    FlowLayout flowLayoutZeroHVgap = new FlowLayout(FlowLayout.LEADING, 0, 0);
    keyChoice =
      new KeyChoice(activeFeature.getEntry().getEntryInformation(),
                    activeFeature.getKey());
    keyChoice.setLayout(flowLayoutZeroHVgap);
    
    JPanel keyAndLocation = new JPanel(new BorderLayout());
    keyAndLocation.add(keyChoice, BorderLayout.WEST);
    keyAndLocation.add(locationText, BorderLayout.CENTER);
    
    ButtonPanel buttonPanel = new ButtonPanel(this, entry_group);
    keyAndLocation.add(buttonPanel, BorderLayout.SOUTH);
    
    northPanel.add(keyAndLocation, BorderLayout.SOUTH);
    tabMapViewer.addTab("Gene Map", northPanel);
    
    //
    // protein map
    List<uk.ac.sanger.artemis.io.Feature> proteins = 
      ProteinMapPanel.getProteinsWithProteinMapElement(
        (GFFStreamFeature) getFeature().getEmblFeature());  
    final BasicProteinMapPanel proteinMap;
    if(proteins != null)
    {
      proteinMap = 
        new BasicProteinMapPanel((GFFStreamFeature) proteins.get(0), 
            chadoGene, selection, this);
      
      final JScrollPane jsp_ppviewer = new JScrollPane(proteinMap);
      jsp_ppviewer.getViewport().setBackground(Color.white);
      jsp_ppviewer.setPreferredSize( proteinMap.getPreferredSize() );
      tabMapViewer.addTab("Protein Map", jsp_ppviewer);
    }
    else
      proteinMap = null;
    
    mainPanel.add(tabMapViewer, BorderLayout.NORTH);
    mainPanel.add(tabPane, BorderLayout.CENTER);
    
    for(int i=0; i<transcripts.size(); i++)
    {
      GFFStreamFeature transcript = (GFFStreamFeature) transcripts.get(i);
      tabPane.insertTab(GeneUtils.getUniqueName(transcript), 
          null, new JPanel(new BorderLayout()), null, i);
    }

    tabPane.insertTab("+", 
        null, new JPanel(new BorderLayout()), "Add a transcript", transcripts.size());
    
    addComponentToTab(tabPane, chadoGene, entry_group, null);
    updateLocation();

    tabPane.addChangeListener(new ChangeListener()
    {
      // This method is called whenever the selected tab changes
      public void stateChanged(ChangeEvent evt) 
      {
        int sel = tabPane.getSelectedIndex();
        if(tabPane.getTitleAt(sel).equals("+"))
        {
          uk.ac.sanger.artemis.io.Feature transcript = transcripts.get(lastTabSelected);
          int status = JOptionPane.showConfirmDialog(BasicGeneBuilderFrame.this,
              "Make a duplicate transcript feature of "+
              GeneUtils.getUniqueName(transcript), 
              "Duplicate "+GeneUtils.getUniqueName(transcript), 
              JOptionPane.OK_CANCEL_OPTION);
          if(status != JOptionPane.OK_OPTION)
          {
            tabPane.setSelectedIndex(lastTabSelected);
            return;
          }
          
          String transcriptName = duplicateTranscript(transcript);
          tabPane.setTitleAt(sel, transcriptName);
        }
        
        if(((JComponent)tabPane.getSelectedComponent()).getComponentCount() < 1)
        {
          addComponentToTab(tabPane, chadoGene, entry_group, this);
        }
        
        setActiveFeature(getSelectedTranscriptFeature());
        nucViewer.repaint();
        if(proteinMap != null)
          proteinMap.repaint();
        lastTabSelected = sel;
      }
    });

    final FlowLayout flow_layout =
      new FlowLayout(FlowLayout.CENTER, 18, 1);
    final JPanel ok_cancel_update_panel = new JPanel(flow_layout);
    
    final JCheckBox oneView = new JCheckBox("Overview", true);
    oneView.addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent e)
      {
        try
        {
          stopListeningAll();
        }
        catch (InvalidRelationException e1){}
        dispose();
        new GeneBuilderFrame(getFeature(), entry_group,
            selection, null);
        System.setProperty("basic", "false");
      }
    });
    ok_cancel_update_panel.add(oneView);
    
    if (!getFeature().isReadOnly())
    {
      final JButton okButton = new JButton("OK");
      okButton.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          try
          {
            stopListeningAll();
          }
          catch (InvalidRelationException e1){}

          Key key = keyChoice.getSelectedItem(); 
          final Location location = getLoc();
          if(location == null)
            return;

          if(setQualifiers() && 
             setActiveFeatureKeyAndLocation(key, location)) 
          {
            stopListening();
            propertiesPanel.updateObsoleteSettings(); 
            dispose();
          }
        }
      });
      ok_cancel_update_panel.add(okButton);
      
      final JButton applyButton = new JButton("Apply");
      applyButton.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          try
          {
            stopListeningAll();
          }
          catch (InvalidRelationException e1){}
          Key key = keyChoice.getSelectedItem(); 
          final Location location = getLoc();

          setQualifiers();
          setActiveFeatureKeyAndLocation(key, location);
          try
          {
            addListeners(chadoGene);
          }
          catch (InvalidRelationException e1){}
        }
      });
      ok_cancel_update_panel.add(applyButton);
      
      final JButton cancelButton = new JButton("Cancel");
      cancelButton.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e) 
        {
          stopListening();
          dispose();
        }
      });
      ok_cancel_update_panel.add(cancelButton);
    }

    mainPanel.add(ok_cancel_update_panel, BorderLayout.SOUTH);

    tabPane.setBorder(empty);
    setTitle(title);
    pack();

    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    int height = getPreferredSize().height;
    if(height > (screen.height*0.9))   
    {
      setSize(getPreferredSize().width+jsp.getVerticalScrollBar().getWidth(), 
              (int)(screen.height*0.9));
    }
    Utilities.centreFrame(this);
    setVisible(true);
    try
    {
      addListeners(chadoGene);
    }
    catch (InvalidRelationException e)
    {
      e.printStackTrace();
    }
    
    lastTabSelected = tabPane.getSelectedIndex();
  }
  
  private Location getLoc()
  {
    try 
    {
      return new Location(locationText.getText());
    }
    catch(LocationParseException exception) 
    {
      final String error_string = exception.getMessage ();
      System.out.println(error_string);
      new MessageDialog(null,
          "Cannot apply changes because of location error: " +
                        error_string);
    }
    return null;
  }
  
  public void dispose()
  {
    dispose(false);
  }
  
  /**
   * Override to ensure the GeneBuilderFrame removes all listeners
   */
  protected void dispose(final boolean reopen)
  {
    try
    {
      stopListeningAll();
    }
    catch(InvalidRelationException e)
    {
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
        ChadoTransactionManager.commit(dbDoc, false,chadoTransactionManager);
      }
    }
    
    if(reopen)
    {
      new BasicGeneBuilderFrame((Feature) chadoGene.getGene().getUserData(), entry_group, 
                                selection, chadoTransactionManager);
    }
    else if(GeneUtils.isBoundaryOK(chadoGene) > 0)
    {
      int result = JOptionPane.showConfirmDialog(this, 
          "Gene model boundary needs fixing.\nFix this now?", 
          "Gene Boundary", JOptionPane.YES_NO_OPTION);
      if(result == JOptionPane.YES_OPTION)
        GeneUtils.checkGeneBoundary(chadoGene);
    }
    
    super.dispose();
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
    
    List<uk.ac.sanger.artemis.io.Feature> transcripts = chado_gene.getTranscripts();
    for(int i=0; i<transcripts.size(); i++)
    {
      uk.ac.sanger.artemis.io.Feature transcript = transcripts.get(i);
      Feature trans = (Feature)transcript.getUserData();
       
       if(trans == null)
         trans = new Feature(transcript);
       
       trans.addFeatureChangeListener(this);
       
       if(trans.getEntry() != null)
         trans.getEntry().addEntryChangeListener(this);
       
       List<uk.ac.sanger.artemis.io.Feature> exons = chado_gene.getSpliceSitesOfTranscript(
           (String)trans.getQualifierByName("ID").getValues().get(0), 
           DatabaseDocument.EXONMODEL);
       
       if(exons == null || exons.size() < 1)
         continue;
       
       if(exons.get(0) instanceof org.gmod.schema.sequence.Feature)
         return;
       
       for(int j=0; j<exons.size(); j++)
       {
         uk.ac.sanger.artemis.io.Feature embl_exon = exons.get(j);

         Feature exon = (Feature)embl_exon.getUserData();
         
         if(exon == null)
           exon = new Feature(embl_exon);
         exon.addFeatureChangeListener(this);
         
         if(exon.getEntry() != null)
           exon.getEntry().addEntryChangeListener(this);
       }
    }
  }
  
  private int getFontHeight()
  {
    final FontMetrics fm = getFontMetrics(Options.getOptions().getFont());
    return fm.getHeight();  
  }
  
  private JLabel createStatusBar()
  {
    JLabel statusLine = new JLabel();
    statusLine.setMinimumSize(new Dimension(100, getFontHeight()));
    statusLine.setPreferredSize(new Dimension(100, getFontHeight()));

    Border loweredbevel = BorderFactory.createLoweredBevelBorder();
    Border raisedbevel = BorderFactory.createRaisedBevelBorder();
    Border compound = BorderFactory.createCompoundBorder(raisedbevel,loweredbevel);
    statusLine.setBorder(compound);
    return statusLine;
  }
  
  private void addToPanel(JComponent section, JPanel container, 
		                  String sectionName, String description, boolean isOpen)
  {
    section.setBackground(Color.WHITE);
    GeneEditorPanel.addDarkSeparator(container);
    OpenSectionButton butt =
      GeneEditorPanel.addOpenClosePanel(sectionName, section, container,
    		description);
    butt.setOpen(isOpen);
    container.add(section);
  }
  
  public Feature getSelectedTranscriptFeature()
  {
    int sel = tabPane.getSelectedIndex();
   
    /*System.out.println("Selected transcript "+
        GeneUtils.getUniqueName(chadoGene.getTranscripts().get(sel))+" :: "+
        tabPane.getTitleAt(sel));*/
    if(sel > chadoGene.getTranscripts().size()-1)
      return null;
    return (Feature) chadoGene.getTranscripts().get(sel).getUserData();
  }
  
  public List<uk.ac.sanger.artemis.io.Feature> getCDSOfSelectedTranscriptFeature()
  {
	Feature transcript = getSelectedTranscriptFeature();
    String transcriptName = GeneUtils.getUniqueName(transcript.getEmblFeature());
	return chadoGene.getSpliceSitesOfTranscript(transcriptName, DatabaseDocument.EXONMODEL);
  }
  
  private void addComponentToTab(final JTabbedPane tabPane, 
                                 final ChadoCanonicalGene chadoGene,
                                 final EntryGroup entry_group,
                                 final ChangeListener changeListener)
  {
    Feature transcript = getSelectedTranscriptFeature();
    String transcriptName = GeneUtils.getUniqueName(transcript.getEmblFeature());
      
    JPanel panel = new JPanel();
    jsp.setViewportView(panel);
    jsp.setBorder(empty);
    panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
    panel.setBackground(Color.WHITE);

    propertiesPanel = new BasicPropertiesPanel(chadoGene, this);
    addToPanel(propertiesPanel, panel, "Properties",
        null,true);

    Feature protein =
        (Feature) chadoGene.getProteinOfTranscript(transcriptName).getUserData();
      
    GeneUtils.addLazyQualifiers((GFFStreamFeature)protein.getEmblFeature());

    //
    // literature
    refPanel = new ReferencesPanel(protein);
    addToPanel(refPanel, panel, "Literature/Dbxref", null,
        !refPanel.isEmpty());
    
    //
    // core text
    qualifier_text_area = new QualifierTextArea();
    cvPanel = new CVPanel(protein);
    matchForm = new MatchPanel(protein, 
          (DocumentEntry)getFeature().getEmblFeature().getEntry());
    
    qualifier_text_area.getDocument().addDocumentListener(
        new TextAreaDocumentListener(qualifier_text_area));
    qualifier_text_area.setText(getQualifierString(entry_group, protein));

    addToPanel(qualifier_text_area, panel, "Core", null,
        !qualifier_text_area.getText().equals(""));

    addToPanel(cvPanel, panel, "Controlled Vocabulary", 
        CVPanel.getDescription(), !cvPanel.isEmpty());
    matchForm.updateFromFeature(protein);
    addToPanel(matchForm, panel, "Match", 
        MatchPanel.getDescription(), !matchForm.isEmpty());

    tabPane.removeChangeListener(changeListener);
    ((JComponent)tabPane.getSelectedComponent()).add(jsp);
    tabPane.addChangeListener(changeListener);
  }

  /**
   *  Return a string containing one qualifier per line.  These are the
   *  original qualifiers, not the qualifiers from the qualifier_text_area.
   **/
  private String getQualifierString(EntryGroup entry_group,
                                    Feature feature) 
  {
    final StringBuffer buffer = new StringBuffer();
    final QualifierVector qualifiers = feature.getQualifiers();       
    
    for(int qualifier_index = 0; qualifier_index < qualifiers.size();
        ++qualifier_index) 
    {
      final Qualifier this_qualifier = (Qualifier)qualifiers.elementAt(qualifier_index);

      //
      // strip out CV / Properties / Match qualifiers
      if( (CVPanel.isCvTag(this_qualifier)) ||
          (PropertiesPanel.isPropertiesTag(this_qualifier, feature)) ||
          (MatchPanel.isMatchTag(this_qualifier)) ||
          (ReferencesPanel.isReferenceTag(this_qualifier)) ||
          (ProteinMapPanel.isProteinMapElement(this_qualifier)) )
        continue;
      
      if(this_qualifier instanceof QualifierLazyLoading)
        ((QualifierLazyLoading)this_qualifier).setForceLoad(true);

      final QualifierInfo qualifier_info =
        getFeature().getEntry().getEntryInformation().getQualifierInfo(this_qualifier.getName());

      final StringVector qualifier_strings =
                       StreamQualifier.toStringVector(qualifier_info, this_qualifier);

      for(int value_index = 0; value_index < qualifier_strings.size();
          ++value_index)
      {
        final String qualifier_string = (String)qualifier_strings.elementAt(value_index);
        buffer.append(qualifier_string + "\n");
      }
    }

    return buffer.toString();
  }
  
  public Feature getFeature()
  {
    return activeFeature;
  }
  
  
  protected void setObsoleteChanged(boolean obs, FeatureVector features)
  {
    propertiesPanel.setObsoleteChanged(obs, features);
  }
  
  private Entry getEntry()
  {
    return getFeature().getEntry();
  }
  
  /**
   *  Remove this object as a feature and entry change listener.
   **/
  public void stopListening() 
  {
    getEntry().removeEntryChangeListener(this);
    getFeature().removeFeatureChangeListener(this);
    if(cvPanel != null)
      getFeature().removeFeatureChangeListener(cvPanel);
    if(matchForm != null)
      getFeature().removeFeatureChangeListener(matchForm);
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
  
  private void stopListeningAll() throws InvalidRelationException
  {
    uk.ac.sanger.artemis.io.Feature embl_gene = 
      (uk.ac.sanger.artemis.io.Feature)chadoGene.getGene();
    Feature gene = (Feature)embl_gene.getUserData();   
    stopListening(gene);
    
    List<uk.ac.sanger.artemis.io.Feature> transcripts = chadoGene.getTranscripts();
    for(int i=0; i<transcripts.size(); i++)
    {
      uk.ac.sanger.artemis.io.Feature transcript = transcripts.get(i);
      Feature trans = (Feature)transcript.getUserData();
       
      stopListening(trans);
      
      List<uk.ac.sanger.artemis.io.Feature> exons = chadoGene.getSpliceSitesOfTranscript(
          (String)trans.getQualifierByName("ID").getValues().get(0), 
          DatabaseDocument.EXONMODEL);
      
      if(exons == null)
        continue;
      
      for(int j=0; j<exons.size(); j++)
      {
        uk.ac.sanger.artemis.io.Feature embl_exon = exons.get(j);
        Feature exon = (Feature)embl_exon.getUserData();
        stopListening(exon);      
      }
    }
  }
  
  protected void setActiveFeature(final Feature feature)
  {
    this.activeFeature = feature;
    updateLocation();
    updateKey();
  }
  
  /**
   *  Set the key, location and qualifiers of the feature to be the same as
   *  what values currently shown in the components.
   *  @return true if and only if action succeeds.  It may fail because of an
   *    illegal location or qualifier, in which case a message will be
   *    displayed before returning.
   **/
  private boolean setProteinFeature(uk.ac.sanger.artemis.io.Feature prot) 
  {
    if(prot == null)
      return true;
    
    Feature protein = (Feature) prot.getUserData();
    
    final Key key;
    final Location location;
    if(isActiveFeature(protein))
    {
      key = keyChoice.getSelectedItem();
      location = getLoc();
    }
    else
    {
      key = protein.getKey();
      location = protein.getLocation();
    }
    
    final QualifierVector qualifiers;

    try 
    {
      qualifiers =
        qualifier_text_area.getParsedQualifiers(getEntry().getEntryInformation ());
 
      QualifierVector oldQualifiers = protein.getQualifiers();
      
      // preserve protein properties
      qualifiers.addAll(propertiesPanel.getProteinProperties(protein));
      
      for(int i=0; i<oldQualifiers.size(); i++)
      {
        Qualifier qualifier = (Qualifier) oldQualifiers.get(i);
        if(BasicProteinMapPanel.isProteinMapElement(qualifier))
        {
          qualifiers.addQualifierValues(qualifier);
        }
      }

      // if using controlled vocab form
      if(cvPanel != null)
      {
        QualifierVector cvQualifiers = cvPanel.getCvQualifiers();
        if(cvQualifiers != null && cvQualifiers.size() > 0)
          qualifiers.addAll(cvQualifiers);
      }
      
      if(matchForm != null)
      {
        QualifierVector orthologQualifiers = matchForm.getMatchQualifiers();
        if(orthologQualifiers != null && orthologQualifiers.size() > 0)
          qualifiers.addAll(orthologQualifiers);
      }
      
      if(refPanel != null)
      {
        QualifierVector referenceQualifiers = refPanel.getQualifiers();
        if(referenceQualifiers != null && referenceQualifiers.size() > 0)
          qualifiers.addAll(referenceQualifiers);
      }
    }
    catch(QualifierParseException exception) 
    {
      final String error_string = exception.getMessage();
      System.out.println(error_string);
      new MessageDialog(this, "Cannot apply changes because of a qualifier " +
                        "error: " + error_string);
      return false;
    }

    try 
    {
      entry_group.getActionController().startAction();

      try 
      {
        protein.set(key, location, qualifiers);
      }
      catch(java.lang.Error err)
      {
        err.printStackTrace();

        if(err.getMessage().indexOf("InvalidRelationException")>-1)
        {
          JScrollPane jsp = new JScrollPane(new JLabel(err.getMessage()));
          jsp.setPreferredSize(new Dimension(200,100));
          JOptionPane.showMessageDialog(null, jsp, 
                      "Error", JOptionPane.ERROR_MESSAGE);
        }
        return false;
      }
    } 
    catch(EntryInformationException e) 
    {
      final String error_string = e.getMessage();
      new MessageDialog(this, "Cannot apply changes: " + error_string);

      return false;
    } 
    catch(OutOfRangeException e) 
    {
      new MessageDialog(this, "Cannot apply changes - the location is out of " +
                        "range for this sequence");
      return false;
    } 
    catch(ReadOnlyException e) 
    {
      new MessageDialog(this, "Cannot apply changes - the feature is " +
                        "read only");
      return false;
    } 
    finally
    {
      entry_group.getActionController ().endAction ();
    }

    dribble();
    return true;
  }
  
  private boolean setFeaturePropertiesByFeature(Feature f, QualifierVector qualifiers) 
  {
    final Key key;
    final Location location;
    if(isActiveFeature(f))
    {
      key = keyChoice.getSelectedItem();
      location = getLoc();
    }
    else
    {
      key = f.getKey();
      location = f.getLocation();
    }
    
    QualifierVector oldQualifiers = f.getQualifiers();
    
    for(int i=0; i<oldQualifiers.size(); i++)
    {
      Qualifier qualifier = (Qualifier) oldQualifiers.get(i);
      if(!PropertiesPanel.isPropertiesTag(qualifier, f))
        qualifiers.addQualifierValues(qualifier);
    }
    
    try 
    {
      entry_group.getActionController().startAction();

      try 
      {
        f.set(key, location, qualifiers);
      }
      catch(java.lang.Error err)
      {
        err.printStackTrace();

        if(err.getMessage().indexOf("InvalidRelationException")>-1)
        {
          JScrollPane jsp = new JScrollPane(new JLabel(err.getMessage()));
          jsp.setPreferredSize(new Dimension(200,100));
          JOptionPane.showMessageDialog(null, jsp, 
                      "Error", JOptionPane.ERROR_MESSAGE);
        }
        return false;
      }
    } 
    catch(EntryInformationException e) 
    {
      final String error_string = e.getMessage();
      new MessageDialog(this, "Cannot apply changes: " + error_string);

      return false;
    } 
    catch(OutOfRangeException e) 
    {
      new MessageDialog(this, "Cannot apply changes - the location is out of " +
                        "range for this sequence");
      return false;
    } 
    catch(ReadOnlyException e) 
    {
      new MessageDialog(this, "Cannot apply changes - the feature is " +
                        "read only");
      return false;
    } 
    finally
    {
      entry_group.getActionController ().endAction ();
    }

    return true;
  }
  
  private boolean isActiveFeature(Feature f)
  {
    return activeFeature.getIDString().equals(f.getIDString());
  }
  
  /**
   * Duplicate a given transcript feature and its children.
   * @param transcripts
   * @return
   */
  private String duplicateTranscript(uk.ac.sanger.artemis.io.Feature lastTranscript)
  {
    Feature transcript =
      GeneViewerPanel.createTranscript(chadoGene, entry_group, 
          lastTranscript.getLocation());
    
    String transcriptName = GeneUtils.getUniqueName(transcript.getEmblFeature());
    Iterator<uk.ac.sanger.artemis.io.Feature> it = 
      chadoGene.getChildren(lastTranscript).iterator();
    
    while(it.hasNext())
    {
      uk.ac.sanger.artemis.io.Feature f = it.next();
      if(f.getKey().getKeyString().equals(DatabaseDocument.EXONMODEL) ||
         f.getKey().getKeyString().equals("pseudogenic_exon") )
      {
        GFFStreamFeature gff_exon = null;
        RangeVector ranges = f.getLocation().getRanges();
        
        for(int i=0; i<ranges.size(); i++)
          gff_exon = GeneViewerPanel.addExonFeature(chadoGene, entry_group, gff_exon, 
              (Range)ranges.get(i), transcriptName, selection, f.getKey(), null);
      }
      else if(!f.getKey().equals("polypeptide"))
      {
        final String name;
        if(f.getKey().getKeyString().equals("five_prime_UTR"))
          name = "5UTR";
        else if(f.getKey().getKeyString().equals("three_prime_UTR"))
          name = "3UTR";
        else
          name = f.getKey().getKeyString();
        try
        {
          GeneViewerPanel.addFeature(
              f.getLocation().getTotalRange(), transcriptName, null,
              transcript.getLocation().isComplement(), 
              true, chadoGene, entry_group, f.getKey(), name);
        }
        catch (OutOfRangeException e)
        {
          e.printStackTrace();
        }
      }
    }
    
    GeneViewerPanel.addProteinFeature(chadoGene, entry_group, transcriptName, transcript);
    return transcriptName;
  }
  
  /**
   *  On Unix machines this method will append the text of the feature to a
   *  file in a current directory called .dribble + <the entry name>
   **/
  private void dribble()
  {
    if(!Options.isUnixHost()) 
      return;

    final String dribble_file_name;

    if(getEntry().getName() != null) 
      dribble_file_name = ".dribble." + getEntry().getName();
    else 
      dribble_file_name = ".dribble.no_name";

    try  
    {
      final Writer writer = new FileWriter(dribble_file_name, true);
      getFeature().writeNative(writer);
      writer.flush();
      writer.close();
    } 
    catch(IOException e) 
    {
      System.err.println("IO exception while accessing " + dribble_file_name +
                         ": " + e.getMessage());
    }
  }
  
  private void updateKey()
  {
    keyChoice.setKey(getFeature().getKey());
  }
  
  private void updateLocation()
  {
    locationText.setText(getFeature().getLocation().toStringShort());
  }
  
  private boolean setActiveFeatureKeyAndLocation(Key key, Location location)
  {
    try
    {
      activeFeature.set(key, location, activeFeature.getQualifiers());
    }
    catch (ReadOnlyException e)
    {
      JOptionPane.showMessageDialog(this, "Key",
          "Read Only", JOptionPane.WARNING_MESSAGE);
      return false;
    }
    catch (EntryInformationException e)
    {
      final String error_string = e.getMessage ();
      JOptionPane.showMessageDialog(this, "Key Error",
          "Cannot apply changes because of key error: " +
           error_string, JOptionPane.WARNING_MESSAGE);
      return false;
    }
    catch (OutOfRangeException e)
    {
      final String error_string = e.getMessage ();
      JOptionPane.showMessageDialog(this, "Location Error",
          "Cannot apply changes because of location error: " +
           error_string, JOptionPane.WARNING_MESSAGE);
      return false;
    }
    return true;
  }
  
  protected boolean setQualifiers()
  {
    Feature gene = (Feature) chadoGene.getGene().getUserData();
    Feature transcript = getSelectedTranscriptFeature();
    
    QualifierVector geneProperties =
      propertiesPanel.getGeneProperties(gene);
    
    if(!setFeaturePropertiesByFeature(gene, geneProperties))
      return false;
    
    QualifierVector transcriptProperties =
      propertiesPanel.getTranscriptProperties(transcript);
    
    String transcriptName = GeneUtils.getUniqueName(transcript.getEmblFeature());
    List<uk.ac.sanger.artemis.io.Feature> exons =
      chadoGene.getSpliceSitesOfTranscript(transcriptName, DatabaseDocument.EXONMODEL);
    
    if(exons != null)
    {
      Feature exon = (Feature) exons.get(0).getUserData();
      QualifierVector exonProperties =
        propertiesPanel.getExonProperties(exon);
      if(!setFeaturePropertiesByFeature(exon, exonProperties))
        return false;
    }

    uk.ac.sanger.artemis.io.Feature protein =
       chadoGene.getProteinOfTranscript(transcriptName);
    
    return (setFeaturePropertiesByFeature(transcript, transcriptProperties) &&
            setProteinFeature(protein));
  }
  
  protected MatchPanel getMatchForm()
  {
    return matchForm;
  }

  protected JTextField getLocationText()
  {
    return locationText;
  }
  
  protected Selection getSelection()
  {
    return selection;
  }
  
  public void entryChanged(EntryChangeEvent event)
  {
    
  }
  
  public void featureChanged(FeatureChangeEvent event)
  {
    updateLocation();
    repaint();
  }
  
}
