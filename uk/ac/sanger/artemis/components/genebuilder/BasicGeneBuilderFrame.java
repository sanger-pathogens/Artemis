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
import java.awt.Toolkit;
import java.util.List;

import javax.swing.BoxLayout;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import uk.ac.sanger.artemis.EntryChangeEvent;
import uk.ac.sanger.artemis.EntryChangeListener;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.GotoEventSource;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.chado.ChadoTransactionManager;
import uk.ac.sanger.artemis.components.QualifierTextArea;
import uk.ac.sanger.artemis.components.Utilities;
import uk.ac.sanger.artemis.components.genebuilder.cv.CVPanel;
import uk.ac.sanger.artemis.components.genebuilder.gff.PropertiesPanel;
import uk.ac.sanger.artemis.components.genebuilder.ortholog.MatchPanel;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.DocumentEntry;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierInfo;
import uk.ac.sanger.artemis.io.QualifierLazyLoading;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.StreamQualifier;
import uk.ac.sanger.artemis.util.StringVector;


public class BasicGeneBuilderFrame extends JFrame
       implements EntryChangeListener, FeatureChangeListener
{
  private static final long serialVersionUID = 1L;
  private Feature feature;
  
  public BasicGeneBuilderFrame(Feature feature,
                          final EntryGroup entry_group,
                          final Selection selection,
                          final GotoEventSource goto_event_source)
  {
    this(feature, entry_group, selection, goto_event_source, null);
  }
  
  public BasicGeneBuilderFrame(Feature feature,
      final EntryGroup entry_group,
      final Selection selection,
      final GotoEventSource goto_event_source,
      final ChadoTransactionManager chadoTransactionManager)
  {
    super();
    
    this.feature = feature;
    final ChadoCanonicalGene chadoGene = 
      ((GFFStreamFeature)feature.getEmblFeature()).getChadoGene();
    
    // set title
    final String title = "Artemis Gene Builder: " + 
          ((Feature)chadoGene.getGene().getUserData()).getIDString() +
           (feature.isReadOnly() ?
          "  -  (read only)" : "");

    final List<uk.ac.sanger.artemis.io.Feature> transcripts =
      chadoGene.getTranscripts();
    
    JPanel mainPanel = (JPanel) getContentPane();
    final JTabbedPane tabPane = new JTabbedPane();
    mainPanel.add(tabPane);
    
    for(int i=0; i<transcripts.size(); i++)
    {
      Feature transcript = (Feature) transcripts.get(i).getUserData();
      tabPane.addTab(transcript.getIDString(), new JPanel(new BorderLayout()));
    }
    
    addComonentToTab(tabPane, chadoGene, transcripts, entry_group, null);
    tabPane.addChangeListener(new ChangeListener()
    {
      // This method is called whenever the selected tab changes
      public void stateChanged(ChangeEvent evt) 
      {
        if(((JComponent)tabPane.getSelectedComponent()).getComponentCount() < 1)
          addComonentToTab(tabPane, chadoGene, transcripts, entry_group, this);
      }
    });
    
    setTitle(title);
    pack();
    
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    int height = getPreferredSize().height;
    if(height > (screen.height*0.9))   
      setSize(getPreferredSize().width, (int)(screen.height*0.9));
    Utilities.centreFrame(this);
    
    setVisible(true);
  }
  
  private void addToPanel(JComponent section, JPanel container, String sectionName)
  {
    section.setBackground(Color.WHITE);
    GeneEditorPanel.addDarkSeparator(container);
    GeneEditorPanel.addOpenClosePanel(sectionName, section, container,
        MatchPanel.getDescription());
    container.add(section);
  }
  
  private void addComonentToTab(final JTabbedPane tabPane, 
                                final ChadoCanonicalGene chadoGene,
                                final List<uk.ac.sanger.artemis.io.Feature> transcripts,
                                final EntryGroup entry_group,
                                final ChangeListener changeListener)
  {
    int sel = tabPane.getSelectedIndex();
    
    Feature transcript = (Feature) transcripts.get(sel).getUserData();
    String transcriptName = GeneUtils.getUniqueName(transcript.getEmblFeature());
      
    JPanel panel = new JPanel();
    JScrollPane jsp = new JScrollPane(panel);
    panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
    panel.setBackground(Color.WHITE);

    PropertiesPanel propertiesPanel = new PropertiesPanel(transcript);
    addToPanel(propertiesPanel, panel, "Properties");

    Feature protein =
        (Feature) chadoGene.getProteinOfTranscript(transcriptName).getUserData();
      
    GeneUtils.addLazyQualifiers((GFFStreamFeature)protein.getEmblFeature());
  
    QualifierTextArea qualifier_text_area = new QualifierTextArea();
    CVPanel cvPanel = new CVPanel(protein);
    MatchPanel matchForm = new MatchPanel(protein, 
          (DocumentEntry)getFeature().getEmblFeature().getEntry());
      
    qualifier_text_area.setText(getQualifierString(entry_group, protein));
    addToPanel(qualifier_text_area, panel, "Core");

    addToPanel(cvPanel, panel, "Controlled Vocabulary");
    matchForm.updateFromFeature(protein);
    addToPanel(matchForm, panel, "Match");

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
      // strip out CV qualifiers
      //
      if( (CVPanel.isCvTag(this_qualifier)) ||
          (PropertiesPanel.isPropertiesTag(this_qualifier, feature)) ||
          (MatchPanel.isMatchTag(this_qualifier)) ||
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
    return feature;
  }

  public void entryChanged(EntryChangeEvent event)
  {
    // TODO Auto-generated method stub
    
  }

  public void featureChanged(FeatureChangeEvent event)
  {
    // TODO Auto-generated method stub
    
  }
  
}
