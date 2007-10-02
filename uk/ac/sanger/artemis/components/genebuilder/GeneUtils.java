/* GeneUtils.java
 *
 * This file is part of Artemis
 * Copyright (C) 2007  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
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
 **/
package uk.ac.sanger.artemis.components.genebuilder;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import uk.ac.sanger.artemis.components.EditMenu;
import uk.ac.sanger.artemis.components.MessageDialog;
import uk.ac.sanger.artemis.components.SelectionMenu;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.Feature;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.sequence.MarkerRange;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.EntryVector;
import uk.ac.sanger.artemis.FeatureKeyQualifierPredicate;
import uk.ac.sanger.artemis.FeaturePredicate;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.GotoEventSource;
import uk.ac.sanger.artemis.Selection;


public class GeneUtils
{
  private static final long serialVersionUID = 1L;
  private static Vector showFeatures = new Vector();
  private static Vector hideFeatures = new Vector();
  
  static 
  {
    showFeatures.add("gene");
    showFeatures.add("pseudogene");
    showFeatures.add("exon-model");
    showFeatures.add("pseudogenic_exon");
  }
  
  static
  {
    hideFeatures.add("polypeptide");
    hideFeatures.add("mRNA");
    hideFeatures.add("pseudogenic_transcript");
  }
  
  /**
   * Given a collection of features, determine if these should be
   * shown or hidden in the Artemis display
   * @param features
   */
  public static void defineShowHideGeneFeatures(final FeatureVector features)
  {
    final DefaultListModel showListModel = new DefaultListModel();
    for(int i=0; i<showFeatures.size(); i++)
      showListModel.addElement(showFeatures.get(i));
    final JList displayList = new JList(showListModel);
    
    
    final DefaultListModel hideListModel = new DefaultListModel();
    for(int i=0; i<hideFeatures.size(); i++)
      hideListModel.addElement(hideFeatures.get(i));
    final JList hideList = new JList(hideListModel);
    
    
    final JButton hide_butt = new JButton("HIDE");
    hide_butt.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        while(!displayList.isSelectionEmpty())
        {
          final String hideKey = (String)displayList.getSelectedValue();
          hideListModel.addElement(hideKey);
          showListModel.removeElement(hideKey);
          
          hideFeatures.add(hideKey);
          if(showFeatures.contains(hideKey))
            showFeatures.remove(hideKey);
        }
      }
    });

    Box bdown = Box.createVerticalBox();
    bdown.add(new JLabel("Features Displayed:"));
    bdown.add(new JScrollPane(displayList));
    bdown.add(hide_butt);
    
    final JPanel hideShowPanel = new JPanel(new BorderLayout());
    hideShowPanel.add(bdown, BorderLayout.CENTER);

    
    final JButton show_butt = new JButton("SHOW");
    show_butt.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        while(!hideList.isSelectionEmpty())
        {
          final String showKey = (String)hideList.getSelectedValue();
          showListModel.addElement(showKey);
          hideListModel.removeElement(showKey);
          
          if(hideFeatures.contains(showKey))
            hideFeatures.remove(showKey);
          showFeatures.add(showKey);
        }
      }
    });

    bdown = Box.createVerticalBox();
    bdown.add(Box.createVerticalGlue());
    bdown.add(new JLabel("Features Hidden:"));
    bdown.add(new JScrollPane(hideList));
    bdown.add(show_butt);
    hideShowPanel.add(bdown, BorderLayout.EAST);

    int select = JOptionPane.showConfirmDialog(null, hideShowPanel,
                            "Gene Model Features Displayed...",
                             JOptionPane.OK_CANCEL_OPTION,
                             JOptionPane.QUESTION_MESSAGE);

    if(select == JOptionPane.CANCEL_OPTION)
      return;
    
    showHideGeneFeatures(features);
  }
  
  /**
   * Based on the hidenFeatures and showFeatures set the GFFStreamFeatures
   * visibility
   * @param features
   */
  public static void showHideGeneFeatures(final FeatureVector features)
  {
    for(int i=0; i<features.size(); i++)
    {
      final Feature feature = features.elementAt(i).getEmblFeature();
      
      if(feature instanceof GFFStreamFeature)
      {
        final String key = feature.getKey().getKeyString();
        if(hideFeatures.contains(key))
          ((GFFStreamFeature)feature).setVisible(false);
        else if(showFeatures.contains(key))
          ((GFFStreamFeature)feature).setVisible(true);
      }
    }
  }
  
  /**
   * Determine based on the feature given if a feature is hidden
   * @param key
   * @return
   */
  public static boolean isHiddenFeature(final String key)
  {
    return hideFeatures.contains(key);
  }
  
  /**
   * 
   * @param frame
   * @param selection
   * @param entry_group
   * @param goto_event_source
   */
  public static void createGeneModel(final JFrame frame,
      final Selection selection,
      final EntryGroup entry_group,
      final GotoEventSource goto_event_source)
  {
    if(!SelectionMenu.checkForSelectionRange(frame, selection))
      return;
    final MarkerRange range = selection.getMarkerRange ();
    final Entry default_entry = entry_group.getDefaultEntry ();

    if (default_entry == null) 
    {
      new MessageDialog (frame, "There is no default entry");
      return;
    }
    
    final QualifierVector qualifiers = new QualifierVector();
    final String uniquename = promptForUniquename(entry_group, 
                                   range.isForwardMarker());
    final Qualifier qualifier = new Qualifier("ID", uniquename);
    qualifiers.add(qualifier);
    
    try
    {
      final Location new_location = range.createLocation ();
      final Key key = new Key("gene");
      final uk.ac.sanger.artemis.Feature geneFeature = 
          default_entry.createFeature(key, new_location, qualifiers);

      final ChadoCanonicalGene chadoGene = new ChadoCanonicalGene();
      chadoGene.setGene(geneFeature.getEmblFeature());
      ((uk.ac.sanger.artemis.io.GFFStreamFeature) 
          (geneFeature.getEmblFeature())).setChadoGene(chadoGene);
      
      // create transcript
      uk.ac.sanger.artemis.Feature transcript = 
        GeneViewerPanel.createTranscript(chadoGene, entry_group);
      ((uk.ac.sanger.artemis.io.GFFStreamFeature)
          (transcript.getEmblFeature())).setChadoGene(chadoGene);
      final String transcriptId = 
         (String)transcript.getQualifierByName("ID").getValues().get(0);
      
      // add exon
      GeneViewerPanel.addExonFeature(chadoGene, entry_group, 
          null, range.getRange(), transcriptId, selection, 
          new Key(DatabaseDocument.EXONMODEL), null);
      
      // add protein
      uk.ac.sanger.artemis.Feature polypep =
        GeneViewerPanel.addProteinFeature(chadoGene, entry_group, transcriptId, transcript);
      
      selection.clear();
      selection.add(polypep);
      
      EditMenu.editSelectedFeatures(entry_group, selection, 
          goto_event_source, polypep, null, null);
    }
    catch(ReadOnlyException e)
    {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    catch(EntryInformationException e)
    {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    catch(OutOfRangeException e)
    {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
  
  /**
   * Prompt the user for an ID
   * @return
   */
  public static String promptForUniquename(final EntryGroup entry_group,
                                            final boolean is_forward)
  {
    final Entry default_entry = entry_group.getDefaultEntry ();
    String id = null;
    
    if(default_entry.getEMBLEntry() instanceof 
        uk.ac.sanger.artemis.io.DatabaseDocumentEntry)
    {  
      while(id == null ||
            id.equals("") ||
            id.equals("to_be_set"))
      {
        String msg = "Provide a unique ID ";
        
        if(!is_forward)
          msg = msg + "for reverse strand : ";
        else
          msg = msg + ": ";
        
        id = JOptionPane.showInputDialog(null,
                           msg,
                           "ID missing ",
                           JOptionPane.QUESTION_MESSAGE).trim();
        
        if(!isUniqueID(entry_group, id))
        {
          JOptionPane.showMessageDialog(null, 
              "ID "+id+" not unique.\nEnter a unique ID.", 
              "ID Not Unique", 
              JOptionPane.WARNING_MESSAGE);
          id = null;
        }
      }
    }
    return id;
  }
  
  /**
   * Test to ensure ID (chado uniquename) is unique.
   * @param entry_group
   * @param id
   * @return
   */
  private static boolean isUniqueID(final EntryGroup entry_group,
                                    final String id)
  {
    final FeaturePredicate predicate =
      new FeatureKeyQualifierPredicate(null, "ID", id, 
                                       false, true);
    final FeatureVector features = entry_group.getAllFeatures();
    for(int i=0; i<features.size(); i++)
    {
      uk.ac.sanger.artemis.Feature feature = features.elementAt(i);
      if(predicate.testPredicate(feature))
        return false;
      
    }
    return true;
  }
  
  /**
   * Given an group of entries determine if they contain a database entry
   * @param entryGroup
   * @return
   */
  public static boolean isDatabaseEntry(final EntryGroup entryGroup)
  {
    final EntryVector entries = entryGroup.getActiveEntries();
    
    for(int i=0; i<entries.size(); i++)
    {
      if( entries.elementAt(i).getEMBLEntry() instanceof DatabaseDocumentEntry )
        return true;
    }
    return false;
  }
  
  public static void main(String args[])
  {
    GeneUtils.defineShowHideGeneFeatures(new FeatureVector());
  }
}
