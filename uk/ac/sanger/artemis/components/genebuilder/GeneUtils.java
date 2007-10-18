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
import java.util.Iterator;
import java.util.List;
import java.util.Set;
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
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.RangeVector;
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
   * Used to reverse complement all the gene model features
   * @param chadoGene
   */
  public static void complementGeneModel(final ChadoCanonicalGene chadoGene)
  {
    if(chadoGene == null)
      return;
    try
    {
      final Feature gene = chadoGene.getGene();
      final boolean complement = gene.getLocation().isComplement();
      gene.setLocation(gene.getLocation().getComplement());
      final Set kids = chadoGene.getChildren(gene);
      final Iterator it = kids.iterator();
      while(it.hasNext())
      {
        final Feature f = (Feature)it.next();
        final RangeVector rv = f.getLocation().getRanges();
        rv.reverse();
        f.setLocation(new Location(rv, !complement));
      }
    }
    catch(ReadOnlyException e)
    {
      e.printStackTrace();
    }
    catch(OutOfRangeException e)
    {
      e.printStackTrace();
    }
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
   */
  public static void duplicateGeneModel(final JFrame frame,
      final FeatureVector features_to_duplicate,
      final EntryGroup entry_group)
  {
    final Vector duplicatedGenes = new Vector();
    for (int i = 0 ; i < features_to_duplicate.size () ; ++i) 
    {
      final uk.ac.sanger.artemis.Feature this_feature = 
                           features_to_duplicate.elementAt(i);
      final GFFStreamFeature gffFeature = (GFFStreamFeature)this_feature.getEmblFeature();
      if(duplicatedGenes.contains(gffFeature.getChadoGene()))
        continue;
      
      duplicatedGenes.add(gffFeature.getChadoGene());
      
      try 
      {
        GFFStreamFeature gene = (GFFStreamFeature)gffFeature.getChadoGene().getGene();
        uk.ac.sanger.artemis.Feature newGeneFeature = ((uk.ac.sanger.artemis.Feature)
            gene.getUserData()).duplicate (true);
        
        final ChadoCanonicalGene chadoGene = gffFeature.getChadoGene();
        final ChadoCanonicalGene newchadoGene = new ChadoCanonicalGene();
        ((GFFStreamFeature)newGeneFeature.getEmblFeature()).setChadoGene(newchadoGene);
        newchadoGene.setGene(newGeneFeature.getEmblFeature());
        
        final List transcripts = chadoGene.getTranscripts();
        for(int j=0; j<transcripts.size(); j++)
        {
          final GFFStreamFeature transcript = (GFFStreamFeature)transcripts.get(j);
          final String transcriptName =
            (String)transcript.getQualifierByName("ID").getValues().get(0);
          
          uk.ac.sanger.artemis.Feature newTranscriptFeature = 
            duplicateFeature(transcript, newchadoGene);
          newchadoGene.addTranscript(newTranscriptFeature.getEmblFeature());
          final String newTranscriptName =
            (String)newTranscriptFeature.getQualifierByName("ID").getValues().get(0);
          
          List newFeatures;
          
          newFeatures= duplicateFeatures(chadoGene.get3UtrOfTranscript(transcriptName), newchadoGene);
          for(int k=0; k<newFeatures.size(); k++)
            newchadoGene.add3PrimeUtr(newTranscriptName, (Feature)newFeatures.get(k));
          
          newFeatures = duplicateFeatures(chadoGene.get5UtrOfTranscript(transcriptName), newchadoGene);
          for(int k=0; k<newFeatures.size(); k++)
            newchadoGene.add5PrimeUtr(newTranscriptName, (Feature)newFeatures.get(k));
          
          newFeatures = duplicateFeatures(chadoGene.getOtherFeaturesOfTranscript(transcriptName), newchadoGene);
          for(int k=0; k<newFeatures.size(); k++)
            newchadoGene.addOtherFeatures(newTranscriptName, (Feature)newFeatures.get(k));

          newFeatures = duplicateFeatures(chadoGene.getSplicedFeaturesOfTranscript(transcriptName), newchadoGene);
          for(int k=0; k<newFeatures.size(); k++)
          {
            uk.ac.sanger.artemis.Feature splicedFeature = 
              (uk.ac.sanger.artemis.Feature)newFeatures.get(k);
            newchadoGene.addSplicedFeatures(newTranscriptName, splicedFeature.getEmblFeature());
          }
          
          uk.ac.sanger.artemis.Feature newProtein = 
            duplicateFeature(chadoGene.getProteinOfTranscript(transcriptName), newchadoGene);
          if(newProtein != null)
            newchadoGene.addProtein(newTranscriptName, newProtein.getEmblFeature());
          
        }
      } 
      catch (ReadOnlyException e) {}
      catch(InvalidRelationException e) {}
    }
    
    duplicatedGenes.clear();
  }
  
  
  private static List duplicateFeatures(final List featuresOfTranscript,
                                 final ChadoCanonicalGene chadoGene) 
          throws ReadOnlyException
  {
    final List newFeatures = new Vector();
    
    if(featuresOfTranscript == null)
      return newFeatures;
    
    for(int i=0; i<featuresOfTranscript.size(); i++)
      newFeatures.add(duplicateFeature(
          (GFFStreamFeature)featuresOfTranscript.get(i), chadoGene));
  
    return newFeatures;
  }
  
  private static uk.ac.sanger.artemis.Feature duplicateFeature(
          final Feature feature, final ChadoCanonicalGene chadoGene) 
          throws ReadOnlyException
  {
    if(feature == null)
      return null;
    uk.ac.sanger.artemis.Feature newFeature = 
      ((uk.ac.sanger.artemis.Feature)feature.getUserData()).duplicate(true);
    ((GFFStreamFeature)newFeature.getEmblFeature()).setChadoGene(chadoGene);
    if(isHiddenFeature(newFeature.getKey().getKeyString()))
      ((GFFStreamFeature)newFeature.getEmblFeature()).setVisible(false);
    /*
    try
    {
      final QualifierVector qv = newFeature.getQualifiers().copy();
      
      for(int i=0; i<qv.size(); i++)
      {
        final Qualifier qualifier = (Qualifier)qv.elementAt(i);
        if(!qualifier.getName().equals("ID") &&
           !qualifier.getName().equals("Parent") &&
           !qualifier.getName().equals("Derives_from") &&
           ChadoTransactionManager.isSpecialTag(qualifier.getName()))
          newFeature.getQualifiers().removeQualifierByName(qualifier.getName());
      }

      newFeature.set(newFeature.getKey(), newFeature.getLocation(), qv);
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
    */
    return newFeature;
  }
  
  /**
   * Create gene model from base selection
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
      final FeatureVector newFeatures = new FeatureVector();
      final Location new_location = range.createLocation ();
      final Key key = new Key("gene");
      final uk.ac.sanger.artemis.Feature geneFeature = 
          default_entry.createFeature(key, new_location, qualifiers);
      newFeatures.add(geneFeature);
      
      final ChadoCanonicalGene chadoGene = new ChadoCanonicalGene();
      chadoGene.setGene(geneFeature.getEmblFeature());
      ((uk.ac.sanger.artemis.io.GFFStreamFeature) 
          (geneFeature.getEmblFeature())).setChadoGene(chadoGene);
      
      // create transcript
      uk.ac.sanger.artemis.Feature transcript = 
        GeneViewerPanel.createTranscript(chadoGene, entry_group);
      newFeatures.add(transcript);
      ((uk.ac.sanger.artemis.io.GFFStreamFeature)
          (transcript.getEmblFeature())).setChadoGene(chadoGene);
      final String transcriptId = 
         (String)transcript.getQualifierByName("ID").getValues().get(0);
      
      // add exon
      GeneViewerPanel.addExonFeature(chadoGene, entry_group, 
          null, new_location.getTotalRange(), transcriptId, selection, 
          new Key(DatabaseDocument.EXONMODEL), null);
      
      // add protein
      uk.ac.sanger.artemis.Feature polypep =
        GeneViewerPanel.addProteinFeature(chadoGene, entry_group, transcriptId, transcript);
      newFeatures.add(polypep);
      
      showHideGeneFeatures(newFeatures);
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
  

  private static void deleteFeature(uk.ac.sanger.artemis.Feature feature)
      throws ReadOnlyException
  {
    if(feature != null && feature.getEntry() != null)
      feature.removeFromEntry();
  }

  /**
   * Delete feature and children in a chado gene model
   * @param feature
   * @param chado_gene
   * @throws ReadOnlyException
   */
  public static void deleteAllFeature(uk.ac.sanger.artemis.Feature feature,
      final ChadoCanonicalGene chado_gene) throws ReadOnlyException
  {
    Set children = chado_gene.getChildren(feature.getEmblFeature());
    deleteFeature(feature);
    chado_gene.deleteFeature(feature.getEmblFeature());

    Feature embl_feature;
    Iterator it = children.iterator();

    while(it.hasNext())
    {
      embl_feature = (Feature) it.next();
      deleteFeature((uk.ac.sanger.artemis.Feature) embl_feature.getUserData());
      chado_gene.deleteFeature(embl_feature);
    }
  }
  
  /**
   * Adjust transcript and gene boundaries
   * @param chado_gene
   */
  public static void checkGeneBoundary(final ChadoCanonicalGene chado_gene)
  {
    final List transcripts = chado_gene.getTranscripts();
    int gene_start = Integer.MAX_VALUE;
    int gene_end = -1;
    
    Range range;
    for(int i=0; i<transcripts.size(); i++)
    {
      final Feature transcript = (Feature)transcripts.get(i);
      range = checkTranscriptBoundary(
          (uk.ac.sanger.artemis.Feature)transcript.getUserData(), chado_gene);
      if(range != null && range.getStart() < gene_start)
        gene_start = range.getStart();
      if(range != null && range.getEnd() > gene_end)
        gene_end = range.getEnd();
    }
    
    if(gene_end == -1 && gene_start == Integer.MAX_VALUE)
      return;
    
    setLocation(chado_gene.getGene(), gene_start, gene_end);
  }
  
  /**
   * Check and adjust transcript boundary
   * @param transcript
   * @param chado_gene
   */
  public static Range checkTranscriptBoundary(
      final uk.ac.sanger.artemis.Feature transcript,
      final ChadoCanonicalGene chado_gene)
  {
    final List transcripts = chado_gene.getTranscripts();

    if(transcripts.contains(transcript.getEmblFeature()))
    {
      checkProteinBoundary(transcript.getEmblFeature(), chado_gene);
      
      final Set children = chado_gene.getChildren(transcript.getEmblFeature());
      int transcript_start = Integer.MAX_VALUE;
      int transcript_end = -1;

      final Iterator it = children.iterator();
      while(it.hasNext())
      {
        final Feature feature = (Feature) it.next();
        final Range range = feature.getLocation().getTotalRange();
        if(range.getStart() < transcript_start)
          transcript_start = range.getStart();
        if(range.getEnd() > transcript_end)
          transcript_end = range.getEnd();
      }

      if(transcript_start == Integer.MAX_VALUE ||
         transcript_end == -1)
        return null;
      
      return setLocation(transcript.getEmblFeature(), 
                   transcript_start, transcript_end);
    }
    else
      JOptionPane.showMessageDialog(null,
          "Select a single transcript and try again.", "Transcript Selection",
          JOptionPane.ERROR_MESSAGE);
    return null;
  }
  
  public static void checkProteinBoundary(final Feature transcript,
                                          final ChadoCanonicalGene chado_gene)
  {
    final String transcriptName = getUniqueName(transcript);
    final Feature protein = chado_gene.getProteinOfTranscript(transcriptName);
    if(protein == null)
      return;
    
    int pp_start = Integer.MAX_VALUE;
    int pp_end = -1;
    
    final List dnaFeatures = new Vector();
    if(chado_gene.get3UtrOfTranscript(transcriptName) != null)
      dnaFeatures.addAll(chado_gene.get3UtrOfTranscript(transcriptName));
    if(chado_gene.get5UtrOfTranscript(transcriptName) != null)
      dnaFeatures.addAll(chado_gene.get5UtrOfTranscript(transcriptName)); 
    
    List exons = chado_gene.getSpliceSitesOfTranscript(transcriptName, DatabaseDocument.EXONMODEL);
    if(exons != null)
      dnaFeatures.addAll(exons);
    
    for(int i=0; i<dnaFeatures.size(); i++)
    {
      Feature dnaFeature = (Feature)dnaFeatures.get(i);
      final Range range = dnaFeature.getLocation().getTotalRange();
      if(range.getStart() < pp_start)
        pp_start = range.getStart();
      if(range.getEnd() > pp_end)
        pp_end = range.getEnd();
    }
    
    if(pp_start == Integer.MAX_VALUE || pp_end == -1)
       return;
    setLocation(protein, pp_start, pp_end);
  }
  
  private static Range setLocation(final Feature f, 
                                   final int start, final int end)
  {
    try
    {
      final RangeVector ranges = new RangeVector();
      final Range range = new Range(start, end);
      ranges.add(range);

      final Location new_location = new Location(ranges, 
                        f.getLocation().isComplement());
      f.setLocation(new_location);
      return range;
    }
    catch(OutOfRangeException e)
    {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    catch(ReadOnlyException e)
    {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    return null;
  }
  
  public static String getUniqueName(final Feature feature)
  {
    try
    {
      return (String)feature.getQualifierByName("ID").getValues().get(0);
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
    return null;
  }
  
  public static void main(String args[])
  {
    GeneUtils.defineShowHideGeneFeatures(new FeatureVector());
  }
}
