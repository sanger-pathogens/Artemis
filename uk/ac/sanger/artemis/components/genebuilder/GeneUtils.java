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
import java.awt.Cursor;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Collection;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import org.gmod.schema.general.DbXRef;
import org.gmod.schema.sequence.FeatureCvTerm;
import org.gmod.schema.sequence.FeatureDbXRef;
import org.gmod.schema.sequence.FeaturePub;
import org.gmod.schema.sequence.FeatureSynonym;

import uk.ac.sanger.artemis.chado.ClusterLazyQualifierValue;
import uk.ac.sanger.artemis.chado.FeatureForUpdatingResidues;
import uk.ac.sanger.artemis.chado.FeatureLocLazyQualifierValue;
import uk.ac.sanger.artemis.components.EditMenu;
import uk.ac.sanger.artemis.components.MessageDialog;
import uk.ac.sanger.artemis.components.SelectionMenu;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.DatabaseInferredFeature;
import uk.ac.sanger.artemis.io.DocumentEntry;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.Feature;
import uk.ac.sanger.artemis.io.GFFDocumentEntry;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.KeyVector;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierLazyLoading;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.RangeVector;
import uk.ac.sanger.artemis.sequence.MarkerRange;
import uk.ac.sanger.artemis.util.ByteBuffer;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.EntryVector;
import uk.ac.sanger.artemis.FeatureKeyQualifierPredicate;
import uk.ac.sanger.artemis.FeaturePredicate;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.GotoEventSource;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.Selection;


public class GeneUtils
{
  private static Vector<String> hideFeatures = new Vector<String>();
  private static JCheckBox showObsolete = new JCheckBox("Show Obsolete Features",false);
  private static String nonCodingTranscripts[] =
                                { "tRNA", "rRNA", "snRNA", "snoRNA", "ncRNA", "scRNA" };
  private static StringVector featuresToUpdateResidues = 
    Options.getOptions().getOptionValues("sequence_update_features");
  
  static
  {
    hideFeatures.add("polypeptide");
    hideFeatures.add(DatabaseDocument.TRANSCRIPT);
    hideFeatures.add("pseudogenic_transcript");
  }
  
  /**
   * Used when a whole sequence is loaded in and the features are loaded
   * lazily
   * @param feature
   */
  public static void addLazyQualifiers(final GFFStreamFeature feature)
  {
    if(feature.isLazyLoaded() || feature.getChadoLazyFeature() == null)
      return;
    
    // synonyms
    final Collection featureSynonyms = feature.getChadoLazyFeature().getFeatureSynonyms();
    
    final Iterator it = featureSynonyms.iterator();
    while(it.hasNext())
    {
      final FeatureSynonym featureSynonym = (FeatureSynonym) it.next();
      final String name = featureSynonym.getSynonym().getCvTerm().getName();
      String value = featureSynonym.getSynonym().getName();
      
      if(!featureSynonym.isCurrent())
        value.concat(GFFStreamFeature.encode(";current=false"));      
      
      Qualifier qualifier = feature.getQualifiers().getQualifierByName(name);
      if(qualifier == null)
        qualifier = new Qualifier(name, value);
      else
        qualifier.addValue(value);
      
      feature.getQualifiers().setQualifier(qualifier);
    }
    
    // dbxrefs
    if(feature.getQualifierByName("Dbxref") == null)
    {
      DbXRef dbxref = feature.getChadoLazyFeature().getDbXRef();
      if(dbxref != null)
      {
        String value = dbxref.getDb().getName() + ":" + 
                       dbxref.getAccession();
        feature.getQualifiers().setQualifier(new Qualifier("Dbxref", value));
        
        if(feature.isReadOnly() && feature.getKey().equals("polypeptide_domain"))
        {
          value= "protein motif:"+value;
          feature.getQualifiers().setQualifier(new Qualifier("inference", value));
        }
      }
    }
   
    
    final Collection featureDbXRefs = feature.getChadoLazyFeature().getFeatureDbXRefs();
    final Iterator it2 = featureDbXRefs.iterator();
    while(it2.hasNext())
    {
      final FeatureDbXRef featureDbXRef = (FeatureDbXRef) it2.next();
      String value = featureDbXRef.getDbXRef().getDb().getName() + ":" + 
                     featureDbXRef.getDbXRef().getAccession();
      
      Qualifier qualifier = feature.getQualifiers().getQualifierByName("Dbxref");
      if(qualifier == null)
        qualifier = new Qualifier("Dbxref", value);
      else
        qualifier.addValue(value);
      feature.getQualifiers().setQualifier(qualifier);
    }
    
    
    // feature cvterms (GO, product....)
    final Collection featureCvTerms = feature.getChadoLazyFeature().getFeatureCvTerms();
    if(featureCvTerms != null)
    {
      final Iterator it3 = featureCvTerms.iterator();
      while(it3.hasNext())
      {
        FeatureCvTerm featureCvTerm = (FeatureCvTerm)it3.next();
        List featureCvTermDbXRefList = null;
       
        if(featureCvTerm.getFeatureCvTermDbXRefs() != null)
          featureCvTermDbXRefList = new Vector(featureCvTerm.getFeatureCvTermDbXRefs());
      
        List featureCvTermPubList = null;
      
        if(featureCvTerm.getFeatureCvTermPubs() != null)
          featureCvTermPubList = new Vector(featureCvTerm.getFeatureCvTermPubs());
       
        ByteBuffer this_buff = new ByteBuffer();
        DatabaseDocument.appendControlledVocabulary(this_buff, null, featureCvTerm,
                                 featureCvTermDbXRefList,featureCvTermPubList, null, false);
      
        final String qualifierString = new String(this_buff.getBytes());
        int ind = qualifierString.indexOf('=');
        final String name  = qualifierString.substring(0, ind);
        final String value = GFFStreamFeature.decode(
          qualifierString.substring(ind+1, qualifierString.length()-1));
      
        Qualifier qualifier = feature.getQualifiers().getQualifierByName(name);
        if(qualifier == null)
          qualifier = new Qualifier(name, value);
        else
          qualifier.addValue(value);
        feature.getQualifiers().setQualifier(qualifier);
      }
    }
    // feature pubs - literature
    final Collection featurePubs = feature.getChadoLazyFeature().getFeaturePubs();
    
    if(featurePubs != null)
    {
      final Iterator it4 = featurePubs.iterator();
      while(it4.hasNext())
      {
        FeaturePub featurePub = (FeaturePub) it4.next();

        Qualifier qualifier = feature.getQualifiers().getQualifierByName(
            "literature");
        if(qualifier == null)
          qualifier = new Qualifier("literature", featurePub.getPub()
              .getUniqueName());
        else
          qualifier.addValue(featurePub.getPub().getUniqueName());
        feature.getQualifiers().setQualifier(qualifier);
      }
    }
    
    feature.setLazyLoaded(true);
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
  
  public static void addSegment(final GFFStreamFeature feature,
                                final RangeVector rangesToAdd,
                                final String transcriptName) 
         throws ReadOnlyException, EntryInformationException
  {
    // add new ID
    final Hashtable id_store = feature.getSegmentRangeStore();
    String prefix[] = null;
    Enumeration enum_ids = id_store.keys();
    while(enum_ids.hasMoreElements())
    {
      String id = (String) enum_ids.nextElement();
      prefix = feature.getPrefix(id, ':');
      if(prefix[0] != null)
        break;
    }

    // USE PREFIX TO CREATE NEW ID
    RangeVector rv = (RangeVector)feature.getLocation().getRanges().clone();
    for(int i=0; i<rangesToAdd.size(); i++)
    {
      final Range range = (Range) rangesToAdd.elementAt(i);
      final String ID;
      if(prefix[0] != null)
      {
        int auto_num = feature.getAutoNumber(prefix[0], ':');
        ID = prefix[0] + ":" + auto_num;
        feature.getSegmentRangeStore().put(ID, range);
      }
      else
      {
        String key = feature.getKey().toString();
        ID = transcriptName + ":" + key + ":1";
        feature.getSegmentRangeStore().put(ID, range);
      }

      if(!rv.containsRange(range))
        rv.add(range);
    }
    
    feature.setQualifier(new Qualifier("ID", feature.getSegmentID( rv )));
  }
  
  /**
   * Used when writing the database entry to a file. This routine
   * forces lazy-loading qualifier values to be read in full.
   * @param entry
   * @param parent
   */
  public static void lazyLoadAll(final Entry entry, final JFrame parent)
  {
    final List lazySimilarityValues = new Vector();
    final List lazyClusterValues = new Vector();
    final FeatureVector features = entry.getAllFeatures();
    // find any lazy values to be loaded
     
    for(int i=0; i<features.size(); i++)
    {
      QualifierVector qualifiers = features.elementAt(i).getQualifiers();
      for(int j=0; j<qualifiers.size(); j++)
      {
        Qualifier qualifier = (Qualifier)qualifiers.get(j);
        if(qualifier instanceof QualifierLazyLoading &&
           !((QualifierLazyLoading)qualifier).isAllLazyValuesLoaded())
        {
          if( ((QualifierLazyLoading)qualifier).getValue(0) instanceof FeatureLocLazyQualifierValue )
            lazySimilarityValues.addAll( ((QualifierLazyLoading)qualifier).getLazyValues() );
          else if( ((QualifierLazyLoading)qualifier).getValue(0) instanceof ClusterLazyQualifierValue )
          {
            List lazyValues = ((QualifierLazyLoading)qualifier).getLazyValues();
            lazyClusterValues.addAll(lazyValues);
          }
          
          ((QualifierLazyLoading)qualifier).setForceLoad(true);
        }
      }
    }
    
    if(lazySimilarityValues.size() > 0 || lazyClusterValues.size() > 0)
    {
      int n = JOptionPane.YES_OPTION;  
      if(parent != null)
        n =  JOptionPane.showConfirmDialog(parent,
          "Load and write to file all qualifers from the database?"+
          "\nThis may take a few minutes.",
          "Load All Data",
          JOptionPane.YES_NO_OPTION);
      
      if(n == JOptionPane.YES_OPTION)
      {
        if(parent != null)
          parent.setCursor(new Cursor(Cursor.WAIT_CURSOR));
        final DatabaseDocument document = 
          (DatabaseDocument)((DocumentEntry)entry.getEMBLEntry()).getDocument();
        
        if(lazySimilarityValues.size() > 0)
          FeatureLocLazyQualifierValue.bulkRetrieve(lazySimilarityValues,document);
        
        if(lazyClusterValues.size() > 0)
          ClusterLazyQualifierValue.setClusterFromValueList(lazyClusterValues, document);

        for(int i=0; i<features.size(); i++)
        {
          GFFStreamFeature feature = (GFFStreamFeature)(features.elementAt(i).getEmblFeature());
          if(feature.isReadOnly() && 
             feature.getKey().equals("polypeptide_domain") &&
             feature.getChadoLazyFeature() != null)
          {
            // load dbxrefs for domains
            if(feature.getQualifierByName("Dbxref") == null)
            {
              DbXRef dbxref = feature.getChadoLazyFeature().getDbXRef();
              if(dbxref != null)
              {
                String value = dbxref.getDb().getName() + ":" + 
                               dbxref.getAccession();
                feature.getQualifiers().setQualifier(new Qualifier("Dbxref", value));
                
                value= "protein motif:"+value;
                feature.getQualifiers().setQualifier(new Qualifier("inference", value));
              }
            }
          }
        }
        if(parent != null)
          parent.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
      }
    }
  }
  
  /**
   * Sorts the elements of the vector using a simple O(n^2) selection
   * sort.
   */
  private static void sort(Vector v)
  {
    int smallest;

    for (int i = 0; i < v.size (); ++i)
    {
      //find smallest remaining element
      smallest = i;
      for(int j = i + 1 ; j < v.size () ; ++j)
      {
        if(((String)v.get(j)).compareTo( (String)v.get(smallest)) < 0)
          smallest = j;
      }
      //exchange smallest and i
      if (smallest != i)
      {
        final String tmp = (String)v.get(i);
        v.setElementAt (v.get(smallest), i);
        v.setElementAt (tmp, smallest);
      }
    }
  }

  
  /**
   * Given a collection of features, determine if these should be
   * shown or hidden in the Artemis display
   * @param features
   */
  public static void defineShowHideGeneFeatures(final FeatureVector features)
  {
    final KeyVector keys = 
      features.elementAt(0).getEntry().getEntryInformation().getSortedValidKeys ();
    final Vector showFeatures = new Vector();
    for(int i=0; i<keys.size(); i++)
    {
      String keyStr = ((Key)keys.get(i)).getKeyString();
      if( !hideFeatures.contains(keyStr) )
        showFeatures.add(keyStr);
    }
    
    sort(hideFeatures);
    
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
          
          hideFeatures.add(hideKey);
          if(showFeatures.contains(hideKey))
            showFeatures.remove(hideKey);
          
          sort(hideFeatures);
          hideListModel.add(hideFeatures.indexOf(hideKey), hideKey);
          showListModel.removeElement(hideKey);
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
          
          if(hideFeatures.contains(showKey))
            hideFeatures.remove(showKey);
          showFeatures.add(showKey);
          sort(showFeatures);
          
          showListModel.add(showFeatures.indexOf(showKey), showKey);
          hideListModel.removeElement(showKey);
        }
      }
    });

    bdown = Box.createVerticalBox();
    bdown.add(Box.createVerticalGlue());
    bdown.add(new JLabel("Features Hidden:"));
    bdown.add(new JScrollPane(hideList));
    bdown.add(show_butt);
    hideShowPanel.add(bdown, BorderLayout.EAST);
    
    hideShowPanel.add(showObsolete, BorderLayout.SOUTH);

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
        if(isObsolete((GFFStreamFeature)feature))
        {
          if(!showObsolete.isSelected())
          {
            ((GFFStreamFeature)feature).setVisible(false);
            continue;
          }
        }
        
        final String key = feature.getKey().getKeyString();
        if(hideFeatures.contains(key))
          ((GFFStreamFeature)feature).setVisible(false);
        else 
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
   * Test if this feature is obsolete
   * @param feature
   * @return
   */
  public static boolean isObsolete(final uk.ac.sanger.artemis.io.GFFStreamFeature feature)
  {
    Qualifier qualifier = feature.getQualifierByName("isObsolete");
    if(qualifier == null)
      return false;
    
    if( ((String)qualifier.getValues().get(0)).equals("true") )
      return true;
    
    return false;
  }
  
  /**
   * 
   */
  public static Vector<ChadoCanonicalGene> duplicateGeneModel(final JFrame frame,
      final FeatureVector features_to_duplicate,
      final EntryGroup entry_group)
  {
    final Vector<ChadoCanonicalGene> duplicatedGenes = new Vector<ChadoCanonicalGene>();
    final Vector<ChadoCanonicalGene> newGenes = new Vector<ChadoCanonicalGene>();
    for (int i = 0 ; i < features_to_duplicate.size () ; ++i) 
    {
      final GFFStreamFeature gffFeature = 
           (GFFStreamFeature)features_to_duplicate.elementAt(i).getEmblFeature();
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
        newGenes.add(newchadoGene);
        ((GFFStreamFeature)newGeneFeature.getEmblFeature()).setChadoGene(newchadoGene);
        newchadoGene.setGene(newGeneFeature.getEmblFeature());
        
        final List<Feature> transcripts = chadoGene.getTranscripts();
        for(int j=0; j<transcripts.size(); j++) // duplicate transcripts and children
        {
          final GFFStreamFeature transcript = (GFFStreamFeature)transcripts.get(j);
          final String transcriptName = getUniqueName(transcript);
          
          uk.ac.sanger.artemis.Feature newTranscriptFeature = duplicateFeature(transcript, newchadoGene);
          newchadoGene.addTranscript(newTranscriptFeature.getEmblFeature());
          final String newTranscriptName = getUniqueName(newTranscriptFeature.getEmblFeature());

          List<uk.ac.sanger.artemis.Feature> newFeatures= 
              duplicateFeatures(chadoGene.get3UtrOfTranscript(transcriptName), newchadoGene);
          for(uk.ac.sanger.artemis.Feature utrFeature: newFeatures)
            newchadoGene.add3PrimeUtr(newTranscriptName, utrFeature.getEmblFeature());
          
          newFeatures = duplicateFeatures(chadoGene.get5UtrOfTranscript(transcriptName), newchadoGene);
          for(uk.ac.sanger.artemis.Feature utrFeature: newFeatures)
            newchadoGene.add5PrimeUtr(newTranscriptName, utrFeature.getEmblFeature());
          
          newFeatures = duplicateFeatures(chadoGene.getOtherFeaturesOfTranscript(transcriptName), newchadoGene);
          for(uk.ac.sanger.artemis.Feature otherFeature: newFeatures)
            newchadoGene.addOtherFeatures(newTranscriptName, otherFeature.getEmblFeature());

          newFeatures = duplicateFeatures(chadoGene.getSplicedFeaturesOfTranscript(transcriptName), newchadoGene);
          for(uk.ac.sanger.artemis.Feature splicedFeature: newFeatures)
            newchadoGene.addSplicedFeatures(newTranscriptName, splicedFeature.getEmblFeature());
          
          uk.ac.sanger.artemis.Feature newProtein = 
            duplicateFeature(chadoGene.getProteinOfTranscript(transcriptName), newchadoGene);
          if(newProtein != null)
            newchadoGene.addProtein(newTranscriptName, newProtein.getEmblFeature());
        }
      } 
      catch (ReadOnlyException e) {}
    }
    
    duplicatedGenes.clear();
    return newGenes;
  }
  
  
  private static List<uk.ac.sanger.artemis.Feature> duplicateFeatures(
                                 final List<Feature> featuresOfTranscript,
                                 final ChadoCanonicalGene chadoGene) 
          throws ReadOnlyException
  {
    final List<uk.ac.sanger.artemis.Feature> newFeatures = 
        new Vector<uk.ac.sanger.artemis.Feature>();
    
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
    
    QualifierVector qualifiers = new QualifierVector();
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
      GFFStreamFeature exonFeature =
         GeneViewerPanel.addExonFeature(chadoGene, entry_group, 
          null, new_location.getTotalRange(), transcriptId, selection, 
          new Key(DatabaseDocument.EXONMODEL), null);
      
      // add protein
      uk.ac.sanger.artemis.Feature polypep =
        GeneViewerPanel.addProteinFeature(chadoGene, entry_group, transcriptId, transcript);
      newFeatures.add(polypep);
      
      // add inferred CDS
      if(DatabaseDocument.CHADO_INFER_CDS)
        DatabaseInferredFeature.createFeature(transcriptId, exonFeature, 
                              chadoGene, entry_group.getDefaultEntry());
      
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
        uk.ac.sanger.artemis.io.DatabaseDocumentEntry ||
       default_entry.getEMBLEntry() instanceof 
        uk.ac.sanger.artemis.io.GFFDocumentEntry)
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
   * Prompt the user for an ID and provide a default automated ID
   * based on the range.
   * @param entry_group
   * @param is_forward
   * @param range
   * @return
   */
  public static String promptForUniquename(final EntryGroup entry_group,
                                           final boolean is_forward,
                                           final Range range)
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
        
        
        id = JOptionPane.showInputDialog(null, msg, 
            default_entry.getName()+":"+
            range.getStart()+".."+
            range.getEnd());
        
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
  
  /**
   * Given an group of entries determine if they contain a GFF entry
   * @param entryGroup
   * @return
   */
  public static boolean isGFFEntry(final EntryGroup entryGroup)
  {
    final EntryVector entries = entryGroup.getActiveEntries();
    
    for(int i=0; i<entries.size(); i++)
    {
      if( entries.elementAt(i).getEMBLEntry() instanceof GFFDocumentEntry )
        return true;
    }
    return false;
  }
  
  
  /**
   * Given a feature determine if it belongs to a database entry
   * @param entryGroup
   * @return
   */
  public static boolean isDatabaseEntry(final Feature feature)
  {
    if( feature.getEntry() instanceof DatabaseDocumentEntry &&
        ((GFFStreamFeature)feature).getDocumentEntry().getDocument() instanceof DatabaseDocument)
      return true;

    return false;
  }
  

  private static void deleteFeature(uk.ac.sanger.artemis.Feature feature)
      throws ReadOnlyException
  {
    if(feature != null && feature.getEntry() != null)
      feature.removeFromEntry();
  }

  public static void deleteAllFeature(uk.ac.sanger.artemis.Feature feature,
      final ChadoCanonicalGene chado_gene) throws ReadOnlyException
  {
    deleteAllFeature(feature, chado_gene, true); 
  }
  
  /**
   * Delete feature and children in a chado gene model
   * @param feature
   * @param chado_gene
   * @throws ReadOnlyException
   */
  public static void deleteAllFeature(uk.ac.sanger.artemis.Feature feature,
      final ChadoCanonicalGene chado_gene, final boolean updateChadoCanonicalGene) throws ReadOnlyException
  {
    Set<Feature> children = chado_gene.getChildren(feature.getEmblFeature());
    deleteFeature(feature);
    if(updateChadoCanonicalGene)
      chado_gene.deleteFeature(feature.getEmblFeature());

    Feature embl_feature;
    Iterator<Feature> it = children.iterator();

    while(it.hasNext())
    {
      embl_feature = it.next();
      deleteFeature((uk.ac.sanger.artemis.Feature) embl_feature.getUserData());
      if(updateChadoCanonicalGene)
        chado_gene.deleteFeature(embl_feature);
    }
  }
  
  /**
   * Check gene model strands for any inconsistencies
   * @param chado_gene
   * @return true is gene model ranges are correct
   */
  public static boolean isStrandOK(final ChadoCanonicalGene chado_gene)
  {
    boolean isRev = chado_gene.getGene().getLocation().isComplement();
    final List<Feature> transcripts = chado_gene.getTranscripts();
    
    for(int i=0; i<transcripts.size(); i++)
    {
      final Feature transcript = (Feature)transcripts.get(i);
      
      if(isRev ^ transcript.getLocation().isComplement())
        return false;
      
      final Feature protein = 
        chado_gene.getProteinOfTranscript(GeneUtils.getUniqueName(transcript));
      if(protein != null && (isRev ^ protein.getLocation().isComplement()))
        return false;
      
      final Set<Feature> children = chado_gene.getChildren(transcript);
      final Iterator<Feature> it = children.iterator();
      while(it.hasNext())
      {
        final Feature feature = it.next();
        if(isRev ^ feature.getLocation().isComplement())
          return false;
      }
    }
    return true;
  }

  /**
   * Check gene model boundaries for any inconsistencies
   * @param chado_gene
   * @return 0 - if consisent
   *         1 - if transcript start or end is outside gene range
   *         2 - if child feature of a transcript is outside the transcript range
   *         3 - if the span of the children features does not match start and end of the transcript
   *         4 - if the protein range does not match CDS
   *         5 - if the gene range does not match the largest transcript range
   */
  public static int isBoundaryOK(final ChadoCanonicalGene chado_gene)
  {
    final Range geneRange = chado_gene.getGene().getLocation().getTotalRange();
    final List<Feature> transcripts = chado_gene.getTranscripts();
    int geneStart = Integer.MAX_VALUE;
    int geneEnd = -1;
    
    for(Feature transcript: transcripts)
    {
      final Range transcriptRange = transcript.getLocation().getTotalRange();
      int transcriptStart = Integer.MAX_VALUE;
      int transcriptEnd = -1;
      int ppStart = Integer.MAX_VALUE;
      int ppEnd   = -1;
      
      if(transcriptRange.getStart() < geneRange.getStart() ||
         transcriptRange.getEnd()   > geneRange.getEnd())
        return 1;
      
      if(transcriptRange.getStart() < geneStart)
        geneStart = transcriptRange.getStart();
      if(transcriptRange.getEnd() > geneEnd)
        geneEnd = transcriptRange.getEnd();
      
      final Feature protein = 
        chado_gene.getProteinOfTranscript(GeneUtils.getUniqueName(transcript));
      String proteinName = null;
      if(protein != null)
        proteinName = GeneUtils.getUniqueName(protein);
      
      final Set<Feature> children = chado_gene.getChildren(transcript);
      final Iterator<Feature> it = children.iterator();
      while(it.hasNext())
      {
        final Feature feature = it.next();
        final Range childRange = feature.getLocation().getTotalRange();
        if(childRange.getStart() < transcriptRange.getStart() ||
           childRange.getEnd()   > transcriptRange.getEnd())
          return 2;

        if(proteinName != null &&
           GeneUtils.getUniqueName(feature).equals(proteinName))
          continue;
        
        if(childRange.getStart() < transcriptStart)
          transcriptStart = childRange.getStart();
        if(childRange.getEnd() > transcriptEnd )
          transcriptEnd = childRange.getEnd();
        
        String keyStr = feature.getKey().getKeyString();
        if( (DatabaseDocument.CHADO_INFER_CDS  && keyStr.equals("CDS")) ||
            (!DatabaseDocument.CHADO_INFER_CDS && keyStr.equals(DatabaseDocument.EXONMODEL)) ||
           feature.getKey().equals("pseudogenic_exon"))
        {
          if(childRange.getStart() < ppStart)
            ppStart = childRange.getStart();
          if(childRange.getEnd() > ppEnd )
            ppEnd   = childRange.getEnd();
        }
      }
      
      if((transcriptRange.getStart() != transcriptStart && transcriptStart < Integer.MAX_VALUE) ||
         (transcriptRange.getEnd()   != transcriptEnd   && transcriptEnd   > -1))
        return 3;

      if(protein != null)
      {
        final Range proteinRange = protein.getLocation().getTotalRange();
        if((proteinRange.getStart() != ppStart && ppStart < Integer.MAX_VALUE) ||
           (proteinRange.getEnd()   != ppEnd   && ppEnd   > -1))
          return 4;
      }
    }
    
    // check gene range
    if((geneRange.getStart() != geneStart && geneStart < Integer.MAX_VALUE) ||
       (geneRange.getEnd()   != geneEnd   && geneEnd   > -1))
      return 5;
    
    return 0;
  }
  
  public static void checkGeneBoundary(final ChadoCanonicalGene chado_gene)
  {
    checkGeneBoundary(chado_gene, true);
  }
  
  protected static Range checkTranscriptBoundary(
      final uk.ac.sanger.artemis.Feature transcript,
      final ChadoCanonicalGene chado_gene)
  {
    return checkTranscriptBoundary(transcript, chado_gene, true);
  }
  
  /**
   * Adjust transcript and gene boundaries
   * @param chado_gene
   */
  public static void checkGeneBoundary(final ChadoCanonicalGene chado_gene, final boolean changeEmblFeature)
  {
    final List<Feature> transcripts = chado_gene.getTranscripts();
    int gene_start = Integer.MAX_VALUE;
    int gene_end = -1;
    
    Range range;
    for(Feature transcript: transcripts)
    {
      range = checkTranscriptBoundary(
          (uk.ac.sanger.artemis.Feature)transcript.getUserData(), chado_gene, changeEmblFeature);
      if(range != null && range.getStart() < gene_start)
        gene_start = range.getStart();
      if(range != null && range.getEnd() > gene_end)
        gene_end = range.getEnd();
    }
    
    if(gene_end == -1 && gene_start == Integer.MAX_VALUE)
      return;
    
    setLocation(chado_gene.getGene(), gene_start, gene_end, changeEmblFeature);
  }
  
  /**
   * Check and adjust transcript boundary
   * @param transcript
   * @param chado_gene
   * @param changeEmblFeature
   * @return
   */
  protected static Range checkTranscriptBoundary(
      final uk.ac.sanger.artemis.Feature transcript,
      final ChadoCanonicalGene chado_gene,
      final boolean changeEmblFeature)
  {
    final List transcripts = chado_gene.getTranscripts();

    if(transcripts.contains(transcript.getEmblFeature()))
    {
      checkProteinBoundary(transcript.getEmblFeature(), chado_gene, changeEmblFeature);
      
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
            transcript_start, transcript_end, changeEmblFeature);
    }
    else
      JOptionPane.showMessageDialog(null,
          "Select a single transcript and try again.", "Transcript Selection",
          JOptionPane.ERROR_MESSAGE);
    return null;
  }
  
  /**
   * Check and adjust protein boundary
   * @param transcript
   * @param chado_gene
   */
  private static void checkProteinBoundary(final Feature transcript,
                                          final ChadoCanonicalGene chado_gene,
                                          final boolean changeEmblFeature)
  {
    final String transcriptName = getUniqueName(transcript);
    final Feature protein = chado_gene.getProteinOfTranscript(transcriptName);
    if(protein == null)
      return;
    
    int pp_start = Integer.MAX_VALUE;
    int pp_end = -1;
    
    final List<Feature> dnaFeatures = new Vector<Feature>();
/*    if(chado_gene.get3UtrOfTranscript(transcriptName) != null)
      dnaFeatures.addAll(chado_gene.get3UtrOfTranscript(transcriptName));
    if(chado_gene.get5UtrOfTranscript(transcriptName) != null)
      dnaFeatures.addAll(chado_gene.get5UtrOfTranscript(transcriptName));*/ 
    
    List<Feature> exons;
    if(DatabaseDocument.CHADO_INFER_CDS)
      exons = chado_gene.getSpliceSitesOfTranscript(transcriptName, "CDS");
    else
      exons = chado_gene.getSpliceSitesOfTranscript(transcriptName, DatabaseDocument.EXONMODEL);
    if(exons != null)
      dnaFeatures.addAll(exons);
    
    exons = chado_gene.getSpliceSitesOfTranscript(transcriptName, "pseudogenic_exon");
    if(exons != null)
      dnaFeatures.addAll(exons);
    
    for(Feature dnaFeature: dnaFeatures)
    {
      final Range range = dnaFeature.getLocation().getTotalRange();
      if(range.getStart() < pp_start)
        pp_start = range.getStart();
      if(range.getEnd() > pp_end)
        pp_end = range.getEnd();
    }
    
    if(pp_start == Integer.MAX_VALUE || pp_end == -1)
       return;
    
    setLocation(protein, pp_start, pp_end, changeEmblFeature);
  }
  
  /**
   * For a feature propagate the uniquename as the prefix for
   * the associated children.
   * @param gene
   * @param newName
   * @param children
   */
  public static void propagateId(final GFFStreamFeature feature,
                                 final String newName, 
                                 final Set children)
  {
    final ChadoCanonicalGene gene = feature.getChadoGene();
    final Iterator it = children.iterator();
    while(it.hasNext())
    {
      final GFFStreamFeature child = (GFFStreamFeature)it.next();
      final Hashtable segmentHash  = child.getSegmentRangeStore();
      
      final String oldId = getUniqueName(child);
      
      final Set childrenOfChild = gene.getChildren(child);
      int index = oldId.lastIndexOf('.');
      
      if(index == -1)
        index = oldId.indexOf(':');
      
      if(index > -1)
      {
        final String newId;
        if(  segmentHash != null && 
           ( segmentHash.size() > 1 ||
             child.getKey().getKeyString().equals(DatabaseDocument.EXONMODEL)))
        {
          final Set idKeys = segmentHash.keySet();
          final Hashtable newSegmentHash = new Hashtable(idKeys.size());
          final Iterator itKeys = idKeys.iterator();
          final Hashtable newIdMapToOldId = new Hashtable(idKeys.size());
          while(itKeys.hasNext())
          {
            String oldKey = (String)itKeys.next();
            index = oldKey.lastIndexOf('.');
            if(index == -1)
              index = oldKey.indexOf(':');
            final String newKey = newName + oldKey.substring(index);
            Object range = segmentHash.get(oldKey);

            newSegmentHash.put(newKey, range);
            newIdMapToOldId.put(newKey, oldKey);
          }
          child.setSegmentRangeStore(newSegmentHash);
          child.setNewIdMapToOldId(newIdMapToOldId);
          newId = child.getSegmentID(child.getLocation().getRanges());
        }
        else
          newId = newName + oldId.substring(index);
        try
        {
          ((uk.ac.sanger.artemis.Feature) child.getUserData())
              .setQualifier(new Qualifier("ID", newId));
          gene.updateUniqueName(oldId, newId, childrenOfChild);
        }
        catch(ReadOnlyException e)
        {
          e.printStackTrace();
        }
        catch(EntryInformationException e)
        {
          e.printStackTrace();
        }
      }
    }
  }
  
  /**
   * Fix parent (Parent and Derives_from) qualifiers. Used when the
   * uniqueName of the parent has been changed.
   * @param oldName
   * @param newName
   * @param children
   */
  public static void fixParentQualifier(final String oldName,
                                        final String newName, 
                                        final Set children)
  {
    final Iterator it = children.iterator();
    while(it.hasNext())
    {
      Feature child = (Feature)it.next();
      QualifierVector qualifiers = child.getQualifiers();
      if( qualifiers.getQualifierByName("Parent") != null &&
          ((String)(qualifiers.getQualifierByName("Parent").getValues().get(0))).equals(oldName) )
      {
        qualifiers.removeQualifierByName("Parent");
        qualifiers.setQualifier(new Qualifier("Parent", newName));
      }
      else if( qualifiers.getQualifierByName("Derives_from") != null &&
          ((String)(qualifiers.getQualifierByName("Derives_from").getValues().get(0))).equals(oldName) )
      {
        qualifiers.removeQualifierByName("Derives_from");
        qualifiers.setQualifier(new Qualifier("Derives_from", newName));
      }
    }
  }
  
  /**
   * Converts a gene to pseudogene or vice-versa
   * @param chado_gene
   * @throws EntryInformationException 
   * @throws ReadOnlyException 
   * @throws OutOfRangeException 
   */
  public static void convertPseudo(final ChadoCanonicalGene chado_gene) 
         throws ReadOnlyException, EntryInformationException, OutOfRangeException
  {
    final Key geneKey = chado_gene.getGene().getKey();
    
    boolean convertToPseudogene = false;
    
    uk.ac.sanger.artemis.Feature gene = (uk.ac.sanger.artemis.Feature)chado_gene.getGene().getUserData();
    if(geneKey.equals("gene"))
    {
      gene.set(new Key("pseudogene"), gene.getLocation(), gene.getQualifiers());
      convertToPseudogene = true;
    }
    else if(geneKey.equals("pseudogene"))
      gene.set(new Key("gene"), gene.getLocation(), gene.getQualifiers());
    else
      return;
    
    final List transcripts = chado_gene.getTranscripts();
    for(int i=0; i<transcripts.size(); i++)
    {
      final uk.ac.sanger.artemis.Feature transcript = (
          uk.ac.sanger.artemis.Feature)((Feature)transcripts.get(i)).getUserData();
      final String transcriptName = getUniqueName(transcript.getEmblFeature());
      final List exons;
      if(convertToPseudogene)
      {
        exons = chado_gene.getSpliceSitesOfTranscript(transcriptName, 
                                         DatabaseDocument.EXONMODEL);
        transcript.set(new Key("pseudogenic_transcript"), transcript.getLocation(), 
            transcript.getQualifiers());
      }
      else
      {
        exons = chado_gene.getSpliceSitesOfTranscript(transcriptName, 
                                                 "pseudogenic_exon");
        transcript.set(new Key(DatabaseDocument.TRANSCRIPT), transcript.getLocation(), transcript.getQualifiers());
      }
      
      if(exons == null)
        continue;
      
      for(int j=0; j<exons.size(); j++)
      {
        final uk.ac.sanger.artemis.Feature exon = (uk.ac.sanger.artemis.Feature)((Feature)exons.get(j)).getUserData();
        exon.resetColour();
        if(convertToPseudogene)
          exon.set(new Key("pseudogenic_exon"), exon.getLocation(), exon.getQualifiers());
        else
          exon.set(new Key(DatabaseDocument.EXONMODEL), exon.getLocation(), exon.getQualifiers());
      }
    }
  }

  private static Range setLocation(final Feature f, 
                                   final int start,
                                   final int end,
                                   final boolean changeEmblFeature)
  {
    try
    {
      final RangeVector ranges = new RangeVector();
      final Range range = new Range(start, end);
      ranges.add(range);

      final Location new_location = new Location(ranges, 
          f.getLocation().isComplement());
      
      if(changeEmblFeature)
        f.setLocation(new_location);
      else
        ((uk.ac.sanger.artemis.Feature)f.getUserData()).setLocation(new_location);
      return range;
    }
    catch (OutOfRangeException e)
    {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    catch (ReadOnlyException e)
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
  
  /**
   * Return an array on non-coding transcripts
   * @return
   */
  public static String[] getNonCodingTranscripts()
  {
    return nonCodingTranscripts;
  }
  
  /**
   * Test if the key given is a non-coding transcript key
   * @param key
   * @return
   */
  public static boolean isNonCodingTranscripts(final Key key)
  {
    for(int i=0; i<nonCodingTranscripts.length; i++)
      if(nonCodingTranscripts[i].equals(key.getKeyString()))
        return true;
    return false;
  }
  
  public static String deriveResidues(final GFFStreamFeature gffFeature)
  {
    final ChadoCanonicalGene chadoGene = gffFeature.getChadoGene();

    boolean isProteinFeature = 
      ((uk.ac.sanger.artemis.Feature)gffFeature.getUserData()).isProteinFeature();
       
    String residues = null;
    try
    {
      String transcriptName =
        chadoGene.getTranscriptFromName(GeneUtils.getUniqueName(gffFeature));

      List<Feature> splicedFeatures = 
        chadoGene.getSplicedFeaturesOfTranscript(transcriptName);
      for (Feature emblFeature: splicedFeatures)
      {
        if (emblFeature.getKey().getKeyString().equals(
            DatabaseDocument.EXONMODEL) ||
            emblFeature.getKey().getKeyString().equals("pseudogenic_exon"))
        {
          uk.ac.sanger.artemis.Feature f = 
            (uk.ac.sanger.artemis.Feature) emblFeature.getUserData();
          if (!isProteinFeature)
            residues = f.getTranslationBases();
          else
            residues = f.getTranslation().toString();
        }
      }
    }
    catch (Exception e){  }

    return residues;
  }
  
  public static FeatureForUpdatingResidues getFeatureForUpdatingResidues(
      final GFFStreamFeature gffFeature)
  {
    if(!isFeatureToUpdateResidues(gffFeature.getKey().getKeyString()))
      return null;
    String residues = deriveResidues(gffFeature);
    if(residues == null)
      return null;
    final FeatureForUpdatingResidues chadoFeature = new FeatureForUpdatingResidues();
    chadoFeature.setStartBase(0);
    chadoFeature.setLength(residues.length());
    chadoFeature.setNewSubSequence(residues);
    chadoFeature.setResidueUpdate(true);
    
    if(gffFeature.getQualifierByName("feature_id") != null)
      chadoFeature.setFeatureId( Integer.parseInt( (String)
          gffFeature.getQualifierByName("feature_id").getValues().get(0)) );
    else
    {
      chadoFeature.setFeatureId(-1);
      chadoFeature.setUniqueName(getUniqueName(gffFeature));
    }
    chadoFeature.setSeqLen(residues.length());
    return chadoFeature;
  }
  
  /**
   * Look at the sequence_update_features property in the options
   * file to see if this is a feature to update residues for.
   * @param keyStr
   * @return
   */
  public static boolean isFeatureToUpdateResidues(final String keyStr)
  {
    if(featuresToUpdateResidues == null)
      return false;
    
    return featuresToUpdateResidues.contains(keyStr);
  }
  
  public static void main(String args[])
  {
    GeneUtils.defineShowHideGeneFeatures(new FeatureVector());
  }
}
