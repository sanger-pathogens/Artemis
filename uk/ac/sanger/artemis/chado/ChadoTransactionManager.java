/* ChadoTransactionManager.java
 *
 * created: July 2005
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2005  Genome Research Limited
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
 */

package uk.ac.sanger.artemis.chado;

import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureSegment;
import uk.ac.sanger.artemis.FeatureSegmentVector;
import uk.ac.sanger.artemis.sequence.SequenceChangeListener;
import uk.ac.sanger.artemis.sequence.SequenceChangeEvent;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.RangeVector;
import uk.ac.sanger.artemis.io.StreamQualifier;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.EntryChangeListener;
import uk.ac.sanger.artemis.EntryChangeEvent;
import uk.ac.sanger.artemis.components.Splash;

import java.util.Collection;
import java.util.Vector;
import java.util.Hashtable;
import java.util.List;
import java.util.Enumeration;
import javax.swing.JOptionPane;

import org.gmod.schema.sequence.FeatureCvTermProp;
import org.gmod.schema.sequence.FeatureLoc;
import org.gmod.schema.sequence.FeatureProp;
import org.gmod.schema.sequence.FeatureDbXRef;
import org.gmod.schema.sequence.FeatureRelationship;
import org.gmod.schema.sequence.FeatureSynonym;
import org.gmod.schema.sequence.FeatureCvTerm;
import org.gmod.schema.sequence.Synonym;

import org.gmod.schema.cv.*;
import org.gmod.schema.pub.Pub;
import org.gmod.schema.general.Db;
import org.gmod.schema.general.DbXRef;
import org.gmod.schema.sequence.FeatureCvTermDbXRef;

/**
 *
 * Chado transaction manager listens for feature, entry and sequence changes.
 * <code>ChadoTransactionManager</code> creates and tracks the feature insertions,
 * deletions, and updates to commit back to the database.
 *
 **/
public class ChadoTransactionManager
       implements FeatureChangeListener, EntryChangeListener, SequenceChangeListener 
{

  private Vector sql = new Vector();
  
  /** GFF3 predefined tags */
  private String reserved_tags[] = 
          {   "ID",
              "Name",
              "Alias",
              "Parent",
              "Target",
              "Gap",
              "Derives_from",
              "Dbxref",
              "Ontology_term",
              "score",
              "codon_start",
              "gff_source",      // program or database
              "gff_seqname" };   // seqID of coord system
           
  //controlled vocab tags
  public static String cv_tags[] =
          {   "GO",
              "controlled_curation",
              "product" };
  
  //synonym tags from cv
  private String synonym_tags[] = null;
  private static String SYNONYM_TAG_CVNAME = "genedb_synonym_type";


  public ChadoTransactionManager()
  {
    
  }
  
  /**
   *  Implementation of the FeatureChangeListener interface.  We listen for
   *  changes in every feature of every entry in this group.
   **/ 
  public void featureChanged(FeatureChangeEvent event)
  {
    if(event.featureHasChanged())
    {
      final GFFStreamFeature feature = 
        (GFFStreamFeature)event.getFeature().getEmblFeature();

      if(event.getType() == FeatureChangeEvent.SEGMENT_CHANGED)
      {
        RangeVector rv_new = event.getNewLocation().getRanges();
        RangeVector rv_old = event.getOldLocation().getRanges();

        Splash.logger4j.debug("SEGMENT_CHANGED");
         
        if(rv_old.size() > rv_new.size()) // segment deleted
        {
          // delete segment
          int ideleted;
          Vector deleted = new Vector();
          boolean found;
          for(ideleted=0; ideleted<rv_old.size(); ideleted++)
          {   
            Range range = (Range)rv_old.get(ideleted);
            found = false;
            for(int j=0; j<rv_new.size(); j++)
            {
              if( ((Range)rv_new.get(j)).equals(range) )
                found = true;     
            }
            
            if(!found)
              deleted.add(new Integer(ideleted));
          }
   
          for(int i=0; i<deleted.size();i++)
          {
            ideleted = ((Integer)deleted.elementAt(i)).intValue();
            Range range_old = (Range)rv_old.elementAt(ideleted);
            String seg_id   = feature.getSegmentID(range_old);
            deleteFeature(seg_id);
            feature.getSegmentRangeStore().remove(seg_id);
          }
          
          String new_id = feature.getSegmentID(rv_new);
          Qualifier qualifier = new Qualifier("ID", new_id);
          try
          {
            feature.setQualifier(qualifier);
          }
          catch(ReadOnlyException e)
          {
            e.printStackTrace();
          }
          catch(EntryInformationException e)
          {
            e.printStackTrace();
          }
          
          // update feature_relationship.rank
          updateFeatureRelationshipRank(feature, rv_new);
        }
        else if(rv_old.size() < rv_new.size()) // feature segment added
        {
          boolean found;
          FeatureSegmentVector segments =
            ((uk.ac.sanger.artemis.Feature)feature.getUserData()).getSegments();
          for(int iadd=0; iadd<segments.size(); iadd++)
          {   
            FeatureSegment segment = segments.elementAt(iadd);
            Range range = segment.getRawRange();
            found = false;
            for(int j=0; j<rv_old.size(); j++)
            {
              if( ((Range)rv_old.get(j)).equals(range) )
                found = true;     
            }
            
            String transcript = 
              (String)feature.getQualifierByName("Parent").getValues().get(0);
            if(!found)
            {
              // find prefix fo ID
              Hashtable id_store = feature.getSegmentRangeStore();
              String prefix[] = null;
              Enumeration enum_ids = id_store.keys();
              while(enum_ids.hasMoreElements())
              {
                String id = (String)enum_ids.nextElement();
                prefix = feature.getPrefix(id, ':');
                if(prefix[0] != null)
                  break;
              }
              
              // USE PREFIX TO CREATE NEW ID
              if(prefix[0] != null)
              {
                int auto_num = feature.getAutoNumber(prefix[0], ':');
                feature.getSegmentRangeStore().put(prefix[0]+":"+auto_num, 
                                                   range);
              }
              else
              {
                String key = feature.getKey().toString();
                feature.getSegmentRangeStore().put(transcript+":"+ key +":1", 
                                                   range);
              }
       
              
              String new_id = feature.getSegmentID(rv_new);
              Qualifier qualifier = new Qualifier("ID", new_id);
              try
              {
                feature.setQualifier(qualifier);
              }
              catch(ReadOnlyException e)
              {
                e.printStackTrace();
              }
              catch(EntryInformationException e)
              {
                e.printStackTrace();
              }
              
              String segment_uniquename = feature.getSegmentID(range);
              insertFeatureSegment(segment,segment_uniquename);
              // update feature_relationship.rank
              updateFeatureRelationshipRank(feature, rv_new);
            }
          }
        }
      }
      else if(event.getType() == FeatureChangeEvent.LOCATION_CHANGED)
      {
        RangeVector rv_new = event.getNewLocation().getRanges();
        RangeVector rv_old = event.getOldLocation().getRanges();

        Splash.logger4j.debug("LOCATION_CHANGED "+feature.getFirstBase()+".."+feature.getLastBase()+
                              "   new="+rv_new.size()+" old="+rv_old.size());
        if(rv_new.size() != rv_old.size())
          return;
        
        ChadoTransaction tsn;
        int ichanged;
        Vector changes = new Vector();
        for(ichanged=0; ichanged<rv_old.size(); ichanged++)
        {
          Range rnew = (Range)rv_new.elementAt(ichanged);
          Range rold = (Range)rv_old.elementAt(ichanged);
   
          if(rnew.getStart() != rold.getStart() ||
             rnew.getEnd()   != rold.getEnd())
            changes.add(new Integer(ichanged));
        }
 
        for(int i=0; i<changes.size();i++)
        {
          ichanged = ((Integer)changes.elementAt(i)).intValue();
          
          Range range_new = (Range)rv_new.elementAt(ichanged);
          Range range_old = (Range)rv_old.elementAt(ichanged);
          String seg_id   = feature.getSegmentID(range_new);
          
          if(seg_id == null)
            seg_id   = feature.getSegmentID(range_old);
          
          if(feature.getSegmentRangeStore() != null)
            feature.getSegmentRangeStore().put(seg_id, range_new);
          
          FeatureLoc featureloc = getFeatureLoc(feature, seg_id, range_new);
          
          tsn = new ChadoTransaction(ChadoTransaction.UPDATE,
                                     featureloc,
                                     feature.getLastModified(), feature);

          sql.add(tsn);
        }
      }
      else if(event.getType() == FeatureChangeEvent.QUALIFIER_CHANGED)
      {
        Splash.logger4j.debug("QUALIFIER_CHANGED "
            +event.getOldQualifiers().getQualifierByName("ID").getValues().get(0));
        
        editKeyAndQualifiers(event.getOldQualifiers(),event.getNewQualifiers(),
            event.getOldKey(), event.getNewKey(),
            feature, FeatureChangeEvent.QUALIFIER_CHANGED);
      }
      else if(event.getType() == FeatureChangeEvent.ALL_CHANGED)
      {
        Splash.logger4j.debug("ALL_CHANGED "+event.getOldKey().toString()+"  "+
                                          event.getNewKey().toString());
        
        editKeyAndQualifiers(event.getOldQualifiers(),event.getNewQualifiers(),
                             event.getOldKey(), event.getNewKey(),
                             feature, FeatureChangeEvent.ALL_CHANGED);
        
        if(event.getOldKey().compareTo( event.getNewKey() ) != 0 &&
           event.getNewKey().toString().equals("gene") &&
           feature.getChadoGene() == null)
        {
          ChadoCanonicalGene chado_gene = new ChadoCanonicalGene();
          chado_gene.setGene(feature);
          feature.setChadoGene(chado_gene);
        }
      }
    }
  }
 
  /**
   *  Invoked when an Entry is changed.
   **/
  public void entryChanged(EntryChangeEvent event)
  {
    if(event.getType() == EntryChangeEvent.FEATURE_ADDED)
    { 
      // if this is a duplicate feature then ignore
      if(event.isDuplicate())
      {
        Feature feature = event.getFeature();
        Qualifier qualifier_uniquename;
        try
        {
          qualifier_uniquename = feature.getQualifierByName("ID");
          Splash.logger4j.debug("FEATURE_ADDED ------> DUPLICATE "+
              (String)(qualifier_uniquename.getValues()).elementAt(0));
        }
        catch(InvalidRelationException e)
        {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }

        return;
      }
     
      Splash.logger4j.debug("FEATURE_ADDED");
      Feature feature = event.getFeature();
      insertFeature(feature);
    }
    else if(event.getType() == EntryChangeEvent.FEATURE_DELETED)
    { 
      try
      {
        Qualifier qualifier_uniquename = event.getFeature().getQualifierByName("ID");
        String feature_uniquename = 
                             (String)(qualifier_uniquename.getValues()).elementAt(0);
        Splash.logger4j.debug("FEATURE_DELETED "+feature_uniquename);
        
        if(event.getFeature().getSegments().size() > 1)
        {
          GFFStreamFeature gff_feature =
            (GFFStreamFeature)event.getFeature().getEmblFeature();
          RangeVector ranges = gff_feature.getLocation().getRanges();
          for(int i=0; i<ranges.size(); i++)
          {
            Range range = (Range)ranges.get(i);
            feature_uniquename = gff_feature.getSegmentID(range);
            deleteFeature(feature_uniquename);
          }    
        }
        else
          deleteFeature(feature_uniquename);
      }
      catch(InvalidRelationException e)
      {
        e.printStackTrace();
      }
    }

//  System.out.println(event.getEntry().getName());
  }
  
  /**
   * Update spliced features rank
   * @param feature
   * @param rv_new
   */ 
  private void updateFeatureRelationshipRank(final GFFStreamFeature feature,
                                             final RangeVector rv_new)
  {   
    // update feature_relationship.rank
    ChadoTransaction tsn;
    Hashtable feature_relationship_rank_store = new Hashtable();
    Qualifier qualifier_relation = feature.getQualifierByName("Parent");
    
    for(int rank=0; rank<rv_new.size(); rank++)
    {
      Range range   = (Range)rv_new.elementAt(rank);
      String seq_id = feature.getSegmentID(range);
      
      org.gmod.schema.sequence.Feature chado_feature = 
               new org.gmod.schema.sequence.Feature();
      chado_feature.setUniqueName(seq_id);

      List featureRelationshipsForSubjectId = null;
      if(qualifier_relation != null)
      {
        StringVector parents = qualifier_relation.getValues();
        if(parents.size() > 0)
          featureRelationshipsForSubjectId = new Vector();
        
        for(int i=0; i<parents.size(); i++)
        {
          org.gmod.schema.sequence.Feature parent =
              new org.gmod.schema.sequence.Feature();
          parent.setUniqueName((String)parents.get(i));
          FeatureRelationship feature_relationship =
              new FeatureRelationship();
          CvTerm cvterm = new CvTerm();
          cvterm.setCvTermId(DatabaseDocument.getCvtermID("part_of").intValue());
          
          feature_relationship.setFeatureByObjectId(parent);
          feature_relationship.setFeatureBySubjectId(chado_feature);
          feature_relationship.setCvTerm(cvterm);
          feature_relationship.setRank(rank);
          featureRelationshipsForSubjectId.add(feature_relationship);
          
          tsn = new ChadoTransaction(ChadoTransaction.UPDATE,
              feature_relationship,
              feature.getLastModified(), feature);
          sql.add(tsn);
        }
      }

      feature_relationship_rank_store.put(seq_id, new Integer(rank));
    }
    feature.setFeature_relationship_rank_store(feature_relationship_rank_store);  
  }
  
  /**
   * Create transaction for inserting a feature.
   * @param feature
   */
  private void insertFeature(final Feature feature)
  {
    String feature_uniquename = null;
    try
    {
      Qualifier qualifier_uniquename = feature.getQualifierByName("ID");

      if(qualifier_uniquename != null)
      {
        feature_uniquename = (String)(qualifier_uniquename.getValues()).elementAt(0);
        Splash.logger4j.debug("FEATURE_ADDED "+feature_uniquename);
      }
      
      while(feature_uniquename == null ||
            feature_uniquename.equals("") ||
            feature_uniquename.equals("to_be_set"))
      {
        feature_uniquename = JOptionPane.showInputDialog(null,
                             "Provide a systematic_id : ",
                             "systematic_id missing in "+
                             feature.getIDString(),
                             JOptionPane.QUESTION_MESSAGE).trim();
      }
      feature.setQualifier(new Qualifier("ID", feature_uniquename));
    }
    catch(EntryInformationException eie)
    {
      eie.printStackTrace();
    }
    catch(ReadOnlyException roe)
    {
      roe.printStackTrace();
    }

    
    FeatureLoc featureloc = getFeatureLoc(
                   (GFFStreamFeature)feature.getEmblFeature(),
                   feature_uniquename, 
                  feature.getLocation().getTotalRange());
    org.gmod.schema.sequence.Feature chado_feature = 
              new org.gmod.schema.sequence.Feature();
    chado_feature.setFeatureLoc(featureloc);
    
    try
    { 
      // relationship attributes
      Qualifier qualifier_relation = feature.getQualifierByName("Parent");
      List featureRelationshipsForSubjectId = null;
      if(qualifier_relation != null)
      {
        StringVector parents = qualifier_relation.getValues();
        if(parents.size() > 0)
          featureRelationshipsForSubjectId = new Vector();
        
        for(int i=0; i<parents.size(); i++)
        {
          org.gmod.schema.sequence.Feature parent =
              new org.gmod.schema.sequence.Feature();
          parent.setUniqueName((String)parents.get(i));
          FeatureRelationship feature_relationship =
              new FeatureRelationship();
          CvTerm cvterm = new CvTerm();
          cvterm.setCvTermId(DatabaseDocument.getCvtermID("part_of").intValue());
          
          feature_relationship.setFeatureByObjectId(parent);
          feature_relationship.setFeatureBySubjectId(chado_feature);
          feature_relationship.setCvTerm(cvterm);
          featureRelationshipsForSubjectId.add(feature_relationship);
        }
      }
      
      qualifier_relation = feature.getQualifierByName("Derives_from");
      if(qualifier_relation != null)
      {
        StringVector derives = qualifier_relation.getValues();
        if(derives.size() > 0 && featureRelationshipsForSubjectId == null)
          featureRelationshipsForSubjectId = new Vector();
        
        for(int i=0; i<derives.size(); i++)
        {
          org.gmod.schema.sequence.Feature parent =
                                      new org.gmod.schema.sequence.Feature();
          parent.setUniqueName((String) derives.get(i));
          FeatureRelationship feature_relationship = new FeatureRelationship();
          CvTerm cvterm = new CvTerm();
          cvterm.setCvTermId(DatabaseDocument.getCvtermID("derives_from")
              .intValue());

          feature_relationship.setFeatureByObjectId(parent);
          feature_relationship.setFeatureBySubjectId(chado_feature);
          feature_relationship.setCvTerm(cvterm);
          featureRelationshipsForSubjectId.add(feature_relationship);
        }
      }
      chado_feature.setFeatureRelationshipsForSubjectId(
                       featureRelationshipsForSubjectId);
    }
    catch(InvalidRelationException ire){}
    
    chado_feature.setUniqueName(feature_uniquename);
    chado_feature.setName(feature_uniquename);

    String key = feature.getKey().toString();
    
    CvTerm cvterm = new CvTerm();
    cvterm.setCvTermId(DatabaseDocument.getCvtermID(key).intValue());
    chado_feature.setCvTerm(cvterm);

    addQualifiers(feature.getQualifiers(), chado_feature);
    // create transaction object
    
    ChadoTransaction tsn = new ChadoTransaction(ChadoTransaction.INSERT,
                               chado_feature,
                               null, (GFFStreamFeature)null);
    sql.add(tsn);  
  }
  
  /**
   * Create transaction for inserting a feature.
   * @param feature
   */
  private void insertFeatureSegment(final FeatureSegment segment,
                                    final String segment_uniquename)
  {
    org.gmod.schema.sequence.Feature chado_feature = 
      new org.gmod.schema.sequence.Feature();
    FeatureLoc featureloc = new FeatureLoc();
    chado_feature.setFeatureLoc(featureloc);
    
    if(segment.isForwardSegment())
      featureloc.setStrand(new Short((short)1));
    else
      featureloc.setStrand(new Short((short)-1));
    
    // codon_start attribute
    Feature feature = segment.getFeature();
    try
    {
      Qualifier qualifier_phase = feature.getQualifierByName("codon_start");
      if(qualifier_phase != null)
      {
        String phase = (String)(qualifier_phase.getValues()).elementAt(0);

        if(phase.equals ("1"))
          featureloc.setPhase(new Integer(0));
        else if(phase.equals("2"))
          featureloc.setPhase(new Integer(1));
        else if(phase.equals("3")) 
          featureloc.setPhase(new Integer(2));
      }
      else
        featureloc.setPhase(null);
      
      // relationship attributes
      Qualifier qualifier_relation = feature.getQualifierByName("Parent");
      List featureRelationshipsForSubjectId = null;
      if(qualifier_relation != null)
      {
        StringVector parents = qualifier_relation.getValues();
        if(parents.size() > 0)
          featureRelationshipsForSubjectId = new Vector();
        
        for(int i=0; i<parents.size(); i++)
        {
          org.gmod.schema.sequence.Feature parent =
              new org.gmod.schema.sequence.Feature();
          parent.setUniqueName((String)parents.get(i));
          FeatureRelationship feature_relationship =
              new FeatureRelationship();
          CvTerm cvterm = new CvTerm();
          cvterm.setCvTermId(DatabaseDocument.getCvtermID("part_of").intValue());
          
          feature_relationship.setFeatureByObjectId(parent);
          feature_relationship.setFeatureBySubjectId(chado_feature);
          feature_relationship.setCvTerm(cvterm);
          featureRelationshipsForSubjectId.add(feature_relationship);
        }
      }
      
      qualifier_relation = feature.getQualifierByName("Derives_from");
      if(qualifier_relation != null)
      {
        StringVector derives = qualifier_relation.getValues();
        if(derives.size() > 0 && featureRelationshipsForSubjectId == null)
          featureRelationshipsForSubjectId = new Vector();
        
        for(int i=0; i<derives.size(); i++)
        {
          org.gmod.schema.sequence.Feature parent =
                                      new org.gmod.schema.sequence.Feature();
          parent.setUniqueName((String) derives.get(i));
          FeatureRelationship feature_relationship = new FeatureRelationship();
          CvTerm cvterm = new CvTerm();
          cvterm.setCvTermId(DatabaseDocument.getCvtermID("derives_from")
              .intValue());

          feature_relationship.setFeatureByObjectId(parent);
          feature_relationship.setFeatureBySubjectId(chado_feature);
          feature_relationship.setCvTerm(cvterm);
          featureRelationshipsForSubjectId.add(feature_relationship);
        }
      }
      chado_feature.setFeatureRelationshipsForSubjectId(
                       featureRelationshipsForSubjectId);
    }
    catch(InvalidRelationException ire){}

    featureloc.setFmin(new Integer(segment.getRawRange().getStart()-1));
    featureloc.setFmax(new Integer(segment.getRawRange().getEnd()));
    chado_feature.setUniqueName(segment_uniquename);
    chado_feature.setName(segment_uniquename);

    String key = feature.getKey().toString();
    
    CvTerm cvterm = new CvTerm();
    cvterm.setCvTermId(DatabaseDocument.getCvtermID(key).intValue());
    chado_feature.setCvTerm(cvterm);

    //addQualifiers(feature.getQualifiers(), chado_feature);
    // create transaction object
    
    ChadoTransaction tsn = new ChadoTransaction(ChadoTransaction.INSERT,
        chado_feature,
        null, (GFFStreamFeature)null);
   
    sql.add(tsn);  
  }
  
  /**
   * Set the transaction for deleting a feature.
   */
  private void deleteFeature(final String uniquename)
  {
    org.gmod.schema.sequence.Feature chado_feature = 
      new org.gmod.schema.sequence.Feature();
    chado_feature.setUniqueName(uniquename);
    
    ChadoTransaction tsn = new ChadoTransaction(ChadoTransaction.DELETE,
        chado_feature,
        null, (GFFStreamFeature)null);

    sql.add(tsn); 
  }

  /**
   * Add qualifiers that are in a <code>QualifierVector</code> to a 
   * <code>ChadoFeature</code>.
   * @param qualifiers		the <code>QualifierVector</code>
   * @param chado_feature	the <code>ChadoFeature</code>
   */
  private void addQualifiers(final QualifierVector qualifiers,
                             final org.gmod.schema.sequence.Feature chado_feature)
  {
    // add qualifiers/attributes
    for(int qualifier_index = 0; qualifier_index < qualifiers.size();
      ++qualifier_index)
    {
      final Qualifier this_qualifier = (Qualifier)qualifiers.elementAt(qualifier_index);
      final String name = this_qualifier.getName();

      // ignore reserved tags
      if(isReservedTag(name) || isSynonymTag(name))
        continue;

      final StringVector qualifier_values = this_qualifier.getValues();
      
      try
      {
        int type_id = DatabaseDocument.getCvtermID( name ).intValue();
        for(int value_index = 0; value_index < qualifier_values.size();
          ++value_index)
        {
 //         if(qualifier_values.elementAt(value_index) instanceof ChadoFeatureProp)
 //         {
 //           System.out.println("HERE ChadoFeatureProp");
 //           chado_feature.addQualifier(type_id, 
 //                        (ChadoFeatureProp)qualifier_values.elementAt(value_index));
 //         }
 //         else
 //         {
            // happens when duplicating features 
            FeatureProp featureprop = new FeatureProp();
            featureprop.setValue((String)qualifier_values.elementAt(value_index));
            CvTerm cvTerm = new CvTerm();
            cvTerm.setCvTermId(type_id);
            featureprop.setCvTerm(cvTerm);
            chado_feature.addFeatureProp(featureprop);
 //           chado_feature.addQualifier(type_id, featureprop);
 //         }
        }
      }
      catch(NullPointerException npe)
      {
        JOptionPane.showMessageDialog(null,
            name+" is not a valid qualifier!",
            "Invalid Qualifier",
            JOptionPane.WARNING_MESSAGE);
      }
    } 
  }

  /**
   * Determine if this is a GFF3 predefined tag.
   * @param tag
   * @return  true if the tag is a GFF3 predefined tag
   */
  private boolean isReservedTag(final String tag)
  {
    for(int i=0; i<reserved_tags.length; i++)
      if(tag.equals(reserved_tags[i]))
        return true;
    return false;
  }
  
  /**
   * Determine if this is a controlled vocabulary tag, e.g GO.
   * @param tag
   * @return  true if the tag is a CV tag
   */
  private boolean isCvTag(final String tag)
  {
    for(int i=0; i<cv_tags.length; i++)
      if(tag.equals(cv_tags[i]))
        return true;
    return false;
  }
  
  /**
   * Determine if this is a GFF3 predefined tag.
   * @param tag
   * @return  true if the tag is a GFF3 predefined tag
   */
  private boolean isSynonymTag(final String tag)
  {
    if(synonym_tags == null)
    {
      synonym_tags = DatabaseDocument.getSynonymTypeNames(SYNONYM_TAG_CVNAME);
      if(synonym_tags == null || synonym_tags.length < 1)
      {
        Splash.logger4j.debug("Using default synonym names");
        synonym_tags = new String[6];
        synonym_tags[0] = "synonym";
        synonym_tags[1] = "gene";
        synonym_tags[2] = "systematic_id";
        synonym_tags[3] = "primary_name";
        synonym_tags[4] = "reserved_name";
        synonym_tags[5] = "obsolete_name";
      }
    }
    
    for(int i=0; i<synonym_tags.length; i++)
      if(tag.equals(synonym_tags[i]))
        return true;
    return false;
  }
  
  /**
   * Compare the old and new keys and qualifiers and find the qualifiers 
   * that have changed or been added and UPDATE, INSERT or DELETE accordingly.
   * @param qualifiers_old	old qualifiers
   * @param qualifiers_new	new qualifiers
   * @param feature		GFF feature that has been changed
   */
  private void editKeyAndQualifiers(final QualifierVector qualifiers_old, 
                                    final QualifierVector qualifiers_new, 
                                    final Key old_key, 
                                    final Key new_key,
                                    final GFFStreamFeature feature,
                                    final int event_type)
  {
    String uniquename = (String)(feature.getQualifierByName("ID").getValues()).elementAt(0);
    ChadoTransaction tsn;

    // updating the key unless just a qualifier changed
    if(event_type != FeatureChangeEvent.QUALIFIER_CHANGED && 
       !new_key.equals(old_key))
    {
      Integer lcvterm_id = DatabaseDocument.getCvtermID(new_key.getKeyString());
      if(lcvterm_id == null)   // chado doesn't recognise this
      {
        JOptionPane.showMessageDialog(null, 
                  new_key.getKeyString()+" is not a valid key!\n"+
                  "There is no CV term set for this key.",
                  "Invalid Feature Key",
                  JOptionPane.WARNING_MESSAGE);
      }
      else
      {  
        org.gmod.schema.sequence.Feature chado_feature =
           new org.gmod.schema.sequence.Feature();
        CvTerm cvterm = new CvTerm();
        cvterm.setCvTermId( lcvterm_id.intValue() );
        
        chado_feature.setCvTerm(cvterm);
        chado_feature.setUniqueName(uniquename);
        
        tsn = new ChadoTransaction(ChadoTransaction.UPDATE,
            chado_feature,
            feature.getLastModified(), feature);

        sql.add(tsn);
      }
    }
    
    // look for qualifiers to DELETE
    for(int qualifier_index = 0; qualifier_index < qualifiers_old.size();
        ++qualifier_index)
    {
      final Qualifier this_qualifier = (Qualifier)qualifiers_old.elementAt(qualifier_index);
      String name = this_qualifier.getName();
      
      if(!qualifiers_new.contains(name))
      {
        if(isReservedTag(name) || isSynonymTag(name) || isCvTag(name))
        {
          handleReservedTags(feature, uniquename, 
              null,
              this_qualifier);
          continue;
        }
        
        // get the cvterm_id for this featureprop/qualifier
        Integer lcvterm_id = DatabaseDocument.getCvtermID(name);
        if(lcvterm_id == null)   // chado doesn't recognise this
        {
          JOptionPane.showMessageDialog(null,
                      name+" is not a valid qualifier!",
                      "Invalid Qualifier",
                      JOptionPane.WARNING_MESSAGE);
          continue;
        }

        processFeatureProp(feature, null, -1, ChadoTransaction.DELETE, 
                           uniquename, lcvterm_id);
      }
    }

    // look for qualifiers to INSERT/UPDATE
    for(int qualifier_index = 0; qualifier_index < qualifiers_new.size();
        ++qualifier_index)
    {
      final Qualifier this_qualifier = (Qualifier)qualifiers_new.elementAt(qualifier_index);
      String name = this_qualifier.getName();          
      int old_index = qualifiers_old.indexOfQualifierWithName(name);

      Qualifier this_old_qualifier = null;
      StringVector old_qualifier_strings = null;
      final StringVector new_qualifier_strings =
                   StreamQualifier.toStringVector(null, this_qualifier);

      if(old_index> -1)  // update qualifier
      {
        this_old_qualifier = (Qualifier)qualifiers_old.elementAt(old_index);

        old_qualifier_strings =
                   StreamQualifier.toStringVector(null, this_old_qualifier);

        // check if anything has changed for this qualifier name
        boolean need_to_update = false;
        for(int value_index = 0; value_index < new_qualifier_strings.size();
            ++value_index)
        {
          String qualifier_string = (String)new_qualifier_strings.elementAt(value_index);
          if(!old_qualifier_strings.contains(qualifier_string))
            need_to_update = true;
        }
        if(!need_to_update &&
            new_qualifier_strings.size() == old_qualifier_strings.size())
          continue;
      }

      if(isReservedTag(name) || isSynonymTag(name) || isCvTag(name))
      {
        handleReservedTags(feature, uniquename, 
                           this_qualifier,
                           qualifiers_old.getQualifierByName(name));
        continue;
      }
      
      // get the cvterm_id for this featureprop/qualifier
      Integer lcvterm_id = null;
      if(!name.equals("timelastmodified"))
      {
        lcvterm_id = DatabaseDocument.getCvtermID(name);

        if(lcvterm_id == null)   // chado doesn't recognise this
        {
          JOptionPane.showMessageDialog(null, 
                    name+" is not a valid qualifier!\n"+
                    "There is no CV term set for this qualifier.",
                    "Invalid Qualifier",
                    JOptionPane.WARNING_MESSAGE);
          continue;
        }
      }

      if(old_index > -1 &&
         new_qualifier_strings.size() == old_qualifier_strings.size())
      {
        //  
        // UPDATE existing featureprop's
        //
        for(int rank = 0; rank < new_qualifier_strings.size();
            ++rank)
        {         
          String qualifier_string = (String)new_qualifier_strings.elementAt(rank);
          int index = qualifier_string.indexOf("=");

          if(index > -1)
            qualifier_string = qualifier_string.substring(index+1);
          
          processFeatureProp(feature, qualifier_string, rank, ChadoTransaction.UPDATE, 
                             uniquename, lcvterm_id);
        }

      }
      else
      {       
        //
        // DELETE any existing featureprops
        //
        if(old_index > -1)
        {
          processFeatureProp(feature, null, -1, ChadoTransaction.DELETE, 
                             uniquename, lcvterm_id);
        }
          
        //
        // INSERT new featureprops
        //
        for(int rank = 0; rank < new_qualifier_strings.size();
            ++rank)
        {
          String qualifier_string = (String)new_qualifier_strings.elementAt(rank);
          int index = qualifier_string.indexOf("=");
          if(index > -1)
            qualifier_string = qualifier_string.substring(index+1);
         
          processFeatureProp(feature, qualifier_string, rank, ChadoTransaction.INSERT, 
                             uniquename, lcvterm_id);
        }

      }
    }

  }
  
  /**
   * 
   * @param feature
   * @param type        ChadoTransaction type DELETE/UPDATE/INSERT
   * @param uniquename
   * @param lcvterm_id
   */
  private void processFeatureProp(final GFFStreamFeature feature,
      final String qualifier_string, final int rank,
      final int type, String uniquename, Integer lcvterm_id)
  {
    ChadoTransaction tsn;
    if(feature.getFeature_relationship_rank_store() != null)
    {
      Hashtable rank_hash = feature.getFeature_relationship_rank_store();
      Enumeration id_keys= rank_hash.keys();
      while(id_keys.hasMoreElements())
      {
        uniquename = (String)id_keys.nextElement();
        FeatureProp featureprop = getFeatureProp(uniquename, qualifier_string,
                                                 lcvterm_id, rank);
        tsn = new ChadoTransaction(type,
            featureprop,
            feature.getLastModified(), feature);
        
        tsn.setUniquename(uniquename);
        sql.add(tsn);
      }
    }
    else
    {
      FeatureProp featureprop = getFeatureProp(uniquename,
                         qualifier_string, lcvterm_id, rank);
    
      tsn = new ChadoTransaction(type,
          featureprop,
          feature.getLastModified(), feature);       
      sql.add(tsn);
    }
  }
  
  /**
   * Handle database transactions involving the GFF3 reserved tags.
   * @param feature         the <code>GFFStreamFeature</code>
   * @param type            the transaction type (INSERT/UPDATE/DELETE)
   * @param new_qualifier   the new qualifier
   * @param old_qualifier   the old qualifier
   */
  private void handleReservedTags(final GFFStreamFeature feature,
                                  String uniquename,
                                  final Qualifier new_qualifier,
                                  final Qualifier old_qualifier)
  {  
    StringVector new_qualifier_strings = null;
    
    if(new_qualifier != null)
      new_qualifier_strings = StreamQualifier.toStringVector(null, new_qualifier);
    
    StringVector old_qualifier_strings;
    
    if(old_qualifier != null)
      old_qualifier_strings = StreamQualifier.toStringVector(null, old_qualifier);
    else
      old_qualifier_strings = new StringVector();
    
    final String qualifier_name;
    
    if(old_qualifier != null)
      qualifier_name = old_qualifier.getName();
    else
      qualifier_name = new_qualifier.getName();
    
    if(qualifier_name.equals("ID"))
    { 
      // this shouldn't be possible
      if(new_qualifier.getValues() == null)
        return;
      
      org.gmod.schema.sequence.Feature chado_feature =
        new org.gmod.schema.sequence.Feature();
     
      chado_feature.setUniqueName((String)new_qualifier.getValues().get(0));
     
      Splash.logger4j.debug(uniquename+"  in handleReservedTags() NEW="+
          (String)new_qualifier.getValues().get(0)+" OLD="+
          (String)old_qualifier.getValues().get(0));
      ChadoTransaction tsn = new ChadoTransaction(ChadoTransaction.UPDATE,
                chado_feature,
                feature.getLastModified(), feature);
      tsn.setOldUniquename( (String)old_qualifier.getValues().get(0) );
     
      sql.add(tsn);
      return;
    }
    
    ChadoTransaction tsn = null;
    // find tags that have been deleted
    for(int i = 0; i < old_qualifier_strings.size(); ++i)
    {
      String qualifier_string = (String)old_qualifier_strings.elementAt(i);
      
      if( new_qualifier_strings == null ||
         !new_qualifier_strings.contains(qualifier_string) )
      {
         int index = qualifier_string.indexOf("=");
         qualifier_string = qualifier_string.substring(index+1);
         
         if(qualifier_name.equals("Dbxref"))
         {
           Splash.logger4j.debug(uniquename+"  in handleReservedTags() DELETE db="+
               qualifier_string.substring(0,index)+" acc="+qualifier_string.substring(index+1));
         
           FeatureDbXRef old_dbxref = getFeatureDbXRef(qualifier_string,
                                                       uniquename);
           
           tsn = new ChadoTransaction(ChadoTransaction.DELETE,
               old_dbxref,
               feature.getLastModified(), feature);
           sql.add(tsn);
         }
         else if(qualifier_name.equals("codon_start"))
         {
           FeatureLoc featureloc = getFeatureLoc(feature, uniquename, 
               feature.getLocation().getTotalRange());
           
           tsn = new ChadoTransaction(ChadoTransaction.UPDATE,
                                      featureloc,
                                      feature.getLastModified(), feature);
           sql.add(tsn);
         }
         else if(isCvTag(qualifier_name))
         {
           /*Qualifier qual = feature.getQualifierByName(qualifier_name);
           StringVector values = qual.getValues();
           int beginIndex = qualifier_string.indexOf("term=")+5;
           int endIndex   = qualifier_string.indexOf(";",beginIndex);
           final String thisTerm;
           if(endIndex > -1)
             thisTerm = qualifier_string.substring(beginIndex, endIndex);
           else
             thisTerm = qualifier_string.substring(beginIndex);*/
           
           Splash.logger4j.debug(uniquename+"  in handleReservedTags() DELETE "+
               qualifier_name+" "+qualifier_string);
           
           FeatureCvTerm feature_cvterm = getFeatureCvTerm(qualifier_name, qualifier_string, 
                                                           uniquename);
           /*Vector rankables = null;
           int rank = 0;
           for(int j=0; j<values.size(); j++)
           {
             String val = (String)values.get(j);
             if(val.indexOf(thisTerm) > -1 &&
                !val.equals(qualifier_string))
             {
               if(rankables == null)
                 rankables = new Vector();
               
               FeatureCvTerm fc = getFeatureCvTerm(qualifier_name, val, 
                                                   uniquename);
               fc.setRank(rank);
               rankables.add(fc);
               rank++;
               System.out.println("UPDATE ------> "+values.get(j));
             }
           }*/
           
           tsn = new ChadoTransaction(ChadoTransaction.DELETE,
                                      feature_cvterm,
                                      feature.getLastModified(), feature);
           sql.add(tsn);
           
           /*if(rankables != null)
           {
             for(int j=0; j<rankables.size(); j++)
             {
               feature_cvterm = (FeatureCvTerm)rankables.get(j);
               tsn = new ChadoTransaction(ChadoTransaction.UPDATE,
                 feature_cvterm,
                 feature.getLastModified(), feature);
             }
           }*/
         }
         else if(isSynonymTag(qualifier_name))
         {
           Splash.logger4j.debug(uniquename+"  in handleReservedTags() DELETE "+qualifier_name+" "+
                              qualifier_string);
           
           FeatureSynonym feature_synonym = getFeatureSynonym(qualifier_name,
                                               qualifier_string, uniquename);
          
           tsn = new ChadoTransaction(ChadoTransaction.DELETE,
               feature_synonym,
               feature.getLastModified(), feature);
           sql.add(tsn);
         }
         
      }
    }
    
    if(new_qualifier_strings == null)
      return;
    
    // find tags that have been inserted
    for(int i = 0; i < new_qualifier_strings.size(); ++i)
    {
      String qualifier_string = (String)new_qualifier_strings.elementAt(i);
      if(!old_qualifier_strings.contains(qualifier_string))
      {
         int index = qualifier_string.indexOf("=");
         qualifier_string = qualifier_string.substring(index+1);
         
         if(qualifier_name.equals("Dbxref"))
         {   
           Splash.logger4j.debug(uniquename+"  in handleReservedTags() INSERT db="+
             qualifier_string.substring(0,index)+" acc="+qualifier_string.substring(index+1));
           FeatureDbXRef new_dbxref = getFeatureDbXRef(qualifier_string,
                                                       uniquename);
           
           tsn = new ChadoTransaction(ChadoTransaction.INSERT,
               new_dbxref,
               feature.getLastModified(), feature);
         }
         else if(qualifier_name.equals("codon_start"))
         {
           FeatureLoc featureloc = getFeatureLoc(feature, uniquename, 
               feature.getLocation().getTotalRange());
           
           tsn = new ChadoTransaction(ChadoTransaction.UPDATE,
                                      featureloc,
                                      feature.getLastModified(), feature);
           //sql.add(tsn);
         }
         else if(isCvTag(qualifier_name))
         {
           Splash.logger4j.debug(uniquename+"  in handleReservedTags() INSERT "+
                                 qualifier_name+" "+qualifier_string);
           FeatureCvTerm feature_cvterm = getFeatureCvTerm(qualifier_name,
              qualifier_string, uniquename);
           tsn = new ChadoTransaction(ChadoTransaction.INSERT, 
                      feature_cvterm,
                      feature.getLastModified(), feature);
         }
         else if(isSynonymTag(qualifier_name))
         {
           Splash.logger4j.debug(uniquename+"  in handleReservedTags() INSERT "+
                                 qualifier_name+" "+qualifier_string);

           FeatureSynonym feature_synonym = getFeatureSynonym(qualifier_name,
                                               qualifier_string, uniquename);

           tsn = new ChadoTransaction(ChadoTransaction.INSERT,
               feature_synonym,
               feature.getLastModified(), feature);
         }
         sql.add(tsn);
      }
    }  
    
  }
  
  
  /**
   * Strip out double quotes around a string.
   * @param s a <code>String</code> to strip quotes
   * @return the resulting <code>String</code>
   */
  private String stripQuotes(String s)
  {
    if(s == null)
      return null;
    
    if(s.startsWith("\"") && s.endsWith("\""))
      s = s.substring(1,s.length()-1);
    
    return s;
  }

  /**
   * Get the FeatureLoc object
   * @param feature
   * @param seg_id
   * @param range_new
   * @return
   */
  private FeatureLoc getFeatureLoc(GFFStreamFeature feature,
                             String seg_id, Range range_new)
  {
    FeatureLoc featureloc = new FeatureLoc();
    org.gmod.schema.sequence.Feature chado_feature = 
      new org.gmod.schema.sequence.Feature();
    chado_feature.setUniqueName(seg_id);
    
    featureloc.setFeatureByFeatureId(chado_feature);
    featureloc.setFmax(new Integer(range_new.getEnd()));
    featureloc.setFmin(new Integer(range_new.getStart()-1));
    
    if(feature.getFeature_relationship_rank_store() != null)
    {
      Hashtable rank_hash = feature.getFeature_relationship_rank_store();
      if(rank_hash.containsKey(seg_id))
      {
        Integer rank = (Integer)rank_hash.get(seg_id);
        featureloc.setRank(rank.intValue());
      }
      else
        featureloc.setRank(0);
    }
    else
      featureloc.setRank(0);
    
    boolean is_complement = feature.getLocation().isComplement();
    if(is_complement)
      featureloc.setStrand(new Short((short) -1));
    else
      featureloc.setStrand(new Short((short) 1));
    
    Qualifier qualifier_phase = feature.getQualifierByName("codon_start");
    if(qualifier_phase != null)
    {
      String phase = (String)(qualifier_phase.getValues()).elementAt(0);

      if(phase.equals ("1"))
        featureloc.setPhase(new Integer(0));
      else if(phase.equals("2"))
        featureloc.setPhase(new Integer(1));
      else if(phase.equals("3")) 
        featureloc.setPhase(new Integer(2));
    }
    else
      featureloc.setPhase(null);
    
    return featureloc;
  }
  
  /**
   * Create the <code>FeatureProp</code> object
   * @param uniquename
   * @param qualifier_string
   * @param lcvterm_id
   * @param rank
   * @return
   */
  private FeatureProp getFeatureProp(final String uniquename,
                                     final String qualifier_string,
                                     final Integer lcvterm_id, 
                                     final int rank)
  {
    FeatureProp featureprop = new FeatureProp();
    org.gmod.schema.sequence.Feature chado_feature   = 
      new org.gmod.schema.sequence.Feature();
    chado_feature.setUniqueName(uniquename);
    CvTerm cvterm = new CvTerm();
    cvterm.setCvTermId(lcvterm_id.intValue());
    featureprop.setValue(stripQuotes(qualifier_string));
    featureprop.setRank(rank);
    featureprop.setCvTerm(cvterm);
    featureprop.setFeature(chado_feature);
    return featureprop;
  }
  
  /**
   * Create the <code>FeatureSynonym</code> object
   * @param qualifier_name
   * @param qualifier_string
   * @param uniqueName
   * @return
   */
  private FeatureSynonym getFeatureSynonym(final String qualifier_name,
                                           final String qualifier_string,
                                           final String uniqueName)
  {
    Integer lcvterm_id = DatabaseDocument.getCvtermID(qualifier_name);
    FeatureSynonym feature_synonym = new FeatureSynonym();
    org.gmod.schema.sequence.Feature chado_feature = 
      new org.gmod.schema.sequence.Feature();
    chado_feature.setUniqueName(uniqueName);
    
    Synonym synonym = new Synonym();
    CvTerm cvterm = new CvTerm();
    cvterm.setCvTermId(lcvterm_id.intValue());
    synonym.setName(qualifier_string);
    synonym.setCvTerm(cvterm);
    
    feature_synonym.setSynonym(synonym);
    feature_synonym.setFeature(chado_feature); 
    return feature_synonym;
  }
  
  /**
   * Create the <code>FeatureDbXRef</code> object
   * @param qualifier_string
   * @param uniqueName
   * @return
   */
  private FeatureDbXRef getFeatureDbXRef(final String qualifier_string,
                                         final String uniqueName)
  {
    int index = qualifier_string.lastIndexOf(":");
    FeatureDbXRef feature_dbxref = new FeatureDbXRef();
    DbXRef dbxref = new DbXRef();
    Db db = new Db();
    db.setName(qualifier_string.substring(0,index));
    dbxref.setDb(db);
    dbxref.setAccession(qualifier_string.substring(index+1));
    feature_dbxref.setDbXRef(dbxref);
    org.gmod.schema.sequence.Feature feat = 
      new org.gmod.schema.sequence.Feature();
    feat.setUniqueName(uniqueName);
    feature_dbxref.setFeature(feat);
    return feature_dbxref;
  }
  
  /**
   * Create the <code>FeatureCvTerm</code> object
   * @param qualifier_name
   * @param qualifier_string
   * @param uniqueName
   * @return
   */
  private FeatureCvTerm getFeatureCvTerm(final String qualifier_name,
                                         String qualifier_string,
                                         final String uniqueName)
  {
    Splash.logger4j.debug("Build FeatureCvTerm for "+uniqueName);
    
    if(qualifier_string.startsWith("\""))
      qualifier_string = qualifier_string.substring(1,qualifier_string.length()-1);
    
    FeatureCvTerm feature_cvterm = new FeatureCvTerm();
    org.gmod.schema.sequence.Feature chado_feature =
      new org.gmod.schema.sequence.Feature();
    chado_feature.setUniqueName(uniqueName);
    feature_cvterm.setFeature(chado_feature);
    
    if(qualifier_name.toLowerCase().equals("product"))
    {
      CvTerm cvTerm = getCvTerm(qualifier_string);

      feature_cvterm.setCvTerm(cvTerm);
      Splash.logger4j.debug("Finished building FeatureCvTerm for "+uniqueName);
      return feature_cvterm;
    }
    
    List featureCvTermProps = new Vector();
    StringVector strings = StringVector.getStrings(qualifier_string, ";");
    for(int i=0; i<strings.size(); i++)
    {    
      String this_qualifier_part = ((String)strings.get(i)).trim();
      String this_qualifier_part_lowercase = this_qualifier_part.toLowerCase();
      
      if(this_qualifier_part_lowercase.startsWith("term="))
      {
        String cvTermName = this_qualifier_part.substring(5);
        CvTerm cvTerm = getCvTerm(cvTermName);

        feature_cvterm.setCvTerm(cvTerm);
        continue;
      }

      if(this_qualifier_part_lowercase.startsWith("cv="))
        continue;
        
      // the WITH column is associated with one or more FeatureCvTermDbXRef
      if(this_qualifier_part_lowercase.startsWith("with="))
      {
        String withStr = this_qualifier_part.substring(5);
        loadDbXRefsAndPubs(withStr, feature_cvterm);
        continue;
      }

      if(this_qualifier_part_lowercase.equals("qualifier=not"))
      {
        feature_cvterm.setNot(true);
        continue;
      }
      
      // N.B. 
      // 1) for /GO the db_xref is a Pub (for primary pubs) 
      //    or FeatureCvTermPub (for others) in /GO
      // 2) for /controlled_curation the db_xref is a FeatureCvTermDbXRef 
      //    or a Pub
      if(this_qualifier_part_lowercase.startsWith("db_xref="))
      {
        String dbxrefStr = this_qualifier_part.substring(8);
        loadDbXRefsAndPubs(dbxrefStr, feature_cvterm);
        continue;
      }
      
      // feature_cvterm_prop's  
      
      if(!this_qualifier_part_lowercase.startsWith("goid=") &&
         !this_qualifier_part_lowercase.startsWith("aspect="))
      {
        int index = this_qualifier_part_lowercase.indexOf('=');
        String prop = this_qualifier_part.substring(index+1);
        
        CvTerm cvTerm = getCvTerm(this_qualifier_part.substring(0,index));
        
        Splash.logger4j.debug("FeatureCvTermProp = "+this_qualifier_part_lowercase);
        FeatureCvTermProp featureCvTermProp = new FeatureCvTermProp();
        featureCvTermProp.setValue(prop);
        featureCvTermProp.setCvTerm(cvTerm);   
        featureCvTermProp.setRank(featureCvTermProps.size());
        
        featureCvTermProps.add(featureCvTermProp);
        
        continue;
      }
    }
    
    feature_cvterm.setFeatureCvTermProps(featureCvTermProps);
    
    Splash.logger4j.debug("Finished building FeatureCvTerm for "+uniqueName);
    return feature_cvterm;
  }
  
  private CvTerm getCvTerm(String cvTermName)
  {
    if(cvTermName.startsWith("\""))
      cvTermName = cvTermName.substring(1, cvTermName.length()-1);
    
    CvTerm cvTerm = DatabaseDocument.getCvTermByCvTermName(cvTermName);

    if(cvTerm != null)
    {
      Splash.logger4j.debug("USE CvTerm from cache, CvTermId="
          + cvTermName + "  -> " + cvTerm.getCvTermId()+  " " +
          cvTerm.getName()+":"+cvTerm.getCv().getName()); 
    }
    else
    {
      Splash.logger4j.warn("CvTerm not found in cache = " + cvTermName);
      cvTerm = new CvTerm();
      cvTerm.setName(cvTermName);
    }
    return cvTerm;
  }
  
  /**
   * Use to load feature_cvterm_dbxref's and feature_cvterm_pub's into a
   * feature_cvterm.
   * Note:
   * 1) for /GO the db_xref is a Pub (for primary pubs) 
   *    or FeatureCvTermPub (for others) in /GO
   * 2) for /controlled_curation the db_xref is a FeatureCvTermDbXRef 
   *    or a Pub
   * @param searchStr
   * @param feature_cvterm
   */
  private void loadDbXRefsAndPubs(final String searchStr,
                                  final FeatureCvTerm feature_cvterm)
  {
    List featureCvTermDbXRefs = null;
    
    StringVector strings = StringVector.getStrings(searchStr, "|");

    for(int i=0; i<strings.size(); i++)
    {
      String this_part = (String)strings.get(i);
      int ind = this_part.indexOf(':');
      
      final String dbName = this_part.substring(0, ind);
      final String accession = this_part.substring(ind+1);

      if(dbName.equals("PMID"))
      {
        // enter as Pub if primary
        
        Pub pub = new Pub();
        pub.setUniqueName(dbName + ":" + accession);
        
        if(feature_cvterm.getPub() == null)
        {
          Splash.logger4j.debug("Set primary Pub for " + 
              dbName + ":" + accession);
          feature_cvterm.setPub(pub);
        }
        else
        {
          // secondary pub
          Splash.logger4j.debug("Set secondary Pub for " + 
              dbName + ":" + accession);
          Collection featureCvTermPubs = feature_cvterm.getFeatureCvTermPubs();
          if(featureCvTermPubs == null ||
             featureCvTermPubs.size() < 1)
          {
            featureCvTermPubs = new Vector();
            feature_cvterm.setFeatureCvTermPubs(featureCvTermPubs);
          }
          featureCvTermPubs.add(pub);
        }
      }
      else  
      {
        // enter as feature_cvterm_dbxref
        Splash.logger4j.debug("CREATE FeatureCvTermDbXRef for " + 
            dbName + ":" + accession);
        
        DbXRef dbxref = new DbXRef();
        dbxref.setAccession(accession);
        Db db = new Db();
        db.setName(dbName);
        dbxref.setDb(db);

        FeatureCvTermDbXRef featureCvTermDbXRef = new FeatureCvTermDbXRef();
        featureCvTermDbXRef.setDbXRef(dbxref);
        
        if(featureCvTermDbXRefs == null)
          featureCvTermDbXRefs = new Vector();
        featureCvTermDbXRefs.add(featureCvTermDbXRef);
      }
    }
    
    if(featureCvTermDbXRefs != null)
      feature_cvterm.setFeatureCvTermDbXRefs(featureCvTermDbXRefs);
  }
  
  /**
   *  Invoked when a deletion or insertion occurs in a Bases object.
   **/
  public void sequenceChanged(final SequenceChangeEvent event)
  {
  }

  public Vector getFeatureInsertUpdate()
  {
    Vector features = null;
    
    for(int i=0; i<sql.size(); i++)
    {
      ChadoTransaction transaction = (ChadoTransaction)sql.get(i);
      if(transaction.getType() == ChadoTransaction.INSERT ||
         transaction.getType() == ChadoTransaction.UPDATE)
      {
        if(transaction.getFeatureObject() instanceof 
           org.gmod.schema.sequence.Feature)
        {
          if(features == null)
            features = new Vector();
          org.gmod.schema.sequence.Feature feature =
            (org.gmod.schema.sequence.Feature)transaction.getFeatureObject();
          features.add( feature.getUniqueName() );
        }
      }
    }
    return features;
  }
  
  /**
   * Commit the transactions back to the database.  
   *
   **/
  public void commit(DatabaseDocument dbDoc)
  {
    int retVal = dbDoc.commit(sql);
    if(retVal > 0)
      sql = new Vector();
  }
  
  
}

