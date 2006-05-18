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
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.EntryChangeListener;
import uk.ac.sanger.artemis.EntryChangeEvent;
import java.util.Vector;
import javax.swing.JOptionPane;


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
              "gff_source",      // program or database
              "gff_seqname" };   // seqID of coord system
           
  private String synonym_tags[] =
          {   "synonym",         // the rest are PSU synonyms
              "gene",
              "primary_name",
              "reserved_name",
              "obsolete_name" };


  /**
   *  Implementation of the FeatureChangeListener interface.  We listen for
   *  changes in every feature of every entry in this group.
   **/ 
  public void featureChanged(FeatureChangeEvent event)
  {
    if(event.featureHasChanged())
    {
      final GFFStreamFeature feature = (GFFStreamFeature)event.getFeature().getEmblFeature();

      if(event.getType() == FeatureChangeEvent.LOCATION_CHANGED)
      {
        RangeVector rv_new = event.getNewLocation().getRanges();
        RangeVector rv_old = event.getOldLocation().getRanges();

        System.out.println("LOCATION_CHANGED "+feature.getFirstBase()+".."+feature.getLastBase()+
                           "   new="+rv_new.size()+" old="+rv_old.size());

        if(rv_new.size() != rv_old.size())
          return;
        
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
 
        ChadoTransaction tsn;
        for(int i=0; i<changes.size();i++)
        {
          ichanged = ((Integer)changes.elementAt(i)).intValue();

          Range range_new = (Range)rv_new.elementAt(ichanged);
          Range range_old = (Range)rv_old.elementAt(ichanged);
          String seg_id   = feature.getSegmentID(range_new);
 
          tsn = new ChadoTransaction(ChadoTransaction.UPDATE, 
                                     seg_id, "featureloc", 
                                     feature.getLastModified(),
                                     feature);

          if(range_new.getStart() != range_old.getStart())
            tsn.addProperty("fmin", Integer.toString(range_new.getStart()-1));
          if(range_new.getEnd() != range_old.getEnd())
            tsn.addProperty("fmax", Integer.toString(range_new.getEnd()));

          sql.add(tsn);
        }
      }
      else if(event.getType() == FeatureChangeEvent.QUALIFIER_CHANGED)
      {
        System.out.println("QUALIFIER_CHANGED ");
      }
      else if(event.getType() == FeatureChangeEvent.ALL_CHANGED)
      {
        System.out.println("ALL_CHANGED");
        
        editKeyAndQualifiers(event.getOldQualifiers(),event.getNewQualifiers(),
                             event.getOldKey(), event.getNewKey(),
                             feature);
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
      Feature feature = event.getFeature();
      String feature_uniquename = null;

      try
      {
        Qualifier qualifier_uniquename = feature.getQualifierByName("ID");

        if(qualifier_uniquename != null)
          feature_uniquename = (String)(qualifier_uniquename.getValues()).elementAt(0);
        
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

      ChadoFeature chado_feature = new ChadoFeature();
      ChadoFeatureLoc featureloc = new ChadoFeatureLoc();
      chado_feature.setFeatureloc(featureloc);
      
      if(feature.isForwardFeature())
        featureloc.setStrand(1);
      else
        featureloc.setStrand(-1);

      // codon_start attribute
      try
      {
        Qualifier qualifier_phase = feature.getQualifierByName("codon_start");
        if(qualifier_phase != null)
        {
          String phase = (String)(qualifier_phase.getValues()).elementAt(0);

          if(phase.equals ("1"))
            featureloc.setPhase(0);
          else if(phase.equals("2"))
            featureloc.setPhase(1);
          else if(phase.equals("3")) 
            featureloc.setPhase(2);
        }
      }
      catch(InvalidRelationException ire){}

      if(feature.isForwardFeature())
        featureloc.setStrand(1);
      else
        featureloc.setStrand(-1);

      featureloc.setFmin(feature.getRawFirstBase()-1);
      featureloc.setFmax(feature.getRawLastBase());
      chado_feature.setUniquename(feature_uniquename);
      chado_feature.setName(feature_uniquename);

      String key = feature.getKey().toString();
      
      ChadoCvterm cvterm = new ChadoCvterm();
      cvterm.setId(DatabaseDocument.getCvtermID(key).longValue());
      chado_feature.setCvterm(cvterm);

      addQualifiers(feature.getQualifiers(), chado_feature);
      // create transaction object
      ChadoTransaction tsn = new ChadoTransaction(ChadoTransaction.INSERT_FEATURE,
                                                  chado_feature);
      sql.add(tsn);
    }
    else if(event.getType() == EntryChangeEvent.FEATURE_DELETED)
    {
      Feature feature = event.getFeature();
      GFFStreamFeature feature_gff = (GFFStreamFeature)event.getFeature().getEmblFeature();
      
      try
      {
        Qualifier qualifier_uniquename = feature.getQualifierByName("ID");
        String feature_uniquename = (String)(qualifier_uniquename.getValues()).elementAt(0);
        ChadoTransaction tsn = new ChadoTransaction(ChadoTransaction.DELETE_FEATURE,
                                                    feature_uniquename, "feature",
                                                    feature_gff.getLastModified(),
                                                    feature_gff);
        sql.add(tsn); 
      }
      catch(InvalidRelationException ire)
      {
        ire.printStackTrace();
      }
    }

//  System.out.println(event.getEntry().getName());
  }

  /**
   * Add qualifiers that are in a <code>QualifierVector</code> to a 
   * <code>ChadoFeature</code>.
   * @param qualifiers		the <code>QualifierVector</code>
   * @param chado_feature	the <code>ChadoFeature</code>
   */
  private void addQualifiers(final QualifierVector qualifiers,
                             final ChadoFeature chado_feature)
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
        long type_id = DatabaseDocument.getCvtermID( name ).longValue();
        for(int value_index = 0; value_index < qualifier_values.size();
          ++value_index)
        {
          chado_feature.addQualifier(type_id, 
                         (ChadoFeatureProp)qualifier_values.elementAt(value_index));
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
   * Determine if this is a GFF3 predefined tag.
   * @param tag
   * @return  true if the tag is a GFF3 predefined tag
   */
  private boolean isSynonymTag(final String tag)
  {
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
                                    final GFFStreamFeature feature)
  {
    final String uniquename = (String)(feature.getQualifierByName("ID").getValues()).elementAt(0);
    ChadoTransaction tsn;

    // updating the key
    if(!new_key.equals(old_key))
    {
      Long lcvterm_id = DatabaseDocument.getCvtermID(new_key.getKeyString());
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
        tsn = new ChadoTransaction(ChadoTransaction.UPDATE,
                                   uniquename, "feature",
                                   feature.getLastModified(),
                                   feature);
        tsn.addProperty("type_id", "'"+ lcvterm_id.toString() +"'");
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
        if(isReservedTag(name) || isSynonymTag(name))
        {
          handleReservedTags(feature, uniquename, 
              null,
              this_qualifier);
          continue;
        }
        
        // get the cvterm_id for this featureprop/qualifier
        Long lcvterm_id = DatabaseDocument.getCvtermID(name);
        if(lcvterm_id == null)   // chado doesn't recognise this
        {
          JOptionPane.showMessageDialog(null,
                      name+" is not a valid qualifier!",
                      "Invalid Qualifier",
                      JOptionPane.WARNING_MESSAGE);
          continue;
        }

        String cvterm_id = lcvterm_id.toString();
        tsn = new ChadoTransaction(ChadoTransaction.DELETE,
                                   uniquename, "featureprop",
                                   feature.getLastModified(),
                                   feature);
        tsn.setConstraint("type_id", cvterm_id);
        sql.add(tsn);
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

      if(isReservedTag(name) || isSynonymTag(name))
      {
        handleReservedTags(feature, uniquename, 
                           this_qualifier,
                           qualifiers_old.getQualifierByName(name));
        continue;
      }
      
      // get the cvterm_id for this featureprop/qualifier
      String cvterm_id = null;
      if(!name.equals("timelastmodified"))
      {
        Long lcvterm_id = DatabaseDocument.getCvtermID(name);

        if(lcvterm_id == null)   // chado doesn't recognise this
        {
          JOptionPane.showMessageDialog(null, 
                    name+" is not a valid qualifier!\n"+
                    "There is no CV term set for this qualifier.",
                    "Invalid Qualifier",
                    JOptionPane.WARNING_MESSAGE);
          continue;
        }
        cvterm_id = lcvterm_id.toString();
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
          
          tsn = new ChadoTransaction(ChadoTransaction.UPDATE,
                                     uniquename, "featureprop",
                                     feature.getLastModified(),
                                     feature);

          tsn.addProperty("value", "'"+ stripQuotes(qualifier_string) +"'");
          tsn.setConstraint("featureprop.type_id", "'"+cvterm_id+"'");
          tsn.setConstraint("rank", Integer.toString(rank));
          sql.add(tsn);
        }

      }
      else
      {       
        //
        // DELETE any existing featureprops
        //
        if(old_index > -1)
        {
          tsn = new ChadoTransaction(ChadoTransaction.DELETE,
                                     uniquename, "featureprop",
                                     feature.getLastModified(),
                                     feature);

          tsn.setConstraint("type_id", cvterm_id);
          sql.add(tsn);
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
         
          tsn = new ChadoTransaction(ChadoTransaction.INSERT,
                                     uniquename, "featureprop",
                                     feature.getLastModified(),
                                     feature);

          tsn.addProperty("value", "'"+ stripQuotes(qualifier_string) +"'");
          tsn.addProperty("type_id", "'"+cvterm_id+"'");
          tsn.addProperty("rank", Integer.toString(rank));

          if(tsn !=null)
            sql.add(tsn);
        }

      }
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
                                  final String uniquename,
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
           index = qualifier_string.lastIndexOf(":");
         
           System.out.println(uniquename+"  in handleReservedTags() DELETE db="+
               qualifier_string.substring(0,index)+" acc="+qualifier_string.substring(index+1));
         
           ChadoDbxref old_dbxref = new ChadoDbxref();
           old_dbxref.setName(qualifier_string.substring(0,index));
           old_dbxref.setAccession(qualifier_string.substring(index+1));

           tsn = new ChadoTransaction(ChadoTransaction.DELETE_DBXREF, 
                                      uniquename, old_dbxref,
                                      feature.getLastModified(),
                                      feature);
         }
         else if(isSynonymTag(qualifier_name))
         {
           System.out.println(uniquename+"  in handleReservedTags() DELETE "+qualifier_name+" "+
                              qualifier_string);
           
           ChadoFeatureSynonym alias = new ChadoFeatureSynonym();
           ChadoSynonym synonym = new ChadoSynonym();
           synonym.setName(qualifier_string);
           alias.setSynonym(synonym);
           alias.setUniquename(uniquename);
          
           tsn = new ChadoTransaction(ChadoTransaction.DELETE_ALIAS, 
                                      alias, 
                                      feature);
           //tsn.setConstraint("synonym.name", "'"+qualifier_string+"'");      
         }
         sql.add(tsn);
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
           index = qualifier_string.lastIndexOf(":");
         
           System.out.println(uniquename+"  in handleReservedTags() INSERT db="+
             qualifier_string.substring(0,index)+" acc="+qualifier_string.substring(index+1));
           ChadoDbxref new_dbxref = new ChadoDbxref();
           new_dbxref.setName(qualifier_string.substring(0,index));
           new_dbxref.setAccession(qualifier_string.substring(index+1));

           tsn = new ChadoTransaction(ChadoTransaction.INSERT_DBXREF, 
                                      uniquename, new_dbxref, 
                                      feature.getLastModified(),
                                      feature);
         }
         else if(isSynonymTag(qualifier_name))
         {
           System.out.println(uniquename+"  in handleReservedTags() INSERT "+qualifier_name+" "+
               qualifier_string);
           Long lcvterm_id = DatabaseDocument.getCvtermID(qualifier_name);
           ChadoFeatureSynonym alias = new ChadoFeatureSynonym();
           ChadoSynonym synonym = new ChadoSynonym();
           ChadoCvterm cvterm = new ChadoCvterm();
           cvterm.setId(lcvterm_id.longValue());
           synonym.setName(qualifier_string);
           synonym.setCvterm(cvterm);
           
           alias.setSynonym(synonym);
           alias.setUniquename(uniquename);

           
           tsn = new ChadoTransaction(ChadoTransaction.INSERT_ALIAS, 
                                      alias, 
                                      feature);
           //tsn.setConstraint("synonym.name", "'"+qualifier_string+"'");
           //tsn.addProperty("type_id", lcvterm_id.toString());
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
    if(s.startsWith("\"") && s.endsWith("\""))
      s = s.substring(1,s.length()-1);
    
    return s;
  }


  /**
   *  Invoked when a deletion or insertion occurs in a Bases object.
   **/
  public void sequenceChanged(final SequenceChangeEvent event)
  {
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

