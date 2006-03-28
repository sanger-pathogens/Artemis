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
import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.sequence.SequenceChangeListener;
import uk.ac.sanger.artemis.sequence.SequenceChangeEvent;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.RangeVector;
import uk.ac.sanger.artemis.io.StreamQualifier;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.EntryChangeListener;
import uk.ac.sanger.artemis.EntryChangeEvent;
import java.util.Vector;
import javax.swing.JOptionPane;
import javax.swing.Box;
import javax.swing.JLabel;
import javax.swing.JTextField;

/**
 *
 * Chado transaction manager listens for feature, entry and sequence changes.
 * <code>ChadoTransactionManager</code> creates and tracks the feature insertions,
 * deletions, and changes to commit back to the database.
 *
 **/
public class ChadoTransactionManager
       implements FeatureChangeListener, EntryChangeListener, SequenceChangeListener 
{

  private Vector sql = new Vector();

  /**
   *
   *  Implementation of the FeatureChangeListener interface.  We listen for
   *  changes in every feature of every entry in this group.
   *
   **/ 
  public void featureChanged(FeatureChangeEvent event)
  {
    if(event.featureHasChanged())
    {
      final GFFStreamFeature feature = (GFFStreamFeature)event.getFeature().getEmblFeature();

      if(event.getType() == event.LOCATION_CHANGED)
      {
        System.out.println("LOCATION_CHANGED ");

        RangeVector rv_new = event.getNewLocation().getRanges();
        RangeVector rv_old = event.getOldLocation().getRanges();

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
                                     seg_id, "featureloc");

          if(range_new.getStart() != range_old.getStart())
            tsn.addProperty("fmin", Integer.toString(range_new.getStart()-1));
          if(range_new.getEnd() != range_old.getEnd())
            tsn.addProperty("fmax", Integer.toString(range_new.getEnd()));

          sql.add(tsn);
        }
      }
      else if(event.getType() == event.QUALIFIER_CHANGED)
      {
        System.out.println("QUALIFIER_CHANGED ");
      }
      else if(event.getType() == event.ALL_CHANGED)
      {
        System.out.println("ALL_CHANGED");
        
        editQualifiers(event.getOldQualifiers(),event.getNewQualifiers(),feature);
      }
    }
  }
 
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

      if(feature.isForwardFeature())
        chado_feature.setStrand(1);
      else
        chado_feature.setStrand(-1);

      // codon_start attribute
      try
      {
        Qualifier qualifier_phase = feature.getQualifierByName("codon_start");
        if(qualifier_phase != null)
        {
          String phase = (String)(qualifier_phase.getValues()).elementAt(0);

          if(phase.equals ("1"))
            chado_feature.setPhase(0);
          else if(phase.equals("2"))
            chado_feature.setPhase(1);
          else if(phase.equals("3")) 
            chado_feature.setPhase(2);
        }
      }
      catch(InvalidRelationException ire){}

      if(feature.isForwardFeature())
        chado_feature.setStrand(1);
      else
        chado_feature.setStrand(-1);

      chado_feature.setFmin(feature.getRawFirstBase()-1);
      chado_feature.setFmax(feature.getRawLastBase());
      chado_feature.setUniquename(feature_uniquename);
      chado_feature.setName(feature_uniquename);

      String key = feature.getKey().toString();
      chado_feature.setType_id(DatabaseDocument.getCvtermID(key).longValue());

      addQualifiers(feature.getQualifiers(), chado_feature);
      // create transaction object
      ChadoTransaction tsn = new ChadoTransaction(ChadoTransaction.INSERT_FEATURE,
                                                  chado_feature);
      System.out.println("---->HERE FEATURE_ADDED "+tsn.getChadoFeature().getUniquename());

      sql.add(tsn);
    }
    else if(event.getType() == EntryChangeEvent.FEATURE_DELETED)
    {
      Feature feature = event.getFeature();
      try
      {
        Qualifier qualifier_uniquename = feature.getQualifierByName("ID");
        String feature_uniquename = (String)(qualifier_uniquename.getValues()).elementAt(0);
        ChadoTransaction tsn = new ChadoTransaction(ChadoTransaction.DELETE_FEATURE,
                                                    feature_uniquename, "feature");
        sql.add(tsn); 
      }
      catch(InvalidRelationException ire)
      {
        ire.printStackTrace();
      }
      System.out.println("HERE FEATURE_DELETED");
    }

//  System.out.println(event.getEntry().getName());
  }

  /**
   *
   * Add qualifiers in a <code>QualifierVector</code> to a <code>ChadoFeature</code>.
   * @param qualifiers		the <code>QualifierVector</code>
   * @param chado_feature	the <code>ChadoFeature</code>
   *
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

      if(name.equals("ID") ||
         name.equals("score") ||
         name.equals("gff_source") ||
         name.equals("gff_seqname") )
        continue;

      final StringVector qualifier_values = this_qualifier.getValues();
     
      long type_id = DatabaseDocument.getCvtermID( name ).longValue();
      for(int value_index = 0; value_index < qualifier_values.size();
        ++value_index)
      {
        chado_feature.addQualifier(type_id, 
                                   (String)qualifier_values.elementAt(value_index));
      }
    } 
  }

  /**
   *
   *  Return a string containing one qualifier per line.  These are the
   *  original qualifiers, not the qualifiers from the qualifier_text_area.
   *  @param qualifiers	the <code>QualifierVector</code>
   *  @return	the <code>String</code> representation of the qualifiers
   *
   **/
  private String getQualifierString(QualifierVector qualifiers)
  {
    final StringBuffer buffer = new StringBuffer();

    for(int qualifier_index = 0; qualifier_index < qualifiers.size();
        ++qualifier_index)
    {
      final Qualifier this_qualifier = (Qualifier)qualifiers.elementAt(qualifier_index);

      final StringVector qualifier_strings =
                       StreamQualifier.toStringVector(null, this_qualifier);

      for(int value_index = 0; value_index < qualifier_strings.size();
          ++value_index)
      {
        final String qualifier_string = (String)qualifier_strings.elementAt(value_index);
        buffer.append(qualifier_string + "\n");
      }
    }

    return buffer.toString();
  }

  /**
   *
   * Find the qualifiers that have changed or been added and UPDATE
   * INSERT or DELETE accordingly.
   *
   * @param qualifiers_old	old qualifiers
   * @param qualifiers_new	new qualifiers
   * @param feature		GFF feature that has been changed
   *
   */
  private void editQualifiers(QualifierVector qualifiers_old, 
                              QualifierVector qualifiers_new, 
                              GFFStreamFeature feature)
  {
  
    String feature_id = (String)(feature.getQualifierByName("ID").getValues()).elementAt(0);
    ChadoTransaction tsn;

    // look for qualifiers to DELETE
    for(int qualifier_index = 0; qualifier_index < qualifiers_old.size();
        ++qualifier_index)
    {
      final Qualifier this_qualifier = (Qualifier)qualifiers_old.elementAt(qualifier_index);
      String name = this_qualifier.getName();

      if(!qualifiers_new.contains(name))
      {
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
                                   feature_id, "featureprop");

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

      if(name.equals("ID") ||
         name.equals("Parent") ||
         name.equals("gff_seqname") ||
         name.equals("gff_source") ||
         name.equals("timelastmodified"))
        continue;

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
        if(!need_to_update)
          continue;
      }

      // get the cvterm_id for this featureprop/qualifier
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
  
      String cvterm_id = lcvterm_id.toString();

      if(old_index > -1 &&
         new_qualifier_strings.size() == old_qualifier_strings.size())
      {
        //  
        // UPDATE existing featureprop's
        //
        for(int value_index = 0; value_index < new_qualifier_strings.size();
            ++value_index)
        {
          String qualifier_string = (String)new_qualifier_strings.elementAt(value_index);
//        if(old_qualifier_strings.contains(qualifier_string))
//          continue;

          tsn = new ChadoTransaction(ChadoTransaction.UPDATE,
                                     feature_id, "featureprop");
          int index = qualifier_string.indexOf("=");

          if(index > -1)
            qualifier_string = qualifier_string.substring(index+1);

          tsn.addProperty("value", "'"+ stripQuotes(qualifier_string) +"'");
          tsn.setConstraint("featureprop.type_id", "'"+cvterm_id+"'");
          tsn.setConstraint("rank", Integer.toString(value_index));
          sql.add(tsn);
//          String[] sql_array = tsn.getSqlQuery();
//          for(int i=0; i<sql_array.length; i++)
//            System.out.println(sql_array[i]);
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
                                     feature_id, "featureprop");

          tsn.setConstraint("type_id", cvterm_id);
          sql.add(tsn);
        }
          
        //
        // INSERT new featureprops
        //
        for(int value_index = 0; value_index < new_qualifier_strings.size();
            ++value_index)
        {
          tsn = new ChadoTransaction(ChadoTransaction.INSERT,
                                     feature_id, "featureprop");
          String qualifier_string = (String)new_qualifier_strings.elementAt(value_index);
          int index = qualifier_string.indexOf("=");
          if(index > -1)
            qualifier_string = qualifier_string.substring(index+1);

          tsn.addProperty("value", "'"+ stripQuotes(qualifier_string) +"'");
          tsn.addProperty("type_id", "'"+cvterm_id+"'");
          tsn.addProperty("rank", Integer.toString(value_index));
          sql.add(tsn);
        }

//      System.out.println("******** "+DatabaseDocument.getCvtermID(name));
      }
    }

  }

  /**
   *
   * Strip out quotes around a string.
   *
   */
  private String stripQuotes(String s)
  {
    if(s.startsWith("\"") && s.endsWith("\""))
      s = s.substring(1,s.length()-1);
    
    return s;
  }


  /**
   *
   *  This method fixes up the location of this Feature when a sequence
   *  changes.
   *
   **/ 
  public void sequenceChanged(final SequenceChangeEvent event)
  {
  }

  /**
   *
   * Commit the transactions back to the database.  
   *
   **/
  public void commit(DatabaseDocument dbDoc)
  {
    dbDoc.commit(sql);
    sql = new Vector();
  }
}

