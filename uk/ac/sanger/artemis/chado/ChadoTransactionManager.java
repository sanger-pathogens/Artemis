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
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.EntryChangeListener;
import uk.ac.sanger.artemis.EntryChangeEvent;
import java.util.Vector;
import javax.swing.JOptionPane;

/**
*
* Class responsible for tracking the feature insertions/deletions/changes
* to be able to commit these back to the database.
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
    System.out.println("HERE");
    if(event.getType() == EntryChangeEvent.FEATURE_ADDED)
      System.out.println("HERE FEATURE_ADDED");

    Feature feature = event.getFeature();

    if(feature == null)
      System.out.println("HERE feature == null");
    else
      System.out.println("HERE feature != null");

    System.out.println(event.getEntry().getName());
  }

  /**
  *
  *  Return a string containing one qualifier per line.  These are the
  *  original qualifiers, not the qualifiers from the qualifier_text_area.
  *
  **/
  private String getQualifierString(QualifierVector qualifiers)
  {
    final StringBuffer buffer = new StringBuffer();

    ChadoTransaction tsn;
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
    for(int qualifier_index = 0; qualifier_index < qualifiers_new.size();
        ++qualifier_index)
    {
      final Qualifier this_qualifier = (Qualifier)qualifiers_new.elementAt(qualifier_index);

      if(!qualifiers_old.contains(this_qualifier))
      {
        String name = this_qualifier.getName();
        int old_index = qualifiers_old.indexOfQualifierWithName(name);

        Qualifier this_old_qualifier = null;
        StringVector old_qualifier_strings = null;
        if(old_index> -1)  // update qualifier
        {
          this_old_qualifier = (Qualifier)qualifiers_old.elementAt(old_index);

          old_qualifier_strings =
                     StreamQualifier.toStringVector(null, this_old_qualifier);
        }

        final StringVector new_qualifier_strings =
                     StreamQualifier.toStringVector(null, this_qualifier);

        // get the cvterm_id for this featureprop/qualifier
        String cvterm_id = DatabaseDocument.getCvtermID(name).toString();

        if(cvterm_id == null)   // chado doesn't recognise this
        {
          JOptionPane.showMessageDialog(null, 
                      name+" is not a valid qualifier!",
                      "Invalid Qualifier",
                      JOptionPane.WARNING_MESSAGE);
          continue;
        }
 
        if(old_index > -1 &&
           new_qualifier_strings.size() == old_qualifier_strings.size())
        {
          //  
          // UPDATE existing featureprop's
          //
          for(int value_index = 0; value_index < new_qualifier_strings.size();
              ++value_index)
          {
            tsn = new ChadoTransaction(ChadoTransaction.UPDATE,
                                       feature_id, "featureprop");
            String qualifier_string = (String)new_qualifier_strings.elementAt(value_index);
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
          // DELETE existing featureprops
          //
          if(old_index > -1)
          {
            tsn = new ChadoTransaction(ChadoTransaction.DELETE,
                                       feature_id, "featureprop");

            tsn.setConstraint("type_id", cvterm_id);
            sql.add(tsn);
//          String[] sql_array = tsn.getSqlQuery();
//          for(int i=0; i<sql_array.length; i++)
//            System.out.println(sql_array[i]);
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
//          String[] sql_array = tsn.getSqlQuery();
//          for(int i=0; i<sql_array.length; i++)
//            System.out.println(sql_array[i]);
          }

//        System.out.println("******** "+DatabaseDocument.getCvtermID(name));
        }
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
  }
}

