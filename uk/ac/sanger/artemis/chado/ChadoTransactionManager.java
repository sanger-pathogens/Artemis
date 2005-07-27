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

import java.util.Vector;


/**
*
* Class responsible for tracking the feature insertions/deletions/changes
* to be able to commit these back to the database.
*
**/
public class ChadoTransactionManager
       implements FeatureChangeListener, SequenceChangeListener 
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
      System.out.println("HERE featureChanged() "+event.getType()+" "+event.getFeature().getIDString());

      ChadoTransaction tsn = null;
      if(event.getType() == event.LOCATION_CHANGED)
      {
        System.out.println("LOCATION_CHANGED ");

        RangeVector rv_new = event.getNewLocation().getRanges();
        RangeVector rv_old = event.getOldLocation().getRanges();

        int ichanged;
        for(ichanged=0; ichanged<rv_old.size(); ichanged++)
        {
          if(!rv_new.contains(rv_old.elementAt(ichanged)))
            break; 
        }
 
        Range range_new = (Range)rv_new.elementAt(ichanged);
        String seg_id   = feature.getSegmentID(range_new);
 
        tsn = new ChadoTransaction(ChadoTransaction.UPDATE, 
                                   seg_id, "featureloc");
        tsn.addProperty("fmin", Integer.toString(range_new.getStart()-1));
        tsn.addProperty("fmax", Integer.toString(range_new.getEnd()));
      }
      else if(event.getType() == event.QUALIFIER_CHANGED)
      {
        System.out.println("QUALIFIER_CHANGED ");
        System.out.println(getQualifierString(event.getOldQualifiers()));
        System.out.println(getQualifierString(event.getNewQualifiers()));
      }
      else if(event.getType() == event.ALL_CHANGED)
        System.out.println("ALL_CHANGED");
 
      if(tsn != null)  
        sql.add(tsn);
    }
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

    for(int qualifier_index = 0; qualifier_index < qualifiers.size();
        ++qualifier_index)
    {
      final Qualifier this_qualifier = (Qualifier)qualifiers.elementAt(qualifier_index);

      final StringVector qualifier_strings =
                       StreamQualifier.toStringVector(null, this_qualifier);

      for(int value_index = 0; value_index < qualifier_strings.size();
          ++value_index)
      {
        final String qualifier_string = qualifier_strings.elementAt(value_index);
        buffer.append(qualifier_string + "\n");
      }
    }

    return buffer.toString();
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

