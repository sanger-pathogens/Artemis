/* DatabaseDocumentEntry.java
 *
 * created: Mar 2005
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2005  Genome Research Limited
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

package uk.ac.sanger.artemis.io;

import java.sql.*;
import uk.ac.sanger.artemis.util.*;

import java.io.*;

/**
 *  A DatabaseDocumentEntry that can read an Entry from a Document containing
 *  relational database entry
 **/

public class DatabaseDocumentEntry extends GFFDocumentEntry 
    implements DocumentEntry 
{

  public DatabaseDocumentEntry(final Document doc)
         throws EntryInformationException, IOException
  {
    super(doc,null);
  }

  /**
   *  Returns true if and only if this entry is read only.  For now this
   *  always returns true - BlastDocumentEntry objects can't be changed.
   **/
  public boolean isReadOnly() 
  {
    return false;
  }

  /**
   *  If the given Sequence can be added directly to this Entry, then return a
   *  copy of it, otherwise create and return a new feature of the appropriate
   *  type for this Entry.
   **/
  protected StreamSequence makeNativeSequence(final Sequence sequence)
  {
    return new FastaStreamSequence(sequence);
  }


  /**
   *  If the given feature can be added directly to this Entry, then return
   *  it, otherwise create and return a new feature of the appropriate type.
   *  @param copy if true then always new a new copy of the Feature.
   **/
  protected SimpleDocumentFeature makeNativeFeature(final Feature feature,
                                                    final boolean copy)
  {
    if (!copy)
    {
      if(feature instanceof DatabaseStreamFeature)
        return (DatabaseStreamFeature)feature;
      else
        return (GFFStreamFeature)feature;
    }
    else
      return new DatabaseStreamFeature(feature);      
  }

}
