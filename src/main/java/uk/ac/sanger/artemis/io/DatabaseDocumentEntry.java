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

import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.chado.ChadoTransactionManager;
import uk.ac.sanger.artemis.util.*;

import java.io.*;

/**
 *  A DatabaseDocumentEntry that can read an Entry from a Document containing
 *  relational database entry
 **/

public class DatabaseDocumentEntry extends GFFDocumentEntry 
    implements DocumentEntry 
{
  private PartialSequence sequence;
  private ChadoTransactionManager ctm = new ChadoTransactionManager();
  
  public DatabaseDocumentEntry(final Document doc, 
                               final ReadListener listener)
         throws EntryInformationException, IOException
  {
    super(doc,listener);
  }

  public DatabaseDocumentEntry()
  {
    super(new GFFEntryInformation());
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
  protected Object makeNativeFeature(final Feature feature,
                                     final boolean copy)
  {
    if (!copy)
    {
      if(feature instanceof DatabaseStreamFeature)
        return (DatabaseStreamFeature)feature;
      else if(feature instanceof GFFStreamFeature)
        return (GFFStreamFeature)feature;
      else
      {
        QualifierVector qualifiers = feature.getQualifiers().copy();
        if(qualifiers.getQualifierByName("ID") == null)
        {
          String idStr = null;
          StringVector v = Options.getOptions().getSystematicQualifierNames();
          for(int i=0; i<v.size(); i++)
          {
            final String sysName = (String)v.get(i);
            if(qualifiers.getQualifierByName(sysName) != null)
            {
              idStr = (String)qualifiers.getQualifierByName(sysName).getValues().get(0);
              break;
            }
          }
          // autogenerate ID
          if(idStr == null)
            idStr = feature.getKey().getKeyString()+":"+feature.getLocation().toString();
          qualifiers.add(new Qualifier("ID", idStr));
        }
        
        return new DatabaseStreamFeature(feature.getKey(), 
                   feature.getLocation(), qualifiers);
      }
    }
    else
      return new DatabaseStreamFeature(feature);      
  }
  
  public void setPartialSequence(final PartialSequence sequence)
  {
    this.sequence = sequence;
  }
  
  public Sequence getSequence() 
  {
    if(sequence == null)
      return super.getSequence();
    return sequence;
  }
  
  /**
   *  Returns true if and only if there have been some changes to this Entry
   *  since the last save.
   **/
  public boolean hasUnsavedChanges() 
  {
    return ctm.hasTransactions();
  }
  
  public ChadoTransactionManager getChadoTransactionManager()
  {
    return ctm;
  }
}
