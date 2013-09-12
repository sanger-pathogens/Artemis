/* 
 * This file is part of Artemis
 *
 * Copyright (C) 2013  Genome Research Limited
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
package uk.ac.sanger.artemis.io;

import static org.junit.Assert.assertTrue;
import junit.framework.Assert;

import org.junit.Test;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.io.Entry;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.sequence.NoSequenceException;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.StringVector;

public class GFFTest
{
  @Test
  /**
   * For a GFF with multiple sequences check the offset position 
   * of a gene is correctly set.
   */
  public void testGFFMultipleFastaOffset()
  {
    try
    {
      final Entry entry = Utils.getEntry("/data/Pf3D7_01_02_v3.gff.gz");
      final EntryGroup egrp = new SimpleEntryGroup();
      egrp.add(new uk.ac.sanger.artemis.Entry(entry));
      final FeatureVector features = egrp.getAllFeatures();

      // change the translation table did cause a problem
      // with the offset that has now been fixed in GFFDocumentEntry
      Utils.changeTranslationTable("11");

      for (int i=0; i<features.size(); i++)
      {
        Feature f = features.elementAt(i);
        try
        {
          final Qualifier q = f.getQualifierByName("ID");
          if (q != null)
          {
            final String id = q.getValues().get(0);
            if (id.equals("PF3D7_0200100"))
            {
              assertTrue("Offset check first base: " + id + 
                  " " +f.getFirstBase() + " != 666083",
                  f.getFirstBase() == 666083);

              assertTrue("Offset check location: " + id+ 
                  " " +f.getEmblFeature().getLocation().getFirstBase() + " != 666083", 
                  f.getEmblFeature().getLocation().getFirstBase() == 666083);
            }
          }
        }
        catch(InvalidRelationException e)
        {
          Assert.fail(e.getMessage());
        }
      }
    }
    catch (OutOfRangeException e)
    {
      Assert.fail(e.getMessage());
    }
    catch (NoSequenceException e)
    {
      Assert.fail(e.getMessage());
    }
  }

}