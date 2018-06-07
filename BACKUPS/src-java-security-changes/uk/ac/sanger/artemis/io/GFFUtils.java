/* GFFUtils.java
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

import java.util.Vector;

import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;


public class GFFUtils
{
  /**
   * Test if this is feature is marked as having a stop codon
   * redefined as selenocysteine
   * @param f
   * @return
   */
  public static boolean isSelenocysteine(Feature f)
  {
    if(!(f instanceof GFFStreamFeature))
      return false;
    try
    {
      ChadoCanonicalGene gffGene = ((GFFStreamFeature)f).getChadoGene();
      if(gffGene == null)
        return false;
      String transcript = gffGene.getTranscriptFromName(
          GeneUtils.getUniqueName(f));
      if(transcript == null)
        return false;
      Feature pep = gffGene.getProteinOfTranscript(transcript);
      if(pep == null)
        return false;
      if(pep.getQualifierByName("stop_codon_redefined_as_selenocysteine") != null)
        return true;
    }
    catch (Exception e){}

    return false;
  }
  
  /**
   * Update the segment range store for GFFStreamFeature with a new location
   * @param gff
   * @param oldLocation
   * @param newLocation
   */
  public static void updateSegmentRangeStore(final GFFStreamFeature gff,
                                             final Location oldLocation, 
                                             final Location newLocation)
  {
    final RangeVector rv_new = newLocation.getRanges();
    final RangeVector rv_old = oldLocation.getRanges();
    if(rv_new.size() != rv_old.size())
    {
      final RangeVector rangesToAdd = new RangeVector();
      for(Range r: rv_new)
        if(!rv_old.containsRange(r))
          rangesToAdd.add(r);
      
      final Vector<Integer> deleted = new Vector<Integer>();
      for(int ideleted = 0; ideleted < rv_old.size(); ideleted++)
        if(!rv_new.containsRange(rv_old.get(ideleted)))
          deleted.add(ideleted);

      try
      {
        if(gff.getQualifierByName("Parent") != null)
          GeneUtils.addSegment(gff, rangesToAdd, 
             gff.getQualifierByName("Parent").getValues().get(0));
      }
      catch(Exception e)
      {
        e.printStackTrace();
      }
      
      for(Integer d: deleted)
        gff.getSegmentRangeStore().remove(
            gff.getSegmentID(rv_old.elementAt(d)));
    }
    else if(gff.getSegmentRangeStore() != null)
    {
      Vector<Integer> changes = new Vector<Integer>();
      for(int i=0; i<rv_old.size(); i++)
      {
        Range rnew = rv_new.elementAt(i);
        Range rold = rv_old.elementAt(i);

        if(rnew.getStart() != rold.getStart() ||
           rnew.getEnd()   != rold.getEnd() ||
           (oldLocation.isComplement(rold) !=
            newLocation.isComplement(rnew)))
          changes.add(i);
      }

      for(Integer c: changes)
      {
        Range rnew = rv_new.elementAt(c);
        String segId = gff.getSegmentID(rnew);
        if(segId == null)
          segId = gff.getSegmentID(rv_old.elementAt(c));
        gff.getSegmentRangeStore().put(segId, rnew);
      }
    }
  }
  
}