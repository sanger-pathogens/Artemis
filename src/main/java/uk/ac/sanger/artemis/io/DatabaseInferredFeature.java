/* DatabaseInferredFeature.java
 *
 * created: 2009
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2009  Genome Research Limited
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

import java.util.List;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.FeatureSegmentVector;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.ReadOnlyException;

/**
 * A feature (CDS and UTRs) that is inferred from a chado
 * database hierarchy (gene, transcript, exon, and polypeptide)
 * defined by the <code>ChadoCanonicalGene</code>.
 */
public class DatabaseInferredFeature
       extends GFFStreamFeature implements FeatureChangeListener
{
  
  public DatabaseInferredFeature(Key key, 
                                 Location location,
                                 QualifierVector qualifiers,
                                 ChadoCanonicalGene gene)
  {
    super(key, location, qualifiers);
    setChadoGene(gene);
  }

  public void addFeatureListeners()
  {
    if(!getKey().equals(Key.CDS))
      return;
    String parent = (String) getQualifierByName("Parent").getValues().get(0);
    if(getChadoGene().getProteinOfTranscript(parent) == null)
      return;

    Feature protein = getChadoGene().getProteinOfTranscript(parent);  
    ((uk.ac.sanger.artemis.Feature)protein.getUserData()).addFeatureChangeListener(this);
    
    List list = getChadoGene().getSpliceSitesOfTranscript(parent, "exon");
    for(int i=0; i<list.size(); i++)
    {
      uk.ac.sanger.artemis.Feature f = 
          ((uk.ac.sanger.artemis.Feature)((Feature)list.get(i)).getUserData());
      f.addFeatureChangeListener(this);
    }

    updateInferredLocations();
  }
  
  private void removeFeatureListeners()
  {
    if(!getKey().equals(Key.CDS))
      return;
    String parent = (String) getQualifierByName("Parent").getValues().get(0);
    if(getChadoGene().getProteinOfTranscript(parent) == null)
      return;

    Feature protein = getChadoGene().getProteinOfTranscript(parent);  
    ((uk.ac.sanger.artemis.Feature)protein.getUserData()).removeFeatureChangeListener(this);
    
    List list = getChadoGene().getSpliceSitesOfTranscript(parent, "exon");
    for(int i=0; i<list.size(); i++)
    {
      uk.ac.sanger.artemis.Feature f = 
          ((uk.ac.sanger.artemis.Feature)((Feature)list.get(i)).getUserData());
      f.removeFeatureChangeListener(this);
    }
  }
  
  public void featureChanged(FeatureChangeEvent event)
  {
    if(event.getType() == FeatureChangeEvent.LOCATION_CHANGED ||
       event.getType() == FeatureChangeEvent.SEGMENT_CHANGED)
    {   
      if(getKey().equals(Key.CDS))
        updateInferredLocations();
    }
  }
  
  public void updateInferredLocations()
  {
    String transcriptID = (String) getQualifierByName("Parent").getValues().get(0);
    if(getChadoGene().getProteinOfTranscript(transcriptID) == null)
      return;
    
    Range proteinRange =
      getChadoGene().getProteinOfTranscript(transcriptID).getLocation().getTotalRange();
    
    Feature exonFeature = 
      ((Feature)(getChadoGene().getSpliceSitesOfTranscript(transcriptID, "exon").get(0)));
    
    RangeVector r_old = exonFeature.getLocation().getRanges();
    RangeVector r_cds = new RangeVector();
    RangeVector r_diff_5 = new RangeVector();
    RangeVector r_diff_3 = new RangeVector();
    boolean isComplement = exonFeature.getLocation().isComplement();
    
    for(int i=0; i<r_old.size(); i++)
    {
      Range exonRange = (Range) r_old.get(i);
      findRanges(proteinRange, exonRange, isComplement,
                 r_diff_5, r_diff_3, r_cds);
    }
    
    //
    // update CDS location
    try
    {
      if(!containsAll(r_cds, getLocation().getRanges()))
        super.setLocation(new Location(r_cds, getLocation().isComplement()));
    }
    catch (ReadOnlyException e)
    {
      e.printStackTrace();
    }
    catch (OutOfRangeException e)
    {
      e.printStackTrace();
    }

    //
    // update UTR locations
    Entry entry = 
       ((uk.ac.sanger.artemis.Feature) getUserData()).getEntry();  
    try
    {
      if (r_diff_5.size() > 0)
      {
        List listUTR = getChadoGene().get5UtrOfTranscript(transcriptID);
        boolean isFound = containsAll(listUTR,  r_diff_5);
          
        if(!isFound)
        {
          Feature utrFeature = addUTR(transcriptID, "five_prime_UTR", r_diff_5, 
              getChadoGene().get5UtrOfTranscript(transcriptID), entry);
          getChadoGene().add5PrimeUtr(transcriptID, utrFeature);
        }
      }
      else
        removeUTR(getChadoGene().get5UtrOfTranscript(transcriptID), entry);
      
      if(r_diff_3.size() > 0)
      {
        List listUTR = getChadoGene().get3UtrOfTranscript(transcriptID);
        boolean isFound = containsAll(listUTR,  r_diff_3);
          
        if(!isFound)
        { 
          Feature utrFeature = addUTR(transcriptID, "three_prime_UTR", r_diff_3, 
              getChadoGene().get3UtrOfTranscript(transcriptID), entry);
          getChadoGene().add3PrimeUtr(transcriptID, utrFeature);
        }
      }
      else
        removeUTR(getChadoGene().get3UtrOfTranscript(transcriptID), entry);
    }
    catch (ReadOnlyException e)
    {
      e.printStackTrace();
    }
    catch (EntryInformationException e)
    {
      e.printStackTrace();
    }
  }
  
  /**
   * Classify the ranges based on the exon and protein ranges. 
   * The range is divided into 5'UTR, 3'UTR and CDS ranges.
   * @param proteinRange
   * @param exonRange
   * @param isComplement
   * @param r_diff_5
   * @param r_diff_3
   * @param r_cds
   */
  private void findRanges(Range proteinRange, 
                          Range exonRange,
                          boolean isComplement,
                          RangeVector r_diff_5,
                          RangeVector r_diff_3,
                          RangeVector r_cds)
  {
    if( !proteinRange.contains(exonRange) )
    {
      if(exonRange.getStart () < proteinRange.getStart () &&
         exonRange.getEnd ()   > proteinRange.getEnd ())
      {
        try
        {
          if(isComplement)
          {
            r_diff_3.add(new Range(exonRange.getStart(), proteinRange.getStart()-1));
            r_diff_5.add(new Range(proteinRange.getEnd()+1, exonRange.getEnd()));
          }
          else
          {
            r_diff_5.add(new Range(exonRange.getStart(), proteinRange.getStart()-1));
            r_diff_3.add(new Range(proteinRange.getEnd()+1, exonRange.getEnd()));
          }
          r_cds.add(new Range(proteinRange.getStart(), proteinRange.getEnd()));
        }
        catch (OutOfRangeException e){ e.printStackTrace(); }
      }
      else if(exonRange.getStart() >= proteinRange.getStart() &&
              exonRange.getStart() < proteinRange.getEnd())
      {   
        // exon is partially CDS and UTR
        try
        {
          r_cds.add(new Range(exonRange.getStart(), proteinRange.getEnd()));
          
          if(isComplement)
            r_diff_5.add(new Range(proteinRange.getEnd()+1, exonRange.getEnd()));
          else
            r_diff_3.add(new Range(proteinRange.getEnd()+1, exonRange.getEnd()));
        }
        catch (OutOfRangeException e){ e.printStackTrace(); }
      }
      else if(exonRange.getEnd() > proteinRange.getStart() &&
              exonRange.getEnd() <= proteinRange.getEnd())
      {
        try
        {
          r_cds.add(new Range(proteinRange.getStart(), exonRange.getEnd()));
          
          if(isComplement)
            r_diff_3.add(new Range(exonRange.getStart(), proteinRange.getStart()-1));
          else
            r_diff_5.add(new Range(exonRange.getStart(), proteinRange.getStart()-1));
        }
        catch (OutOfRangeException e){ e.printStackTrace(); }
      }
      else if(exonRange.getStart() < proteinRange.getStart())
        r_diff_5.add(exonRange);
      else if(exonRange.getStart() > proteinRange.getStart())
        r_diff_3.add(exonRange);
    }
    else
      r_cds.add(exonRange);
  }
  
  /**
   * Override to ensure that associated features get updated.
   */
  public void setLocation (final Location location)
     throws ReadOnlyException, OutOfRangeException 
  {
    if(getQualifierByName("Parent") != null)
      removeFeatureListeners();
    super.setLocation(location);

    if(getKey().equals(Key.CDS) && getQualifierByName("Parent") != null)
    {
      // update the exon feature
      String transcriptID = (String) getQualifierByName("Parent").getValues().get(0);
      Feature exonFeature = 
        ((Feature)(getChadoGene().getSpliceSitesOfTranscript(transcriptID, "exon").get(0)));
      
      RangeVector exonRanges = exonFeature.getLocation().getRanges();
      RangeVector cdsRanges = getLocation().getRanges();
      FeatureSegmentVector segments =
        ((uk.ac.sanger.artemis.Feature)exonFeature.getUserData()).getSegments();
      
      for(int i=0; i<cdsRanges.size(); i++)
      {
        Range cdsRange = (Range) cdsRanges.get(i);
        if(!exonRanges.containsRange(cdsRange))
        {
          for(int j=0; j<segments.size(); j++)
          {
            Range exonRange = segments.elementAt(j).getRawRange();
            if( exonRange.overlaps(cdsRange) &&
               (exonRange.getStart() == cdsRange.getStart() || 
                exonRange.getEnd()   == cdsRange.getEnd()) )
            {
              segments.elementAt(j).setRange(cdsRange);
            }
          }
        }
      }
      
    }
    else if(getKey().getKeyString().equals("five_prime_UTR"))
    {
      
    }
    
    if(getQualifierByName("Parent") != null)
      addFeatureListeners();
  }

  /**
   * Add a UTR to an entry.
   * @param transcriptName
   * @param keyName
   * @param ranges
   * @param existingUTRs
   * @param entry
   * @return
   * @throws ReadOnlyException
   * @throws EntryInformationException
   */
  private Feature addUTR(String transcriptName, 
                         String keyName,
                         RangeVector ranges,
                         List existingUTRs,
                         Entry entry)
          throws ReadOnlyException, EntryInformationException
  {
    removeUTR(existingUTRs, entry);   

    QualifierVector qualifiers = new QualifierVector();
    qualifiers.add(new Qualifier("ID", transcriptName+
        ( keyName.equals("five_prime_UTR") ? ":5UTR" : ":3UTR" ) ));
    qualifiers.add(new Qualifier("Parent", transcriptName));
    
    Location l = new Location(ranges, getLocation().isComplement());
    DatabaseInferredFeature utrFeature = new DatabaseInferredFeature(
        new Key(keyName), l, qualifiers, getChadoGene());
     
    uk.ac.sanger.artemis.Feature f = new uk.ac.sanger.artemis.Feature(utrFeature);
    f.setEntry(entry);
    
    return entry.getEMBLEntry().add(utrFeature);
  }
  
  /**
   * Remove a list of UTRs from an entry.
   * @param existingUTRs
   * @param entry
   * @throws ReadOnlyException
   */
  private void removeUTR(List existingUTRs,
                         Entry entry) throws ReadOnlyException
  {
    if (existingUTRs == null)
      return;
    
    for (int i = 0; i < existingUTRs.size(); i++)
      entry.getEMBLEntry().remove((Feature) existingUTRs.get(i));
    existingUTRs.clear();
  }
  
  /**
   * Check that a list of UTRs contains all the ranges in a 
   * new <code>RangeVector</code>.
   * @param listUTR
   * @param newRanges
   * @return
   */
  private boolean containsAll(List listUTR, 
                              RangeVector newRanges)
  {
    RangeVector oldRanges = new RangeVector();
    if(listUTR != null)
    {
      for(int i=0;i<listUTR.size(); i++)
      {
        Feature featureUTR = (Feature) listUTR.get(i);
        oldRanges.addAll(featureUTR.getLocation().getRanges());
      }
    }
    return containsAll(newRanges, oldRanges);
  }
  
  /**
   * Check if a <code>RangeVector</code> contains all the ranges
   * in a second <code>RangeVector</code>.
   * @param newRanges
   * @param oldRanges
   * @return
   */
  private boolean containsAll(RangeVector newRanges, 
                              RangeVector oldRanges)
  {
    if(newRanges.size() != oldRanges.size())
      return false;
    for(int i=0; i<oldRanges.size(); i++)
    {
      if(!newRanges.containsRange( (Range)oldRanges.get(i) ))
        return false;
    }
    return true;
  }
  
  /**
   * Create a CDS
   * @param transcriptId
   * @param exonFeature
   * @param chadoGene
   * @param entry
   */
  public static void createFeature(final String transcriptId,
                                   final Feature exonFeature,
                                   final ChadoCanonicalGene chadoGene,
                                   final Entry entry)
  {
    QualifierVector qualifiers = new QualifierVector();
    qualifiers.add(new Qualifier("ID", transcriptId+":CDS"));
    qualifiers.add(new Qualifier("Parent", transcriptId));
    
    DatabaseInferredFeature cdsFeature = new DatabaseInferredFeature(
        Key.CDS, exonFeature.getLocation(), qualifiers, chadoGene);
    uk.ac.sanger.artemis.Feature cds = new uk.ac.sanger.artemis.Feature(cdsFeature);

    try
    {
      chadoGene.addSplicedFeatures(transcriptId, cdsFeature);
      entry.add(cds, true);
      cdsFeature.addFeatureListeners();
    }
    catch (ReadOnlyException e)
    {
      e.printStackTrace();
    }
    catch (EntryInformationException e)
    {
      e.printStackTrace();
    }
    catch (OutOfRangeException e)
    {
      e.printStackTrace();
    } 
  }
  
  /**
   * 
   * @param entryGroup
   */
  public static void addListenersToEntryGroup(EntryGroup entryGroup)
  {
    FeatureVector features = entryGroup.getAllFeatures();
    for(int i=0; i<features.size(); i++)
      if(features.elementAt(i).getEmblFeature() instanceof DatabaseInferredFeature)
        ((DatabaseInferredFeature)features.elementAt(i).getEmblFeature()).addFeatureListeners();
  }
}