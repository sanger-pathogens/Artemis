/*
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
package uk.ac.sanger.artemis.circular.digest;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.circular.DNADraw;
import uk.ac.sanger.artemis.circular.Track;
import uk.ac.sanger.artemis.circular.TrackManager;
import uk.ac.sanger.artemis.circular.Wizard;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.sequence.MarkerRange;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;

import javax.swing.JOptionPane;

public class Utils
{
	//private static final Logger logger = Logger.getLogger(Utils.class);

  protected static ReportDetails findCutSitesFromEmbossReport(Reader reader)
  {
    EmbossTableParser etp = new EmbossTableParser();
    ReportDetails ret = new ReportDetails();
    try
    {
      BufferedReader br = new BufferedReader(reader);
      ret.cutSites = etp.parse(br);
      ret.length = etp.getLength();
    }
    catch (IOException exp)
    {
      throw new RuntimeException("Couldn't read, or parse results");
    }
    return ret;
  }
	
  protected static List<FragmentBand> findCutSitesFromExperiment(Reader reader) 
            throws IOException
  {
    BufferedReader br = new BufferedReader(reader);
    String line;
    
    List<FragmentBand> list = new ArrayList<FragmentBand>();
    while ((line = br.readLine()) != null)
    { 
      FragmentBand band = new FragmentBand();
      band.genomeFragmentLength = Integer.parseInt(line.trim());
      
      list.add(band);
    }
    return list;
  }

  /**
   * Creates the DNAPlotter components and adds a forward and reverse track for
   * the restriction enzyme results.
   * 
   * @param rd
   * @param entryGroup
   * @return
   */
  protected static DNADraw createDNADrawFromReportDetails(
      final ReportDetails rd, final EntryGroup entryGroup)
  {
    DNADraw dna = new DNADraw();
    dna.setArtemisEntryGroup(entryGroup);

    int sequenceLength = rd.length;
    Hashtable<String, Object> lineAttr = new Hashtable<String, Object>();
    lineAttr.put("lsize", new Integer(1));
    lineAttr.put("circular", new Boolean(true));
    lineAttr.put("start", new Integer(0));
    lineAttr.put("end", new Integer(sequenceLength));
    dna.setLineAttributes(lineAttr);

    // set ticks
    int div;
    if (sequenceLength < 1000)
      div = 100;
    else if (sequenceLength < 10000)
      div = 1000;
    else if (sequenceLength < 100000)
      div = 10000;
    else
      div = 100000;
    int tickNum = sequenceLength / div;
    int tick = tickNum * (div / 10);
    
    if(tick < 1)
      tick = 2;
    while ((sequenceLength % tick) < (div / 10))
    {
      tickNum++;
      tick = tickNum * (div / 10);
    }
    dna.setMinorTickInterval(tick);
    dna.setTickInterval(tick);
    dna.setOpaque(true);

    final Entry newEntry = dna.getArtemisEntryGroup().createEntry("RE");
    newEntry.getEMBLEntry().setName("Restriction Sites");
    Track forward = new Track(0.9, "misc_feature", true, false, newEntry);
    forward.setQualifier("note");
    forward.setQualifierValue("plus");
    Track reverse = new Track(0.85, "misc_feature", true, false, newEntry);
    reverse.setQualifier("note");
    reverse.setQualifierValue("minus");
    Track fileTrack = new Track(0.8, "CDS", true, true, entryGroup
        .getSequenceEntry());

    Wizard.tracks = new Track[3];
    Wizard.tracks[0] = forward;
    Wizard.tracks[1] = reverse;
    Wizard.tracks[2] = fileTrack;

    int counter = 0;
    /*if (rd.cutSites.size() == 1)
    {
      CutSite cutSite = rd.cutSites.get(0);
      // dna.addFeatureToTrack(createFeature(cutSite.getEnd(), rd.length
      // + cutSite.getStart(), 1, sequenceLength), forward, false);
      addFeature(cutSite.getFivePrime(), rd.length, counter, newEntry, dna,
          cutSite.isForward());
      
      dna.setGeneticMarker(new Vector());
      TrackManager trackManager = dna.getTrackManager();
      if (trackManager == null)
      {
        trackManager = new TrackManager(dna);
        dna.setTrackManager(trackManager);
      }
      trackManager.update(Wizard.tracks);

      return dna;
    }*/

    Iterator<CutSite> it = rd.cutSites.iterator();
    CutSite firstCutSite = null;
    int lastCutPos = 1;

    while (it.hasNext())
    {
      CutSite cutSite = it.next();

      if (counter == 0)
      {
        firstCutSite = cutSite;
        lastCutPos = cutSite.getThreePrime();
        counter++;
        continue;
      }

      addFeature(lastCutPos, cutSite.getFivePrime(), counter, newEntry, dna,
          cutSite.isForward());
      lastCutPos = cutSite.getFivePrime();
      counter++;
    }

    if (firstCutSite != null)
    {
      addFeature(lastCutPos, rd.length, -1, newEntry, dna, firstCutSite
          .isForward());
      addFeature(1, firstCutSite.getFivePrime(), -1, newEntry, dna,
          firstCutSite.isForward());
    }

    dna.setGeneticMarker(new Vector());
    TrackManager trackManager = dna.getTrackManager();
    if (trackManager == null)
    {
      trackManager = new TrackManager(dna);
      dna.setTrackManager(trackManager);
    }
    trackManager.update(Wizard.tracks);

    return dna;
  }

  /**
   * Add a new feature to the entry with the given coordinates
   * 
   * @param coord1
   * @param coord2
   * @param counter
   * @param newEntry
   * @param dna
   * @param isForward
   */
  private static void addFeature(int coord1, int coord2, int counter,
      Entry entry, DNADraw dna, boolean isForward)
  {
    String colour;
    if (counter < 0)
      colour = "7";
    else if (counter % 2 == 0)
      colour = "2";
    else
      colour = "5";

    try
    {
      QualifierVector qualifiers = new QualifierVector();
      qualifiers.add(new Qualifier("colour", colour));
      final MarkerRange r;

      if (isForward)
        qualifiers.add(new Qualifier("note", "plus"));
      else
        qualifiers.add(new Qualifier("note", "minus"));
      // if(isForward)
      r = new MarkerRange(dna.getArtemisEntryGroup().getSequenceEntry()
          .getBases().getForwardStrand(), coord1, coord2);
      /*
       * else r = new MarkerRange(
       * dna.getArtemisEntryGroup().getSequenceEntry().
       * getBases().getReverseStrand(), coord1,coord2);
       */

      // int len = r.getCount();
      uk.ac.sanger.artemis.io.Feature f = new uk.ac.sanger.artemis.io.EmblStreamFeature(
          new Key("misc_feature"), r.createLocation(), qualifiers);
      entry.add(new uk.ac.sanger.artemis.Feature(f), false);

      // if(!isForward)
      // System.out.println(counter+" "+coord1+".."+coord2+"   "+f.getLocation().toStringShort());
    }
    catch (Exception e)
    {
      e.printStackTrace();
    }
  }

}

class ReportDetails
{
  int length;
  List<CutSite> cutSites;
}
