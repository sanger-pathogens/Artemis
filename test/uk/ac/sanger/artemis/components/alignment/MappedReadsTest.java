/* 
 * This file is part of Artemis
 *
 * Copyright (C) 2014  Genome Research Limited
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
package uk.ac.sanger.artemis.components.alignment;

import static org.junit.Assert.assertTrue;

import java.awt.GraphicsEnvironment;
import java.net.URL;
import java.util.Hashtable;
import java.util.List;

import javax.swing.JFrame;

import junit.framework.Assert;
import uk.ac.sanger.artemis.io.Utils;
import org.junit.BeforeClass;
import org.junit.Test;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.components.EntryEdit;
import uk.ac.sanger.artemis.components.alignment.BamUtils;
import uk.ac.sanger.artemis.components.alignment.BamView;
import uk.ac.sanger.artemis.components.alignment.MappedReads;
import uk.ac.sanger.artemis.components.alignment.ReadCount;

public class MappedReadsTest 
{
  private static BamView bv;
  private static FeatureVector fv;
  
  @BeforeClass
  public static void setUp() 
  {
    // ignore if in headless mode with no x11
    if(GraphicsEnvironment.getLocalGraphicsEnvironment().isHeadless())
      return;
    URL entryFile = MappedReadsTest.class.getResource("/data/MAL_8h.bam");
    System.setProperty("bam", entryFile.getFile());
    final EntryGroup egrp = Utils.getEntryGroup("/data/MAL1.embl.gz");
    final EntryEdit ee = new EntryEdit(egrp);
    ee.setVisible(true);

    while( (bv = ee.getJamView()) == null) 
    {
      // wait for BamView to be constructed
      try {
        Thread.sleep(100);
      } catch(Exception e){};
    }
    
    // get a gene feature
    fv = new FeatureVector();
    final FeatureVector features = egrp.getAllFeatures();
    for(int i=0; i<features.size(); i++) 
    {
      Feature f = features.elementAt(i);
      if(f.getSystematicName().equals("PFA0110w"))
        fv.add(f);
    }
  }
  
  @Test
  /**
   * Test the read count for a gene
   */
  public void readCounts()
  {
    // ignore if in headless mode with no x11
    if(GraphicsEnvironment.getLocalGraphicsEnvironment().isHeadless())
      return;

    final Hashtable<String, List<ReadCount>> featureReadCount =
        BamUtils.calculateMappedReads(bv, fv, false, true, false, null, null);
    final List<ReadCount> cnts = featureReadCount.get("PFA0110w");
    
    ReadCount c = cnts.get(0);
    assertTrue(1495.f == c.senseCnt);
    assertTrue(998.f == c.antiCnt);
  }
  
  @Test
  /**
   * Read count for a gene excluding the intron
   */
  public void readCountsExcludeIntron()
  {
    // ignore if in headless mode with no x11
    if(GraphicsEnvironment.getLocalGraphicsEnvironment().isHeadless())
      return;

    final Hashtable<String, List<ReadCount>> featureReadCount =
        BamUtils.calculateMappedReads(bv, fv, false, false, false, null, null);
    final List<ReadCount> cnts = featureReadCount.get("PFA0110w");
    
    ReadCount c = cnts.get(0);
    assertTrue(1494.f == c.senseCnt);
    assertTrue(997.f == c.antiCnt);
  }
  
  @Test
  /**
   * Tes the read count for a gene not including those in the intron
   */
  public void rpkm()
  {
    // ignore if in headless mode with no x11
    if(GraphicsEnvironment.getLocalGraphicsEnvironment().isHeadless())
      return;

    String refName = (String) bv.getCombo().getSelectedItem();
    int thisLength = bv.getSeqLengths().get(refName);
    int mappedReads[] = BamUtils.calc(bv, refName, thisLength, 
        false, null);

    Hashtable<String, List<ReadCount>> featureReadCount =
        BamUtils.calculateMappedReads(bv, fv, false, true, false, mappedReads, null);
    final List<ReadCount> cnts = featureReadCount.get("PFA0110w");
    
    ReadCount c = cnts.get(0);
    assertTrue(183768.719f == c.senseCnt);
    assertTrue(122676.375f == c.antiCnt);
  }
}
