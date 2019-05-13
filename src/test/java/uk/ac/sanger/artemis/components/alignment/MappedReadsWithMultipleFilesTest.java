/* 
 * This file is part of Artemis
 *
 * Copyright (C) 2018  Genome Research Limited
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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.awt.GraphicsEnvironment;
import java.net.URL;
import java.util.Hashtable;
import java.util.List;

import uk.ac.sanger.artemis.io.Utils;

import org.junit.After;
import org.junit.Test;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.components.EntryEdit;
import uk.ac.sanger.artemis.components.alignment.BamUtils;
import uk.ac.sanger.artemis.components.alignment.BamView;
import uk.ac.sanger.artemis.components.alignment.ReadCount;

/**
 * Unit test to exercise RPKM counts,
 * including use of multiple hand-crafted bam files.
 * 
 * @author kp11
 *
 */
public class MappedReadsWithMultipleFilesTest 
{
  private BamView bv;
  private FeatureVector fv;
  
  private EntryEdit ee; 
  
  /** 
   * Gene EBA181 (reverse feature).
   * Base range: 110984..111057,111196..111284,111401..111479,111572..116033 
   */
  private static final String TEST_GENE_ID_EBA181 = "PFA0125c";
  
  /** 
   * Gene EBA181 (forward feature).
   * Base range: 154410..156754,156891..159200,159375..159468,159553..159665,159790..159879,159990..160170,160322..160660 
   */
  private static final String TEST_GENE_ID_PFA0175W = "PFA0175w";

  /** Number of bases in EBA181. */
  private static final int TEST_GENE_EBA181_LEN = (111057 - 110984 + 1) + 
													(111284 - 111196 + 1) + 
													(111479 - 111401 + 1) + 
													(116033 - 111572 + 1);
  
  /** Number of bases in PFA0175W. */
  private static final int TEST_GENE_PFA0175W_LEN = (156754 - 154410 + 1) + 
													(159200 - 156891 + 1) + 
													(159468 - 159375 + 1) + 
													(159665 - 159553 + 1) +
													(159879 - 159790 + 1) + 
													(160170 - 159990 + 1) +
													(160660 - 160322 + 1);
  
  private static final int FILE1_TEST_GENE_ID_EBA181_NUM_MAPPED_READS = 3;
  private static final int FILE2_TEST_GENE_ID_EBA181_NUM_MAPPED_READS = 8;
  
  private static final int FILE1_TEST_GENE_ID_PFA0175W_NUM_MAPPED_READS = 4;
  private static final int FILE2_TEST_GENE_ID_PFA0175W_NUM_MAPPED_READS = 1;
  
  private static final int FILE1_TOTAL_NUM_MAPPED_READS = 
		FILE1_TEST_GENE_ID_EBA181_NUM_MAPPED_READS + FILE1_TEST_GENE_ID_PFA0175W_NUM_MAPPED_READS;
  
  private static final int FILE2_TOTAL_NUM_MAPPED_READS = 
  		FILE2_TEST_GENE_ID_EBA181_NUM_MAPPED_READS + FILE2_TEST_GENE_ID_PFA0175W_NUM_MAPPED_READS;
  
  
  /**
   * Load a Bamview with test bam files.
   * 
   * @param useBothFeatures boolean - whether or not to use just one or both test genes
   * @param bamFiles String varargs - one or more bam files
   */
  private void createBamView(boolean useBothFeatures, String ... bamFiles) 
  {
    // ignore if in headless mode with no x11
    if(GraphicsEnvironment.isHeadless())
      return;
    
    int idx = 0;
    StringBuilder buf = new StringBuilder(100);
    for (String bamFile : bamFiles)
    {
    	URL entryFile = MappedReadsWithMultipleFilesTest.class.getResource(bamFile);
    	buf.append(entryFile.getFile());
    	++idx;
    	
    	if (idx < bamFiles.length)
    	{
    		buf.append(",");
    	}
    }
    
    System.setProperty("bam", buf.toString());
    
    final EntryGroup egrp = Utils.getEntryGroup("/data/embl/MAL1.embl.gz");
    ee = new EntryEdit(egrp);
    ee.setVisible(true);

    while( (bv = ee.getJamView()) == null) 
    {
      // wait for BamView to be constructed
      try {
        Thread.sleep(1000);
      } catch(Exception e){};
    }
    
    /*
     * Add gene features.
     */
    fv = new FeatureVector();
    final FeatureVector features = egrp.getAllFeatures();
    for(int i=0; i<features.size(); i++) 
    {
      Feature f = features.elementAt(i);

      // Add the test gene(s)
      
      if(f.getSystematicName().equals(TEST_GENE_ID_EBA181))
    	  fv.add(f);

      if(useBothFeatures && f.getSystematicName().equals(TEST_GENE_ID_PFA0175W))
    	  fv.add(f);
    }
    
  }
  
  /**
   * Cleanup by closing the GUI window
   * after each test.
   */
  @After
  public void cleanUp()
  {
	  if (ee != null)
	  {
		  ee.dispose();     
	  }
	  
	  System.clearProperty("bam");
	  
	  if (fv != null)
	  {
		  fv.removeAllElements();
		  fv = null;
	  }
  }
  
  /**
   * Just a double-check method to test the data set-up.
   */
  private void confirmFeatureTypes()
  {
	  Feature f = fv.elementAt(0);
	  if (f.getSystematicName().equals(TEST_GENE_ID_EBA181))
	  {
		  assertTrue(!f.isForwardFeature());
	  }
	  else
	  {
		  assertTrue(f.isForwardFeature());
	  }
	  
	  if (fv.size() == 2)
	  {
		  f = fv.elementAt(1);
		  if (f.getSystematicName().equals(TEST_GENE_ID_PFA0175W))
		  {
			  assertTrue(f.isForwardFeature());
		  }
		  else
		  {
			  assertTrue(!f.isForwardFeature());
		  }
	  }
	  
  }
  
  /**
   * Test the RPKMs for: single gene feature (EBA181), using bam file 1.
   */
  @Test
  public void testRpkmSingleGeneForBam1()
  {
    // ignore if in headless mode with no x11
    if(GraphicsEnvironment.isHeadless())
      return;
    
    createBamView(false, "/data/rpkm-tests/RPKM_sorted_test_1.bam");

    confirmFeatureTypes();
    
    String refName = (String) bv.getCombo().getSelectedItem();
    int thisLength = bv.getSeqLengths().get(refName);
    
    // Get the total number of mapped reads for the alignment file
    int mappedReads[] = BamUtils.calc(bv, refName, thisLength, false, null);
    assertEquals(1, mappedReads.length);
    assertTrue("Check number of mapped reads", FILE1_TOTAL_NUM_MAPPED_READS == mappedReads[0]);
    
 	// Get the total number of mapped reads for the alignment file
    int mappedReadsByStrand[] = BamUtils.calc(bv, refName, thisLength, true, null);
    assertEquals(1, mappedReadsByStrand.length);
    assertTrue("Check number of mapped reads by strand", 
    		FILE1_TOTAL_NUM_MAPPED_READS == mappedReadsByStrand[0]);
    
    /* ============ Calculate RPKM with default options, for EBA181 ============= */
    
    Hashtable<String, List<ReadCount>> featureRpkms =
        BamUtils.calculateMappedReads(bv, fv, false, true, false, mappedReads, null);
 
    // We're only examining one feature in this test...
    assertEquals("Check no. of feature counts", 1, featureRpkms.size()); 
    
    List<ReadCount> cnts = featureRpkms.get(TEST_GENE_ID_EBA181);
    assertEquals("Check no. of rpkm counts", 1, cnts.size());
    
    ReadCount rpkm = cnts.get(0);
    
    assertTrue("Check rpkm value for single file, introns included (antisense)", 
    		Float.compare( 
    			(3 / (((float) FILE1_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    			rpkm.antiCnt) == 0);
    assertTrue("Check rpkm value for single file, introns included (sense)", 
    		Float.compare(0.0f, rpkm.senseCnt)  == 0);
    
    rpkm = null;
    featureRpkms = null;
    cnts = null;
    
    
    /* ============ Calculate RPKM with intron reads not included, for EBA181 ============= */
    
    featureRpkms = BamUtils.calculateMappedReads(bv, fv, false, false, false, mappedReads, null);
    
    cnts = featureRpkms.get(TEST_GENE_ID_EBA181);
    assertEquals("Check no. of rpkm counts", 1, cnts.size());
    
    rpkm = cnts.get(0);
    
    assertTrue("Check rpkm value for single file, introns excluded (antisense)", 
    		Float.compare(
    				(3 / (((float) FILE1_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.antiCnt) == 0);  
    assertTrue("Check rpkm value for single file, introns excluded (sense)", 
    		Float.compare(0.0f, rpkm.senseCnt) == 0);
    
    rpkm = null;
    featureRpkms = null;
    cnts = null;
    
    
    /* ============ Calculate RPKM [useStrandTag=true], with intron reads not included, for EBA181 ============= */
    
    // Set useStrandTag=true (picks up XS attributes)...
    featureRpkms = BamUtils.calculateMappedReads(bv, fv, false, false, true, mappedReads, null);
    
    cnts = featureRpkms.get(TEST_GENE_ID_EBA181);
    
    rpkm = cnts.get(0);
    assertEquals("Check no. of rpkm counts", 1, cnts.size());
 
    // EBA181 is a reverse feature so...
    assertTrue("Check rpkm value for single file feature, by strand, introns excluded (sense)",
    		Float.compare(
    				(2 / (((float) FILE1_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.senseCnt) == 0);  
    assertTrue("Check rpkm value for single file feature, by strand, introns excluded (antisense)", 
    		Float.compare(
    				(1 / (((float) FILE1_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.antiCnt) == 0);

  }
  
  /**
   * Test RPKMs for: single gene feature (EBA181), using bam file 2.
   */
  @Test
  public void testRpkmSingleGeneForBam2()
  {
	// ignore if in headless mode with no x11
    if(GraphicsEnvironment.isHeadless())
      return;
    
    createBamView(false, "/data/rpkm-tests/RPKM_sorted_test_2.bam");

    String refName = (String) bv.getCombo().getSelectedItem();
    int thisLength = bv.getSeqLengths().get(refName);
    
    // Get the total number of mapped reads for the alignment file
    int mappedReads[] = BamUtils.calc(bv, refName, thisLength, false, null);
    assertEquals(1, mappedReads.length);
    assertTrue("Check number of mapped reads", FILE2_TOTAL_NUM_MAPPED_READS == mappedReads[0]);
    
    
    /* ============ Calculate RPKM with default options, for EBA181 ============= */
    
    Hashtable<String, List<ReadCount>> featureRpkms =
        BamUtils.calculateMappedReads(bv, fv, false, true, false, mappedReads, null);
    
    // We're only examining one feature in this test...
    assertEquals("Check no. of feature counts", 1, featureRpkms.size());
    
    List<ReadCount> cnts = featureRpkms.get(TEST_GENE_ID_EBA181);
    assertEquals("Check no. of rpkm counts", 1, cnts.size());
    
    ReadCount rpkm = cnts.get(0);
  
    assertTrue("Check rpkm value for single file, introns included (antisense)", 
    		Float.compare(
    				(8 / (((float) FILE2_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.antiCnt) == 0);
    assertTrue("Check rpkm value for single file, introns included (sense)", 
    		Float.compare(0.0f, rpkm.senseCnt) == 0);
    
    rpkm = null;
    featureRpkms = null;
    cnts = null;
    
    
    /* ============ Calculate RPKM with intron reads not included, for EBA181 ============= */
    
    featureRpkms = BamUtils.calculateMappedReads(bv, fv, false, false, false, mappedReads, null);
     
    cnts = featureRpkms.get(TEST_GENE_ID_EBA181);
    assertEquals("Check no. of rpkm counts", 1, cnts.size());
    
    rpkm = cnts.get(0);
      
    assertTrue("Check rpkm value for single file, introns excluded (antisense)", 
    		Float.compare(
    				(5 / (((float) FILE2_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.antiCnt) == 0);
    assertTrue("Check rpkm value for single file, introns included (sense)", 
    		Float.compare(0.0f, rpkm.senseCnt) == 0);
 
  }
  
  /**
   * Test RPKMs for: single gene feature (EBA181), two bam files.
   */
  @Test
  public void testRpkmSingleGeneForBam1AndBam2()
  {
	// ignore if in headless mode with no x11
    if(GraphicsEnvironment.isHeadless())
      return;
    
    createBamView(false, "/data/rpkm-tests/RPKM_sorted_test_1.bam", "/data/rpkm-tests/RPKM_sorted_test_2.bam");

    String refName = (String) bv.getCombo().getSelectedItem();
    int thisLength = bv.getSeqLengths().get(refName);
    
    // Get the total number of mapped reads for the alignment file
    int mappedReads[] = BamUtils.calc(bv, refName, thisLength, false, null);
    assertEquals(2, mappedReads.length);
    
    // mappedReads is indexed by file number.
    assertTrue("Check number of mapped reads for bam 1", FILE1_TOTAL_NUM_MAPPED_READS == mappedReads[0]);
    assertTrue("Check number of mapped reads for bam 2", FILE2_TOTAL_NUM_MAPPED_READS == mappedReads[1]);
    
    
    /* ============ Calculate RPKM with default options, for EBA181 ============= */
    
    Hashtable<String, List<ReadCount>> featureRpkms =
        BamUtils.calculateMappedReads(bv, fv, false, true, false, mappedReads, null);
    
    List<ReadCount> cnts = featureRpkms.get(TEST_GENE_ID_EBA181);
    assertEquals("Check no. of rpkm counts", 2, cnts.size());
    
    ReadCount rpkm = cnts.get(0);
  
    assertTrue("Check rpkm value for two bams, introns included, file 1 (antisense)", 
    		Float.compare(
    				(3 / (((float) FILE1_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.antiCnt) == 0);
    assertTrue("Check rpkm value for two bams, introns included, file 1 (sense)", 
    		Float.compare(0.0f, rpkm.senseCnt) == 0);
    
	rpkm = cnts.get(1);

    assertTrue("Check rpkm value for two bams, introns included, file 2 (antisense)", 
    		Float.compare(
    				(8 / (((float) FILE2_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.antiCnt) == 0);
    assertTrue("Check rpkm value for two bams, introns included, file 2 (sense)", 
    		Float.compare(0.0f, rpkm.senseCnt) == 0);
    
    rpkm = null;
    featureRpkms = null;
    cnts = null;
    
    
    /* ============ Calculate RPKM with intron reads not included, for EBA181 ============= */
    
    featureRpkms = BamUtils.calculateMappedReads(bv, fv, false, false, false, mappedReads, null);
     
    cnts = featureRpkms.get(TEST_GENE_ID_EBA181);
    assertEquals("Check no. of rpkm counts", 2, cnts.size());
    
    rpkm = cnts.get(0);
    
    assertTrue("Check rpkm value for two bams, introns excluded, file 1 (antisense)", 
    		Float.compare(
    				(3 / (((float) FILE1_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.antiCnt) == 0);
    assertTrue("Check rpkm value for two bams, introns included, file 1 (sense)", 
    		Float.compare(0.0f, rpkm.senseCnt) == 0);
    
    rpkm = cnts.get(1);
 
    assertTrue("Check rpkm value for two bams, introns excluded, file 2 (antisense)", 
    		Float.compare(
    				(5 / (((float) FILE2_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.antiCnt) == 0);
    assertTrue("Check rpkm value for two bams, introns included, file 2 (sense)", 
    		Float.compare(0.0f, rpkm.senseCnt) == 0);
    
    rpkm = null;
    featureRpkms = null;
    cnts = null;
    
    
    /* ============ Calculate RPKM [useStrandTag=true], with intron reads not included, for EBA181 ============= */
    
    // Set useStrandTag=true (picks up XS attributes)...
    featureRpkms = BamUtils.calculateMappedReads(bv, fv, false, false, true, mappedReads, null);
    
    cnts = featureRpkms.get(TEST_GENE_ID_EBA181);
    assertEquals("Check no. of rpkm counts", 2, cnts.size());
    
    rpkm = cnts.get(0);
      
    // EBA181 is a reverse feature so...
    assertTrue("Check rpkm value for two bams, by strand, introns excluded (sense)", 
    		Float.compare(
    				(2 / (((float) FILE1_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.senseCnt) == 0);  
    assertTrue("Check rpkm value for two bams, by strand, introns excluded (antisense)", 
    		Float.compare(
    				(1 / (((float) FILE1_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.antiCnt) == 0);
    
    rpkm = cnts.get(1);
   
    // EBA181 is a reverse feature so...
    assertTrue("Check rpkm value for two bams, by strand, introns excluded (sense)", 
    		Float.compare(
    				(1 / (((float) FILE2_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.senseCnt) == 0);  
    assertTrue("Check rpkm value for two bams, by strand, introns excluded (antisense)", 
    		Float.compare(
    				(4 / (((float) FILE2_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.antiCnt) == 0);
  }
  
  /**
   * Test RPKMs for: two gene features (EBA181 and PFA0175W), two bam files. 
   */
  @Test
  public void testRpkmMultiGenesForBam1AndBam2()
  {
	// ignore if in headless mode with no x11
    if(GraphicsEnvironment.isHeadless())
      return;
    
    // Use both features
    createBamView(true, "/data/rpkm-tests/RPKM_sorted_test_1.bam", "/data/rpkm-tests/RPKM_sorted_test_2.bam");

    String refName = (String) bv.getCombo().getSelectedItem();
    int thisLength = bv.getSeqLengths().get(refName);
    
    // Get the total number of mapped reads for the alignment file
    int mappedReads[] = BamUtils.calc(bv, refName, thisLength, false, null);
    assertEquals(2, mappedReads.length);
    
    // mappedReads is indexed by file number.
    assertTrue("Check number of mapped reads for bam 1", FILE1_TOTAL_NUM_MAPPED_READS == mappedReads[0]);
    assertTrue("Check number of mapped reads for bam 2", FILE2_TOTAL_NUM_MAPPED_READS == mappedReads[1]);
    
    
    /* ============ Calculate RPKM with default options, for EBA181 ============= */
    
    Hashtable<String, List<ReadCount>> featureRpkms =
        BamUtils.calculateMappedReads(bv, fv, false, true, false, mappedReads, null);
    
    List<ReadCount> cnts = featureRpkms.get(TEST_GENE_ID_EBA181);
    assertEquals("Check no. of rpkm counts", 2, cnts.size());
    
    ReadCount rpkm = cnts.get(0);
   
    assertTrue("Check rpkm value for two bams, introns included, file 1 (antisense)", 
    		Float.compare(
    				(3 / (((float) FILE1_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.antiCnt) == 0);
    assertTrue("Check rpkm value for two bams, introns included, file 1 (sense)", 
    		Float.compare(0.0f, rpkm.senseCnt) == 0);
    
	rpkm = cnts.get(1);
     
    assertTrue("Check rpkm value for two bams, introns included, file 2 (antisense)", 
    		Float.compare(
    				(8 / (((float) FILE2_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.antiCnt) == 0);
    assertTrue("Check rpkm value for two bams, introns included, file 2 (sense)", 
    		Float.compare(0.0f, rpkm.senseCnt) == 0);
    
    rpkm = null;
    featureRpkms = null;
    cnts = null;
    
    
    /* ============ Calculate RPKM with intron reads not included, for EBA181 ============= */
    
    featureRpkms = BamUtils.calculateMappedReads(bv, fv, false, false, false, mappedReads, null);
     
    cnts = featureRpkms.get(TEST_GENE_ID_EBA181);
    assertEquals("Check no. of rpkm counts", 2, cnts.size());
    
    rpkm = cnts.get(0);
 
    assertTrue("Check rpkm value for two bams, introns excluded, file 1 (antisense)", 
    		Float.compare(
    				(3 / (((float) FILE1_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.antiCnt) == 0);
    assertTrue("Check rpkm value for two bams, introns included, file 1 (sense)", 
    		Float.compare(0.0f, rpkm.senseCnt) == 0);
    
    rpkm = cnts.get(1);
    
    assertTrue("Check rpkm value for two bams, introns excluded, file 2 (antisense)", 
    		Float.compare(
    				(5 / (((float) FILE2_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.antiCnt) == 0);
    assertTrue("Check rpkm value for two bams, introns included, file 2 (sense)", 
    		Float.compare(0.0f, rpkm.senseCnt) == 0);
    
    rpkm = null;
    featureRpkms = null;
    cnts = null;
    
    /* ============ Calculate RPKM with default options, for PFA0175w ============= */
    
    featureRpkms =
        BamUtils.calculateMappedReads(bv, fv, false, true, false, mappedReads, null);
    
    cnts = featureRpkms.get(TEST_GENE_ID_PFA0175W);
    assertEquals("Check no. of rpkm counts", 2, cnts.size());
    
    rpkm = cnts.get(0);
     
    // Reverse read specified by sam flag value...
    assertTrue("Check rpkm value for two bams, introns included, file 1 (antisense)", 
    		Float.compare(
    				(1 / (((float) FILE1_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_PFA0175W_LEN / 1000.0f))),
    				rpkm.antiCnt) == 0);
    
    assertTrue("Check rpkm value for two bams, introns included, file 1 (sense)", 
    		Float.compare(
    				(3 / (((float) FILE1_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_PFA0175W_LEN / 1000.0f))),
    				rpkm.senseCnt) == 0);
    
	rpkm = cnts.get(1);
    
    assertTrue("Check rpkm value for two bams, introns included, file 2 (antisense)", 
    		Float.compare(
    				(1 / (((float) FILE2_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_PFA0175W_LEN / 1000.0f))),
    				rpkm.antiCnt) == 0);
    assertTrue("Check rpkm value for two bams, introns included, file 2 (sense)", 
    		Float.compare(0.0f, rpkm.senseCnt) == 0);
    
    rpkm = null;
    featureRpkms = null;
    cnts = null;
    
    
    /* ============ Calculate RPKM with intron reads not included, for PFA0175w ============= */
    
    featureRpkms = BamUtils.calculateMappedReads(bv, fv, false, false, false, mappedReads, null);
     
    cnts = featureRpkms.get(TEST_GENE_ID_PFA0175W);
    assertEquals("Check no. of rpkm counts", 2, cnts.size());
    
    rpkm = cnts.get(0);
    
    assertTrue("Check rpkm value for two bams, introns excluded, file 1 (antisense)", 
    		Float.compare(
    				(1 / (((float) FILE1_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_PFA0175W_LEN / 1000.0f))),
    				rpkm.antiCnt) == 0);
    assertTrue("Check rpkm value for two bams, introns excluded, file 1 (sense)", 
    		Float.compare(
    				(3 / (((float) FILE1_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_PFA0175W_LEN / 1000.0f))),
    				rpkm.senseCnt) == 0);
    
    rpkm = cnts.get(1);
    
    // No mapped reads that fall in exon regions.
    assertTrue("Check rpkm value for two bams, introns excluded, file 2 (sense)", 
    		Float.compare(0.0f, rpkm.senseCnt) == 0);
    assertTrue("Check rpkm value for two bams, introns excluded, file 2 (antisense)", 
    		Float.compare(0.0f, rpkm.antiCnt) == 0);
    
    rpkm = null;
    featureRpkms = null;
    cnts = null;
    
    
    /* ============ Calculate RPKM (reads entirely contained) with intron reads not included, for EBA181 ============= */
    
    /* 
     * Set the contained flag to true in the following call. This instructs htsjdk to only return
     * reads that fall entirely within the feature region (i.e start base to last base).
     * This is used in the "no overlaps" Artemis option.
     */
    featureRpkms = BamUtils.calculateMappedReads(bv, fv, true, false, false, mappedReads, null);
     
    cnts = featureRpkms.get(TEST_GENE_ID_EBA181);
    assertEquals("Check no. of rpkm counts", 2, cnts.size());
    
    rpkm = cnts.get(0);

    assertTrue("Check rpkm value for two bams, introns excluded, no overlap, file 1 (antisense)", 
    		Float.compare(
    				(3 / (((float) FILE1_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.antiCnt) == 0);
    assertTrue("Check rpkm value for two bams, introns included, no overlap, file 1 (sense)", 
    		Float.compare(0.0f, rpkm.senseCnt) == 0);
    
    rpkm = cnts.get(1);
    
    // File 2 contains one read that overlaps the end of the gene.
    assertTrue("Check rpkm value for two bams, introns excluded, no overlap, file 2 (antisense)", 
    		Float.compare(
    				(4 / (((float) FILE2_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.antiCnt) == 0);
    assertTrue("Check rpkm value for two bams, introns included, no overlap, file 2 (sense)", 
    		Float.compare(0.0f, rpkm.senseCnt) == 0);
    
    rpkm = null;
    featureRpkms = null;
    cnts = null;
    
    
    /* ============ Calculate RPKM [useStrandTag=true], with intron reads not included, for EBA181 ============= */
    
    // Set useStrandTag=true (picks up XS attributes)...
    featureRpkms = BamUtils.calculateMappedReads(bv, fv, false, false, true, mappedReads, null);
    
    cnts = featureRpkms.get(TEST_GENE_ID_EBA181);
    assertEquals("Check no. of rpkm counts", 2, cnts.size());
    
    rpkm = cnts.get(0);
        
    // EBA181 is a reverse feature so...
    assertTrue("Check rpkm value for two bams, by strand, introns excluded (sense)", 
    		Float.compare(
    				(2 / (((float) FILE1_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.senseCnt) == 0);  
    assertTrue("Check rpkm value for two bams, by strand, introns excluded (antisense)", 
    		Float.compare(
    				(1 / (((float) FILE1_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.antiCnt) == 0);
    
    rpkm = cnts.get(1);
      
    // EBA181 is a reverse feature so...
    assertTrue("Check rpkm value for two bams, by strand, introns excluded (sense)", 
    		Float.compare(
    				(1 / (((float) FILE2_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.senseCnt) == 0);  
    assertTrue("Check rpkm value for two bams, by strand, introns excluded (antisense)", 
    		Float.compare(
    				(4 / (((float) FILE2_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_EBA181_LEN / 1000.0f))),
    				rpkm.antiCnt) == 0);
    
    rpkm = null;
    featureRpkms = null;
    cnts = null;
    
    
    /* ============ Calculate RPKM [useStrandTag=true], with intron reads not included, for PFA0175W ============= */
    
    // Set useStrandTag=true (picks up XS attributes)...
    featureRpkms = BamUtils.calculateMappedReads(bv, fv, false, false, true, mappedReads, null);
    
    cnts = featureRpkms.get(TEST_GENE_ID_PFA0175W);
    assertEquals("Check no. of rpkm counts", 2, cnts.size());
    
    rpkm = cnts.get(0);
    
    // PFA0175W is a forward feature so...
    assertTrue("Check rpkm value for two bams, by strand, introns excluded (sense)", 
    		Float.compare(
    				(3 / (((float) FILE1_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_PFA0175W_LEN / 1000.0f))),
    				rpkm.senseCnt) == 0);  
    assertTrue("Check rpkm value for two bams, by strand, introns excluded (antisense)", 
    		Float.compare(
    				(1 / (((float) FILE1_TOTAL_NUM_MAPPED_READS / 1000000.0f) * (TEST_GENE_PFA0175W_LEN / 1000.0f))),
    				rpkm.antiCnt) == 0);
    
    rpkm = cnts.get(1);
    
    // PFA0175W is a forward feature so. We expect zero values as the only read is in an intron region.
    assertTrue("Check rpkm value for two bams, by strand, introns excluded (sense)", 
    		Float.compare(0.0f, rpkm.senseCnt) == 0);  
    assertTrue("Check rpkm value for two bams, by strand, introns excluded (antisense)", 
    		Float.compare(0.0f, rpkm.antiCnt) == 0);
    
  }
  
  /**
   * And do some testing of pure read counts for: 
   * two gene features (EBA181 and PFA0175W), two bam files.
   */
  @Test
  public void testCountsMultiGenesForBam1AndBam2()
  {
	// ignore if in headless mode with no x11
    if(GraphicsEnvironment.isHeadless())
      return;
    
    // Use both features
    createBamView(true, "/data/rpkm-tests/RPKM_sorted_test_1.bam", "/data/rpkm-tests/RPKM_sorted_test_2.bam");
 
    
    /* ============ Calculate the counts with default options, for EBA181 ============= */
    
    Hashtable<String, List<ReadCount>> featureCounts =
        BamUtils.calculateMappedReads(bv, fv, false, true, false, null, null);
    
    List<ReadCount> cnts = featureCounts.get(TEST_GENE_ID_EBA181);
    assertEquals("Check no. of counts", 2, cnts.size());
    
    ReadCount readCounts = cnts.get(0);
   
    assertTrue("Check counts for two bams, introns included, file 1 (antisense)", 
    		Float.compare(3.0f, readCounts.antiCnt) == 0);
    assertTrue("Check counts for two bams, introns included, file 1 (sense)", 
    		Float.compare(0.0f, readCounts.senseCnt) == 0);
    
    readCounts = cnts.get(1);
     
    assertTrue("Check counts for two bams, introns included, file 2 (antisense)", 
    		Float.compare(8.0f, readCounts.antiCnt) == 0);
    assertTrue("Check counts for two bams, introns included, file 2 (sense)", 
    		Float.compare(0.0f, readCounts.senseCnt) == 0);
    
    readCounts = null;
    featureCounts = null;
    cnts = null;
    
    
    /* ============ Calculate the counts with intron reads not included, for EBA181 ============= */
    
    featureCounts = BamUtils.calculateMappedReads(bv, fv, false, false, false, null, null);
     
    cnts = featureCounts.get(TEST_GENE_ID_EBA181);
    assertEquals("Check no. of counts", 2, cnts.size());
    
    readCounts = cnts.get(0);
 
    assertTrue("Check counts for two bams, introns excluded, file 1 (antisense)", 
    		Float.compare(3.0f, readCounts.antiCnt) == 0);
    assertTrue("Check counts for two bams, introns included, file 1 (sense)", 
    		Float.compare(0.0f, readCounts.senseCnt) == 0);
    
    readCounts = cnts.get(1);
    
    assertTrue("Check counts for two bams, introns excluded, file 2 (antisense)", 
    		Float.compare(5.0f, readCounts.antiCnt) == 0);
    assertTrue("Check counts for two bams, introns included, file 2 (sense)", 
    		Float.compare(0.0f, readCounts.senseCnt) == 0);
    
    readCounts = null;
    featureCounts = null;
    cnts = null;
    
    /* ============ Calculate the counts with default options, for PFA0175w ============= */
    
    featureCounts =
        BamUtils.calculateMappedReads(bv, fv, false, true, false, null, null);
    
    cnts = featureCounts.get(TEST_GENE_ID_PFA0175W);
    assertEquals("Check no. of counts", 2, cnts.size());
    
    readCounts = cnts.get(0);
     
    // Reverse read specified by sam flag value...
    assertTrue("Check counts for two bams, introns included, file 1 (antisense)", 
    		Float.compare(1.0f, readCounts.antiCnt) == 0);
    
    assertTrue("Check counts for two bams, introns included, file 1 (sense)", 
    		Float.compare(3.0f, readCounts.senseCnt) == 0);
    
    readCounts = cnts.get(1);
    
    assertTrue("Check counts for two bams, introns included, file 2 (antisense)", 
    		Float.compare(1.0f, readCounts.antiCnt) == 0);
    assertTrue("Check counts for two bams, introns included, file 2 (sense)", 
    		Float.compare(0.0f, readCounts.senseCnt) == 0);
    
    readCounts = null;
    featureCounts = null;
    cnts = null;
    
    
    /* ============ Calculate the counts with intron reads not included, for PFA0175w ============= */
    
    featureCounts = BamUtils.calculateMappedReads(bv, fv, false, false, false, null, null);
     
    cnts = featureCounts.get(TEST_GENE_ID_PFA0175W);
    assertEquals("Check no. of counts", 2, cnts.size());
    
    readCounts = cnts.get(0);
    
    assertTrue("Check counts for two bams, introns excluded, file 1 (antisense)", 
    		Float.compare(1.0f, readCounts.antiCnt) == 0);
    assertTrue("Check counts for two bams, introns excluded, file 1 (sense)", 
    		Float.compare(3.0f, readCounts.senseCnt) == 0);
    
    readCounts = cnts.get(1);
    
    // No mapped reads that fall in exon regions.
    assertTrue("Check counts for two bams, introns excluded, file 2 (sense)", 
    		Float.compare(0.0f, readCounts.senseCnt) == 0);
    assertTrue("Check counts for two bams, introns excluded, file 2 (antisense)", 
    		Float.compare(0.0f, readCounts.antiCnt) == 0);
    
    readCounts = null;
    featureCounts = null;
    cnts = null;
    
    
    /* ============ Calculate the counts (reads entirely contained) with intron reads not included, for EBA181 ============= */
    
    /* 
     * Set the contained flag to true in the following call. This instructs htsjdk to only return
     * reads that fall entirely within the feature region (i.e start base to last base).
     * This is used in the "no overlaps" Artemis option.
     */
    featureCounts = BamUtils.calculateMappedReads(bv, fv, true, false, false, null, null);
     
    cnts = featureCounts.get(TEST_GENE_ID_EBA181);
    assertEquals("Check no. of counts", 2, cnts.size());
    
    readCounts = cnts.get(0);

    assertTrue("Check counts for two bams, introns excluded, no overlap, file 1 (antisense)", 
    		Float.compare(3.0f,  readCounts.antiCnt) == 0);
    assertTrue("Check counts for two bams, introns included, no overlap, file 1 (sense)", 
    		Float.compare(0.0f, readCounts.senseCnt) == 0);
    
    readCounts = cnts.get(1);
    
    // File 2 contains one read that overlaps the end of the gene.
    assertTrue("Check counts for two bams, introns excluded, no overlap, file 2 (antisense)", 
    		Float.compare(4.0f, readCounts.antiCnt) == 0);
    assertTrue("Check counts for two bams, introns included, no overlap, file 2 (sense)", 
    		Float.compare(0.0f, readCounts.senseCnt) == 0);
    
    readCounts = null;
    featureCounts = null;
    cnts = null;
    
    
    /* ============ Calculate the counts [useStrandTag=true], with intron reads not included, for EBA181 ============= */
    
    // Set useStrandTag=true (picks up XS attributes)...
    featureCounts = BamUtils.calculateMappedReads(bv, fv, false, false, true, null, null);
    
    cnts = featureCounts.get(TEST_GENE_ID_EBA181);
    assertEquals("Check no. of counts", 2, cnts.size());
    
    readCounts = cnts.get(0);
        
    // EBA181 is a reverse feature so...
    assertTrue("Check counts for two bams, by strand, introns excluded (sense)", 
    		Float.compare(2.0f, readCounts.senseCnt) == 0);  
    assertTrue("Check counts for two bams, by strand, introns excluded (antisense)", 
    		Float.compare(1.0f, readCounts.antiCnt) == 0);
    
    readCounts = cnts.get(1);
      
    // EBA181 is a reverse feature so...
    assertTrue("Check counts for two bams, by strand, introns excluded (sense)", 
    		Float.compare(1.0f, readCounts.senseCnt) == 0);  
    assertTrue("Check counts for two bams, by strand, introns excluded (antisense)", 
    		Float.compare(4.0f, readCounts.antiCnt) == 0);
    
    readCounts = null;
    featureCounts = null;
    cnts = null;
    
    
    /* ============ Calculate the counts [useStrandTag=true], with intron reads not included, for PFA0175W ============= */
    
    // Set useStrandTag=true (picks up XS attributes)...
    featureCounts = BamUtils.calculateMappedReads(bv, fv, false, false, true, null, null);
    
    cnts = featureCounts.get(TEST_GENE_ID_PFA0175W);
    assertEquals("Check no. of counts", 2, cnts.size());
    
    readCounts = cnts.get(0);
    
    // PFA0175W is a forward feature so...
    assertTrue("Check counts for two bams, by strand, introns excluded (sense)", 
    		Float.compare(3.0f, readCounts.senseCnt) == 0);  
    assertTrue("Check counts for two bams, by strand, introns excluded (antisense)", 
    		Float.compare(1.0f, readCounts.antiCnt) == 0);
    
    readCounts = cnts.get(1);
    
    // PFA0175W is a forward feature so. We expect zero values as the only read is in an intron region.
    assertTrue("Check counts for two bams, by strand, introns excluded (sense)", 
    		Float.compare(0.0f, readCounts.senseCnt) == 0);  
    assertTrue("Check counts for two bams, by strand, introns excluded (antisense)", 
    		Float.compare(0.0f, readCounts.antiCnt) == 0);
    
  }

}
