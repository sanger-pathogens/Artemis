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

import static org.junit.Assert.*;

import org.junit.Test;

import java.io.File;
import java.net.URL;

import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.cram.CRAIIndex;
import junit.framework.AssertionFailedError;


/**
 * JUnit tests for the BamUtils class.
 * 
 * @author kp11
 *
 */
public class BamUtilsTest {

	@Test
	public void testIsCramFile() 
	{
		assertTrue(BamUtils.isCramFile("/a/file/path/file.cram"));
		assertTrue(BamUtils.isCramFile("/a/file/path/file.cram.cram"));
		assertTrue(BamUtils.isCramFile("file.cram"));
		assertTrue(BamUtils.isCramFile("http://mywebsite/file.cram"));
		assertTrue(BamUtils.isCramFile("https://mywebsite/file.cram"));
		assertTrue(BamUtils.isCramFile("ftp://myftpsite/file.cram"));
		
		assertFalse(BamUtils.isCramFile("/a/file/path/file.bam"));
		assertFalse(BamUtils.isCramFile("/a/file/path/file"));
		assertFalse(BamUtils.isCramFile("http://mywebsite/file"));
		assertFalse(BamUtils.isCramFile("ftp://myftpsite/file"));
	}
	
	@Test
	public void testIsBamFile() 
	{
		assertTrue(BamUtils.isBamFile("/a/file/path/file.bam"));
		assertTrue(BamUtils.isBamFile("/a/file/path/file.cram.bam"));
		assertTrue(BamUtils.isBamFile("file.bam"));
		assertTrue(BamUtils.isBamFile("http://mywebsite/file.bam"));
		assertTrue(BamUtils.isBamFile("https://mywebsite/file.bam"));
		assertTrue(BamUtils.isBamFile("ftp://myftpsite/file.bam"));
		
		assertFalse(BamUtils.isBamFile("/a/file/path/file.cram"));
		assertFalse(BamUtils.isBamFile("/a/file/path/file"));
		assertFalse(BamUtils.isBamFile("http://mywebsite/file"));
		assertFalse(BamUtils.isBamFile("ftp://myftpsite/file"));
	}
	
	@Test
	public void testConstructIndexFilename() 
	{
		assertEquals("/a/file/path/file.bam.bai", BamUtils.constructIndexFilename("/a/file/path/file.bam"));
		assertEquals("http://mywebsite/file.bam.bai", BamUtils.constructIndexFilename("http://mywebsite/file.bam"));
		assertEquals("https://mywebsite/file.bam.bai", BamUtils.constructIndexFilename("https://mywebsite/file.bam"));
		assertEquals("ftp://myftpsite/file.bam.bai", BamUtils.constructIndexFilename("ftp://myftpsite/file.bam"));
		
		assertEquals("/a/file/path/file.cram.crai", BamUtils.constructIndexFilename("/a/file/path/file.cram"));
		assertEquals("http://mywebsite/file.cram.crai", BamUtils.constructIndexFilename("http://mywebsite/file.cram"));
		assertEquals("https://mywebsite/file.cram.crai", BamUtils.constructIndexFilename("https://mywebsite/file.cram"));
		assertEquals("ftp://myftpsite/file.cram.crai", BamUtils.constructIndexFilename("ftp://myftpsite/file.cram"));
		
		assertNull(BamUtils.constructIndexFilename("ftp://myftpsite/file.jpeg"));
		
	}
	
	@Test
	public void testCreateTempIndexFile() throws Exception
	{
		File indexFile = null;
				
		indexFile = BamUtils.createTempIndexFile(new URL("file:/a/file/path/file.bam.bai"));
		assertTrue(indexFile.exists());
		
		indexFile = BamUtils.createTempIndexFile(new URL("http://mywebsite/file.bam.bai"));
		assertTrue(indexFile.exists());
		
		indexFile = BamUtils.createTempIndexFile(new URL("ftp://myftpsite/file.bam.bai"));
		assertTrue(indexFile.exists());
		
		indexFile = BamUtils.createTempIndexFile(new URL("file:/a/file/path/file.cram.crai"));
		assertTrue(indexFile.exists());
		
		indexFile = BamUtils.createTempIndexFile(new URL("http://mywebsite/file.cram.crai"));
		assertTrue(indexFile.exists());
		
		indexFile = BamUtils.createTempIndexFile(new URL("ftp://myftpsite/file.cram.crai"));
		assertTrue(indexFile.exists());
	}
	
	@Test
	public void testCreateIndexFileFromScratch() throws Exception
	{
		File fileBai = File.createTempFile("unit_test_bam_index_file1", ".bai");
		fileBai.deleteOnExit();
		assertTrue(fileBai.exists());
		assertTrue(fileBai.length()==0);
		BamUtils.createIndexFileFromScratch("data/MAL_8h.bam", fileBai);
		assertTrue(fileBai.length()>0);
		
		File fileCrai = File.createTempFile("unit_test_bam_index_file1", ".crai");
		fileCrai.deleteOnExit();
		assertTrue(fileCrai.exists());
		assertTrue(fileCrai.length()==0);
		BamUtils.createIndexFileFromScratch("data/NV.cram", fileCrai);
		assertTrue(fileCrai.length()>0);
	}
	
	@Test
	public void testSamRecordEqualityCheck()
	{
		byte [] bases = new String("ATGACAAGGGTCAGTTTTAGAGA").getBytes();

		SAMRecord rec1 = new SAMRecord(new SAMFileHeader());
		rec1.setReadName("rec1");
		rec1.setAlignmentStart(10);
		rec1.setReadBases(bases);
		rec1.setFlags(8);

		SAMRecord rec2 = new SAMRecord(new SAMFileHeader());
		rec2.setReadName("rec1");
		rec2.setAlignmentStart(10);
		rec2.setReadBases(bases);
		rec2.setFlags(8);
		
		assertTrue(BamUtils.samRecordEqualityCheck(rec1, rec1));
		assertTrue(BamUtils.samRecordEqualityCheck(rec1, rec2));

		// Check different read name
		rec1.setReadName("another_read");
		assertFalse(BamUtils.samRecordEqualityCheck(rec1, rec2));
		rec1.setReadName("rec1");
		
		// Check different start alignment
		rec1.setAlignmentStart(12);
		assertFalse(BamUtils.samRecordEqualityCheck(rec1, rec2));
		rec1.setAlignmentStart(10);
		
		// Check different flags
		rec1.setFlags(256);
		assertFalse(BamUtils.samRecordEqualityCheck(rec1, rec2));
		rec1.setFlags(8);
	
	}
	
	@Test
	public void testGetIndexFileForFile() throws Exception
	{
		// File ================================================
		
		// Check getting an already existing BAM index file
		//
		File alreadyExistingBamIndexFile = BamUtils.getIndexFile("data/MAL_8h.bam");
		assertNotNull(alreadyExistingBamIndexFile);
		assertTrue(alreadyExistingBamIndexFile.exists());
		assertEquals("MAL_8h.bam.bai", alreadyExistingBamIndexFile.getName());
		
		// Check getting an index file for a BAM file with no index.
		// Should create one.
		//
		File checkFile = new File("data/MAL_8h_noindex.bam.bai");
		if (checkFile.exists())
		{
			// Remove any leftover index from a previous run...
			checkFile.delete();
		}
		try 
		{
			File createdBamIndexFile = BamUtils.getIndexFile("data/MAL_8h_noindex.bam");
			assertNotNull(createdBamIndexFile);
			assertTrue(createdBamIndexFile.exists());
			assertEquals("MAL_8h_noindex.bam.bai", createdBamIndexFile.getName());
			assertTrue(createdBamIndexFile.length()>0);
		} 
		finally
		{
			if (checkFile.exists())
			{
				// Remove the created index file
				checkFile.delete();
			}
		}
		
		// Check getting an index file for a CRAM file with no index.
		// Should create one.
		//
		File createdCramIndexFile = BamUtils.getIndexFile("data/NV.cram");
		assertNotNull(createdCramIndexFile);
		assertTrue(createdCramIndexFile.exists());
		assertEquals("NV.cram.crai", createdCramIndexFile.getName());
		assertTrue(createdCramIndexFile.length()>0);
				
		// Check getting an already existing CRAM index file
		//
		File alreadyExistingCramIndexFile = BamUtils.getIndexFile("data/NV.cram");
		assertNotNull(alreadyExistingCramIndexFile);
		assertTrue(alreadyExistingCramIndexFile.exists());
		assertEquals("NV.cram.crai", alreadyExistingCramIndexFile.getName());
	}
	
	@Test
	public void testGetIndexFileForFtpBam() throws Exception
	{
		System.out.println("testGetIndexFileForFtpBam(): This test creates index files and may be slow to complete");
		
		// FTP Seekable streams ================================================
		// These will block so commented out for normal test runs.
		
		// Check creating an index file for a BAM file with no index for FTP
		// Should create local copy.
		//
		File createdBamIndexFtpFile = null;
		try
		{
			createdBamIndexFtpFile = BamUtils.getIndexFile("ftp://ftp.sanger.ac.uk/pub/project/pathogens/kp11/NV_noindex.bam");
			throw new AssertionFailedError("Expected an exception to be thrown");
		}
		catch (SAMException e)
		{
			assertTrue(e.getMessage().contains("Failed to find an index file"));
		}
		
		assertNull(createdBamIndexFtpFile);
		
		// Check getting the index file for a BAM file with an index for FTP
		// Should create local copy.
		//
		File existingBamIndexFtpFile = BamUtils.getIndexFile("ftp://ftp.sanger.ac.uk/pub/project/pathogens/kp11/NV.bam");
		assertNotNull(existingBamIndexFtpFile);
		assertTrue(existingBamIndexFtpFile.exists());
		assertTrue(existingBamIndexFtpFile.getName().endsWith(BAMIndex.BAMIndexSuffix));
		assertTrue(existingBamIndexFtpFile.length()>0);
	}
	
	@Test
	public void testCreateIndexFileForFtpCram() throws Exception
	{
		// Check creating an index file for a CRAM file with no index for FTP
		// Should not be allowed.
		//
		File createdCramIndexFtpFile = null;
		try 
		{
			createdCramIndexFtpFile = BamUtils.getIndexFile("ftp://ftp.sanger.ac.uk/pub/project/pathogens/kp11/NV_noindex.cram");
			throw new AssertionFailedError("Expected an exception to be thrown");
		} 
		catch (SAMException e)
		{
			assertTrue(e.getMessage().contains("Failed to find an index file"));
		}
		
		assertNull(createdCramIndexFtpFile);
		
	}
	
	@Test
	public void testGetExistingIndexFileForFtpCram() throws Exception
	{
		System.out.println("testGetExistingIndexFileForFtpCram(): This test creates index files and may be slow to complete");
		
		// Check getting the index file for a CRAM file with an index for FTP
		// Should create local copy.
		//
		File existingCramIndexFtpFile = BamUtils.getIndexFile("ftp://ftp.sanger.ac.uk/pub/project/pathogens/kp11/NV.cram");
		assertNotNull(existingCramIndexFtpFile);
		assertTrue(existingCramIndexFtpFile.exists());
		assertTrue(existingCramIndexFtpFile.getName().endsWith(CRAIIndex.CRAI_INDEX_SUFFIX));
		assertTrue(existingCramIndexFtpFile.length()>0);
	}
}
