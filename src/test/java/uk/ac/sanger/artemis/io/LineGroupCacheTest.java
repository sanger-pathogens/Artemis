/* 
 * This file is part of Artemis
 *
 * Copyright (C) 2019  Genome Research Limited
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

import static org.junit.Assert.*;

import java.io.StringReader;
import java.io.StringWriter;
import java.util.List;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;

import uk.ac.sanger.artemis.util.LinePushBackReader;

/**
 * JUnit test for the LineGroupCache class.
 * 
 * @author kp11
 *
 */
public class LineGroupCacheTest
{

	private static final String EMBL_HDR = "FH   Key             Location/Qualifiers\n";
	
	private static final String EMBL_MISC1 = "ID   1    standard; DNA; HTG; 200111453 BP.\n";
	private static final String EMBL_MISC2 = "AC   chromosome:WSB_EiJ_v1:1:1:200111453:1\n";
	
	// Hacked sequence
	private static final String EMBL_SEQUENCE = 
			"SQ   Sequence  10 BP;   11 A;   12 C;   13 G;   14 T;   15 other;\n" + 
			"     actgaaacag ttcctgcctg gctttgctga aaggggtaat aaaatcttct gtgtaaaggc        60\n";
	
	private static final String EMBL_FEATURE = 
			"FT   CDS             1..10\n" + 
			"FT                   /note=test\n";
	
	private static final int NUM_EMBL_HDRS_FOR_TESTS = 4;
	
	private static final int NUM_EMBL_MISC_FEATURES_FOR_TESTS = 2;
	
	private LineGroupCache groups = new LineGroupCache();
	
	// Test line groups
	private LineGroup hdr1;
	private LineGroup hdr2;
	private LineGroup hdr3;
	private LineGroup hdr4;
	private EmblMisc misc1;
	private EmblMisc misc2;
	private EmblStreamFeature feature;
	private FeatureTable ft;
	private StreamSequence sequence;
	
	@Mock 
	private PublicDBDocumentEntry entry;
	
	
	/**
	 * Add some test data to our line group cache.
	 * 
	 * @throws Exception
	 */
	@Before
	public void setUp() throws Exception 
	{
		MockitoAnnotations.initMocks(this);
		
		// Add some headers
		hdr1 = new FeatureHeader(new LinePushBackReader(new StringReader(EMBL_HDR)));
		groups.add(hdr1);
		hdr2 = new FeatureHeader(new LinePushBackReader(new StringReader(EMBL_HDR)));
		groups.add(hdr2);
		hdr3 = new FeatureHeader(new LinePushBackReader(new StringReader(EMBL_HDR)));
		groups.add(hdr3);
		hdr4 = new FeatureHeader(new LinePushBackReader(new StringReader(EMBL_HDR)));
		groups.add(hdr4);
		
		// Add some misc lien groups
		misc1 = new EmblMisc(new LinePushBackReader(new StringReader(EMBL_MISC1)));
		groups.add(misc1);
		misc2 = new EmblMisc(new LinePushBackReader(new StringReader(EMBL_MISC2)));
		groups.add(misc2);
		
		
		// Create a feature table and add a feature
		QualifierVector qualifiers = new QualifierVector();
		ft = new StreamFeatureTable();
		qualifiers.add(new Qualifier("note", "test"));
		feature = new EmblStreamFeature(new Key("CDS"), new Location("1..10"), qualifiers);
		ft.add(feature);
		groups.add(ft);
		
		// Add a sequence
		sequence = new EmblStreamSequence(new LinePushBackReader(new StringReader(EMBL_SEQUENCE)));
		groups.add(sequence);
	}
	
	/**
	 * Clear up at end of each test.
	 */
	@After
	public void tearDown()
	{
		groups = null;
		hdr1 = hdr2 = hdr3 = hdr4 = null;
		ft = null;
		sequence = null;
		
	}
	
	/**
	 * Test the findRelevantList() method.
	 * Also implicitly uses the add() method.
	 * @throws Exception
	 */
	@Test
	public void testFindRelevantList() throws Exception
	{
		// When
		
		List<LineGroup> result = groups.findRelevantList(hdr1);
		
		// Then
		
		assertNotNull( result );
		assertEquals( NUM_EMBL_HDRS_FOR_TESTS, result.size() );
		
		// When
		
		result = groups.findRelevantList(misc1);
		
		// Then
		
		assertNotNull( result );
		assertEquals( NUM_EMBL_MISC_FEATURES_FOR_TESTS, result.size() );
		
	}
	
	/**
	 * Test the getFeatureTable() method.
	 * Also implicitly uses the add() method.
	 */
	@Test
	public void testGetFeatureTable()
	{
		// When
		
		FeatureTable result = groups.getFeatureTable();
		
		// Then
		
		assertNotNull( result );
		assertTrue( result instanceof FeatureTable );
		assertEquals( ft, (FeatureTable)result );
	}
	
	/**
	 * Test the getSequence() method.
	 * Also implicitly uses the add() method.
	 * @throws Exception
	 */
	@Test
	public void testGetSequence() throws Exception
	{
		// When
		
		LineGroup result = groups.getSequence();
		
		// Then
		
		assertNotNull( result );
		assertTrue( result instanceof StreamSequence );
		assertEquals( sequence, (StreamSequence)result );
	}
	
	/**
	 * Additional tests for the add() method.
	 * @throws Exception
	 */
	@Test
	public void testAdd() throws Exception
	{
		// When 
		
		groups.clearMiscLineGroups();
		// This should not get added
		groups.add(new GFFMisc(new LinePushBackReader(new StringReader(LineGroupCache.FORWARD_REF_DELIM))));
		
		List<LineGroup> lineGrpList = groups.getAllLineGroups();
		
		// Then
		
		// Misc list should be empty
		assertEquals( NUM_EMBL_HDRS_FOR_TESTS + 2, lineGrpList.size() ); // +3 for feature table and sequence
	}
	
	/**
	 * Test the testRemoveFeatureTable() method.
	 * @throws Exception
	 */
	@Test
	public void testRemoveFeatureTable() throws Exception
	{					
		// When
		
		LineGroup removedFT = groups.removeFeatureTable(ft);
		
		// Then
		
		assertNotNull( removedFT );
		assertEquals( ft, removedFT );
		assertNull( groups.getFeatureTable() );
	}
	
	/**
	 * Test the testRemoveSequence() method.
	 * @throws Exception
	 */
	@Test
	public void testRemoveSequence() throws Exception
	{
		// When
		
		LineGroup removedSeq = groups.removeSequence(sequence);
		
		// Then
		
		assertNotNull( removedSeq );
		assertEquals( sequence, removedSeq );
		assertNull( groups.getSequence() );
	}
	
	/**
	 * Test the getHeadersAsText() method.
	 */
	@Test
	public void testGetHeadersAsText()
	{
		// When
		
		String result = groups.getHeadersAsText();
		
		// Then
		
		assertEquals(EMBL_MISC1 + EMBL_MISC2 + EMBL_HDR + EMBL_HDR + EMBL_HDR + EMBL_HDR, result);
		
		// When
		
		groups.clear();
		result = groups.getHeadersAsText();
		
		// Then
		
		assertNull(result);
	}
	
	/**
	 * Test the getAllLineGroups() method.
	 */
	@Test
	public void testGetAllLineGroups()
	{
		
		// When
		
		List<LineGroup> lineGrpList = groups.getAllLineGroups();
		
		// Then
		
		assertEquals( NUM_EMBL_MISC_FEATURES_FOR_TESTS + NUM_EMBL_HDRS_FOR_TESTS + 2, lineGrpList.size() ); // +2 for feature table and sequence
		
		// When
		
		groups.clear();
		lineGrpList = groups.getAllLineGroups();
				
		// Then
				
		assertEquals( 0, lineGrpList.size() );
	}
	
	/**
	 * Test the clear() method.
	 */
	@Test
	public void testClear()
	{
		// When
		
		groups.clear();
		List<LineGroup> lineGrpList = groups.getAllLineGroups();
		
		// Then
		
		assertNotNull( lineGrpList );
		assertEquals( 0, lineGrpList.size() );	
		
		// Check multiple calls
		groups.clear();
		groups.clear();
		
		assertNull(groups.getFeatureTable());
		assertNull(groups.getSequence());
	}
	
	/**
	 * Test the clearMiscLineGroups() method.
	 */
	@Test
	public void testClearMiscLineGroups()
	{
		// When
		
		groups.clearMiscLineGroups();
		List<LineGroup> lineGrpList = groups.getAllLineGroups();
		
		// Then
		
		assertEquals( NUM_EMBL_HDRS_FOR_TESTS + 2, lineGrpList.size() ); // +2 for feature table and sequence
		
		// Check multiple calls
		groups.clearMiscLineGroups();
		groups.clearMiscLineGroups();
	}
	
	
	/**
	 * Test the writeToStream() method.
	 */
	@Test
	public void testWriteToStream() throws Exception
	{
		// =============================================
		// Test with headers, feature table and sequence
		// =============================================
		
		// Given 
		
		StringWriter writer = new StringWriter();
		
		// When
		
		groups.writeToStream(entry, writer);
		String result = writer.toString();
		
		// Then
		
		assertNotNull(result);
		
		assertEquals(
				EMBL_MISC1 + 
				EMBL_MISC2 + 
				EMBL_HDR + 
				EMBL_HDR + 
				EMBL_HDR + 
				EMBL_HDR + 
				EMBL_FEATURE +
				EMBL_SEQUENCE +
				"//\n", 
				result);
		
		// =============================================
		// Test with sequence only
		// =============================================
		
		// Given 
		
		writer = new StringWriter();
		
		// When
		
		groups.clear();
		groups.add(sequence);
		groups.writeToStream(entry, writer);
		result = writer.toString();
		
		// Then
		
		assertNotNull(result);
		
		assertEquals(
				EMBL_SEQUENCE, 
				result);
		
		// =============================================
		// Test with feature table only
		// =============================================
		
		// Given 
		
		writer = new StringWriter();
		ft.add(feature);
		
		// When
		
		groups.clear();
		groups.add(ft);
		groups.writeToStream(entry, writer);
		result = writer.toString();

		// Then
		
		assertNotNull(result);
		
		assertEquals(
				EMBL_FEATURE, 
				result);
	}
	
}
