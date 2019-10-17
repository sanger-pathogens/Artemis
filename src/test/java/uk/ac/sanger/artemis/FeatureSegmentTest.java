package uk.ac.sanger.artemis;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;

import static org.mockito.Mockito.*;

import uk.ac.sanger.artemis.io.EmblStreamSequence;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.LocationParseException;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.Strand;
import uk.ac.sanger.artemis.util.OutOfRangeException;

/**
 * Unit tests for the FeatureSegment class.
 * 
 * @author kp11
 *
 */
public class FeatureSegmentTest
{

	private static final int START_POS = 70;
	private static final int END_POS   = 80;
	
	private static final int FEATURE_START_POS = 10;
	private static final int FEATURE_END_POS   = 150;
	
	// NUM_SEQ_BASES bases...
	private static final String SEQUENCE = "ggcgaggcggggaaagcactgcgcgctgacggtggtgctgattgtattttttcagcgtctcagcgcgtcgtgacggcacttagtctgcccgttgaggcgttgtgtgtctgcggggtgttttgtgcggtggtgagcgtgtgaggggggatgacggggtgtaaaaaagccgcccgcaggcggcgatgttcagtcgttgtcagtgtccagtgagtagtttttaaagcggatgacctcctgaccgagccagccgtttatctcgcggatcctgtcctgtaacgggataagctcattgcggacaaagacctttgccactttctcaatatcacccagcgacccgacgttctccggcttgccacccatcaactgaaaggggatgcggtgcgcgtccagcaggtcagcggcgctggcttttttgatattaaaaaaatcgtccttcgtcgccacttcactgagggggataattttaatgccgtcggctttcccctgtggggcatagagaaacaggtttttaaagttgttgcggcctttcgacttgaccatgttttcgcgaagcatttcgatatcgttgcgatcctgcacggcatcggtgacatacatgatgtatccggcatgtgcgccattttcgtaatacttgcggcggaacaacgtggccgactcattcagccaggcagagttaagggcgctgagatattccggcaggccgtacagctcctgattaatatccggctccagcaggtgaaacacggagccgggcgcgaaggctgtcggctcgttgaaggacggcacccaccagtaaacatcctcttccacgccacggcgggtatattttgccggtgaggtttccagtctgatgaccttaccggtggtgctgtaacgcttttccagaaacgcattaccgaacaccagaaaatccagcacaaagcggctgaaatcctgctgggaaagccatggatgcgggataaatgtcgaggccagaatattgcgtttgacgtaaatcggcgagctgtgatgcacggcagcccgcaggctttttgccagaccggtaaagctgaccggtggctcataccatctgccgttactgatgcactcgacgtaatccagaatgtcacggcggtcgagtaccggcaccggctcaccaaaggtgaatgcctccattttcgggccgctggcggtcattgtttttgccgcaggttgcggtgttttcccttttttcttgctcatcagtaaaactccagaatggtggatgtcagcggggtgctgataccggcggtgagtggctcatttaacagggcgtgcatggtcgcccaggcgaggtcggcgtggctggcttcctcgctgcggctggcctcataggtggcgctgcgtccgctgctggtcatggtcttgcggatagccataaacgagctggtgatgtcggtggcgctgacgtcgtattccagacagccacggcggataacgtcttttgccttgagcaccattgcggttttcatttccggcgtgtagcggatatcacgcgcggcgggatagaacgagcgcacgagctggaacacgccgacaccgaggccggtggcatcaataccgatgtattcgacgttgtatttttcggtgagtttgcggatggattccgcctgggtggcaaagtccatgcctttccactggtgacgctcaagtattctgaatttgccaccggccaccaccggcggtgccagtaccacgcatccggcgctgtcgccacggtgtgacgggtcgtaaccaatccataccgggcgggagccgaacggattggcggcaaacggcgcatagtcttcccattcttccagcgtgtcgaccatgcagcgttgcagctcctcgaacgggaacaccgatgccttgtcgtcaacaaattcacacatgaacaggtttttaaaatcgtcggcgctgttttcgcgt";
	
	@Mock
	private Feature feature;
	
	private Range range;
	private Bases bases;
	private Strand fwdStrand;
	
	
	@Before
	public void setupForTest() throws Exception
	{
		MockitoAnnotations.initMocks(this);
		
		range = new Range (START_POS, END_POS);
		bases = new Bases( new EmblStreamSequence(SEQUENCE) );
		fwdStrand = bases.getForwardStrand();
	}
	
	@After
	public void cleanupAfterTest() throws Exception
	{
		range = null;
		bases = null;
		fwdStrand = null;
	}
	
	@Test
	public void testSetStartPosition() throws OutOfRangeException, LocationParseException
	{
		// When

		when(feature.getStrand()).thenReturn(fwdStrand);
		when(feature.getLocation()).thenReturn(new Location(new Range(FEATURE_START_POS, FEATURE_END_POS)));
		when(feature.isForwardFeature()).thenReturn(true);
		
		FeatureSegment segment = new FeatureSegment(
				feature,
                range);
		
		// Sanity check
		assertEquals(START_POS, segment.getStart().getPosition());
		assertEquals(END_POS, segment.getEnd().getPosition());
		
		
		try
		{
			segment.setStartPosition(50);
		}
		catch (OutOfRangeException e)
		{
			fail("setStartPosition threw unexpected OutOfRangeException exeption: " + e.getMessage());
		}
		
		// Then (also intrinsically test updateRange)
		
		assertEquals(50, segment.getStart().getPosition());
		assertEquals(50, segment.getRawRange().getStart());
		assertEquals(END_POS, segment.getRawRange().getEnd());
	}
	
	@Test
	public void testSetEndPosition() throws OutOfRangeException, LocationParseException
	{
		// When

		when(feature.getStrand()).thenReturn(fwdStrand);
		when(feature.getLocation()).thenReturn(new Location(new Range(FEATURE_START_POS, FEATURE_END_POS)));
		when(feature.isForwardFeature()).thenReturn(true);
		
		FeatureSegment segment = new FeatureSegment(
				feature,
                range);
		
		// Sanity check
		assertEquals(START_POS, segment.getStart().getPosition());
		assertEquals(END_POS, segment.getEnd().getPosition());
		
		try
		{
			segment.setEndPosition(100);
		}
		catch (OutOfRangeException e)
		{
			fail("setEndPosition threw unexpected OutOfRangeException exeption: " + e.getMessage());
		}
		
		// Then (also intrinsically test updateRange)
		
		assertEquals(100, segment.getEnd().getPosition());
		assertEquals(START_POS, segment.getRawRange().getStart());
		assertEquals(100, segment.getRawRange().getEnd());
	}
	
	@Test
	public void testGetMarkerRange() throws OutOfRangeException
	{
		// When

		when(feature.getStrand()).thenReturn(fwdStrand);
		when(feature.getLocation()).thenReturn(new Location(new Range(FEATURE_START_POS, FEATURE_END_POS)));
		when(feature.isForwardFeature()).thenReturn(true);
		
		FeatureSegment segment = new FeatureSegment(
				feature,
                range);

		
		// Then 
		
		assertNotNull(segment.getMarkerRange());
		assertEquals(START_POS, segment.getMarkerRange().getStart().getRawPosition());
		assertEquals(END_POS, segment.getMarkerRange().getEnd().getRawPosition());
		assertEquals(START_POS, segment.getMarkerRange().getStart().getPosition());
		assertEquals(END_POS, segment.getMarkerRange().getEnd().getPosition());
	}
	
	@Test
	public void testGetBaseCount() throws OutOfRangeException 
	{
		// When

		when(feature.getStrand()).thenReturn(fwdStrand);
		when(feature.getLocation()).thenReturn(new Location(new Range(FEATURE_START_POS, FEATURE_END_POS)));
		when(feature.isForwardFeature()).thenReturn(true);
		
		FeatureSegment segment = new FeatureSegment(
				feature,
                range);

		
		// Then
		
		assertEquals(END_POS-START_POS+1, segment.getBaseCount());
	}
	
	@Test
	public void testGetBases() throws OutOfRangeException
	{
		// When

		when(feature.getStrand()).thenReturn(fwdStrand);
		when(feature.getLocation()).thenReturn(new Location(new Range(FEATURE_START_POS, FEATURE_END_POS)));
		when(feature.isForwardFeature()).thenReturn(true);
		
		FeatureSegment segment = new FeatureSegment(
				feature,
                range);

		
		// Then

		assertEquals("gtgacggcact", segment.getBases());
		
	}
	
	@Test
	public void testGetFrameID() throws OutOfRangeException
	{
		// Given
		
		FeatureSegmentVector segs = new FeatureSegmentVector();
		
		// When

		when(feature.getStrand()).thenReturn(fwdStrand);
		when(feature.getLocation()).thenReturn(new Location(new Range(FEATURE_START_POS, FEATURE_END_POS)));
		when(feature.isForwardFeature()).thenReturn(true);
		
		FeatureSegment segment = new FeatureSegment(
				feature,
                range);

		segs.add(segment);
		when(feature.getSegments()).thenReturn(segs);
		
		int frameID = segment.getFrameID();
		
		// Then

		assertEquals(2, frameID);
		
	}

}
