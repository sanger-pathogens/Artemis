package uk.ac.sanger.artemis.sequence;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;

import static org.mockito.Mockito.*;

import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.io.EmblStreamSequence;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.LocationParseException;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.Strand;
import uk.ac.sanger.artemis.util.OutOfRangeException;

/**
 * Basic unit tests for the Strand class.
 * 
 * @author kp11
 *
 */
public class StrandTest
{
	private static final int FEATURE_START_POS = 10;
	private static final int FEATURE_END_POS   = 150;
	
	// NUM_SEQ_BASES bases...
	private static final String SEQUENCE = "ggcgaggcggggaaagcactgcgcgctgacggtggtgctgattgtattttttcagcgtctcagcgcgtcgtgacggcacttagtctgcccgttgaggcgttgtgtgtctgcggggtgttttgtgcggtggtgagcgtgtgaggggggatgacggggtgtaaaaaagccgcccgcaggcggcgatgttcagtcgttgtcagtgtccagtgagtagtttttaaagcggatgacctcctgaccgagccagccgtttatctcgcggatcctgtcctgtaacgggataagctcattgcggacaaagacctttgccactttctcaatatcacccagcgacccgacgttctccggcttgccacccatcaactgaaaggggatgcggtgcgcgtccagcaggtcagcggcgctggcttttttgatattaaaaaaatcgtccttcgtcgccacttcactgagggggataattttaatgccgtcggctttcccctgtggggcatagagaaacaggtttttaaagttgttgcggcctttcgacttgaccatgttttcgcgaagcatttcgatatcgttgcgatcctgcacggcatcggtgacatacatgatgtatccggcatgtgcgccattttcgtaatacttgcggcggaacaacgtggccgactcattcagccaggcagagttaagggcgctgagatattccggcaggccgtacagctcctgattaatatccggctccagcaggtgaaacacggagccgggcgcgaaggctgtcggctcgttgaaggacggcacccaccagtaaacatcctcttccacgccacggcgggtatattttgccggtgaggtttccagtctgatgaccttaccggtggtgctgtaacgcttttccagaaacgcattaccgaacaccagaaaatccagcacaaagcggctgaaatcctgctgggaaagccatggatgcgggataaatgtcgaggccagaatattgcgtttgacgtaaatcggcgagctgtgatgcacggcagcccgcaggctttttgccagaccggtaaagctgaccggtggctcataccatctgccgttactgatgcactcgacgtaatccagaatgtcacggcggtcgagtaccggcaccggctcaccaaaggtgaatgcctccattttcgggccgctggcggtcattgtttttgccgcaggttgcggtgttttcccttttttcttgctcatcagtaaaactccagaatggtggatgtcagcggggtgctgataccggcggtgagtggctcatttaacagggcgtgcatggtcgcccaggcgaggtcggcgtggctggcttcctcgctgcggctggcctcataggtggcgctgcgtccgctgctggtcatggtcttgcggatagccataaacgagctggtgatgtcggtggcgctgacgtcgtattccagacagccacggcggataacgtcttttgccttgagcaccattgcggttttcatttccggcgtgtagcggatatcacgcgcggcgggatagaacgagcgcacgagctggaacacgccgacaccgaggccggtggcatcaataccgatgtattcgacgttgtatttttcggtgagtttgcggatggattccgcctgggtggcaaagtccatgcctttccactggtgacgctcaagtattctgaatttgccaccggccaccaccggcggtgccagtaccacgcatccggcgctgtcgccacggtgtgacgggtcgtaaccaatccataccgggcgggagccgaacggattggcggcaaacggcgcatagtcttcccattcttccagcgtgtcgaccatgcagcgttgcagctcctcgaacgggaacaccgatgccttgtcgtcaacaaattcacacatgaacaggtttttaaaatcgtcggcgctgttttcgcgt";
	
	@Mock
	private Feature feature;
	
	private Bases bases;
	private Strand fwdStrand;
	private Strand revStrand;
	
	
	
	@Before
	public void setupForTest() throws Exception
	{
		MockitoAnnotations.initMocks(this);
		
		bases = new Bases( new EmblStreamSequence(SEQUENCE) );
		fwdStrand = bases.getForwardStrand();
		revStrand = bases.getReverseStrand();
	}
	
	@After
	public void cleanupAfterTest() throws Exception
	{
		bases = null;
		fwdStrand = null;
		revStrand = null;
	}
	
	@Test
	public void testGetStrandBases() throws OutOfRangeException, LocationParseException
	{
		// When

		when(feature.getStrand()).thenReturn(fwdStrand);
		when(feature.getLocation()).thenReturn(new Location(new Range(FEATURE_START_POS, FEATURE_END_POS)));
		when(feature.isForwardFeature()).thenReturn(true);
		
		String bases = fwdStrand.getStrandBases();
		String revBases = revStrand.getStrandBases();
		
		// Then
		
		assertEquals(SEQUENCE, bases);
		assertEquals(Bases.reverseComplement(SEQUENCE), revBases);
	}
	
	@Test
	public void testGetDirection()
	{
		
		assertEquals(Strand.FORWARD, fwdStrand.getDirection());
		assertEquals(Strand.REVERSE, revStrand.getDirection());
		
	}
	
	@Test
	public void testGetSequenceLength()
	{
		
		assertEquals(SEQUENCE.length(), fwdStrand.getSequenceLength());
		assertEquals(SEQUENCE.length(), revStrand.getSequenceLength());
		
	}
	
	@Test
	public void testMakeMarker() throws OutOfRangeException
	{
		// When
		
		Marker marker1 = fwdStrand.makeMarker(60);
		Marker marker2 = revStrand.makeMarker(90);
		
		// Then
		
		assertNotNull(marker1);
		assertEquals(60, marker1.getRawPosition());
		assertNotNull(marker2);
		assertEquals(SEQUENCE.length()-90+1, marker2.getRawPosition());
		
	}
	
	@Test
	public void testMakeMarkerFromRawPosition() throws OutOfRangeException
	{
		// When
		
		Marker marker1 = fwdStrand.makeMarkerFromRawPosition(60);
		Marker marker2 = revStrand.makeMarkerFromRawPosition(90);
		
		// Then
		
		assertNotNull(marker1);
		assertEquals(60, marker1.getRawPosition());
		assertNotNull(marker2);
		assertEquals(90, marker2.getRawPosition());
		
	}
	
	@Test
	public void testMakeMarkerRangeFromRawPositions() throws OutOfRangeException
	{
		// When
		
		MarkerRange marker1 = fwdStrand.makeMarkerRangeFromRawPositions(10,30);
		MarkerRange marker2 = revStrand.makeMarkerRangeFromRawPositions(200,300);
		
		// Then
		
		assertNotNull(marker1);
		assertEquals(10, marker1.getRawStart().getPosition());
		assertEquals(30, marker1.getRawEnd().getPosition());
		assertNotNull(marker2);
		assertEquals(SEQUENCE.length()-200+1, marker2.getRawStart().getPosition());
		assertEquals(SEQUENCE.length()-300+1, marker2.getRawEnd().getPosition());
		
	}
	
	@Test
	public void testGetRawPosition()
	{
		// When

		int posFwd = fwdStrand.getRawPosition(50);
		int posRev = revStrand.getRawPosition(60);
		
		// Then
		
		assertEquals(50, posFwd);
		assertEquals(SEQUENCE.length()-60+1, posRev);
	}

	@Test
	public void testGetTranslation() throws OutOfRangeException
	{
		// When

		AminoAcidSequence seqFwd = fwdStrand.getTranslation(
				new Range(10,50),
	            true);
		AminoAcidSequence seqRev = fwdStrand.getTranslation(
				new Range(100,200),
	            true);
	
	
		// Then

		assertEquals("gkalradggadci", seqFwd.toString());
		assertEquals("lcvcgvfcavvsv*ggmtgckkaarrrrcsvvv", seqRev.toString());

	}
	
	@Test
	public void testGetSpacedTranslation() throws OutOfRangeException
	{
		// When

		AminoAcidSequence seqFwd = fwdStrand.getSpacedTranslation(
				new Range(10,50),
	            true);
		AminoAcidSequence seqRev = fwdStrand.getSpacedTranslation(
				new Range(100,200),
	            true);
	
	
		// Then

		assertEquals("g  k  a  l  r  a  d  g  g  a  d  c  i  ", seqFwd.toString());
		assertEquals("l  c  v  c  g  v  f  c  a  v  v  s  v  *  g  g  m  t  g  c  k  k  a  a  r  r  r  r  c  s  v  v  v  ", seqRev.toString());

	}
	
	@Test
	public void testGetACount() throws OutOfRangeException
	{
		// When

		int countFwd = fwdStrand.getACount();
		int countRev = revStrand.getACount();
	
		// Then

		assertEquals(383, countFwd);
		assertEquals(465, countRev);
	}
	
	@Test
	public void testGetCCount() throws OutOfRangeException
	{
		// When

		int countFwd = fwdStrand.getCCount();
		int countRev = revStrand.getCCount();
	
	
		// Then

		assertEquals(506, countFwd);
		assertEquals(575, countRev);
	}
	
	@Test
	public void testGetGCount() throws OutOfRangeException
	{
		// When

		int countFwd = fwdStrand.getGCount();
		int countRev = revStrand.getGCount();
	
	
		// Then

		assertEquals(575, countFwd);
		assertEquals(506, countRev);
	}
	
	@Test
	public void testGetTCount() throws OutOfRangeException
	{
		// When

		int countFwd = fwdStrand.getTCount();
		int countRev = revStrand.getTCount();
	
	
		// Then

		assertEquals(465, countFwd);
		assertEquals(383, countRev);
	}

}
