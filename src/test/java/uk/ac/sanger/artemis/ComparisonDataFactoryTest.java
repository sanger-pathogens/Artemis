package uk.ac.sanger.artemis;

import static org.junit.Assert.*;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.util.List;

import org.junit.Before;
import org.junit.Test;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;

import static org.mockito.Mockito.*;

import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.LinePushBackReader;

/**
 * JUnit test for the ComparisonDataFactory
 * class and associated functionality.
 * 
 * @author kp11
 *
 */
public class ComparisonDataFactoryTest
{
	@Mock
	Document comparisonDoc;
	
	@Before
	public void setUp() throws Exception {

		MockitoAnnotations.initMocks(this);
	}

	/**
	 * Test the peekFirstLine method.
	 * @throws Exception
	 */
	@Test
	public void testPeekFirstLine() throws Exception
	{
		// Test first line retrieval works
		
		String testText = "Line-1\nLine-2\nLine3";
		StringReader reader = new StringReader(testText);
		LinePushBackReader pushbackReader = new LinePushBackReader(reader);
		
		String line = ComparisonDataFactory.peekFirstLine(pushbackReader, "test-filename.txt");
		
		assertEquals("Line-1", line);
		// Check first line is still available and has not been consumed
		assertEquals("Line-1", pushbackReader.readLine());
		
		// Test end of file handling
		
		testText = "";
		line = null;
		
		reader = new StringReader(testText);
		
		try 
		{
			line = ComparisonDataFactory.peekFirstLine(new LinePushBackReader(reader), "test-filename.txt");
			fail("Expected an IOException for end of file");
		} 
		catch (IOException e) 
		{
			// Expected
			assertNotNull(e.getMessage());
			assertTrue(e.getMessage().contains("test-filename.txt"));
		}	
	}
	
	/**
	 * Test the readHeaders method.
	 * @throws Exception
	 */
	@Test
	public void testReadHeaders() throws Exception
	{
		// Test reading of comment headers
		
		String testText = "#Comment-1\n# Comment-2\n#  Comment 3\nData Line\n";
		StringReader reader = new StringReader(testText);
		LinePushBackReader pushbackReader = new LinePushBackReader(reader);
		
		List<String> headers = ComparisonDataFactory.readHeaders(pushbackReader, "test-filename.txt");
		
		assertNotNull(headers);
		assertEquals(3, headers.size());
		assertEquals("#Comment-1", headers.get(0));
		assertEquals("# Comment-2", headers.get(1));
		assertEquals("#  Comment 3", headers.get(2));
		
		assertEquals("Data Line", pushbackReader.readLine());
		
		// Test end of file handling
		
		testText = "";
		headers = null;
		
		reader = new StringReader(testText);
		
		try 
		{
			headers = ComparisonDataFactory.readHeaders(pushbackReader, "test-filename.txt");
			fail("Expected an IOException for end of file");
		} 
		catch (IOException e) 
		{
			// Expected
			assertNotNull(e.getMessage());
			assertTrue(e.getMessage().contains("test-filename.txt"));
		}	
	}
	
	/**
	 * Test loading of Blast web site generated comparison data.
	 * @throws Exception
	 */
	@Test
	public void testReadComparisonDataForBlastWebSiteHitTable() throws Exception
	{
		// TBlastx...
		
		// Given
		InputStream testFile = ComparisonDataFactoryTest.class.getResourceAsStream(
				"/data/act-comparison-files/web-tblastx-comparison-file.txt");
		
		// When
		when( comparisonDoc.getReader() ).thenReturn(new InputStreamReader(testFile));
		
		ComparisonData result = ComparisonDataFactory.readComparisonData (comparisonDoc);
		
		// Then
		assertTrue(result instanceof BlastWebSiteHitTableComparisonData);
		
		// Check data...
		
		AlignMatch matches [] = result.getMatches();
		assertEquals(62, matches.length);
		
		AlignMatch firstDataLine = matches[0];
		
		assertEquals(3678, firstDataLine.getScore());
		assertEquals(100,   firstDataLine.getPercentID());
		assertEquals(2458, firstDataLine.getQuerySequenceStart());
		assertEquals(7239, firstDataLine.getQuerySequenceEnd());
		assertEquals(3358, firstDataLine.getSubjectSequenceStart());
		assertEquals(8139, firstDataLine.getSubjectSequenceEnd());
		
		AlignMatch lastDataLine = matches[61];
		
		assertEquals(20,   lastDataLine.getScore());
		assertEquals(31,   lastDataLine.getPercentID());
		assertEquals(1413,  lastDataLine.getQuerySequenceStart());
		assertEquals(1478,  lastDataLine.getQuerySequenceEnd());
		assertEquals(6687,  lastDataLine.getSubjectSequenceStart());
		assertEquals(6752,  lastDataLine.getSubjectSequenceEnd());
		
		testFile.close();
		
		
		// Blastn...
		
		// Given
		testFile = ComparisonDataFactoryTest.class.getResourceAsStream(
				"/data/act-comparison-files/web-blastn-comparison-file.txt");
		
		// When
		when( comparisonDoc.getReader() ).thenReturn(new InputStreamReader(testFile));
		
		result = ComparisonDataFactory.readComparisonData (comparisonDoc);
		
		// Then
		assertTrue(result instanceof BlastWebSiteHitTableComparisonData);
		
		// Check data...
		
		matches = result.getMatches();
		assertEquals(3, matches.length);
		
		firstDataLine = matches[0];
		
		assertEquals(10072, firstDataLine.getScore());
		assertEquals(100,  firstDataLine.getPercentID());
		assertEquals(2101, firstDataLine.getQuerySequenceStart());
		assertEquals(7554, firstDataLine.getQuerySequenceEnd());
		assertEquals(3001, firstDataLine.getSubjectSequenceStart());
		assertEquals(8454, firstDataLine.getSubjectSequenceEnd());
		
		lastDataLine = matches[2];
		
		assertEquals(1773,   lastDataLine.getScore());
		assertEquals(95,     lastDataLine.getPercentID());
		assertEquals(1141,   lastDataLine.getQuerySequenceStart());
		assertEquals(2100,   lastDataLine.getQuerySequenceEnd());
		assertEquals(1861,   lastDataLine.getSubjectSequenceStart());
		assertEquals(2820,   lastDataLine.getSubjectSequenceEnd());
		
		testFile.close();
		
		
		// Empty Blastn file...
		
		// Given
		testFile = new ByteArrayInputStream(new byte[0]);
		
		// When
		when( comparisonDoc.getReader() ).thenReturn(new InputStreamReader(testFile));
		
		try 
		{
			result = ComparisonDataFactory.readComparisonData (comparisonDoc);
			fail("Expected IO Excption for an end of file");
		}
		catch (IOException e) 
		{
			// Expected
			assertTrue(e.getMessage().contains("End of file"));
		}
		
		testFile.close();
	}
	
	/**
	 * Test loading of SSAHA comparison data.
	 * @throws Exception
	 */
	@Test
	public void testReadComparisonDataForSSAHA() throws Exception
	{
		// Check identification of file type...
		
		// Given
		StringBuffer dataBuf = new StringBuffer();
		dataBuf.append("F\tQueryName1\t10\t20\tSubjectName1\t30\t40\t35\t95\n");
		dataBuf.append("F\tQueryName1\t11\t21\tSubjectName1\t31\t41\t36\t96\n");
		dataBuf.append("F\tQueryName1\t12\t22\tSubjectName1\t32\t42\t37\t97\n");
		
		StringReader reader = new StringReader(dataBuf.toString());
		
		// When
		when( comparisonDoc.getReader() ).thenReturn(reader);
		
		ComparisonData result = ComparisonDataFactory.readComparisonData (comparisonDoc);
		
		// Then
		assertTrue(result instanceof SSAHAComparisonData);
		

		// Check data...
		
		AlignMatch matches [] = result.getMatches();
		assertEquals(3, matches.length);
		
		AlignMatch firstDataLine = matches[0];
		
		assertEquals(35,   firstDataLine.getScore());
		assertEquals(95,   firstDataLine.getPercentID());
		assertEquals(10,   firstDataLine.getQuerySequenceStart());
		assertEquals(20,   firstDataLine.getQuerySequenceEnd());
		assertEquals(30,   firstDataLine.getSubjectSequenceStart());
		assertEquals(40,   firstDataLine.getSubjectSequenceEnd());
		
		reader.close();
		
	}
	
	/**
	 * Test loading of crunch comparison data.
	 * @throws Exception
	 */
	@Test
	public void testReadComparisonDataForCrunch() throws Exception
	{
		// Check identification of file type...
		
		// Given
		InputStream testFile = ComparisonDataFactoryTest.class.getResourceAsStream(
				"/data/act-comparison-files/crunch-comparison-file.crunch");
		
		// When
		when( comparisonDoc.getReader() ).thenReturn(new InputStreamReader(testFile));
		
		ComparisonData result = ComparisonDataFactory.readComparisonData (comparisonDoc);
		
		// Then
		assertTrue(result instanceof MSPcrunchComparisonData);
		

		// Check data...
		
		AlignMatch matches [] = result.getMatches();
		assertEquals(126, matches.length);
		
		AlignMatch firstDataLine = matches[0];
		
		assertEquals(9583, firstDataLine.getScore());
		assertEquals(99,   firstDataLine.getPercentID());
		assertEquals(1474, firstDataLine.getQuerySequenceStart());
		assertEquals(6490, firstDataLine.getQuerySequenceEnd());
		assertEquals(1474, firstDataLine.getSubjectSequenceStart());
		assertEquals(6490, firstDataLine.getSubjectSequenceEnd());
		
		AlignMatch lastDataLine = matches[125];
		
		assertEquals(42,   lastDataLine.getScore());
		assertEquals(87,   lastDataLine.getPercentID());
		assertEquals(521,  lastDataLine.getQuerySequenceStart());
		assertEquals(565,  lastDataLine.getQuerySequenceEnd());
		assertEquals(66,   lastDataLine.getSubjectSequenceStart());
		assertEquals(110,  lastDataLine.getSubjectSequenceEnd());
		
		testFile.close();
		
	}
	
	/**
	 * Test loading of Blast m8 formatted comparison data (with no headers).
	 * @throws Exception
	 */
	@Test
	public void testReadComparisonDataForBlastM8NoHdrs() throws Exception
	{
		// Given
		InputStream testFile = ComparisonDataFactoryTest.class.getResourceAsStream(
				"/data/act-comparison-files/blast-m8-comparison-file.txt");
		
		// When
		when( comparisonDoc.getReader() ).thenReturn(new InputStreamReader(testFile));
		
		ComparisonData result = ComparisonDataFactory.readComparisonData (comparisonDoc);
		
		// Then
		assertTrue(result instanceof BlastM8ComparisonData);
		
		testFile.close();
	}
	
	/**
	 * Test loading of Blast m8 formatted comparison data (with header).
	 * @throws Exception
	 */
	@Test
	public void testReadComparisonDataForBlastM8WithHdrs() throws Exception
	{
		// Given
		InputStream testFile = ComparisonDataFactoryTest.class.getResourceAsStream(
				"/data/act-comparison-files/blast-m8-comparison-file-withhdr.txt");
		
		// When
		when( comparisonDoc.getReader() ).thenReturn(new InputStreamReader(testFile));
		
		ComparisonData result = ComparisonDataFactory.readComparisonData (comparisonDoc);
		
		// Then
		assertTrue(result instanceof BlastM8ComparisonData);
		
		
		// Check data...
		
		AlignMatch matches [] = result.getMatches();
		assertEquals(3, matches.length);
		
		AlignMatch firstDataLine = matches[0];
		
		assertEquals(5454, firstDataLine.getScore());
		assertEquals(100,  firstDataLine.getPercentID());
		assertEquals(2101, firstDataLine.getQuerySequenceStart());
		assertEquals(7554, firstDataLine.getQuerySequenceEnd());
		assertEquals(3001, firstDataLine.getSubjectSequenceStart());
		assertEquals(8454, firstDataLine.getSubjectSequenceEnd());
		
		AlignMatch lastDataLine = matches[2];
		
		assertEquals(960,   lastDataLine.getScore());
		assertEquals(100,   lastDataLine.getPercentID());
		assertEquals(1141,  lastDataLine.getQuerySequenceStart());
		assertEquals(2100,  lastDataLine.getQuerySequenceEnd());
		assertEquals(1861,  lastDataLine.getSubjectSequenceStart());
		assertEquals(2820,  lastDataLine.getSubjectSequenceEnd());
		
		testFile.close();
	}
	
	/**
	 * Test response to unrecognised file format.
	 * @throws Exception
	 */
	@Test
	public void testUnknownFileFormat() throws Exception
	{
		// Given
		String testText = "10\t23\n25\t12\n";
		StringReader reader = new StringReader(testText);
		
		// When
		when( comparisonDoc.getReader() ).thenReturn(reader);
		
		// Then
		try
		{
			@SuppressWarnings("unused")
			ComparisonData result = ComparisonDataFactory.readComparisonData (comparisonDoc);
			fail("Expected exception for an unknown file format");
		} 
		catch (Exception e)
		{
			// Expected
			assertEquals("cannot understand the comparison file format", e.getMessage());
		}
		
		reader.close();
	}
	
	/**
	 * Test response to an incorrect no. of fields in a Blast web site generated file.
	 * @throws Exception
	 */
	@Test
	public void testIncorrectFieldNumForBlastWebSiteComparisonData1() throws Exception
	{
		// Given
		StringBuffer dataBuf1 = new StringBuffer();
		
		dataBuf1.append("# blastn\n");
		dataBuf1.append("# Iteration: 0\n");
		dataBuf1.append("# Query: NC_001954_1_1\n");
		dataBuf1.append("# RID: YUYABC86114\n");
		dataBuf1.append("# Database: n/a\n");
		dataBuf1.append("# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score\n");
		dataBuf1.append("# 3 hits found\n");
		
		//  Not enough fields
		dataBuf1.append("NC_001954_1_1\tNC_001954_1_1\t100.000\t5454\t0\t0\t2101\t7554\t3001\t8454\t0.0\n");
				
		StringReader reader = new StringReader(dataBuf1.toString());
		
		// When
		when( comparisonDoc.getReader() ).thenReturn(reader);
		
		// Then
		try
		{
			@SuppressWarnings("unused")
			ComparisonData result = ComparisonDataFactory.readComparisonData (comparisonDoc);
			fail("Expected exception for incorrect number of fields");
		} 
		catch (Exception e)
		{
			// Expected
			assertTrue(e.getMessage().contains("unexpected number of fields for this line"));
		}
		
		reader.close();
	}
	
	/**
	 * Test response to an incorrect no. of fields in a Blast web site generated file.
	 * @throws Exception
	 */
	@Test
	public void testIncorrectFieldNumForBlastWebSiteComparisonData2() throws Exception
	{
		// Given
		StringBuffer dataBuf1 = new StringBuffer();
		
		dataBuf1.append("# blastn\n");
		dataBuf1.append("# Iteration: 0\n");
		dataBuf1.append("# Query: NC_001954_1_1\n");
		dataBuf1.append("# RID: YUYABC86114\n");
		dataBuf1.append("# Database: n/a\n");
		dataBuf1.append("# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score\n");
		dataBuf1.append("# 3 hits found\n");
		
		//  Too many fields
		dataBuf1.append("NC_001954_1_1\tNC_001954_1_1\t100.000\t5454\t0\t0\t2101\t7554\t3001\t8454\t0.0\34\23\26\132\n");
		
		StringReader reader = new StringReader(dataBuf1.toString());
		
		// When
		when( comparisonDoc.getReader() ).thenReturn(reader);
		
		// Then
		try
		{
			@SuppressWarnings("unused")
			ComparisonData result = ComparisonDataFactory.readComparisonData (comparisonDoc);
			fail("Expected exception for incorrect number of fields");
		} 
		catch (Exception e)
		{
			// Expected
			assertTrue(e.getMessage().contains("unexpected number of fields for this line"));
		}
		
		reader.close();
	}

}
