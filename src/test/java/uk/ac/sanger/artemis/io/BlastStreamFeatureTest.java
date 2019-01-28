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

import java.io.IOException;
import java.io.StringReader;

import org.junit.Test;

import uk.ac.sanger.artemis.util.LinePushBackReader;
import uk.ac.sanger.artemis.io.ReadFormatException;

/**
 * JUnit test class for the BlastStreamFeature class.
 * Expected Blast header format is:<br/>
 * @code query_id subject_id pct_identity aln_length n_of_mismatches gap_openings q_start q_end s_start s_end e_value bit_score
 * 
 * @author kp11
 *
 */
public class BlastStreamFeatureTest
{
	// Qualifier name constants
	private static final String BLAST_SCORE_QUALIFIER_NAME = "blast_score";
	private static final String SCORE_QUALIFIER_NAME = "score";
	private static final String PERCENT_ID_QUALIFIER_NAME = "percent_id";
	private static final String QUERY_ID_QUALIFIER_NAME = "query_id";
	private static final String SUBJECT_START_QUALIFIER_NAME = "subject_start";
	private static final String SUBJECT_END_QUALIFIER_NAME = "subject_end";
	private static final String SUBJECT_ID_QUALIFIER_NAME = "subject_id";
	private static final String NOTE_QUALIFIER_NAME = "note";
	
	/**
	 * Helper method to create a BlastStreamFeature from a Blast file line
	 * @param blastFileLine String
	 * @return BlastStreamFeature
	 * @throws Exception 
	 */
	private BlastStreamFeature getBlastStreamFeature(String blastFileLine) throws Exception
	{
		BlastStreamFeature blastFeature = null;
		
		StringReader reader = new StringReader(blastFileLine);
		
		LinePushBackReader stream = new LinePushBackReader(reader);
		
		try
		{
			blastFeature = BlastStreamFeature.readFromStream(stream);
		} 
		catch (InvalidRelationException | IOException e)
		{
			throw e;
		}
		finally
		{
			try
			{
				stream.close();
				reader.close();
			}
			catch (Exception e)
			{
				fail("Unable to close streams");
			}
		}
		
		return blastFeature;
	}
	
	@Test
	public void testBlastStreamFeature1() throws Exception
	{
		BlastStreamFeature blastFeature = 
				getBlastStreamFeature("H9U8IX202BYF90	233	86.667	15	2	0	127	141	51	65	1.0	17.7");
		
		Location loc1 = blastFeature.getLocation();
		assertEquals(51, loc1.getFirstBase());
		assertEquals(65, loc1.getLastBase());
		assertEquals(51, loc1.getRanges().get(0).getStart());
		assertEquals(65, loc1.getRanges().get(0).getEnd());
		assertFalse(loc1.isComplement());
		
		assertEquals(BlastStreamFeature.BLAST_FEATURE_NAME, blastFeature.getKey().getKeyString());
		
		assertEquals(1, blastFeature.getQualifierByName(BLAST_SCORE_QUALIFIER_NAME).getValues().size());
		assertEquals("17.7", blastFeature.getQualifierByName(BLAST_SCORE_QUALIFIER_NAME).getValues().get(0));
		
		assertEquals(1, blastFeature.getQualifierByName(SCORE_QUALIFIER_NAME).getValues().size());
		assertEquals("86.667", blastFeature.getQualifierByName(SCORE_QUALIFIER_NAME).getValues().get(0));
		
		assertEquals(1, blastFeature.getQualifierByName(PERCENT_ID_QUALIFIER_NAME).getValues().size());
		assertEquals("86.667", blastFeature.getQualifierByName(PERCENT_ID_QUALIFIER_NAME).getValues().get(0));
		
		assertEquals(1, blastFeature.getQualifierByName(QUERY_ID_QUALIFIER_NAME).getValues().size());
		assertEquals("H9U8IX202BYF90", blastFeature.getQualifierByName(QUERY_ID_QUALIFIER_NAME).getValues().get(0));
		
		assertEquals(1, blastFeature.getQualifierByName(SUBJECT_START_QUALIFIER_NAME).getValues().size());
		assertEquals("51", blastFeature.getQualifierByName(SUBJECT_START_QUALIFIER_NAME).getValues().get(0));
		
		assertEquals(1, blastFeature.getQualifierByName(SUBJECT_END_QUALIFIER_NAME).getValues().size());
		assertEquals("65", blastFeature.getQualifierByName(SUBJECT_END_QUALIFIER_NAME).getValues().get(0));
		
		assertEquals(1, blastFeature.getQualifierByName(SUBJECT_ID_QUALIFIER_NAME).getValues().size());
		assertEquals("233", blastFeature.getQualifierByName(SUBJECT_ID_QUALIFIER_NAME).getValues().get(0));
		
		assertEquals(1, blastFeature.getQualifierByName(NOTE_QUALIFIER_NAME).getValues().size());
		assertEquals("hit to 233 51..65  score: 17.7  percent id: 86.667  e-value: 1.0", blastFeature.getQualifierByName(NOTE_QUALIFIER_NAME).getValues().get(0));
		
		
	}
	
	@Test
	public void testBlastStreamFeature2() throws Exception
	{
		BlastStreamFeature blastFeature = 
				getBlastStreamFeature("H9U8IX202BYF90	233	100.000	8	0	0	21	28	138	131	3.8	15.9");
		
		Location loc1 = blastFeature.getLocation();
		assertEquals(131, loc1.getFirstBase());
		assertEquals(138, loc1.getLastBase());
		assertEquals(131, loc1.getRanges().get(0).getStart());
		assertEquals(138, loc1.getRanges().get(0).getEnd());
		assertTrue(loc1.isComplement());
		
		assertEquals(BlastStreamFeature.BLAST_FEATURE_NAME, blastFeature.getKey().getKeyString());
		
		assertEquals("15.9", blastFeature.getQualifierByName(BLAST_SCORE_QUALIFIER_NAME).getValues().get(0));
		assertEquals("100.000", blastFeature.getQualifierByName(SCORE_QUALIFIER_NAME).getValues().get(0));
		assertEquals("100.000", blastFeature.getQualifierByName(PERCENT_ID_QUALIFIER_NAME).getValues().get(0));
		assertEquals("H9U8IX202BYF90", blastFeature.getQualifierByName(QUERY_ID_QUALIFIER_NAME).getValues().get(0));
		assertEquals("138", blastFeature.getQualifierByName(SUBJECT_START_QUALIFIER_NAME).getValues().get(0));
		assertEquals("131", blastFeature.getQualifierByName(SUBJECT_END_QUALIFIER_NAME).getValues().get(0));
		assertEquals("233", blastFeature.getQualifierByName(SUBJECT_ID_QUALIFIER_NAME).getValues().get(0));
	}
	
	@Test
	public void testBlastStreamFeature3() throws Exception
	{
		BlastStreamFeature blastFeature = 
				getBlastStreamFeature("H9U8IX202BYF90	233	100.000	8	0	0	155	162	181	174	3.8	15.9");
		
		Location loc1 = blastFeature.getLocation();
		assertEquals(174, loc1.getFirstBase());
		assertEquals(181, loc1.getLastBase());
		assertEquals(174, loc1.getRanges().get(0).getStart());
		assertEquals(181, loc1.getRanges().get(0).getEnd());
		assertTrue(loc1.isComplement());
		
		assertEquals(BlastStreamFeature.BLAST_FEATURE_NAME, blastFeature.getKey().getKeyString());
		
		assertEquals("15.9", blastFeature.getQualifierByName(BLAST_SCORE_QUALIFIER_NAME).getValues().get(0));
		assertEquals("100.000", blastFeature.getQualifierByName(SCORE_QUALIFIER_NAME).getValues().get(0));
		assertEquals("100.000", blastFeature.getQualifierByName(PERCENT_ID_QUALIFIER_NAME).getValues().get(0));
		assertEquals("H9U8IX202BYF90", blastFeature.getQualifierByName(QUERY_ID_QUALIFIER_NAME).getValues().get(0));
		assertEquals("181", blastFeature.getQualifierByName(SUBJECT_START_QUALIFIER_NAME).getValues().get(0));
		assertEquals("174", blastFeature.getQualifierByName(SUBJECT_END_QUALIFIER_NAME).getValues().get(0));
		assertEquals("233", blastFeature.getQualifierByName(SUBJECT_ID_QUALIFIER_NAME).getValues().get(0));
	}
	
	@Test
	public void testBlastStreamFeature4() throws Exception
	{
		BlastStreamFeature blastFeature = 
				getBlastStreamFeature("Z	819	0.0	8	0	0	198	191	34	41	3.8	0.0");
		
		Location loc1 = blastFeature.getLocation();
		assertEquals(34, loc1.getFirstBase());
		assertEquals(41, loc1.getLastBase());
		assertEquals(34, loc1.getRanges().get(0).getStart());
		assertEquals(41, loc1.getRanges().get(0).getEnd());
		assertTrue(loc1.isComplement());
		
		assertEquals(BlastStreamFeature.BLAST_FEATURE_NAME, blastFeature.getKey().getKeyString());
		
		assertEquals("0.0", blastFeature.getQualifierByName(BLAST_SCORE_QUALIFIER_NAME).getValues().get(0));
		assertEquals("0.0", blastFeature.getQualifierByName(SCORE_QUALIFIER_NAME).getValues().get(0));
		assertEquals("0.0", blastFeature.getQualifierByName(PERCENT_ID_QUALIFIER_NAME).getValues().get(0));
		assertEquals("Z", blastFeature.getQualifierByName(QUERY_ID_QUALIFIER_NAME).getValues().get(0));
		assertEquals("34", blastFeature.getQualifierByName(SUBJECT_START_QUALIFIER_NAME).getValues().get(0));
		assertEquals("41", blastFeature.getQualifierByName(SUBJECT_END_QUALIFIER_NAME).getValues().get(0));
		assertEquals("819", blastFeature.getQualifierByName(SUBJECT_ID_QUALIFIER_NAME).getValues().get(0));
	}
	
	@Test
	public void testInvalidBlastLine() throws Exception
	{
		@SuppressWarnings("unused")
		BlastStreamFeature blastFeature = null;
		
		try 
		{
			blastFeature = 
					getBlastStreamFeature("Z	819	0.0	8	0	0	198	191	34	41	3.8");
			fail("Expected a ReadFormatException to be thrown");
		}
		catch (ReadFormatException e)
		{
			// Expected
		}
		
		blastFeature = getBlastStreamFeature("");
		assertNull("Expected a ReadFormatException to be thrown", blastFeature);
	}
}
