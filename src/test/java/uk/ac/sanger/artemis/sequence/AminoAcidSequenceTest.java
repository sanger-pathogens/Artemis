/* 
 * This file is part of Artemis
 *
 * Copyright (C) 2017  Genome Research Limited
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
package uk.ac.sanger.artemis.sequence;

import static org.junit.Assert.*;

import java.net.URL;

import org.junit.Test;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.io.DocumentEntryFactory;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.DocumentFactory;

import org.junit.Before;
import org.junit.After;
import org.junit.Assert;


/**
 * Unit test for the AminoAcidSequence class.
 * 
 * @author kp11
 *
 */
public class AminoAcidSequenceTest 
{

	private static final String TEST_DNA_SEQUENCE = "aaaaacaagaatacaaccacgactagaagcaggagtataatcatgattcaacaccagcatccacccccgcctcgacgccggcgtctactcctgcttgaagacgaggatgcagccgcggctggaggcgggggtgtagtcgtggtttaatactagtattcatcctcgtcttgatgctggtgtttattcttgtttZZZaaaaacaagaatacaaccacgactagaagcaggagtataatcatgattcaacaccagcatccacccccgcctcgacgccggcgtctactcctgcttgaagacgaggatgcagccgcggctggaggcgggggtgtagtcgtggtttaatactagtattcatcctcgtcttgatgctggtgtttattcttgttt";
	private static final String EXPECTED_AMINO_ACID_RESULT = "knknttttrsrsiimiqhqhpppprrrrllllededaaaaggggvvvv#y+yssss*cwclflf.knknttttrsrsiimiqhqhpppprrrrllllededaaaaggggvvvv#y+yssss*cwclflf";
	private static final String EXPECTED_SPACED_AMINO_ACID_RESULT = "k  n  k  n  t  t  t  t  r  s  r  s  i  i  m  i  q  h  q  h  p  p  p  p  r  r  r  r  l  l  l  l  e  d  e  d  a  a  a  a  g  g  g  g  v  v  v  v  #  y  +  y  s  s  s  s  *  c  w  c  l  f  l  f  .  k  n  k  n  t  t  t  t  r  s  r  s  i  i  m  i  q  h  q  h  p  p  p  p  r  r  r  r  l  l  l  l  e  d  e  d  a  a  a  a  g  g  g  g  v  v  v  v  #  y  +  y  s  s  s  s  *  c  w  c  l  f  l  f  ";
	private static final String SHORT_AMINO_ACID_SEQUENCE = "galmfwkqespvicyhrndt";
	private static final String SHORT_AMINO_ACID_SEQUENCE_WITH_STOP_CODONS = "galm#fwkqesp+vicyhrn*dt";

	
	/*
	 * Helper method for loading sequence data from file.
	 */
	private Bases loadBases() 
	{
		
		URL url = BasesTest.class.getResource("/etc/af063097.embl");

		final EntryInformation artemisEntryInformation = Options.getArtemisEntryInformation();
		final Document entryDocument = DocumentFactory.makeDocument(url.getFile());

		try {
			final uk.ac.sanger.artemis.io.Entry emblEntry = DocumentEntryFactory
					.makeDocumentEntry(artemisEntryInformation, entryDocument, null);
			Entry entry = new Entry(emblEntry);
		      
			return entry.getBases();
		} 
		catch (Exception e) {
			org.junit.Assert.fail(e.getMessage());
		}
		
		return null;
	}

	
	@Test
	public void testGetTranslationForString() 
	{
		
		boolean unknownIsXFlag = false;
		AminoAcidSequence seq = AminoAcidSequence.getTranslation(TEST_DNA_SEQUENCE, unknownIsXFlag);
		
		assertEquals("Check amino acid translation", EXPECTED_AMINO_ACID_RESULT, seq.toString());
		
		unknownIsXFlag = true;
		seq = AminoAcidSequence.getTranslation(TEST_DNA_SEQUENCE, unknownIsXFlag);
		
		assertEquals("Check amino acid translation with uknown=x flag set", EXPECTED_AMINO_ACID_RESULT.replace(".","x"), seq.toString());
	}
	
	@Test
	public void testGetTranslationForArray() 
	{
		
		boolean unknownIsXFlag = false;
		AminoAcidSequence seq = AminoAcidSequence.getTranslation(TEST_DNA_SEQUENCE.toCharArray(), unknownIsXFlag);
		
		assertEquals("Check amino acid translation", EXPECTED_AMINO_ACID_RESULT, seq.toString());
		
		unknownIsXFlag = true;
		seq = AminoAcidSequence.getTranslation(TEST_DNA_SEQUENCE.toCharArray(), unknownIsXFlag);
		
		assertEquals("Check amino acid translation with uknown=x flag set", EXPECTED_AMINO_ACID_RESULT.replace(".","x"), seq.toString());
	}
	
	@Test
	public void testGetSpacedTranslationString() 
	{
		
		boolean unknownIsXFlag = false;
		AminoAcidSequence seq = AminoAcidSequence.getSpacedTranslation(TEST_DNA_SEQUENCE, unknownIsXFlag);
		
		assertEquals("Check amino acid translation", EXPECTED_SPACED_AMINO_ACID_RESULT, seq.toString());
		
		unknownIsXFlag = true;
		seq = AminoAcidSequence.getSpacedTranslation(TEST_DNA_SEQUENCE, unknownIsXFlag);
		
		assertEquals("Check amino acid translation with uknown=x flag set", EXPECTED_SPACED_AMINO_ACID_RESULT.replace(".","x"), seq.toString());
	}
	
	@Test
	public void testGetSpacedTranslationArray() 
	{
		
		boolean unknownIsXFlag = false;
		AminoAcidSequence seq = AminoAcidSequence.getSpacedTranslation(TEST_DNA_SEQUENCE.toCharArray(), unknownIsXFlag);
		
		assertEquals("Check amino acid translation", EXPECTED_SPACED_AMINO_ACID_RESULT, seq.toString());
		
		unknownIsXFlag = true;
		seq = AminoAcidSequence.getSpacedTranslation(TEST_DNA_SEQUENCE.toCharArray(), unknownIsXFlag);
		
		assertEquals("Check amino acid translation with uknown=x flag set", EXPECTED_SPACED_AMINO_ACID_RESULT.replace(".","x"), seq.toString());
	}
	
	@Test
	public void testGetCodonTranslationString() throws Exception 
	{
		
		assertEquals("Check insufficient number of letters", '.', AminoAcidSequence.getCodonTranslation("aa"));
		assertEquals("Check single codon translation", 't', AminoAcidSequence.getCodonTranslation("act"));
		assertEquals("Check single stop codon translation 1", '*', AminoAcidSequence.getCodonTranslation("tga"));
		assertEquals("Check single stop codon translation 2", '#', AminoAcidSequence.getCodonTranslation("taa"));
		assertEquals("Check single stop codon translation 3", '+', AminoAcidSequence.getCodonTranslation("tag"));
		assertEquals("Check single start codon translation", 'm', AminoAcidSequence.getCodonTranslation("atg"));
			
	}
	
	@Test
	public void testGetCodonTranslationLetters() throws Exception 
	{
		
		assertEquals("Check incorrect base", '.', AminoAcidSequence.getCodonTranslation('a','a','Z'));
		assertEquals("Check single codon translation", 't', AminoAcidSequence.getCodonTranslation('a','c','t'));
		assertEquals("Check single stop codon translation 1", '*', AminoAcidSequence.getCodonTranslation('t','g','a'));
		assertEquals("Check single stop codon translation 2", '#', AminoAcidSequence.getCodonTranslation('t','a','a'));
		assertEquals("Check single stop codon translation 3", '+', AminoAcidSequence.getCodonTranslation('t','a','g'));
		assertEquals("Check single start codon translation", 'm', AminoAcidSequence.getCodonTranslation('a','t','g'));
			
	}
	
	@Test 
	public void testLength() throws Exception 
	{
		
		AminoAcidSequence seq = new AminoAcidSequence(SHORT_AMINO_ACID_SEQUENCE);
		
		assertEquals(20, seq.length());
		
	}
	
	@Test 
	public void testElementAt() throws Exception 
	{
		
		AminoAcidSequence seq = new AminoAcidSequence(SHORT_AMINO_ACID_SEQUENCE);
		
		assertEquals('g', seq.elementAt(0));
		assertEquals('w', seq.elementAt(5));
		assertEquals('t', seq.elementAt(SHORT_AMINO_ACID_SEQUENCE.length()-1));
	}
	
	@Test
	public void testGetMolecularWeight() throws Exception 
	{
		
		AminoAcidSequence seq = new AminoAcidSequence(SHORT_AMINO_ACID_SEQUENCE);
		assertEquals(2395.7551, seq.getMolecularWeight(), 0.0003f);
		
		seq = new AminoAcidSequence("");
		assertEquals(0.0f, seq.getMolecularWeight(), 0.1f);
		
	}
	
	@Test
	public void testToString() throws Exception 
	{
		
		AminoAcidSequence seq = new AminoAcidSequence(SHORT_AMINO_ACID_SEQUENCE);
		assertEquals(SHORT_AMINO_ACID_SEQUENCE, seq.toString());
		
	}
	
	@Test
	public void testCheckForMatch() throws Exception 
	{
		
		AminoAcidSequence seq = new AminoAcidSequence(SHORT_AMINO_ACID_SEQUENCE);
		
		assertTrue(seq.checkForMatch(new AminoAcidSequence(SHORT_AMINO_ACID_SEQUENCE)));
		assertTrue(seq.checkForMatch(new AminoAcidSequence("c" + SHORT_AMINO_ACID_SEQUENCE + "a")));
		assertTrue(seq.checkForMatch(new AminoAcidSequence("c" + SHORT_AMINO_ACID_SEQUENCE + "a")));
		assertFalse(seq.checkForMatch(new AminoAcidSequence("ca")));
		
		seq = new AminoAcidSequence("aaaaa");
		assertFalse(seq.checkForMatch(new AminoAcidSequence(SHORT_AMINO_ACID_SEQUENCE)));
		
		seq = new AminoAcidSequence("xxxxx");
		assertTrue("Wildcard match", seq.checkForMatch(new AminoAcidSequence(SHORT_AMINO_ACID_SEQUENCE)));
		
		seq = new AminoAcidSequence("les");
		assertTrue("Wildcard match", seq.checkForMatch(new AminoAcidSequence("xxxxx")));
		
		seq = new AminoAcidSequence("fwxq");
		assertTrue("Wildcard match", seq.checkForMatch(new AminoAcidSequence(SHORT_AMINO_ACID_SEQUENCE)));
		
		seq = new AminoAcidSequence("");
		assertTrue("Empty match", seq.checkForMatch(new AminoAcidSequence("")));
	}
	
	@Test
	public void testFindMatch() throws Exception 
	{
		// BTW length of strand is 33593.
		
		final Bases bases = loadBases();
		Marker marker = new Marker(bases.getForwardStrand(), 1);
		
		// TODO - Why does this method not pick up the matching g AA at position 1 - bug?
		
		// Search forwards on forward strand...
		
		AminoAcidSequence seq = new AminoAcidSequence("g");
		/* Args: final boolean search_backwards, final boolean search_fwd_strand, final boolean search_bwd_strand) */
		MarkerRange markRange = seq.findMatch(bases, marker, false, true, false);
        
		assertNotNull("Check we got a marker range back", markRange);
		assertEquals("Forward strand forward search match start position", 6, markRange.getStart().getPosition());
		assertEquals("Forward strand forward search match end position", 8, markRange.getEnd().getPosition());
		
		seq = new AminoAcidSequence("e");
		markRange = seq.findMatch(bases, marker, false, true, false);
		assertNotNull("Check we got a marker range back", markRange);
		assertEquals("Forward strand forward search match start position", 4, markRange.getStart().getPosition());
		assertEquals("Forward strand forward search match end position", 6, markRange.getEnd().getPosition());
		
		seq = new AminoAcidSequence("+s");
		markRange = seq.findMatch(bases, marker, false, true, false);
		assertNotNull("Check we got a marker range back", markRange);
		assertEquals("Forward strand forward search match start position", 81, markRange.getStart().getPosition());
		assertEquals("Forward strand forward search match end position", 86, markRange.getEnd().getPosition());
		
		
		// Search backwards on forward strand...
		
		marker = new Marker(bases.getForwardStrand(), bases.getForwardStrand().getSequenceLength());
		seq = new AminoAcidSequence("ag");
		markRange = seq.findMatch(bases, marker, true, true, false);
		assertNotNull("Check we got a marker range back", markRange);
		
		assertEquals("Forward strand reverse search match start position", 33581, markRange.getStart().getPosition());
		assertEquals("Forward strand reverse search match end position", 33586, markRange.getEnd().getPosition());
		
		seq = new AminoAcidSequence("*");
		markRange = seq.findMatch(bases, marker, true, true, false);
		assertNotNull("Check we got a marker range back", markRange);
		assertEquals("Forward strand reverse search match start position", 33564, markRange.getStart().getPosition());
		assertEquals("Forward strand reverse search match end position", 33566, markRange.getEnd().getPosition());
		
		
		// Search forwards on reverse strand...
		
		marker = new Marker(bases.getReverseStrand(), 1);
		seq = new AminoAcidSequence("v");
		markRange = seq.findMatch(bases, marker, true, false, true);
		assertNotNull("Check we got a marker range back", markRange);
		assertEquals("Reverse strand forward search match start position", 1, markRange.getStart().getPosition());
		assertEquals("Reverse strand forward search match end position", 3, markRange.getEnd().getPosition());
		
		// Search backwards on reverse strand...
		
		marker = new Marker(bases.getReverseStrand(), bases.getLength());
		seq = new AminoAcidSequence("e");
		markRange = seq.findMatch(bases, marker, false, false, true);
		assertNotNull("Check we got a marker range back", markRange);
		assertEquals("Reverse strand backward search match start position", 33541, markRange.getStart().getPosition());
		assertEquals("Reverse strand backward search match end position", 33543, markRange.getEnd().getPosition());
		
		
		// Test dud search criteria
		marker = new Marker(bases.getForwardStrand(), 1);
		seq = new AminoAcidSequence("ZZZZ");
		markRange = seq.findMatch(bases, marker, false, true, false);
		assertNull("Confirn that we did not get a match", markRange);

	}
	
	@Test
	public void testContainsStopCodon() throws Exception
	{
		
		AminoAcidSequence seq = new AminoAcidSequence(SHORT_AMINO_ACID_SEQUENCE_WITH_STOP_CODONS);
		assertTrue(seq.containsStopCodon());
		
		seq = new AminoAcidSequence(SHORT_AMINO_ACID_SEQUENCE);
		assertFalse(seq.containsStopCodon());
		
		seq = new AminoAcidSequence("");
		assertFalse(seq.containsStopCodon());
		
	}
	
	@Test
	public void testIsStopCodon() throws Exception 
	{
		
		assertTrue(AminoAcidSequence.isStopCodon('*'));
		assertTrue(AminoAcidSequence.isStopCodon('+'));
		assertTrue(AminoAcidSequence.isStopCodon('#'));
		assertFalse(AminoAcidSequence.isStopCodon('a'));
		assertFalse(AminoAcidSequence.isStopCodon('_'));
	}
	
	@Test
	public void testGetOneLetterCode() throws Exception 
	{
		assertEquals('g', AminoAcidSequence.getOneLetterCode("Gly"));
		assertEquals('a', AminoAcidSequence.getOneLetterCode("Ala"));
		assertEquals('l', AminoAcidSequence.getOneLetterCode("Leu"));
		assertEquals('m', AminoAcidSequence.getOneLetterCode("Met"));
		assertEquals('f', AminoAcidSequence.getOneLetterCode("Phe"));
		assertEquals('w', AminoAcidSequence.getOneLetterCode("Trp"));
		assertEquals('k', AminoAcidSequence.getOneLetterCode("Lys"));
		assertEquals('q', AminoAcidSequence.getOneLetterCode("Gln"));
		assertEquals('e', AminoAcidSequence.getOneLetterCode("Glu"));
		assertEquals('s', AminoAcidSequence.getOneLetterCode("Ser"));
		assertEquals('p', AminoAcidSequence.getOneLetterCode("Pro"));
		assertEquals('v', AminoAcidSequence.getOneLetterCode("Val"));
		assertEquals('i', AminoAcidSequence.getOneLetterCode("Ile"));
		assertEquals('c', AminoAcidSequence.getOneLetterCode("Cys"));
		assertEquals('y', AminoAcidSequence.getOneLetterCode("Tyr"));
		assertEquals('h', AminoAcidSequence.getOneLetterCode("His"));
		assertEquals('r', AminoAcidSequence.getOneLetterCode("Arg"));
		assertEquals('n', AminoAcidSequence.getOneLetterCode("Asn"));
		assertEquals('d', AminoAcidSequence.getOneLetterCode("Asp"));
		assertEquals('t', AminoAcidSequence.getOneLetterCode("Thr"));
		
		assertEquals('g', AminoAcidSequence.getOneLetterCode("GLY"));
		assertEquals('a', AminoAcidSequence.getOneLetterCode("ALA"));
	}
	
	@Test
	public void testGetThreeLetterAbbreviation() throws Exception 
	{
		assertEquals("Gly", AminoAcidSequence.getThreeLetterAbbreviation('g'));
		assertEquals("Ala", AminoAcidSequence.getThreeLetterAbbreviation('a'));
		assertEquals("Leu", AminoAcidSequence.getThreeLetterAbbreviation('l'));
		assertEquals("Met", AminoAcidSequence.getThreeLetterAbbreviation('m'));
		assertEquals("Phe", AminoAcidSequence.getThreeLetterAbbreviation('f'));
		assertEquals("Trp", AminoAcidSequence.getThreeLetterAbbreviation('w'));
		assertEquals("Lys", AminoAcidSequence.getThreeLetterAbbreviation('k'));
		assertEquals("Gln", AminoAcidSequence.getThreeLetterAbbreviation('q'));
		assertEquals("Glu", AminoAcidSequence.getThreeLetterAbbreviation('e'));
		assertEquals("Ser", AminoAcidSequence.getThreeLetterAbbreviation('s'));
		assertEquals("Pro", AminoAcidSequence.getThreeLetterAbbreviation('p'));
		assertEquals("Val", AminoAcidSequence.getThreeLetterAbbreviation('v'));
		assertEquals("Ile", AminoAcidSequence.getThreeLetterAbbreviation('i'));
		assertEquals("Cys", AminoAcidSequence.getThreeLetterAbbreviation('c'));
		assertEquals("Tyr", AminoAcidSequence.getThreeLetterAbbreviation('y'));
		assertEquals("His", AminoAcidSequence.getThreeLetterAbbreviation('h'));
		assertEquals("Arg", AminoAcidSequence.getThreeLetterAbbreviation('r'));
		assertEquals("Asn", AminoAcidSequence.getThreeLetterAbbreviation('n'));
		assertEquals("Asp", AminoAcidSequence.getThreeLetterAbbreviation('d'));
		assertEquals("Thr", AminoAcidSequence.getThreeLetterAbbreviation('t'));
		
		assertEquals("Gly", AminoAcidSequence.getThreeLetterAbbreviation('G'));
		assertEquals("Ala", AminoAcidSequence.getThreeLetterAbbreviation('A'));
	}
	
	@Test
	public void testGetSymbolIndex() throws Exception 
	{
		
		assertEquals(7, AminoAcidSequence.getSymbolIndex('g'));
		assertEquals(0, AminoAcidSequence.getSymbolIndex('a'));
		assertEquals(10, AminoAcidSequence.getSymbolIndex('l'));
		assertEquals(12, AminoAcidSequence.getSymbolIndex('m'));
		assertEquals(13, AminoAcidSequence.getSymbolIndex('f'));
		assertEquals(17, AminoAcidSequence.getSymbolIndex('w'));
		assertEquals(11, AminoAcidSequence.getSymbolIndex('k'));
		assertEquals(5, AminoAcidSequence.getSymbolIndex('q'));
		assertEquals(6, AminoAcidSequence.getSymbolIndex('e'));
		assertEquals(15, AminoAcidSequence.getSymbolIndex('s'));
		assertEquals(14, AminoAcidSequence.getSymbolIndex('p'));
		assertEquals(19, AminoAcidSequence.getSymbolIndex('v'));
		assertEquals(9, AminoAcidSequence.getSymbolIndex('i'));
		assertEquals(4, AminoAcidSequence.getSymbolIndex('c'));
		assertEquals(18, AminoAcidSequence.getSymbolIndex('y'));
		assertEquals(8, AminoAcidSequence.getSymbolIndex('h'));
		assertEquals(1, AminoAcidSequence.getSymbolIndex('r'));
		assertEquals(2, AminoAcidSequence.getSymbolIndex('n'));
		assertEquals(3, AminoAcidSequence.getSymbolIndex('d'));
		assertEquals(16, AminoAcidSequence.getSymbolIndex('t'));
		assertEquals(20, AminoAcidSequence.getSymbolIndex('*'));
		assertEquals(21, AminoAcidSequence.getSymbolIndex('#'));
		assertEquals(22, AminoAcidSequence.getSymbolIndex('+'));
		assertEquals(23, AminoAcidSequence.getSymbolIndex('.'));
		assertEquals(23, AminoAcidSequence.getSymbolIndex('x'));
		assertEquals(24, AminoAcidSequence.getSymbolIndex('u'));
		
		assertEquals(7, AminoAcidSequence.getSymbolIndex('G'));
		assertEquals(0, AminoAcidSequence.getSymbolIndex('A'));
	}
	
	@Test
	public void testGetSymbolFromIndex() throws Exception 
	{
		
		assertEquals('g', AminoAcidSequence.getSymbolFromIndex(7));
		assertEquals('a', AminoAcidSequence.getSymbolFromIndex(0));
		assertEquals('l', AminoAcidSequence.getSymbolFromIndex(10));
		assertEquals('m', AminoAcidSequence.getSymbolFromIndex(12));
		assertEquals('f', AminoAcidSequence.getSymbolFromIndex(13));
		assertEquals('w', AminoAcidSequence.getSymbolFromIndex(17));
		assertEquals('k', AminoAcidSequence.getSymbolFromIndex(11));
		assertEquals('q', AminoAcidSequence.getSymbolFromIndex(5));
		assertEquals('e', AminoAcidSequence.getSymbolFromIndex(6));
		assertEquals('s', AminoAcidSequence.getSymbolFromIndex(15));
		assertEquals('p', AminoAcidSequence.getSymbolFromIndex(14));
		assertEquals('v', AminoAcidSequence.getSymbolFromIndex(19));
		assertEquals('i', AminoAcidSequence.getSymbolFromIndex(9));
		assertEquals('c', AminoAcidSequence.getSymbolFromIndex(4));
		assertEquals('y', AminoAcidSequence.getSymbolFromIndex(18));
		assertEquals('h', AminoAcidSequence.getSymbolFromIndex(8));
		assertEquals('r', AminoAcidSequence.getSymbolFromIndex(1));
		assertEquals('n', AminoAcidSequence.getSymbolFromIndex(2));
		assertEquals('d', AminoAcidSequence.getSymbolFromIndex(3));
		assertEquals('t', AminoAcidSequence.getSymbolFromIndex(16));
		assertEquals('*', AminoAcidSequence.getSymbolFromIndex(20));
		assertEquals('#', AminoAcidSequence.getSymbolFromIndex(21));
		assertEquals('+', AminoAcidSequence.getSymbolFromIndex(22));
		assertEquals('.', AminoAcidSequence.getSymbolFromIndex(23));
		assertEquals('u', AminoAcidSequence.getSymbolFromIndex(24));
		
		try 
		{
			AminoAcidSequence.getSymbolFromIndex(25);
			Assert.fail("Expected an Error to be thrown");
		}
		catch (Throwable t) 
		{
			// Expected - Error throw
		}
		
	}
}
