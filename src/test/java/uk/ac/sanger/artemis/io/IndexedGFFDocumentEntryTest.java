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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.mockito.Mockito.when;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.HashMap;

import org.junit.Before;
import org.junit.Test;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;

import htsjdk.samtools.util.BlockCompressedInputStream;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.components.variant.TabixReader;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.util.FileDocument;
import uk.ac.sanger.artemis.util.OutOfRangeException;

/**
 * Unit test for the IndexedGFFDocumentEntry class.
 * 
 * @author kp11
 *
 */
public class IndexedGFFDocumentEntryTest
{

	private final String INDEXED_GFF_FILE = "/data/gff/indexed-gff.gff.gz";
	private final String INDEXED_GFF_FASTA_FILE = "/data/fasta/non-indexed-fasta.fa";
	
	private static final HashMap<Integer, String[]> ENUMERATION_FEATURES = new HashMap<Integer, String[]>();
	static {
		// Check both contigs
		ENUMERATION_FEATURES.put(1, new String[]{ "Pfalciparum_REP_20", "1", "360" });
		ENUMERATION_FEATURES.put(2, new String[]{ "Pfalciparum_REP_15", "361", "1418" });
		ENUMERATION_FEATURES.put(3, new String[]{ "Pfalciparum_REP_35", "2160", "3858" });
		ENUMERATION_FEATURES.put(4, new String[]{ "Pfalciparum_REP_5", "8856", "9021" });
		ENUMERATION_FEATURES.put(5, new String[]{ "Pfalciparum_REP_25", "9313", "9529" });
		ENUMERATION_FEATURES.put(6, new String[]{ "Pfalciparum_REP_55", "11315", "27336" });
		ENUMERATION_FEATURES.put(7, new String[]{ "PF3D7_0100100", "29510", "37126" });
		ENUMERATION_FEATURES.put(8, new String[]{ "PF3D7_0100100.1", "29510", "37126" });
		ENUMERATION_FEATURES.put(9, new String[]{ "PF3D7_0100100.1:pep", "29510", "37126" });
		ENUMERATION_FEATURES.put(10, new String[]{ "PF3D7_0100200", "38982", "40207" });
		ENUMERATION_FEATURES.put(11, new String[]{ "PF3D7_0100200.1", "38982", "40207" });
		ENUMERATION_FEATURES.put(12, new String[]{ "PF3D7_0100200.1:pep", "38982", "40207" });
		ENUMERATION_FEATURES.put(13, new String[]{ "PF3D7_0100300", "42367", "46507" });
		
		ENUMERATION_FEATURES.put(629, new String[]{ "PF3D7_0109600.1:exon{3,2,1}", "377874", "378599" });
		ENUMERATION_FEATURES.put(630, new String[]{ "Pf3D7_02__new:1..1142", "640852", "641993" });
		
		ENUMERATION_FEATURES.put(1579, new String[]{ "PF3D7_0200600.1:exon{2,1}", "686137", "687409" });
		ENUMERATION_FEATURES.put(1580, new String[]{ "PF3D7_0214200.1:exon{2,1}", "1216006", "1216905" });
		ENUMERATION_FEATURES.put(1581, new String[]{ "PF3D7_0216500.1:exon{1,2}", "1319338", "1324554" });
		ENUMERATION_FEATURES.put(1582, new String[]{ "PF3D7_0218000.1:exon{1}", "1379711", "1380703" });
	};
	
	private final static String [][] EXPECTED_CONTIG1_FEATURES = 
		{
				{ "Pfalciparum_REP_15", "361", "1418" },
				{ "Pfalciparum_REP_35", "2160", "3858" },
				{ "Pfalciparum_REP_5", "8856", "9021" },
				{ "Pfalciparum_REP_25", "9313", "9529" },
				{ "Pfalciparum_REP_55", "11315", "27336" },
				{ "PF3D7_0100100", "29510", "37126" },
				{ "PF3D7_0100100.1", "29510", "37126" },
				{ "PF3D7_0100100.1:pep", "29510", "37126" },
				{ "PF3D7_0100100.1:exon{1,2}", "29510", "37126" }
		};
	
	private final static String [][] EXPECTED_CONTIG2_FEATURES = 
		{
				{ "Pf3D7_02__new:1..1142", "640852", "641993" },
				{ "PF3D7_0200100", "666083", "672019" },
				{ "PF3D7_0200100.1", "666083", "672019" },
				{ "PF3D7_0200100.1:pep", "666083", "672019" },
				{ "PF3D7_0200100.1:exon{1,2}", "666083", "672019" }
		};
	
	@Mock
	private FileDocument doc;
	
	@Before
	public void setUp() throws Exception {

		MockitoAnnotations.initMocks(this);
	}
	
	/**
	 * Load a data file.
	 * @param file String
	 * @return Entry
	 */
	protected Entry loadFile(String file)
	{
		return Utils.getEntry(file);
	}
	
	protected File getFile(String filePath)
	{
		URL filePathUrl = IndexedGFFDocumentEntryTest.class.getResource(filePath);
		File file = new File(filePathUrl.getFile());
		
		return file;
	}
	
	/**
	 * Test the query method.
	 * @throws IOException 
	 */
	@Test
	public void testQuery() throws IOException 
	{	
		TabixReader reader = new TabixReader(
				IndexedGFFDocumentEntryTest.class.getResource(INDEXED_GFF_FILE).getFile());

		TabixReader.Iterator tabixIterator = reader.query(-1, 1, 19);
		assertNull(tabixIterator);
	}
	
	/**
	 * General testing of feature enumeration code.
	 * @throws IOException
	 * @throws OutOfRangeException
	 * @throws InvalidRelationException
	 */
	@Test
	public void testFeatureEnumeration() throws IOException, OutOfRangeException, InvalidRelationException
	{
		int numFeatures = 0;
		
		// Given 
		
		// We don't want GUI pop-ups displayed during unit tests...
		DocumentEntryFactory.setDisplayIndexingQuestionForGffs(false);
		
		final File file = getFile(INDEXED_GFF_FILE);
		final BlockCompressedInputStream stream = new BlockCompressedInputStream(file);

		// When
		
		when( doc.getFile() ).thenReturn (file);
		when( doc.getInputStream() ).thenReturn (stream);
		
		final IndexedGFFDocumentEntry indexedGffEntry = new IndexedGFFDocumentEntry(doc);
		final Entry fastaEntry = loadFile(INDEXED_GFF_FASTA_FILE);
		
		final Bases bases = new Bases(fastaEntry.getSequence());
		final EntryGroup group = new SimpleEntryGroup(bases);
		final uk.ac.sanger.artemis.Entry artFastaEntry = new uk.ac.sanger.artemis.Entry(bases, fastaEntry);
		final uk.ac.sanger.artemis.Entry artGffEntry = new uk.ac.sanger.artemis.Entry(bases, indexedGffEntry);
		group.add(artFastaEntry);
		group.add(artGffEntry);
		
		// Then
		
		final uk.ac.sanger.artemis.io.FeatureEnumeration feature_enumerator = indexedGffEntry.features();
		while (feature_enumerator.hasMoreFeatures())
		{
			Feature feature = feature_enumerator.nextFeature();
			
			++numFeatures;
			
			/*System.out.println(
					numFeatures + ": " + 
					feature.getQualifierByName("ID").getValues().elementAt(0) + " " + 
					feature.getFirstBase() + " " + 
					feature.getLastBase());*/
			
			// Check a sample of values...
			
			if (ENUMERATION_FEATURES.containsKey(numFeatures)) 
			{
				assertEquals("Feature index " + numFeatures + " ID", ENUMERATION_FEATURES.get(numFeatures)[0], feature.getQualifierByName("ID").getValues().elementAt(0));
				assertEquals("Feature index " + numFeatures + " start", ENUMERATION_FEATURES.get(numFeatures)[1], Integer.toString(feature.getFirstBase()));
				assertEquals("Feature index " + numFeatures + " end", ENUMERATION_FEATURES.get(numFeatures)[2], Integer.toString(feature.getLastBase()));
			}
		}

		assertEquals(1582, numFeatures);
		
		// Test getListOfContigs while we're at it - 3 contigs (one with no features)
		assertEquals(3, indexedGffEntry.getListOfContigs().size());
		
		//Contig 1
		FeatureVector list = indexedGffEntry.getFeaturesInRange(new Range(361, 34762));
		assertEquals(9, list.size());
		
		int idx = 0;
		for (Feature feature : list)
		{
			//System.out.println("Contig1: " + feature.getQualifierByName("ID").getValues().elementAt(0) + " " + feature.getFirstBase() + " " + feature.getLastBase());
			assertEquals(EXPECTED_CONTIG1_FEATURES[idx][0], feature.getQualifierByName("ID").getValues().elementAt(0));
			assertEquals(EXPECTED_CONTIG1_FEATURES[idx][1], Integer.toString(feature.getFirstBase()));
			assertEquals(EXPECTED_CONTIG1_FEATURES[idx][2], Integer.toString(feature.getLastBase()));
			
			++idx;
		}
		
		// Contig 2 - start features
		list = indexedGffEntry.getFeaturesInRange(new Range(640852, 672020));
		assertEquals(5, list.size());
		
		idx = 0;
		for (Feature feature : list)
		{
			//System.out.println("Contig2: " + feature.getQualifierByName("ID").getValues().elementAt(0) + " " + feature.getFirstBase() + " " + feature.getLastBase());
			assertEquals(EXPECTED_CONTIG2_FEATURES[idx][0], feature.getQualifierByName("ID").getValues().elementAt(0));
			assertEquals(EXPECTED_CONTIG2_FEATURES[idx][1], Integer.toString(feature.getFirstBase()));
			assertEquals(EXPECTED_CONTIG2_FEATURES[idx][2], Integer.toString(feature.getLastBase()));
			
			++idx;
		}
		
	}
}
