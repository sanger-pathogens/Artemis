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
import static org.mockito.Mockito.when;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.URL;
import java.util.zip.GZIPInputStream;

import org.apache.commons.lang3.StringUtils;
import org.junit.Before;
import org.junit.Test;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.util.FileDocument;


/**
 * JUnit test suite for the SimpleDocumentEntry class
 * and underlying functionality (more of an integration test really!).
 * 
 * This test module carries out multiple tests in each test method
 * - this is done for efficiency reasons to avoid having to 
 *   repeatedly load test files.
 * 
 * @author kp11
 *
 */
public class SimpleDocumentEntryTest
{

	// GFF
	
	// Indexed
	private final String INDEXED_GFF_FILE = "/data/gff/indexed-gff.gff.gz";
	private final String INDEXED_GFF_FASTA_FILE = "/data/gff/indexed-gff.fa";
	private final int NUM_INDEXED_GFF_FEATURES = 2238;
	
	// Non-indexed
	private final String GFF_FILE = "/data/gff/notindexed.gff.gz";
	
	private final String GFF_OUTPUT_FILE = "out.gff3";
	private final int NUM_GFF_FEATURES = 1584;
	private final int GFF_SEQ_LEN = 1587953;
	private final String GFF_FEATURE_ID_KEY = "ID";
	
	private final String GFF_SEQ_START = "tgaaccctaaaacctaaaccctaaaccctaaaccctgaaccctaaaccctgaaccctaaaccctaaaccctgaac";
	private final String GFF_SEQ_END = "tagggttcagggtttaggtgtcagggttca";
			
	private final String GFF_HDRS = 
			"##gff-version 3\n" + 
			"##sequence-region Pf3D7_01_v3 1 640851\n" + 
			"##sequence-region Pf3D7_02_v3 1 947102\n";
	
	// GTF
	
	private final String GTF_FILE = "/data/gtf/clostridium_cellulosi.gtf.gz";
	
	private final String GTF_OUTPUT_FILE = "out.gtf";
	private final String GTF_FASTA_FILE = "/data/fasta/clostridium_cellulosi.fa.gz";
	private final int NUM_GTF_FEATURES = 12315;
	private final int GTF_SEQ_LEN = 2229578;
	
	private final String GTF_SEQ_START = "atgcagtccttttcggaagtatggt";
	private final String GTF_SEQ_END = "aaacaaattaaatgttcggaggtttaagagc";
	
	private final String GTF_HDRS = 
			"#!genome-build DG5\n" + 
			"#!genome-version DG5\n" + 
			"#!genome-date 2014-08\n" + 
			"#!genome-build-accession GCA_000953215.1\n" + 
			"#!genebuild-last-updated 2014-08\n";
	
	// EMBL
	
	private final String EMBL_FILE = "/data/embl/MAL1.embl.gz";
	private final String EMBL_OUTPUT_FILE = "out.embl";
	private final int NUM_EMBL_FEATURES = 248;
	private final int EMBL_SEQ_LEN = 643292;
	private final String EMBL_FEATURE_ID_KEY = "locus_tag";
	
	private final String EMBL_SEQ_START = "ctaaacctaaacctaaaccctgaaccctaaaccctaaaccctgaaccctaaaccctgaaccctgaaccctaaac";
	private final String EMBL_SEQ_END = "ttagggttcagggtttagggtttagggtttagggaatagggt";
		
	private final String EMBL_HDRS = 
			"ID   Pf3D7_01.embl; SV ; ; ; ; ; 643292 BP.\n" + 
			"FH   Key             Location/Qualifiers\n" + 
			"FH\n";
	
	// Genbank
	
	private final String GBK_FILE = "/data/genbank/SeaOtter_chrMT.gbk.gz";
	private final String GBK_OUTPUT_FILE = "out.gbk";
	private final int NUM_GBK_FEATURES = 53;
	private final int GBK_SEQ_LEN = 16431;
	
	private final String GBK_SEQ_START = "gttaatgtagcttataaataaagcaaggcactgaaaatgcctagaagagtcataagactc" + 
			"catagacataaaggtttggtcctagccttcctattaattattaacagaattacacatgca";
	private final String GBK_SEQ_END = "acttatactggtgccacgcatgttaatctcacttactaatccattaaaacttcctattcaaaatgaagc" +
			"tatctatagatgtgaattcccacctctatcacccccggactt";
		
	private final String GBK_HDRS = 
			"LOCUS       NC_009692              16431 bp    DNA     circular MAM 15-APR-2009\n" + 
			"DEFINITION  Enhydra lutris mitochondrion, complete genome.\n" + 
			"ACCESSION   NC_009692\n" + 
			"VERSION     NC_009692.1\n" + 
			"DBLINK      Project: 20145\n" + 
			"            BioProject: PRJNA20145\n" + 
			"KEYWORDS    RefSeq.\n" + 
			"SOURCE      mitochondrion Enhydra lutris (sea otter)\n" + 
			"  ORGANISM  Enhydra lutris\n" + 
			"            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;\n" + 
			"            Mammalia; Eutheria; Laurasiatheria; Carnivora; Caniformia;\n" + 
			"            Mustelidae; Lutrinae; Enhydra.\n" + 
			"REFERENCE   1  (bases 1 to 16431)\n" + 
			"  AUTHORS   Yonezawa,T., Nikaido,M., Kohno,N., Fukumoto,Y., Okada,N. and\n" + 
			"            Hasegawa,M.\n" + 
			"  TITLE     Molecular phylogenetic study on the origin and evolution of\n" + 
			"            Mustelidae\n" + 
			"  JOURNAL   Gene 396 (1), 1-12 (2007)\n" + 
			"   PUBMED   17449200\n" + 
			"REFERENCE   2  (bases 1 to 16431)\n" + 
			"  CONSRTM   NCBI Genome Project\n" + 
			"  TITLE     Direct Submission\n" + 
			"  JOURNAL   Submitted (24-JUL-2007) National Center for Biotechnology\n" + 
			"            Information, NIH, Bethesda, MD 20894, USA\n" + 
			"REFERENCE   3  (bases 1 to 16431)\n" + 
			"  AUTHORS   Yonezawa,T., Nikaido,M., Kohno,N., Fukumoto,Y., Okada,N. and\n" + 
			"            Hasegawa,M.\n" + 
			"  TITLE     Direct Submission\n" + 
			"  JOURNAL   Submitted (22-JAN-2007) The Graduate University for Advanced\n" + 
			"            Studies (SOKENDAI), The Department of Biosystems Science; Shonan\n" + 
			"            Village, Hayama, Kanagawa 240-0193, Japan\n" + 
			"COMMENT     REVIEWED REFSEQ: This record has been curated by NCBI staff. The\n" + 
			"            reference sequence was derived from AB291077.\n" + 
			"            COMPLETENESS: full length.\n" + 
			"FEATURES             Location/Qualifiers\n";
	
	// Fasta
	
	private final String FASTA_FILE = "/data/fasta/test.fa.gz";
	private final String FASTA_OUTPUT_FILE = "out.fa";
	private final int FASTA_SEQ_LEN = 240;
	private final String FASTA_TEXT = 
			"atgacaaggcttccattactaaaacgacctcgcagaaaccgaaaaagtgcagccgttcga" + 
			"tctataattcaagaaacccaactctgttctagtgacttgatctggcccatctttcttaaa" + 
			"tctataattcaagaaacccaactctgttctagtgacttgatctggcccatctttcttaaa" + 
			"ttgacaaggcttccattactaaaacgacctcgcagaaaccgaaaaagtgcagccgttcga";
	
	// Blast
	private final String BLAST_FILE = "/data/blast/blast_test_results.blast";
	private final String BLAST_FASTA_FILE = "/data/blast/blast_test.fa";
	private final String BLAST_OUTPUT_FILE = "out.blast";
	
	private final int NUM_BLAST_FEATURES = 46;
	
	public static final boolean LOG = false;
	
	@Mock
	private FileDocument saveToDoc;
	
	@Before
	public void setUp() throws Exception {

		MockitoAnnotations.initMocks(this);
	}
	
	
	// =========================================================
	// Utility methods
	// =========================================================
	
	protected File getFile(String filePath)
	{
		URL filePathUrl = SimpleDocumentEntryTest.class.getResource(filePath);
		File file = new File(filePathUrl.getFile());
		
		return file;
	}
	
	// Basic file comparison
	protected boolean doBasicFileComparison(File file1, File file2, boolean ignoreCase) throws IOException
	{
		boolean result = true;
		
		final boolean isFile1Zipped = file1.getName().endsWith(".gz");
		final boolean isFile2Zipped = file2.getName().endsWith(".gz");
		
		try (
			final BufferedReader fileStream1 = new BufferedReader( 
					new InputStreamReader(
							( isFile1Zipped ?  new GZIPInputStream(new FileInputStream(file1)) : new FileInputStream(file1) ) 
					)
			);
			final BufferedReader fileStream2 = new BufferedReader( 
					new InputStreamReader(
							( isFile2Zipped ?  new GZIPInputStream(new FileInputStream(file2)) : new FileInputStream(file2) ) 
					)
			);
		)
		{
			boolean finished = false;
			while (!finished)
			{
				String line1 = fileStream1.readLine();
				if (line1 != null)
					line1 = line1.replaceAll("\\r\\n|\\r|\\n", "");
				
				String line2 = fileStream2.readLine();
				if (line2 != null)
					line2 = line2.replaceAll("\\r\\n|\\r|\\n", "");
				
				result = ignoreCase ? StringUtils.equalsIgnoreCase(line1, line2) : StringUtils.equals(line1, line2);
				if (!result)
					return result;
				
				if (line1 == null && line2 == null) 
				{
					finished = true;
				}
			}
		}
		
		return result;
	}
	
	protected void printFeatures(String title, Entry entry) throws InvalidRelationException
	{
		final FeatureEnumeration features = entry.features();
		int i = 0;
		
		System.out.println(title + ": ");
				
		while (features.hasMoreFeatures())
		{
			Feature ft = features.nextFeature();
			
			if (ft.getQualifierByName(GFF_FEATURE_ID_KEY) != null)
			{
				log("printFeatures",
						"" + i + ": " + ft.getKey().getKeyString() + " " + 
						ft.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) + " " +
						ft.getFirstBase() + " " +
						ft.getLastBase());
			}
			else
			{
				log("printFeatures", "" + i + ": " + ft.getKey().getKeyString());
			}
			
			++i;
		}
	}
	
	protected void printFeatures(String title, uk.ac.sanger.artemis.Entry entry) throws InvalidRelationException
	{
		final int numFeatures = entry.getFeatureCount();
		
		System.out.println(title + ": ");
				
		for (int i = 0; i < numFeatures; i++)
		{
			uk.ac.sanger.artemis.Feature ft = entry.getFeature(i);
			
			if (ft.getQualifierByName(GFF_FEATURE_ID_KEY) != null)
			{
				log("printFeatures",
						"" + i + ": " + ft.getKey().getKeyString() + " " + 
						ft.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) + " " +
						ft.getFirstBase() + " " +
						ft.getLastBase());
			}
			else
			{
				log("printFeatures", "" + i + ": " + ft.getKey().getKeyString());
			}
		}
	}
	
	/**
	 * Load a data file.
	 * @param file String
	 * @return Entry
	 */
	protected Entry loadFile(String file)
	{
		long startTime = System.currentTimeMillis();
		final Entry entry = Utils.getEntry(file);
		long endTime = System.currentTimeMillis();
		
		log("loadFile", "File load time (" + file + "): " + (endTime - startTime) + " milliseconds");
		
		return entry;
	}
	
	/**
	 * Output some useful info for testing purposes.
	 * @param method String
	 * @param message String
	 */
	public void log(String method, String message)
	{
		if (LOG)
		{
			System.out.println(method + ": " + message);
		}
	}
	
	// =========================================================
	// Tests
	// =========================================================
	
	/**
	 * Exercise Embl file loading/saving.
	 * 
	 * @throws Exception
	 */
	@Test
	public void testEmblFile() throws Exception
	{
		// Given
		
		final File savefile = File.createTempFile(EMBL_OUTPUT_FILE, "");
		final OutputStream outStream = new FileOutputStream(savefile);
		
		log("testEmblFile", "Writing to file: " + savefile.getAbsolutePath());

		// When
		
		when( saveToDoc.getFile() ).thenReturn (savefile);
		when( saveToDoc.getWriter() ).thenReturn (new FileWriter(savefile));
		when( saveToDoc.getOutputStream() ).thenReturn (outStream);
		
		final Entry entry = loadFile(EMBL_FILE);
		
		assertTrue(entry instanceof EmblDocumentEntry); 
		
		((EmblDocumentEntry)entry).save(saveToDoc);
		
		
		// Then
		
		assertEquals( NUM_EMBL_FEATURES, entry.getFeatureCount() );
		assertEquals( EMBL_HDRS, entry.getHeaderText() );
		
		Sequence sequence = entry.getSequence();
		String seqString = String.valueOf( sequence.getCharSubSequence(1, sequence.length()) );
		
		assertTrue(sequence instanceof EmblStreamSequence);
		assertEquals( EMBL_SEQ_LEN, seqString.length() );
		assertEquals( 261151, sequence.getACount() );
		assertEquals( 65283, sequence.getCCount() );
		assertEquals( 66886, sequence.getGCount() );
		assertEquals( 249972, sequence.getTCount() );
		assertTrue(seqString.indexOf(EMBL_SEQ_START) == 0);
		assertTrue(seqString.lastIndexOf(EMBL_SEQ_END) == (seqString.length()-EMBL_SEQ_END.length()));

		//printFeatures("testEmblFile", entry);
		
		// Check the feature table
		assertNotNull(((SimpleDocumentEntry)entry).getFeatureTable());
		assertEquals(NUM_EMBL_FEATURES, ((SimpleDocumentEntry)entry).getFeatureTable().getFeatureCount());
		
		int idx = 0;
		
		// Repeat region
		
		Feature nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "repeat_region", nextFeature.getKey().getKeyString() );
		assertEquals( "Pfalciparum_REP_20", nextFeature.getQualifierByName(EMBL_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 1, nextFeature.getFirstBase() );
		assertEquals( 583, nextFeature.getLastBase() );
				
		idx = 14;
		nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "ncRNA", nextFeature.getKey().getKeyString() );
		assertEquals( "RNAzID:26", nextFeature.getQualifierByName(EMBL_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 40820, nextFeature.getFirstBase() );
		assertEquals( 41021, nextFeature.getLastBase() );
			
		idx = 18;
		nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "CDS", nextFeature.getKey().getKeyString() );
		assertEquals( "PFA0020w", nextFeature.getQualifierByName(EMBL_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 50586, nextFeature.getFirstBase() );
		assertEquals( 51859, nextFeature.getLastBase() );
		
		nextFeature = entry.getFeatureAtIndex(++idx);
		assertEquals( "CDS_motif", nextFeature.getKey().getKeyString() );
		assertEquals( "cdsMotif50846-50860", nextFeature.getQualifierByName(EMBL_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 50847, nextFeature.getFirstBase() );
		assertEquals( 50860, nextFeature.getLastBase() );
		
		idx = 75;
		nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "5'UTR", nextFeature.getKey().getKeyString() );
		assertEquals( "PFA0180w:mRNA", nextFeature.getQualifierByName(EMBL_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 161286, nextFeature.getFirstBase() );
		assertEquals( 161365, nextFeature.getLastBase() );
		
		idx = 83;
		nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "3'UTR", nextFeature.getKey().getKeyString() );
		assertEquals( "PFA0205w:mRNA", nextFeature.getQualifierByName(EMBL_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 181026, nextFeature.getFirstBase() );
		assertEquals( 181049, nextFeature.getLastBase() );
		
		idx = 88;
		nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "CDS", nextFeature.getKey().getKeyString() );
		assertEquals( "PFA0225w", nextFeature.getQualifierByName(EMBL_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 202774, nextFeature.getFirstBase() );
		assertEquals( 204381, nextFeature.getLastBase() );
		
		// Test indexOf and contains
		assertEquals(idx, entry.indexOf(nextFeature));
		assertTrue(entry.contains(nextFeature));
		
		idx = 177;
		nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "rRNA", nextFeature.getKey().getKeyString() );
		assertEquals( "MAL1_28s", nextFeature.getQualifierByName(EMBL_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 478428, nextFeature.getFirstBase() );
		assertEquals( 482531, nextFeature.getLastBase() );
		
		// The very last feature
		
		idx = 247;
		nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "repeat_region", nextFeature.getKey().getKeyString() );
		assertEquals( "Pfalciparum_REP_45", nextFeature.getQualifierByName(EMBL_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 641897, nextFeature.getFirstBase() );
		assertEquals( 642097, nextFeature.getLastBase() );
		
		// Test indexOf and contains
		assertEquals(idx, entry.indexOf(nextFeature));
		assertTrue(entry.contains(nextFeature));
		
	}
	
	/**
	 * Exercise GenBank file loading/saving.
	 * 
	 * @throws Exception
	 */
	@Test
	public void testGenbankFile() throws Exception
	{
		// Given
		
		final File savefile = File.createTempFile(GBK_OUTPUT_FILE, "");
		final OutputStream outStream = new FileOutputStream(savefile);
		
		log("testGenbankFile", "Writing to file: " + savefile.getAbsolutePath());
		
		// When
		
		when( saveToDoc.getFile() ).thenReturn (savefile);
		when( saveToDoc.getWriter() ).thenReturn (new FileWriter(savefile));
		when( saveToDoc.getOutputStream() ).thenReturn (outStream);
		
		final Entry entry = loadFile(GBK_FILE);
		
		assertTrue(entry instanceof GenbankDocumentEntry);
		
		((GenbankDocumentEntry)entry).save(saveToDoc);
		
		
		// Then
		
		assertEquals( NUM_GBK_FEATURES, entry.getFeatureCount() );
		assertEquals( GBK_HDRS, entry.getHeaderText() );
		
		Sequence sequence = entry.getSequence();
		String seqString = String.valueOf( sequence.getCharSubSequence(1, sequence.length()) );
		
		assertTrue(sequence instanceof GenbankStreamSequence);
		assertEquals( GBK_SEQ_LEN, seqString.length() );
		assertEquals( 5340, sequence.getACount() );
		assertEquals( 4418, sequence.getCCount() );
		assertEquals( 2337, sequence.getGCount() );
		assertEquals( 4336, sequence.getTCount() );
		assertTrue(seqString.indexOf(GBK_SEQ_START) == 0);
		assertTrue(seqString.lastIndexOf(GBK_SEQ_END) == (seqString.length()-GBK_SEQ_END.length()));

		//printFeatures("testGenbankFile", entry);
		
		// Check the feature table
		assertNotNull(((SimpleDocumentEntry)entry).getFeatureTable());
		assertEquals(NUM_GBK_FEATURES, ((SimpleDocumentEntry)entry).getFeatureTable().getFeatureCount());
		
		int idx = 0;
		
		// Repeat region
		
		Feature nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "source", nextFeature.getKey().getKeyString() );
		assertEquals( "Enhydra lutris", nextFeature.getQualifierByName("organism").getValues().elementAt(0) );
		assertEquals( 1, nextFeature.getFirstBase() );
		assertEquals( 16431, nextFeature.getLastBase() );
				
		idx = 14;
		nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "tRNA", nextFeature.getKey().getKeyString() );
		assertEquals( "tRNA-Ala", nextFeature.getQualifierByName("product").getValues().elementAt(0) );
		assertEquals( 5032, nextFeature.getFirstBase() );
		assertEquals( 5100, nextFeature.getLastBase() );
		
		idx = 16;
		nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "rep_origin", nextFeature.getKey().getKeyString() );
		assertNull( nextFeature.getQualifierByName("product") );
		assertEquals( 5175, nextFeature.getFirstBase() );
		assertEquals( 5210, nextFeature.getLastBase() );
			
		idx = 19;
		nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "gene", nextFeature.getKey().getKeyString() );
		assertEquals( "COX1", nextFeature.getQualifierByName("gene").getValues().elementAt(0) );
		assertEquals( 5344, nextFeature.getFirstBase() );
		assertEquals( 6888, nextFeature.getLastBase() );
		
		nextFeature = entry.getFeatureAtIndex(++idx);
		assertEquals( "CDS", nextFeature.getKey().getKeyString() );
		assertEquals( "COX1", nextFeature.getQualifierByName("gene").getValues().elementAt(0) );
		assertEquals( 5344, nextFeature.getFirstBase() );
		assertEquals( 6888, nextFeature.getLastBase() );
		
		// Test indexOf and contains
		assertEquals(idx, entry.indexOf(nextFeature));
		assertTrue(entry.contains(nextFeature));
		
		idx = 49;
		nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "CDS", nextFeature.getKey().getKeyString() );
		assertEquals( "CYTB", nextFeature.getQualifierByName("gene").getValues().elementAt(0) );
		assertEquals( "1", nextFeature.getQualifierByName("codon_start").getValues().elementAt(0) );
		assertEquals( "2", nextFeature.getQualifierByName("transl_table").getValues().elementAt(0) );
		assertEquals( "cytochrome b", nextFeature.getQualifierByName("product").getValues().elementAt(0) );
		assertEquals( "YP_001382358.1", nextFeature.getQualifierByName("protein_id").getValues().elementAt(0) );
		assertEquals( "GeneID:5333285", nextFeature.getQualifierByName("db_xref").getValues().elementAt(0) );
		assertEquals( 14175, nextFeature.getFirstBase() );
		assertEquals( 15314, nextFeature.getLastBase() );
		
		// The very last feature
		
		idx = 52;
		nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "D-loop", nextFeature.getKey().getKeyString() );
		assertNull( nextFeature.getQualifierByName("product") );
		assertEquals( 15448, nextFeature.getFirstBase() );
		assertEquals( 16431, nextFeature.getLastBase() );
		
		// Test indexOf and contains
		assertEquals(idx, entry.indexOf(nextFeature));
		assertTrue(entry.contains(nextFeature));
		
	}
	
	/**
	 * Exercise non-indexed GFF file loading/saving.
	 * 
	 * @throws Exception
	 */
	@Test
	public void testNonIndexedGFF() throws Exception
	{
		// Given
		
		final File savefile = File.createTempFile(GFF_OUTPUT_FILE, "");
		final OutputStream outStream = new FileOutputStream(savefile);
		
		log("testNonIndexedGFF", "Writing to file: " + savefile.getAbsolutePath());
		
		// When
		
		when( saveToDoc.getFile() ).thenReturn (savefile);
		when( saveToDoc.getWriter() ).thenReturn (new FileWriter(savefile));
		when( saveToDoc.getOutputStream() ).thenReturn (outStream);
		
		final Entry entry = loadFile(GFF_FILE);
		
		assertTrue(entry instanceof GFFDocumentEntry);
		
		((GFFDocumentEntry)entry).save(saveToDoc);
		
		
		// Then
		
		assertEquals( NUM_GFF_FEATURES, entry.getFeatureCount() );
		assertEquals( GFF_HDRS, entry.getHeaderText() );
		
		Sequence sequence = entry.getSequence();
		String seqString = String.valueOf( sequence.getCharSubSequence(1, sequence.length()) );
		
		assertTrue(sequence instanceof FastaStreamSequence);
		assertEquals( GFF_SEQ_LEN, seqString.length() );
		assertEquals( 639816, sequence.getACount() );
		assertEquals( 159003, sequence.getCCount() );
		assertEquals( 159425, sequence.getGCount() );
		assertEquals( 629709, sequence.getTCount() );
		
		assertTrue(seqString.indexOf(GFF_SEQ_START) == 0);
		assertTrue(seqString.lastIndexOf(GFF_SEQ_END) == (seqString.length()-GFF_SEQ_END.length()));

		//printFeatures("testNonIndexedGFF", entry);
		
		// Check the feature table
		assertNotNull(((SimpleDocumentEntry)entry).getFeatureTable());
		assertEquals(NUM_GFF_FEATURES, ((SimpleDocumentEntry)entry).getFeatureTable().getFeatureCount());
		
		int idx = 0;
		
		// Repeat region
		
		Feature nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "repeat_region", nextFeature.getKey().getKeyString() );
		assertEquals( "Pfalciparum_REP_20", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 1, nextFeature.getFirstBase() );
		assertEquals( 360, nextFeature.getLastBase() );
				
		// Gene
		
		idx = 20;
		nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "gene", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0200300", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 35927, nextFeature.getFirstBase() );
		assertEquals( 37249, nextFeature.getLastBase() );
				
		nextFeature = entry.getFeatureAtIndex(++idx);
		assertEquals( "CDS", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0200300.1:exon{1}", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 35927, nextFeature.getFirstBase() );
		assertEquals( 37249, nextFeature.getLastBase() );
		
		nextFeature = entry.getFeatureAtIndex(++idx);
		assertEquals( "polypeptide", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0200300.1:pep", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 35927, nextFeature.getFirstBase() );
		assertEquals( 37249, nextFeature.getLastBase() );
		
		nextFeature = entry.getFeatureAtIndex(++idx);
		assertEquals( "mRNA", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0200300.1", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 35927, nextFeature.getFirstBase() );
		assertEquals( 37249, nextFeature.getLastBase() );
		
		// Next gene
		
		nextFeature = entry.getFeatureAtIndex(++idx);
		assertEquals( "gene", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0200400", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 38287, nextFeature.getFirstBase() );
		assertEquals( 39303, nextFeature.getLastBase() );
		
		nextFeature = entry.getFeatureAtIndex(++idx);
		assertEquals( "CDS", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0200400.1:exon{2,1}", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 38287, nextFeature.getFirstBase() );
		assertEquals( 39303, nextFeature.getLastBase() );
		
		nextFeature = entry.getFeatureAtIndex(++idx);
		assertEquals( "mRNA", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0200400.1", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 38287, nextFeature.getFirstBase() );
		assertEquals( 39303, nextFeature.getLastBase() );
		
		nextFeature = entry.getFeatureAtIndex(++idx);
		assertEquals( "polypeptide", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0200400.1:pep", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 38287, nextFeature.getFirstBase() );
		assertEquals( 39303, nextFeature.getLastBase() );
		
		// Test indexOf and contains
		assertEquals(idx, entry.indexOf(nextFeature));
		assertTrue(entry.contains(nextFeature));
		
		// Pseudogene
		
		idx = 1535;
		nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "pseudogene", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0222300", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 883452, nextFeature.getFirstBase() );
		assertEquals( 884312, nextFeature.getLastBase() );
		
		nextFeature = entry.getFeatureAtIndex(++idx);
		assertEquals( "polypeptide", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0222300.1:pep", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 883452, nextFeature.getFirstBase() );
		assertEquals( 884312, nextFeature.getLastBase() );
		
		nextFeature = entry.getFeatureAtIndex(++idx);
		assertEquals( "pseudogenic_transcript", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0222300.1", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 883452, nextFeature.getFirstBase() );
		assertEquals( 884312, nextFeature.getLastBase() );
		
		nextFeature = entry.getFeatureAtIndex(++idx);
		assertEquals( "pseudogenic_exon", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0222300.1:exon{1}", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 883452, nextFeature.getFirstBase() );
		assertEquals( 884312, nextFeature.getLastBase() );
		
		// The very last feature
		
		idx = 1583;
		nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "repeat_region", nextFeature.getKey().getKeyString() );
		assertEquals( "Pf3D7_02__new:946548..947102", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 946548, nextFeature.getFirstBase() );
		assertEquals( 947102, nextFeature.getLastBase() );
		
		// Test indexOf and contains
		assertEquals(idx, entry.indexOf(nextFeature));
		assertTrue(entry.contains(nextFeature));
		
	}
	
	/**
	 * Exercise indexed GFF file loading/saving.
	 * 
	 * @throws Exception
	 */
	@Test
	public void testIndexedGFF() throws Exception
	{
		// Given
		
		final File savefile = File.createTempFile(GFF_OUTPUT_FILE, "");
		final OutputStream outStream = new FileOutputStream(savefile);
		
		log("testIndexedGFF", "Writing to file: " + savefile.getAbsolutePath());
		
		// We don't want GUI pop-ups displayed during unit tests...
		DocumentEntryFactory.setDisplayIndexingQuestionForGffs(false);
		
		
		// When
		
		when( saveToDoc.getFile() ).thenReturn (savefile);
		when( saveToDoc.getWriter() ).thenReturn (new FileWriter(savefile));
		when( saveToDoc.getOutputStream() ).thenReturn (outStream);
		
		final Entry entry = loadFile(INDEXED_GFF_FILE);
		final Entry fastaEntry = loadFile(INDEXED_GFF_FASTA_FILE);
		
		final Bases bases = new Bases(fastaEntry.getSequence());
		final EntryGroup group = new SimpleEntryGroup(bases);
		final uk.ac.sanger.artemis.Entry artFastaEntry = new uk.ac.sanger.artemis.Entry(bases, fastaEntry);
		final uk.ac.sanger.artemis.Entry artGffEntry = new uk.ac.sanger.artemis.Entry(bases, entry);
		group.add(artFastaEntry);
		group.add(artGffEntry);
		
		((IndexedGFFDocumentEntry)entry).save(saveToDoc);
	
		
	    // Then
		
		assertTrue(entry instanceof IndexedGFFDocumentEntry);
		assertEquals( NUM_INDEXED_GFF_FEATURES, artGffEntry.getFeatureCount() );
		assertNull( artGffEntry.getHeaderText() ); // Method not implemented for GFF

		
		String seqString = artGffEntry.getBases().toString();
		
		assertEquals( GFF_SEQ_LEN, seqString.length() );
		assertTrue(seqString.indexOf(GFF_SEQ_START) == 0);
		assertTrue(seqString.lastIndexOf(GFF_SEQ_END) == (seqString.length()-GFF_SEQ_END.length()));

		//printFeatures("testIndexedGFF", artGffEntry);
		
		int idx = 0;
		
		// Repeat region
		
		uk.ac.sanger.artemis.Feature nextFeature = artGffEntry.getFeature(idx);
		assertEquals( "repeat_region", nextFeature.getKey().getKeyString() );
		assertEquals( "Pfalciparum_REP_20", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 1, nextFeature.getFirstBase() );
		assertEquals( 360, nextFeature.getLastBase() );
				
		// Gene
		
		idx = 12;
		nextFeature = artGffEntry.getFeature(idx);
		assertEquals( "gene", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0100200", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 1547747, nextFeature.getFirstBase() );
		assertEquals( 1548972, nextFeature.getLastBase() );
		
		nextFeature = artGffEntry.getFeature(++idx);
		assertEquals( "mRNA", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0100200.1", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 1547747, nextFeature.getFirstBase() );
		assertEquals( 1548972, nextFeature.getLastBase() );
		
		nextFeature = artGffEntry.getFeature(++idx);
		assertEquals( "polypeptide", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0100200.1:pep", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 1547747, nextFeature.getFirstBase() );
		assertEquals( 1548972, nextFeature.getLastBase() );
		
		nextFeature = artGffEntry.getFeature(++idx);
		assertEquals( "CDS", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0100200.1:exon:2", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 1547747, nextFeature.getFirstBase() );
		assertEquals( 1547800, nextFeature.getLastBase() );
		
		idx = 22;
		nextFeature = artGffEntry.getFeature(idx);
		assertEquals( "gene", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0100400", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 50363, nextFeature.getFirstBase() );
		assertEquals( 51636, nextFeature.getLastBase() );
		
		nextFeature = artGffEntry.getFeature(++idx);
		assertEquals( "mRNA", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0100400.1", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 50363, nextFeature.getFirstBase() );
		assertEquals( 51636, nextFeature.getLastBase() );
		
		nextFeature = artGffEntry.getFeature(++idx);
		assertEquals( "polypeptide", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0100400.1:pep", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 50363, nextFeature.getFirstBase() );
		assertEquals( 51636, nextFeature.getLastBase() );
		
		nextFeature = artGffEntry.getFeature(++idx);
		assertEquals( "CDS", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0100400.1:exon:2", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 50572, nextFeature.getFirstBase() );
		assertEquals( 51636, nextFeature.getLastBase() );
		
		// Pseudogene
		
		idx = 2069;
		nextFeature = artGffEntry.getFeature(idx);
		assertEquals( "pseudogene", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0220400", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 1461280, nextFeature.getFirstBase() );
		assertEquals( 1463392, nextFeature.getLastBase() );
		
		nextFeature = artGffEntry.getFeature(++idx);
		assertEquals( "pseudogenic_exon", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0220400.1:exon:1", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 1461280, nextFeature.getFirstBase() );
		assertEquals( 1461318, nextFeature.getLastBase() );
		
		nextFeature = artGffEntry.getFeature(++idx);
		assertEquals( "pseudogenic_transcript", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0220400.1", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 1461280, nextFeature.getFirstBase() );
		assertEquals( 1463392, nextFeature.getLastBase() );
		
		nextFeature = artGffEntry.getFeature(++idx);
		assertEquals( "pseudogenic_exon", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0220400.1:exon:2", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 1461420, nextFeature.getFirstBase() );
		assertEquals( 1461977, nextFeature.getLastBase() );
		
		nextFeature = artGffEntry.getFeature(++idx);
		assertEquals( "pseudogenic_exon", nextFeature.getKey().getKeyString() );
		assertEquals( "PF3D7_0220400.1:exon:3", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 1461980, nextFeature.getFirstBase() );
		assertEquals( 1463392, nextFeature.getLastBase() );
		
		// The very last feature
		
		idx = 2237;
		nextFeature = artGffEntry.getFeature(idx);
		assertEquals( "repeat_region", nextFeature.getKey().getKeyString() );
		assertEquals( "Pf3D7_02__new:946548..947102", nextFeature.getQualifierByName(GFF_FEATURE_ID_KEY).getValues().elementAt(0) );
		assertEquals( 1587399, nextFeature.getFirstBase() );
		assertEquals( 1587953, nextFeature.getLastBase() );
		
		// Test indexOf and contains
		Feature testFeature = entry.getFeatureAtIndex(10);
		assertEquals(10, entry.indexOf(testFeature));
		assertTrue(entry.contains(testFeature));
		
		testFeature = entry.getFeatureAtIndex(70);
		assertEquals(70, entry.indexOf(testFeature));
		assertTrue(entry.contains(testFeature));
		
	}
	
	/**
	 * Exercise Fasta file loading/saving.
	 * 
	 * @throws Exception
	 */
	@Test
	public void testNonIndexedFastaFile() throws Exception
	{
		// Given
		
		File savefile = File.createTempFile(FASTA_OUTPUT_FILE, "");
		final OutputStream outStream = new FileOutputStream(savefile);
		
		log("testNonIndexedFastaFile", "Writing to file: " + savefile.getAbsolutePath());
		
		// We don't want GUI pop-ups displayed during unit tests...
		DocumentEntryFactory.setDisplayIndexingQuestionForGffs(false);
		
		
		// When
		
		when( saveToDoc.getFile() ).thenReturn (savefile);
		when( saveToDoc.getWriter() ).thenReturn (new FileWriter(savefile));
		when( saveToDoc.getOutputStream() ).thenReturn (outStream);
		
		final Entry entry = loadFile(FASTA_FILE);
		
		assertTrue(entry instanceof EmblDocumentEntry);
			
		((EmblDocumentEntry)entry).save(saveToDoc);
			
		//printFeatures("testNonIndexedFastaFile", entry);
		
		// Then
		
		final Sequence sequence = entry.getSequence();

		assertTrue(sequence instanceof FastaStreamSequence);
		assertEquals(FASTA_SEQ_LEN, sequence.length());
		assertEquals(FASTA_TEXT, String.valueOf( sequence.getCharSubSequence(1, sequence.length()) ));
		// Check no. of bases
		assertEquals(38, sequence.getGCount());
		assertEquals(77, sequence.getACount());
		assertEquals(63, sequence.getTCount());
		assertEquals(62, sequence.getCCount());
		
		// Check Fasta headers
		String[] hdrs = ((FastaStreamSequence)sequence).getFastaHeaderStrings();
		assertNotNull(hdrs);
		assertEquals(2, hdrs.length);
		assertEquals("AM884176.1", hdrs[0]);
		assertEquals("AM884177.1", hdrs[1]);
		
		// Check the saved file is what we expect...
		assertTrue(doBasicFileComparison(savefile, getFile(FASTA_FILE), true));
		
	}
	
	/**
	 * Exercise GTF file loading/saving.
	 * 
	 * @throws Exception
	 */
	@Test
	public void testGTF() throws Exception
	{
		// Given
		
		final File savefile = File.createTempFile(GTF_OUTPUT_FILE, "");
		final OutputStream outStream = new FileOutputStream(savefile);
		
		log("testGTF", "Writing to file: " + savefile.getAbsolutePath());
		
		// When
		
		when( saveToDoc.getFile() ).thenReturn (savefile);
		when( saveToDoc.getWriter() ).thenReturn (new FileWriter(savefile));
		when( saveToDoc.getOutputStream() ).thenReturn (outStream);
		
		final Entry entry = loadFile(GTF_FILE);
		final Entry fastaEntry = loadFile(GTF_FASTA_FILE);
		
		final Bases bases = new Bases(fastaEntry.getSequence());
		final EntryGroup group = new SimpleEntryGroup(bases);
		final uk.ac.sanger.artemis.Entry artFastaEntry = new uk.ac.sanger.artemis.Entry(bases, fastaEntry);
		final uk.ac.sanger.artemis.Entry artGffEntry = new uk.ac.sanger.artemis.Entry(bases, entry);
		group.add(artFastaEntry);
		group.add(artGffEntry);
		
		((GFFDocumentEntry)entry).save(saveToDoc);
		
		
		// Then
		
		assertEquals( NUM_GTF_FEATURES, entry.getFeatureCount() );
		assertEquals( GTF_HDRS, entry.getHeaderText() );
		
		Sequence sequence = fastaEntry.getSequence();
		String seqString = String.valueOf( sequence.getCharSubSequence(1, sequence.length()) );
		
		assertTrue(sequence instanceof FastaStreamSequence);
		assertEquals( GTF_SEQ_LEN, seqString.length() );
		assertEquals( 620990, sequence.getACount() );
		assertEquals( 492667, sequence.getCCount() );
		assertEquals( 491650, sequence.getGCount() );
		assertEquals( 624271, sequence.getTCount() );
		
		assertTrue(seqString.indexOf(GTF_SEQ_START) == 0);
		assertTrue(seqString.lastIndexOf(GTF_SEQ_END) == (seqString.length()-GTF_SEQ_END.length()));

		// printFeatures("testGTF", entry);
		
		// Check the feature table
		assertNotNull(((SimpleDocumentEntry)entry).getFeatureTable());
		assertEquals(NUM_GTF_FEATURES, ((SimpleDocumentEntry)entry).getFeatureTable().getFeatureCount());
		
		int idx = 0;
		
		// Gene
		
		Feature nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "gene", nextFeature.getKey().getKeyString() );
		assertEquals( "CCDG5_0001", nextFeature.getQualifierByName("gene_id").getValues().elementAt(0) );
		assertEquals( "dnaA", nextFeature.getQualifierByName("gene_name").getValues().elementAt(0) );
		assertEquals( 1, nextFeature.getFirstBase() );
		assertEquals( 1332, nextFeature.getLastBase() );
				
		// Transcript
		
		idx = 14;
		nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "transcript", nextFeature.getKey().getKeyString() );
		assertEquals( "CCDG5_0003", nextFeature.getQualifierByName("gene_id").getValues().elementAt(0) );
		assertEquals( "CDZ23154", nextFeature.getQualifierByName("transcript_id").getValues().elementAt(0) );
		assertEquals( "ena", nextFeature.getQualifierByName("gene_source").getValues().elementAt(0) );
		assertEquals( "protein_coding", nextFeature.getQualifierByName("gene_biotype").getValues().elementAt(0) );
		assertNull( nextFeature.getQualifierByName("transcript_name") );
		assertEquals( "ena", nextFeature.getQualifierByName("transcript_source").getValues().elementAt(0) );
		assertEquals( "protein_coding", nextFeature.getQualifierByName("transcript_biotype").getValues().elementAt(0) );
		assertEquals( 2740, nextFeature.getFirstBase() );
		assertEquals( 2961, nextFeature.getLastBase() );
		
		// Exon
		
		idx = 64;
		nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "exon", nextFeature.getKey().getKeyString() );
		assertEquals( "CCDG5_0011", nextFeature.getQualifierByName("gene_id").getValues().elementAt(0) );
		assertEquals( "CDZ23162", nextFeature.getQualifierByName("transcript_id").getValues().elementAt(0) );
		assertEquals( "ena", nextFeature.getQualifierByName("gene_source").getValues().elementAt(0) );
		assertEquals( "protein_coding", nextFeature.getQualifierByName("gene_biotype").getValues().elementAt(0) );
		assertEquals( "ena", nextFeature.getQualifierByName("transcript_source").getValues().elementAt(0) );
		assertEquals( "protein_coding", nextFeature.getQualifierByName("transcript_biotype").getValues().elementAt(0) );
		assertEquals( "CDZ23162-1", nextFeature.getQualifierByName("exon_id").getValues().elementAt(0) );
		assertEquals( 12722, nextFeature.getFirstBase() );
		assertEquals( 13171, nextFeature.getLastBase() );
		
		// End features
		
		// CDS
		
		idx = 12313;
		nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "CDS", nextFeature.getKey().getKeyString() );
		assertEquals( "CCDG5_2088", nextFeature.getQualifierByName("gene_id").getValues().elementAt(0) );
		assertEquals( "CDZ25168", nextFeature.getQualifierByName("transcript_id").getValues().elementAt(0) );
		assertEquals( "ena", nextFeature.getQualifierByName("gene_source").getValues().elementAt(0) );
		assertEquals( "protein_coding", nextFeature.getQualifierByName("gene_biotype").getValues().elementAt(0) );
		assertEquals( "ena", nextFeature.getQualifierByName("transcript_source").getValues().elementAt(0) );
		assertEquals( "protein_coding", nextFeature.getQualifierByName("transcript_biotype").getValues().elementAt(0) );
		assertEquals( "1", nextFeature.getQualifierByName("exon_number").getValues().elementAt(0) );
		assertEquals( 2226882, nextFeature.getFirstBase() );
		assertEquals( 2228033, nextFeature.getLastBase() );
		
		// start_codon
		
		idx = 12314;
		nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "start_codon", nextFeature.getKey().getKeyString() );
		assertEquals( "CCDG5_2088", nextFeature.getQualifierByName("gene_id").getValues().elementAt(0) );
		assertEquals( "CDZ25168", nextFeature.getQualifierByName("transcript_id").getValues().elementAt(0) );
		assertEquals( "ena", nextFeature.getQualifierByName("gene_source").getValues().elementAt(0) );
		assertEquals( "protein_coding", nextFeature.getQualifierByName("gene_biotype").getValues().elementAt(0) );
		assertEquals( "ena", nextFeature.getQualifierByName("transcript_source").getValues().elementAt(0) );
		assertEquals( "protein_coding", nextFeature.getQualifierByName("transcript_biotype").getValues().elementAt(0) );
		assertEquals( "1", nextFeature.getQualifierByName("exon_number").getValues().elementAt(0) );
		assertEquals( 2228031, nextFeature.getFirstBase() );
		assertEquals( 2228033, nextFeature.getLastBase() );
		
		// stop_codon
		
		idx = 12310;
		nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "stop_codon", nextFeature.getKey().getKeyString() );
		assertEquals( "CCDG5_2088", nextFeature.getQualifierByName("gene_id").getValues().elementAt(0) );
		assertEquals( "CDZ25168", nextFeature.getQualifierByName("transcript_id").getValues().elementAt(0) );
		assertEquals( "ena", nextFeature.getQualifierByName("gene_source").getValues().elementAt(0) );
		assertEquals( "protein_coding", nextFeature.getQualifierByName("gene_biotype").getValues().elementAt(0) );
		assertEquals( "ena", nextFeature.getQualifierByName("transcript_source").getValues().elementAt(0) );
		assertEquals( "protein_coding", nextFeature.getQualifierByName("transcript_biotype").getValues().elementAt(0) );
		assertEquals( "1", nextFeature.getQualifierByName("exon_number").getValues().elementAt(0) );
		assertEquals( 2226879, nextFeature.getFirstBase() );
		assertEquals( 2226881, nextFeature.getLastBase() );
		
		// Test indexOf and contains
		assertEquals(idx, entry.indexOf(nextFeature));
		assertTrue(entry.contains(nextFeature));
		
	}
	
	/**
	 * Exercise Blast file loading/saving.
	 * 
	 * @throws Exception
	 */
	@Test
	public void testBlast() throws Exception
	{
		// Given
		
		final File savefile = File.createTempFile(BLAST_OUTPUT_FILE, "");
		final OutputStream outStream = new FileOutputStream(savefile);
		
		log("testBlast", "Writing to file: " + savefile.getAbsolutePath());
		
		
		// When
		
		when( saveToDoc.getFile() ).thenReturn (savefile);
		when( saveToDoc.getWriter() ).thenReturn (new FileWriter(savefile));
		when( saveToDoc.getOutputStream() ).thenReturn (outStream);
		
		final Entry entry = loadFile(BLAST_FILE);
		final Entry fastaEntry = loadFile(BLAST_FASTA_FILE);
		
		final Bases bases = new Bases(fastaEntry.getSequence());
		final EntryGroup group = new SimpleEntryGroup(bases);
		final uk.ac.sanger.artemis.Entry artFastaEntry = new uk.ac.sanger.artemis.Entry(bases, fastaEntry);
		final uk.ac.sanger.artemis.Entry artBlastEntry = new uk.ac.sanger.artemis.Entry(bases, entry);
		group.add(artFastaEntry);
		group.add(artBlastEntry);
		
		((BlastDocumentEntry)entry).save(saveToDoc);
		
		
		// Then
		
		assertEquals( NUM_BLAST_FEATURES, entry.getFeatureCount() );
		assertNull( entry.getHeaderText() );

		//printFeatures("testBlast", entry);
		
		// Check the feature table
		assertNotNull(((SimpleDocumentEntry)entry).getFeatureTable());
		assertEquals(NUM_BLAST_FEATURES, ((SimpleDocumentEntry)entry).getFeatureTable().getFeatureCount());
		
		// Blast features
	
		// First
		int idx = 0;
		Feature nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "BLASTCDS", nextFeature.getKey().getKeyString() );
		assertEquals( "15.9", nextFeature.getQualifierByName("blast_score").getValues().elementAt(0) );
		assertEquals( "85.714", nextFeature.getQualifierByName("score").getValues().elementAt(0) );
		assertEquals( "85.714", nextFeature.getQualifierByName("percent_id").getValues().elementAt(0) );
		assertEquals( "C7D1KZ100B4HIU", nextFeature.getQualifierByName("query_id").getValues().elementAt(0) );
		assertEquals( "3", nextFeature.getQualifierByName("subject_start").getValues().elementAt(0) );
		assertEquals( "16", nextFeature.getQualifierByName("subject_end").getValues().elementAt(0) );
		assertEquals( "200", nextFeature.getQualifierByName("subject_id").getValues().elementAt(0) );
		assertEquals( "hit to 200 3..16  score: 15.9  percent id: 85.714  e-value: 3.8", nextFeature.getQualifierByName("note").getValues().elementAt(0) );
		assertEquals( 3, nextFeature.getFirstBase() );
		assertEquals( 16, nextFeature.getLastBase() );

		// Last
		idx = NUM_BLAST_FEATURES-1;
		nextFeature = entry.getFeatureAtIndex(idx);
		assertEquals( "BLASTCDS", nextFeature.getKey().getKeyString() );
		assertEquals( "15.9", nextFeature.getQualifierByName("blast_score").getValues().elementAt(0) );
		assertEquals( "100.000", nextFeature.getQualifierByName("score").getValues().elementAt(0) );
		assertEquals( "100.000", nextFeature.getQualifierByName("percent_id").getValues().elementAt(0) );
		assertEquals( "C7D1KZ100BYF90", nextFeature.getQualifierByName("query_id").getValues().elementAt(0) );
		assertEquals( "181", nextFeature.getQualifierByName("subject_start").getValues().elementAt(0) );
		assertEquals( "174", nextFeature.getQualifierByName("subject_end").getValues().elementAt(0) );
		assertEquals( "200", nextFeature.getQualifierByName("subject_id").getValues().elementAt(0) );
		assertEquals( "hit to 200 181..174  score: 15.9  percent id: 100.000  e-value: 3.8", nextFeature.getQualifierByName("note").getValues().elementAt(0) );
		assertEquals( 174, nextFeature.getFirstBase() );
		assertEquals( 181, nextFeature.getLastBase() );
	}
}
