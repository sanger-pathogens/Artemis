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
package uk.ac.sanger.artemis;

import static org.junit.Assert.*;

import java.io.Writer;
import java.util.Calendar;
import java.util.Date;
import java.util.Map;
import java.util.HashMap;
import java.io.StringWriter;
import java.awt.Color;
import java.io.IOException;

import org.junit.Test;
import org.junit.Before;
import org.junit.After;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.Utils;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.sequence.AminoAcidSequence;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.Marker;
import uk.ac.sanger.artemis.sequence.Strand;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.util.OutOfRangeException;

/**
 * Unit test for the Feature class.
 * 
 * @author kp11
 *
 */
public class FeatureTest 
{
	private static final String TEST_FEATURE_ID = "PF3D7_0103500.1";
	private static final String TEST_FEATURE_ID2 = "PF3D7_0100400";
	private static final String POLYPEPTIDE_FEATURE_ID = "PF3D7_0221400.1";
	private static final String GENE_FEATURE_ID = "PF3D7_0108700";
	private static final String EXON_FEATURE_ID = "PF3D7_0100600.1";
	private static final String PSEUDOGENE_FEATURE_ID = "PF3D7_0100500";
	private static final String REPEAT_REGION_FEATURE_ID = "Pf3D7_02__new:946548..947102";
	
	private static final String COLOUR_QUALIFIER_NAME = "colour";
	private static final String PRODUCT_QUALIFIER_NAME = "product";
	private static final String PARENT_QUALIFIER_NAME = "Parent";
	private static final String START_CODON_QUALIFIER_NAME = "codon_start";
	private static final String PARTIAL_QUALIFIER_NAME = "partial";
	private static final String NOTE_QUALIFIER_NAME = "note";
	
	private static final String CDS_FEATURE_ID = "CDS";
	private static final String POLYPEPTIDE_ID = "polypeptide";
	private static final String PSEUDOGENE_ID = "pseudogene";
	private static final String GENE_ID = "gene";
	private static final String REPEAT_REGION_ID = "repeat_region";
	
	private static final String EXPECTED_NATIVE_FORMAT = 
			"Pf3D7_01_v3	chado	CDS	154176	156520	.	+	0	ID=PF3D7_0103500.1:exon:1;Parent=PF3D7_0103500.1;colour=10\n" + 
			"Pf3D7_01_v3	chado	CDS	156657	158966	.	+	0	ID=PF3D7_0103500.1:exon:2;Parent=PF3D7_0103500.1;colour=10\n" + 
			"Pf3D7_01_v3	chado	CDS	159141	159234	.	+	0	ID=PF3D7_0103500.1:exon:3;Parent=PF3D7_0103500.1;colour=10\n" + 
			"Pf3D7_01_v3	chado	CDS	159319	159431	.	+	0	ID=PF3D7_0103500.1:exon:4;Parent=PF3D7_0103500.1;colour=10\n" + 
			"Pf3D7_01_v3	chado	CDS	159556	159645	.	+	0	ID=PF3D7_0103500.1:exon:5;Parent=PF3D7_0103500.1;colour=10\n" + 
			"Pf3D7_01_v3	chado	CDS	159756	159936	.	+	0	ID=PF3D7_0103500.1:exon:6;Parent=PF3D7_0103500.1;colour=10\n" + 
			"Pf3D7_01_v3	chado	CDS	160088	160426	.	+	0	ID=PF3D7_0103500.1:exon:7;Parent=PF3D7_0103500.1;colour=10\n";
	
	private static final String EXPECTED_AMINO_ACIDS_PIRO_OUTPUT = 
			">BL;CDS, Pf3D7_01_02_v3.gff.gz 154176:160426 forward MW:219419\n" + 
			"MTFLCTMNSCEELNITNLIISLIKKNKVHDDLLCVKTLYLFIKTNKVPISLKLYNEKQNIFVSICTIICNKLKIYLYHIL\n" + 
			"CKTGDGEKKNLYLFIKNNLLLLNQDAVLDYDICNNIHKACVTIELFLKIIILCFKHINDNKHNDKHNDKHNDKHNDIHND\n" + 
			"KHNDIHNDKHNDKHNDIHNDIHNDIHNDIHKDVIFLACINNISMLLHILENRNNNMTNDKMCNEHFYLIIIKYIKYLFVH\n" + 
			"IYGNKNILNDKNNKQVLYSCLTCIIKNLLENTIKQNKKIKYKTLKLLYIILQKIQDIYIIDILFNRLSLYLFLLYKTCEP\n" + 
			"SANITNFILKIFLISSNQILLYYYKNNQILDQYKEVHESTHNYEKINVYYKKYKLFLQEHKKKAHKTVLTFDGKHNVDSK\n" + 
			"DSTKHICLHRNTENNQQNEKYDNHKIDDTHIIDNDNNSLTQQFNYVNKLNYSTQNLSNIYFFSYYILTNYDKNIKDIFPI\n" + 
			"CKVIIDHYPFLNKNLVLMSFIFILSEMFYDKNEHWLNHLNVLYGNGLSNYLQSCDISTHDVHSITNLLNNNNNNNSCDDN\n" + 
			"NINHHPYDQTYIPVKNVNIKSPQYDHQKDTNQNILSPLHNGLSFIFKHLIRMLKEDDKNNNHIINNIENYSVESLNYYIE\n" + 
			"NYYFKYMFNIFSFKKQENKYLLKYLKGYIYYRYIFNIIYPNLSPEHFTIFHNILDLYTVSNFNSIQNENLLVSSYTQKNI\n" + 
			"LSENNIIYSLDFFVQENITKDSDNLNKYVQPMNECNQENETASFSFEKFKNIKEAFKNENITLINQVSMFSLLFIERNII\n" + 
			"DNFLDDVIFAAWDNVEKNEKIKIHAQYKIDSSEHMNKGSHDNVTKRWWEQIYYYVDNNDTNKMIMKSEEKNEEKEKYILK\n" + 
			"YKYLHFLNYYLDSLILMTYINLYMNIRGINRDPITLRVYTNSDTENNNKSDGSNKSDGSNKSDGSNKSDGSNKSDGSNKS\n" + 
			"DGSNKSDGRNKSDDNNKSDGSNKSNKNYNRSNKCNSSCNKNNRNYFIKAMQHDIMKSILYEYYPMFNNYFSLAIIIKRIH\n" + 
			"DKINIKNILEKIMPFLSCNKVKKNYNKLTKENIQHNTNDDYMNYHTSLVMICLNKAFYLSYLCMYQKLIEKTFIQNFVYL\n" + 
			"LLKMNNACNYNKDLAKYTIINIYYYIKQINKFKNIKQDTQKNLCFNYNIYDNLHINEKEKYKNRLSQNQPFYIINNIQQH\n" + 
			"QIINIINYHHDVLSSYMYKKIINIDSIKNVTKILKLTQFLVMYNYYTYLYQDVCLNIIKYTKKTHFSYLSKEIKNQLYFL\n" + 
			"ILQMFNYILYLYYKYVNHNRIYICEKNQYFEQIKESILKGSIYVNLHHLSKKIKNQEQRDSSRYTNDEVVKNMNNDNTNK\n" + 
			"RVNNNNNNYHTIKPYIIQTTQKDKTFINQKDKTFIKQKDKTFINQKDKTFIKQKDHLPIFQHNEEEKKKLFLQNYIFFKL\n" + 
			"YKDIFNQNYTDIINNEKIKTEHEHKIKDLINNYINDLNKNINKTKISNLLIDDDFNDQPEEKKNNNLNHVQNLCDNNYFK\n" + 
			"HSNHNDFMSYADIRHTTSHIFHFTKGFLYNENLYLRFSAHLCILRCLYIFSTRLYELYPKIHQIWIYLKINFYKNNYMND\n" + 
			"ILVLKIINYIITIDDKYSVDRILNEIFPQVFDRIKTFETNKEICKQSYEYKFLQNTLLFFLNISRQERYFEKTHSEIFFF\n" + 
			"SLKCLNIIMNEEIKKISLNIICNIYLNNIPKIKKAIDSIICIKDKISFLYDDKKLYNTNEKYIMDNLFIKENVRDIEIMD\n" + 
			"VLNYLSCILNIINLKQIISLIIHIDIISLHFLVFYMSILNKQYKYKCRFNQHTLMVFDHFKDS\n" + 
			"*\n";
	
	private static final String EXPECTED_BASES_RESULT = 
			">PF3D7_0103500.1:exon{1,2,3,4,5,6,7}, 154176:160426 forward\n" + 
			"atgacttttctatgtacaatgaattcctgcgaagaattgaatatcactaatcttataatt\n" + 
			"tccttgattaaaaaaaataaggttcatgatgatttgttatgtgtaaaaactttatatctt\n" + 
			"tttataaaaacaaacaaagttcctatatccttaaaattatataatgagaaacaaaatatc\n" + 
			"tttgtaagcatatgcactattatatgtaataaattaaaaatatacttatatcatatatta\n" + 
			"tgcaaaaccggagatggagagaaaaaaaatttatatctattcataaaaaataatttatta\n" + 
			"ttattaaatcaagatgctgttctcgattatgatatatgtaataatatacacaaggcatgt\n" + 
			"gtaacgatagaattgtttttaaaaattataattttgtgttttaaacatattaatgataat\n" + 
			"aagcataatgataaacataatgataaacataatgataagcataatgatatacataatgat\n" + 
			"aaacataatgatatacataatgataaacataatgataaacataacgatatacataatgat\n" + 
			"atacataatgatatacataatgatatacataaggatgttatattccttgcatgcattaat\n" + 
			"aacatatcaatgttattacatatattagaaaatagaaataataatatgactaatgacaaa\n" + 
			"atgtgtaatgaacatttttatttaataataataaaatatataaaatatttatttgtacat\n" + 
			"atttatggtaacaaaaatatattgaatgataaaaataataaacaagttctttattcttgt\n" + 
			"ctaacatgtataattaaaaatttattagaaaataccataaaacaaaataaaaaaataaaa\n" + 
			"tataaaacattaaagttattatatataattcttcaaaaaattcaagacatatatataata\n" + 
			"gatatattatttaatagactctctttatatttatttcttctatataaaacgtgtgaacca\n" + 
			"tctgcaaatataacaaattttatactcaaaatatttcttatctcatcgaatcaaatactt\n" + 
			"ttatattattataaaaacaatcaaatattagatcaatataaagaagttcatgaaagtaca\n" + 
			"cataattatgagaaaataaatgtatattataaaaaatataaactcttcttacaagaacac\n" + 
			"aaaaagaaggcacacaagactgtattaacattcgatgggaagcataatgtagatagtaaa\n" + 
			"gatagtacaaagcatatatgtttacatagaaatacggaaaataatcaacaaaatgaaaag\n" + 
			"tatgataatcataaaatagatgatacacatattatagataacgataacaacagcttaact\n" + 
			"caacagtttaattatgtaaataaattaaattatagtacacaaaatttatcaaatatttat\n" + 
			"tttttttcttattatatactaaccaattatgataaaaatataaaagacatatttccaata\n" + 
			"tgtaaagttattatagaccattacccttttttaaataaaaatcttgtacttatgagtttc\n" + 
			"atttttattctttccgaaatgttttatgataaaaatgaacattggttaaatcacttgaat\n" + 
			"gtcttatatggaaatgggttatctaactatttacaaagttgtgatataagtacacacgac\n" + 
			"gtgcatagtataacaaatttgctaaataataataataataataatagttgtgatgataat\n" + 
			"aatattaatcatcatccttatgaccagacatatataccggtaaaaaatgtaaatataaaa\n" + 
			"tcaccacaatatgaccaccaaaaggatacgaaccaaaatatattatcacctcttcataat\n" + 
			"ggcttgtcttttatttttaaacatttaataagaatgctaaaagaagatgataaaaataat\n" + 
			"aatcatattattaataatattgaaaactattctgtggaatctttaaattactacattgaa\n" + 
			"aattattacttcaaatatatgtttaatattttttcatttaaaaaacaagaaaataaatat\n" + 
			"cttctcaaatatttaaaaggttatatatattatagatatatattcaacatcatatatcct\n" + 
			"aatttatctcctgaacattttactatcttccataatatattagatttatatacagtctcc\n" + 
			"aattttaattctatacaaaatgaaaatcttcttgtttcatcatatacacagaaaaatatt\n" + 
			"ttaagcgagaataatatcatatattctttagacttttttgttcaagaaaatattacaaaa\n" + 
			"gattctgataatcttaataaatatgtgcaacctatgaatgagtgtaaccaagaaaatgaa\n" + 
			"acggcttccttctcctttgaaaaatttaaaaacataaaggaggcttttaaaaatgaaaac\n" + 
			"atcacgcttattaaccaagtcagcatgttctctcttctttttatagaacgaaatataata\n" + 
			"gataactttctagacgatgttatttttgccgcctgggataatgttgaaaaaaatgaaaaa\n" + 
			"ataaaaattcatgcgcagtataaaatagattcatcagaacatatgaacaaaggatcacat\n" + 
			"gataatgtaaccaaaaggtggtgggaacagatttactactatgtagataataacgacaca\n" + 
			"aataaaatgataatgaaaagtgaagaaaaaaatgaagaaaaagaaaaatatatattaaaa\n" + 
			"tataaatatttgcatttccttaattattacttggactccttaattcttatgacatatatc\n" + 
			"aacctttacatgaatataagaggaattaatagggaccctataacgttaagggtatataca\n" + 
			"aattcggacacggaaaataataataaaagtgatggtagtaataaaagtgatggtagtaat\n" + 
			"aaaagtgatggtagtaataaaagtgatggtagtaataaaagtgatggtagtaataaaagt\n" + 
			"gatggtagtaataaaagtgatggtaggaataaaagtgatgataataataaaagtgatggt\n" + 
			"agtaataaaagtaataaaaactataatagaagcaacaaatgtaatagtagttgtaataaa\n" + 
			"aacaacagaaattattttatcaaagctatgcaacatgatataatgaaaagtatcttatat\n" + 
			"gaatattatcctatgtttaataattatttctccttagccataataatcaaaagaatacac\n" + 
			"gacaaaataaatattaaaaatatcctagaaaaaattatgccttttcttagttgtaataaa\n" + 
			"gtaaaaaaaaattataataaattaacaaaggaaaatatacaacataatacaaatgatgat\n" + 
			"tatatgaattatcatacatcactagtaatgatttgtttaaataaagccttctatttatca\n" + 
			"tatctatgtatgtaccaaaaactaatagaaaaaacctttatccaaaattttgtttattta\n" + 
			"ctcttaaaaatgaataatgcatgcaattataataaagatctagccaaatataccataata\n" + 
			"aatatatattattatataaaacaaataaataaatttaaaaatataaaacaggatacacaa\n" + 
			"aagaatctatgtttcaattataatatttatgacaacttacacataaatgaaaaagaaaaa\n" + 
			"tataaaaatagattatcacaaaatcaacctttttatattattaacaatatccaacaacat\n" + 
			"caaattattaatataattaattatcatcatgatgttttatcttcatatatgtataaaaaa\n" + 
			"attattaacatcgattcaataaaaaacgttacgaaaatattaaagttgacacaattttta\n" + 
			"gttatgtacaattattatacttatttatatcaagatgtatgtctcaatattattaaatat\n" + 
			"accaagaaaactcatttttcttatctaagtaaagaaattaaaaatcaactatatttcctc\n" + 
			"atcctacaaatgtttaattacatactttatttatattataaatatgtaaatcacaacaga\n" + 
			"atatatatatgtgagaaaaatcaatattttgaacaaataaaagaatccatccttaagggt\n" + 
			"agtatctacgtaaatcttcatcatctctcaaaaaaaataaaaaatcaagagcaaagggat\n" + 
			"tcatccagatatacaaatgatgaagttgtgaaaaatatgaataatgataatacaaacaag\n" + 
			"agagttaataataataataataattaccatactattaagccttacattattcagacaact\n" + 
			"caaaaggataaaacatttatcaatcaaaaggataaaacatttatcaaacaaaaggataaa\n" + 
			"acatttatcaatcaaaaggataaaacatttatcaaacaaaaggatcatttacccattttc\n" + 
			"caacataacgaagaagagaaaaaaaaattatttttacaaaactacatcttcttcaaatta\n" + 
			"tataaagatatattcaatcaaaattacacagatattataaataacgaaaagataaaaaca\n" + 
			"gaacatgaacataaaattaaagatcttataaataattacataaatgacctaaacaaaaat\n" + 
			"ataaataaaacaaaaatttcaaatctacttattgatgatgatttcaatgatcaaccagaa\n" + 
			"gaaaaaaaaaacaataacctaaatcatgttcaaaatttatgtgataataattattttaaa\n" + 
			"cattctaatcataacgatttcatgtcgtatgctgacattcgacatactacatctcacatt\n" + 
			"tttcattttacgaaaggttttctttacaacgaaaacttgtatttaagattcagtgctcat\n" + 
			"ctttgcatactgagatgcctgtatattttttctaccagactttatgagctgtatcctaag\n" + 
			"attcatcaaatatggatatatttgaagatcaatttttataaaaataattatatgaatgat\n" + 
			"atccttgttttgaaaataataaattatataataaccatagatgataaatattcagtggat\n" + 
			"agaattttaaatgaaatatttccacaagtattcgatagaataaaaacattcgaaacgaat\n" + 
			"aaagaaatatgtaaacagagttatgaatacaaatttttacagaacactttactatttttt\n" + 
			"ttaaatatttcaagacaggaacgatatttcgaaaaaactcattcggaaatatttttcttc\n" + 
			"tccttaaaatgccttaacataataatgaatgaagaaataaaaaagatatctttaaatata\n" + 
			"atatgtaatatatatcttaataacattcccaagataaagaaagctattgatagtattata\n" + 
			"tgtattaaggataagataagcttcttatacgatgataaaaagttatataacacaaacgaa\n" + 
			"aaatatattatggataatttgtttataaaagaaaatgtacgggatattgaaataatggat\n" + 
			"gtgttgaactatttatcatgtatattaaatataataaatctaaagcaaatcatttcctta\n" + 
			"ataattcatatagatattatttctcttcattttctagtcttttatatgtctatattaaat\n" + 
			"aagcaatataaatataaatgtagattcaatcaacacactcttatggtgtttgatcatttt\n" + 
			"aaggattcataa\n";
	
	private static final String EXPECTED_AMINO_ACIDS_FOR_FEATURE = 
			">PF3D7_0103500.1:exon{1,2,3,4,5,6,7} PF3D7_0103500.1:exon{1,2,3,4,5,6,7} undefined product 154176:160426 forward MW:219419\n" + 
			"MTFLCTMNSCEELNITNLIISLIKKNKVHDDLLCVKTLYLFIKTNKVPISLKLYNEKQNI\n" + 
			"FVSICTIICNKLKIYLYHILCKTGDGEKKNLYLFIKNNLLLLNQDAVLDYDICNNIHKAC\n" + 
			"VTIELFLKIIILCFKHINDNKHNDKHNDKHNDKHNDIHNDKHNDIHNDKHNDKHNDIHND\n" + 
			"IHNDIHNDIHKDVIFLACINNISMLLHILENRNNNMTNDKMCNEHFYLIIIKYIKYLFVH\n" + 
			"IYGNKNILNDKNNKQVLYSCLTCIIKNLLENTIKQNKKIKYKTLKLLYIILQKIQDIYII\n" + 
			"DILFNRLSLYLFLLYKTCEPSANITNFILKIFLISSNQILLYYYKNNQILDQYKEVHEST\n" + 
			"HNYEKINVYYKKYKLFLQEHKKKAHKTVLTFDGKHNVDSKDSTKHICLHRNTENNQQNEK\n" + 
			"YDNHKIDDTHIIDNDNNSLTQQFNYVNKLNYSTQNLSNIYFFSYYILTNYDKNIKDIFPI\n" + 
			"CKVIIDHYPFLNKNLVLMSFIFILSEMFYDKNEHWLNHLNVLYGNGLSNYLQSCDISTHD\n" + 
			"VHSITNLLNNNNNNNSCDDNNINHHPYDQTYIPVKNVNIKSPQYDHQKDTNQNILSPLHN\n" + 
			"GLSFIFKHLIRMLKEDDKNNNHIINNIENYSVESLNYYIENYYFKYMFNIFSFKKQENKY\n" + 
			"LLKYLKGYIYYRYIFNIIYPNLSPEHFTIFHNILDLYTVSNFNSIQNENLLVSSYTQKNI\n" + 
			"LSENNIIYSLDFFVQENITKDSDNLNKYVQPMNECNQENETASFSFEKFKNIKEAFKNEN\n" + 
			"ITLINQVSMFSLLFIERNIIDNFLDDVIFAAWDNVEKNEKIKIHAQYKIDSSEHMNKGSH\n" + 
			"DNVTKRWWEQIYYYVDNNDTNKMIMKSEEKNEEKEKYILKYKYLHFLNYYLDSLILMTYI\n" + 
			"NLYMNIRGINRDPITLRVYTNSDTENNNKSDGSNKSDGSNKSDGSNKSDGSNKSDGSNKS\n" + 
			"DGSNKSDGRNKSDDNNKSDGSNKSNKNYNRSNKCNSSCNKNNRNYFIKAMQHDIMKSILY\n" + 
			"EYYPMFNNYFSLAIIIKRIHDKINIKNILEKIMPFLSCNKVKKNYNKLTKENIQHNTNDD\n" + 
			"YMNYHTSLVMICLNKAFYLSYLCMYQKLIEKTFIQNFVYLLLKMNNACNYNKDLAKYTII\n" + 
			"NIYYYIKQINKFKNIKQDTQKNLCFNYNIYDNLHINEKEKYKNRLSQNQPFYIINNIQQH\n" + 
			"QIINIINYHHDVLSSYMYKKIINIDSIKNVTKILKLTQFLVMYNYYTYLYQDVCLNIIKY\n" + 
			"TKKTHFSYLSKEIKNQLYFLILQMFNYILYLYYKYVNHNRIYICEKNQYFEQIKESILKG\n" + 
			"SIYVNLHHLSKKIKNQEQRDSSRYTNDEVVKNMNNDNTNKRVNNNNNNYHTIKPYIIQTT\n" + 
			"QKDKTFINQKDKTFIKQKDKTFINQKDKTFIKQKDHLPIFQHNEEEKKKLFLQNYIFFKL\n" + 
			"YKDIFNQNYTDIINNEKIKTEHEHKIKDLINNYINDLNKNINKTKISNLLIDDDFNDQPE\n" + 
			"EKKNNNLNHVQNLCDNNYFKHSNHNDFMSYADIRHTTSHIFHFTKGFLYNENLYLRFSAH\n" + 
			"LCILRCLYIFSTRLYELYPKIHQIWIYLKINFYKNNYMNDILVLKIINYIITIDDKYSVD\n" + 
			"RILNEIFPQVFDRIKTFETNKEICKQSYEYKFLQNTLLFFLNISRQERYFEKTHSEIFFF\n" + 
			"SLKCLNIIMNEEIKKISLNIICNIYLNNIPKIKKAIDSIICIKDKISFLYDDKKLYNTNE\n" + 
			"KYIMDNLFIKENVRDIEIMDVLNYLSCILNIINLKQIISLIIHIDIISLHFLVFYMSILN\n" + 
			"KQYKYKCRFNQHTLMVFDHFKDS\n";
	
	// For gene PF3D7_0108700
	private static final String EXPECTED_GENE_BASES_1 = 
			"atgaaattttatagtacttttgttatttgttttattattttaaaagtatgtttaagtaaaaacataaatttaaacgatgaagggaaaaagaaacaaagcggcataataaatgaaaatgaggatggaaaaaataataaaagtaataataataaaaaagttcaacataataaaattaatcgtaaaggaaataatagtgtaactaaaaaaacggatgaacataaaaatgatggtgatgcggaaaataaaaaaaatggtgaacataaaaatgatggtgatgggaaaaataaaaaaaatggtgaacataaaaatgatggtgatgggaaaaataaaaaaatggatgaacatttgtatgatattgatggggaaaataaaaaaaatggtgaacatttgtatgatgttgatggggaaaataaaaaaaatggtgaacatttgtatgatattgatggggaaaataaaaaaaatggtgaacatttgtatgatattgatggggaaaataaaaaaacggatgaacatttgtatgatgttgatgggggaaaaaataatattcttgaatataataacattgatcaacataatgaatttccagaaagcttagaaaacgataattcttatgataaacttttcgatgatgttgaattaagggatattattcataacgagaaattttttgagaatttaaaaaatgtgaataataacgatgtacataattttttattagtcgataaagaaatggaaagaagaaaggaagaggagaaaaaaaagagtagaaaattagaagaggacgaggaagacgatgaagatgaagaggaagacaatgatgagaaatggaaagaagaaaataaaaataatagtagtagtagtagtagtagtagtagtaataattttaataaaaattatgatatttatgataatagagatctttataattatgacaatcttgtggatattacacaagaagaacatataaatgtccctaatataaacgatgaaattaatgaagggcgtaaacatgataaggaagattttttaatatttgagaaaagattaaatgaagctgaagatagaaatgattcgagtgatggggataataaaatggaagatttaataaataaccattttgacgattcaggtgataatggatcttatgaaaatattcaagtaaaaaaagaattatataataaaaatgaaaagattataaaagggagagataataatctgaatgatggagatattcatatgaataaaaataatacaagtaataaattgaataataataataataataataataatgacaataatatgaatagtgagatgaaatataaaacatttaaagaagagtacgaagaacaaatattaagtaaaccagaatcgataaaatatttttttgaagtaataacagaaatattaatagtaataagattaggtttaatatatagatatacgaattttttaataccttttaagaattttcttatttcaaaaactttagaaaatgtcgaaatggtatgtgttacatttttatctatacataaatttggtagagaaaaatatcctgatatatatggtttttttttgataatcctaattttgtatatattaaaatatgtgtttaaaaaaatgatgtataaattatattataataaatatggtatgggtaagaaatattctttgaaaaatcaagaaagtgaaaatatgtacatgcagaatttattaaaaaggattttgtataatgtagaaaaaaagaaggtcttatctaactatgaaagtatattacaagatatattagaaaacacagataacttattagtaaatagtaatataataaataaagaaaataaaaatatttatacagatatgatatcgaatattaataatataggagtctttacatatatatctacacaagcattaaaaaaaattaatcataaatcggatttaataattaatgaattaacaaccaatggttttattgatgaagatacaagtaacaaaatggataaactcaaaaacttcgcttccataaatgaatatccttatgtagatattaaaaagtttaaagataatctcgtcaacgaatatgatgagcactttatattaaacgataaaatgagatacagtggatttcatagtatggaaaatgacaagagtcctccgaatgataagtataacatggatagttataatatggatagttataatatggaaaacatgaacatggaaaacatgaacatggataatatgaacatggaaaacatgaacatggataatagcaacatggataatagcaacatggataatatgaacatggataatagcaacatggataatagccacgtggatatgaataacatgaataataataacatgaataataataacatgggtgtacatataaatatgaataacaattatttatctaataaatataaacacaattctgcatatgatgaaaataacgtgcaggactttttacagaaatcaaaaaagcaagattatggacatgtaacgaatacactgaaagatgggaagcaaatgatttataataacatggatccgaataataaaaaagatgataattattttaagaaaatgatgagagattcagacttatcaaattataacggaaggaaaagtatgtatttaaatgaatcattagaagatgacatgaaaaataatgattatgataataagattaataataaaaattttgattatgcaaataataatagcacatatattgtagaaggtaaagaaattcaacaaggtgaaaaaaataaaacaaatgatatggatgaaaaaaaaaaagatgaatcctttattcattatgacaataaaaaaataataaataatttaaaccccttctcccctttgaatgattataataatatttctaataataaaaagaataccataaacgaagaacaaaataataatgtaaatgaatctgatgttataacatacatagcaaattcttatataaataaaccacctagtaatataaatagatttagtcaaacttctatttataatacacataataaggtaaccaaccccttattaaatacacaaagaaatgtaccaaaagaagtaacacaacgtaaaaacacggaagcatctgtaaataataatgtcgaagaaaattatctgagaaacaaggttacgaaccatgatacggaacaaagcgtggagaacaaattggtaagagtacgatacatgtgaaaaataatttaaaatataaaaattgatatgtagattgtgtaaaaagggctatttctttgtgaaatatgactttaatataatacatatatatatgtatatatgtatatatgtatgtatgtatttatttatttatttatttatttatttatttatgtatgttttatattttgttaggatgaagacaagttaaagatgagtactcacatggtaaataatatatttgaaggtacaaaggaatatataaaaacaaatgaagatatgaataatataagtaaagcaggtgatccctttttgcataatcaagaaggtaaatatataaacaagttaatgatgagaagaaatataagcatactttatatgttataatatatatatatatatatatatatatatatatatatgcatatatttttatttatttatatatatgtgtgtgttaatttttgtagggtctttttcttatgcgccacccccttatcaaggtgtagagaactccagtaataaaaaccaaaagtatctgacccaacgaaagtaataaaaaaaataaatatacacaaataaataaataaataaatatatatatatatatatatatatatatatatataattattattatttttattcttttagggaacgacaacaaattgttgctactaaatctccatttaattaa";
	
	private static final String [][] TEST_GENE_AA_COUNTS = 
	{
		{"a","10"},
	    {"r","28"},
	    {"n","220"},
	    {"d","105"},
	    {"c","13"},
	    {"q","23"},
	    {"e","101"},
	    {"g","45"},
	    {"h","31"},
	    {"i","110"},
	    {"l","83"},
	    {"k","141"},
	    {"m","56"},
	    {"f","51"},
	    {"p","19"},
	    {"s","70"},
	    {"t","41"},
	    {"w","1"},
	    {"y","84"},
	    {"v","51"},
	    {"*","2"},
	    {"#","11"},
	    {"+","3"},
	    {".","0"},
	    {"u","0"}
	};

	private static final String [][] TEST_GENE_BASE_COUNTS = 
	{
		{"a","1810"},
	    {"t","1185"},
	    {"g","601"},
	    {"c","301"}
	};
	
	private static final Map<String, Integer> TEST_GENE_BASE_POSITIONAL_COUNTS = new HashMap<String, Integer>();
	static {
		TEST_GENE_BASE_POSITIONAL_COUNTS.put("a0",633);
		TEST_GENE_BASE_POSITIONAL_COUNTS.put("a1",719);
		TEST_GENE_BASE_POSITIONAL_COUNTS.put("a2",458);
		TEST_GENE_BASE_POSITIONAL_COUNTS.put("t0",261);
		TEST_GENE_BASE_POSITIONAL_COUNTS.put("t1",351);
		TEST_GENE_BASE_POSITIONAL_COUNTS.put("t2",573);
		TEST_GENE_BASE_POSITIONAL_COUNTS.put("g0",312);
		TEST_GENE_BASE_POSITIONAL_COUNTS.put("g1",130);
		TEST_GENE_BASE_POSITIONAL_COUNTS.put("g2",159);
		TEST_GENE_BASE_POSITIONAL_COUNTS.put("c0",93);
		TEST_GENE_BASE_POSITIONAL_COUNTS.put("c1",99);
		TEST_GENE_BASE_POSITIONAL_COUNTS.put("c2",109);
	};
	
	private EntryGroup egrp;
	private FeatureVector features;
	
	/*
	 * Utility method for displaying all qualifiers for the feature.
	 */
	private void displayQualifiers(Feature f)
	{
		for (Qualifier q : f.getQualifiers())
		{
			System.out.println("Qualifier: " + q.getName().toString() + " " + q.getValues().toString());
		}
	}
	
	private void displayRequiredQualifiers(Feature f, String key)
	{
		Entry entry = f.getEntry();
		EntryInformation entryInf = entry.getEntryInformation();
		StringVector vec = entryInf.getRequiredQualifiers(new Key(key));
		
		if (vec != null)
		{
			System.out.println("Required qualifiers: " + vec.toString());
		}
		else 
		{
			System.out.println("No required values found");
		}
	}
	
	private void displayBases(Feature f) 
	{
		System.out.println("Bases: " + f.getBases().toString());
	}
	
	private void displaySegments(Feature f)
	{
		FeatureSegmentVector segments = f.getSegments();
		
		if (segments != null) {
				
			for (int i = 0; i < segments.size(); i++)
			{
				FeatureSegment segment = segments.elementAt(i);
				System.out.println(segment.getStart() + " " + segment.getEnd() + " " + segment.getBases());
			}
		}
	}
	
	@Before
	public void setup() 
	{
		
		egrp = Utils.getEntryGroup("/data/Pf3D7_01_02_v3.gff.gz");

	    features = egrp.getAllFeatures();
	}

	@After
	public void cleanup() 
	{
		
		egrp = null;
		features = null;
	}
	
	@Test
	public void testGetEmblFeature() 
	{

	    final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
	    
	    assertNotNull(f.getEmblFeature());
		
	}
	
	@Test
	public void testWriteNativeEmbl() throws IOException
	{
		final StringWriter writer = new StringWriter();
		
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		f.writeNative(writer);
		
		assertNotNull(writer.getBuffer());
		
		String output = writer.getBuffer().toString();
		assertEquals("Check native format creation", EXPECTED_NATIVE_FORMAT, output);

	}

	@Test
	public void testWritePIROfFeatureEmbl() throws IOException
	{
		final StringWriter writer = new StringWriter();
		
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		f.writePIROfFeature(writer);
		
		assertNotNull(writer.getBuffer());
		
		String output = writer.getBuffer().toString();
		assertEquals("Check PIRO database record creation", EXPECTED_AMINO_ACIDS_PIRO_OUTPUT, output);

	}
	
	@Test
	public void testWriteBasesOfFeatureEmbl() throws IOException
	{
		
		final StringWriter writer = new StringWriter();
		
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		f.writeBasesOfFeature(writer);
		
		assertNotNull(writer.getBuffer());
		
		String output = writer.getBuffer().toString();
		assertEquals("Check bases for feature " + TEST_FEATURE_ID, EXPECTED_BASES_RESULT, output);
	
	}
	
	@Test
	public void testWriteAminoAcidsOfFeatureEmbl() throws IOException
	{
		final StringWriter writer = new StringWriter();
		
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		f.writeAminoAcidsOfFeature(writer);
		
		assertNotNull(writer.getBuffer());
		
		String output = writer.getBuffer().toString();
		assertEquals("Check amino acids for feature " + TEST_FEATURE_ID, EXPECTED_AMINO_ACIDS_FOR_FEATURE, output);
	}
	
	@Test
	public void testGetUpstreamBasesEmbl()
	{
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		final String upstreamBases = f.getUpstreamBases(154176);
		
		assertNotNull(upstreamBases);
		//System.err.println("KEV: " + upstreamBases);
	}
	
	@Test
	public void testGetDownstreamBasesEmbl()
	{
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		final String downstreamBases = f.getDownstreamBases(154176);
		
		assertNotNull(downstreamBases);
	}
	
	@Test
	public void testToStringEmbl()
	{
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		final String featureText = f.toString();
		
		assertNotNull(featureText);
		assertEquals("Check toString", EXPECTED_NATIVE_FORMAT, featureText);
		
	}
	
	@Test
	public void testIsForwardFeature()
	{
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		final String featureText = f.toString();
		
		assertNotNull(featureText);
		assertEquals("Check toString", EXPECTED_NATIVE_FORMAT, featureText);
	}
	
	@Test
	public void testGetKey()
	{
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		final Key key = f.getKey();
		
		assertNotNull(key);
		assertEquals("Check key", CDS_FEATURE_ID, key.getKeyString());
	}
	
	@Test
	public void testIsProteinFeature() throws Exception
	{
		Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		boolean isProtein = f.isProteinFeature();
		assertTrue("Is CDS a protein feature: " + f.getKey().toString(), isProtein);
		
		f = Utils.getFeatureByIdAndKey(POLYPEPTIDE_FEATURE_ID, POLYPEPTIDE_ID, features);
		isProtein = f.isProteinFeature();
		assertTrue("Is a Polypeptide a protein feature: " + f.getKey().toString(), isProtein);
		
		// // Cheat a bit - by creating an exon feature on the fly.......
		f = Utils.getFeatureByIdAndKey(TEST_FEATURE_ID, CDS_FEATURE_ID, features);
		f.set(new Key("exon"), f.getLocation(), f.getQualifiers());
		isProtein = f.isProteinFeature();
		assertTrue("Is an exon a protein feature: " + f.getKey().toString(), isProtein);
		
		f = Utils.getFeatureByIdAndKey(PSEUDOGENE_FEATURE_ID, PSEUDOGENE_ID, features);
		isProtein = f.isProteinFeature();
		assertFalse("Check unknown feature type: " + f.getKey().toString(), isProtein);
		
		// Missing a couple of test cases here.
	}
	
	@Test
	public void testIsCDS()
	{
		Feature f = Utils.getFeatureByIdAndKey(TEST_FEATURE_ID, CDS_FEATURE_ID, features);
		boolean isCds = f.isCDS();
		
		assertTrue("Should be a CDS feature", isCds);
		
		f = Utils.getFeatureByIdAndKey(POLYPEPTIDE_FEATURE_ID, POLYPEPTIDE_ID, features);
		isCds = f.isCDS();
		assertFalse("Shouldn't be a CDS feature: " + f.getKey().toString(), isCds);
	}
	
	@Test
	public void testGetLocation()
	{
		Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		Location loc = f.getLocation();
		
		assertNotNull(loc);
		assertEquals("CDS Location first base index", f.getFirstBase(),  loc.getFirstBase());
		assertEquals("CDS Location last base index", f.getLastBase(),  loc.getLastBase());
		
		
		f = Utils.getFeatureByIdAndKey(POLYPEPTIDE_FEATURE_ID, POLYPEPTIDE_ID, features);
		loc = f.getLocation();
		
		assertNotNull(loc);
		assertEquals("Polypeptide Location first base index", f.getFirstBase(),  loc.getFirstBase());
		assertEquals("Polypeptide Location last base index", f.getLastBase(),  loc.getLastBase());
	
	}
	
	@Test
	public void testGetQualifiers() 
	{
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		final QualifierVector qualVec = f.getQualifiers();
		
		
		// displayQualifiers(f);
		
		int i = 0;
		assertNotNull(qualVec);
		assertEquals("Check number of feature qualifiers", 4, qualVec.size());
		assertEquals("Check qualifier name", qualVec.elementAt(i).getName(), "Parent");
		assertEquals("Check qualifier name", qualVec.elementAt(++i).getName(), COLOUR_QUALIFIER_NAME);
		assertEquals("Check qualifier name", qualVec.elementAt(++i).getName(), "codon_start");
		assertEquals("Check qualifier name", qualVec.elementAt(++i).getName(), "ID");
	

	}
	
	@Test 
	public void testGetEntry()
	{
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		assertNotNull(f.getEntry());
		
	}
	
	@Test
	public void testGetCodonStart() throws Exception
	{
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID2, features);
		
		int codonStartIdx = f.getCodonStart();
		assertEquals("Check non-existent start codon", 1, codonStartIdx);
		
		// Add the specific score qualifier into the feature for purposes of testing
		f.addQualifierValues(new Qualifier(START_CODON_QUALIFIER_NAME, "1"));
		codonStartIdx = f.getCodonStart();
		assertEquals("Check start codon", 1, codonStartIdx);
		
	}
	
	@Test
	public void testGetScore() throws Exception
	{
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID2, features);
		
		int score = f.getScore();
		assertEquals("Check no score", -1, score);
		
		// Add the specific score qualifier into the feature for purposes of testing
		f.addQualifierValues(new Qualifier("score", "35"));
		score = f.getScore();
		assertEquals("Check score", 35, score);	
		f.removeQualifierByName("score");
		
		f.addQualifierValues(new Qualifier("score", "150"));
		score = f.getScore();
		assertEquals("Check score greater than 100", 100, score);
		f.removeQualifierByName("score");
		
		f.addQualifierValues(new Qualifier("score", "-20"));
		score = f.getScore();
		assertEquals("Check score less than 0", 0, score);
		f.removeQualifierByName("score");
	}
	
	@Test
	public void testSetEntry() 
	{
		final Feature f = Utils.getCDSFeatureByIdPrefix("PF3D7_0100400", features);
		
		Entry entry = f.getEntry();
		entry.setName("TestEntry");
		f.setEntry(entry);
		
		assertEquals("Check new entry has been set", "TestEntry", f.getEntry().getName());
	}
	
	@Test
	public void testSetV1() throws Exception
	{
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID2, features);
		
		QualifierVector vec = new QualifierVector();
		vec.add(new Qualifier("product", "a product"));
		
		f.set(
				new Key(CDS_FEATURE_ID),
                	new Location("60000..60089"),
                	vec
        );
		
		assertEquals("Check new key value", CDS_FEATURE_ID, f.getKey().getKeyString());
		assertEquals("Check new Location", 60000, f.getLocation().getFirstBase());
		assertEquals("Check new Location", 60089, f.getLocation().getLastBase());
		assertEquals("Check new QualifierVector size", 1, f.getQualifiers().size());
		assertEquals("Check new Qualifier name", "product", f.getQualifiers().get(0).getName());
		assertEquals("Check new Qualifier value", "a product", f.getQualifiers().get(0).getValues().get(0));
	}
	
	@Test
	public void testSetV2() throws Exception
	{
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID2, features);
		
		QualifierVector vec = new QualifierVector();
		vec.add(new Qualifier("product", "a product"));
		
		Calendar dt = Calendar.getInstance();
		dt.set(2017, 9, 4, 16, 05, 10);
		
		f.set(
				dt.getTime(),
				new Key(CDS_FEATURE_ID),
                	new Location("60000..60089"),
                	vec
        );
		
		assertEquals("Check new key value", CDS_FEATURE_ID, f.getKey().getKeyString());
		assertEquals("Check new Location", 60000, f.getLocation().getFirstBase());
		assertEquals("Check new Location", 60089, f.getLocation().getLastBase());
		assertEquals("Check new QualifierVector size", 1, f.getQualifiers().size());
		assertEquals("Check new Qualifier name", "product", f.getQualifiers().get(0).getName());
		assertEquals("Check new Qualifier value", "a product", f.getQualifiers().get(0).getValues().get(0));
	}
	
	@Test
	public void testAddQualifierValues() throws Exception
	{
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID2, features);
		int numStartingQualifiers = f.getQualifiers().size();
		
		
		Qualifier newQual1 = new Qualifier("new qualifier 1", "1");
		Qualifier qualifier = f.addQualifierValues(newQual1);
		
		assertEquals("new qualifier 1", qualifier.getName());
		assertEquals("1", qualifier.getValues().get(0));
		assertEquals("Check number of qualifiers has increased by 1", numStartingQualifiers+1, f.getQualifiers().size());
		assertEquals("new qualifier 1", f.getQualifiers().getQualifierByName("new qualifier 1").getName());
		assertEquals("1", f.getQualifiers().getQualifierByName("new qualifier 1").getValues().elementAt(0));
		
		assertEquals("Check number of qualifiers has changed for EmblFeature", numStartingQualifiers+1, f.getEmblFeature().getQualifiers().size());
		assertEquals("new qualifier 1", f.getEmblFeature().getQualifierByName("new qualifier 1").getName());
		assertEquals("1", f.getEmblFeature().getQualifierByName("new qualifier 1").getValues().elementAt(0));
		
		
		
		Qualifier newQual2 = new Qualifier("new qualifier 2", "2");
		qualifier = f.addQualifierValues(newQual2);
		
		assertEquals("new qualifier 2", qualifier.getName());
		assertEquals("2", qualifier.getValues().get(0));
		assertEquals("Check number of qualifiers has now increased by 2", numStartingQualifiers+2, f.getQualifiers().size());
		assertEquals("new qualifier 2", f.getQualifiers().getQualifierByName("new qualifier 2").getName());
		assertEquals("2", f.getQualifiers().getQualifierByName("new qualifier 2").getValues().elementAt(0));
		
		assertEquals("Check number of qualifiers has changed for EmblFeature", numStartingQualifiers+2, f.getEmblFeature().getQualifiers().size());
		assertEquals("new qualifier 2", f.getEmblFeature().getQualifierByName("new qualifier 2").getName());
		assertEquals("2", f.getEmblFeature().getQualifierByName("new qualifier 2").getValues().elementAt(0));
		
		
		
		// Clear the qualifiers and add a brand new one from scratch
		f.getQualifiers().clear();
		Qualifier newQual3 = new Qualifier("new qualifier 3", "3");
		qualifier = f.addQualifierValues(newQual3);
		
		assertEquals("new qualifier 3", qualifier.getName());
		assertEquals("3", qualifier.getValues().get(0));
		assertEquals("Check number of qualifiers is now just 1", 1, f.getQualifiers().size());
		assertEquals("new qualifier 3", f.getQualifiers().getQualifierByName("new qualifier 3").getName());
		assertEquals("3", f.getQualifiers().getQualifierByName("new qualifier 3").getValues().elementAt(0));
		
		assertEquals("Check number of qualifiers for EmblFeature", 1, f.getEmblFeature().getQualifiers().size());
		assertEquals("new qualifier 3", f.getEmblFeature().getQualifierByName("new qualifier 3").getName());
		assertEquals("3", f.getEmblFeature().getQualifierByName("new qualifier 3").getValues().elementAt(0));

	}
	
	@Test
	public void testSetQualifier() throws Exception
	{
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		int numStartingQualifiers = f.getQualifiers().size();
		
		Qualifier newQual1 = new Qualifier(COLOUR_QUALIFIER_NAME, "4");
		f.setQualifier(newQual1);
		assertEquals("Check number of qualifiers has not changed", numStartingQualifiers, f.getQualifiers().size());
		assertEquals(COLOUR_QUALIFIER_NAME, f.getQualifiers().getQualifierByName(COLOUR_QUALIFIER_NAME).getName());
		assertEquals("4", f.getQualifiers().getQualifierByName(COLOUR_QUALIFIER_NAME).getValues().elementAt(0));
		
		assertEquals("Check number of qualifiers has not changed", numStartingQualifiers, f.getEmblFeature().getQualifiers().size());
		assertEquals(COLOUR_QUALIFIER_NAME, f.getEmblFeature().getQualifierByName(COLOUR_QUALIFIER_NAME).getName());
		assertEquals("4", f.getEmblFeature().getQualifierByName(COLOUR_QUALIFIER_NAME).getValues().elementAt(0));
		
	}
	
	@Test
	public void testRemoveQualifierByName() throws Exception
	{
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		int numStartingQualifiers = f.getQualifiers().size();
		
		assertTrue("We must have atleast one qualifier for this test", (numStartingQualifiers > 1));
		assertNotNull("Check test qualifier exists", f.getQualifiers().getQualifierByName(COLOUR_QUALIFIER_NAME));
		f.removeQualifierByName(COLOUR_QUALIFIER_NAME);
		
		assertEquals("Check number of qualifiers has reduced by 1", numStartingQualifiers-1, f.getQualifiers().size());
		assertNull("Confirm qualifier was removed", f.getQualifiers().getQualifierByName(COLOUR_QUALIFIER_NAME));
		
		// Try to remove it again
		try 
		{
			f.removeQualifierByName(COLOUR_QUALIFIER_NAME);
		} 
		catch(Exception e) 
		{
			org.junit.Assert.fail("Removing a non-existent qualifier caused an exception: " + e.getMessage());
		}
	}
	
	@Test
	public void testContainsText() 
	{
		// Only basic tests
		
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		
		assertTrue(
			"Check qualifier text search", 
			f.containsText(
				"PF3D7_0103500.1", // Parent
                false,
                false,
                null)
		);
		
		assertTrue(
			"Check qualifier partial text search", 
			f.containsText(
				"pf3D7", 
                true, // partial
                true, // ignore case
                null)
		);
		
		assertFalse(
			"Check non-existent qualifier text",
			f.containsText(
				"dummy-qual",
                false,
                false,
                null)
		);
		
	}
	
	@Test 
	public void testFindOrReplaceText()
	{
		// Only basic tests
		
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		assertTrue(
			"Check qualifier text replacement", 
			f.findOrReplaceText(
				"PF3D7_0103500.1",
                false,
                false,
                false,
                null,
                "new-qualifier")
		);
		
		// Should no longer be found as we have changed its name to new-qualifier
		assertFalse(
			f.containsText(
				"PF3D7_0103500.1",
                false,
                false,
                null));
		
		assertTrue(f.containsText(
				"new-qualifier",
                false,
                false,
                null));
	}
	
	@Test
	public void testHasValidStartCodon() throws Exception
	{
		// TODO This needs further tests for AA translation, force flag etc
		
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		boolean hasStartCodon = f.hasValidStartCodon(false);
		assertTrue("Feature does not have a start codon qualifier", hasStartCodon);
		assertTrue("Feature is a CDS", f.isCDS());
		
		// Add one
		f.addQualifierValues(new Qualifier(START_CODON_QUALIFIER_NAME, "startcodon_value"));
		f.addQualifierValues(new Qualifier(PARTIAL_QUALIFIER_NAME, "partial_value"));
		
		hasStartCodon = f.hasValidStartCodon(false); 
		assertTrue("Feature now has a start codon qualifier", hasStartCodon);
		
	}
	
	@Test 
	public void testHasValidStopCodon() 
	{
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		final boolean hasStopCodon = f.hasValidStopCodon(false);
		
		// TODO Needs other test cases - minimal test
		assertTrue("Feature has valid stop codon", hasStopCodon);
	}
	
	@Test
	public void testHasValidEMBLKey() throws Exception
	{
		Feature f = Utils.getFeatureByIdAndKey(POLYPEPTIDE_FEATURE_ID, POLYPEPTIDE_ID, features);
		boolean valid = f.hasValidEMBLKey();
		
		assertTrue("Feature has valid EMBL Key", valid);
		
		f = Utils.getFeatureByIdAndKey(PSEUDOGENE_FEATURE_ID, PSEUDOGENE_ID, features);
		valid = f.hasValidEMBLKey();
		
		assertTrue("Feature has valid EMBL Key", valid);
		
		f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		valid = f.hasValidEMBLKey();
		
		assertTrue("Feature has valid EMBL Key", valid);
		
		// No way to test negative case as not allowed to set incorrect keys.
	}
	
	@Test
	public void testHasRequiredQualifiers() throws Exception
	{
		//Feature f = Utils.getFeatureByIdAndKey(POLYPEPTIDE_FEATURE_ID, POLYPEPTIDE_ID, features);
		/*egrp = Utils.getEntryGroup("/data/MAL1.embl.gz");
	    features = egrp.getAllFeatures();
	    
		Feature f = Utils.getFeatureByIdAndKey("PFA0410w", CDS_FEATURE_ID, features);
		
		assertNotNull(f);
		displayRequiredQualifiers(f, "PFA0410w");
		
		assertTrue("Feature has valid EMBL qualifiers", f.hasRequiredQualifiers());
		
		f.addQualifierValues(new Qualifier("dummy unknown qualifier", "rubbish value"));
		
		assertFalse("Feature has invalid EMBL qualifiers", f.hasRequiredQualifiers());*/
	}
	
	@Test
	public void testFixStopCodon() throws Exception
	{
		// This feature already has a stop codon
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);

		boolean fixed = f.fixStopCodon();
		assertTrue(fixed);
		
		// Additional test required - Adding of missing stop codon
	}
	
	@Test
	public void testTrimStart() throws Exception 
	{
		// This feature already has a start codon in the right position
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		
		//displayBases(f);
		boolean result = f.trimStart(true, false);
		//displayBases(f);
		
		assertTrue(result);	
		
		// Additional test required - trim start to position start codon at beginning of feature
	}
	
	@Test
	public void testGetLabel() throws Exception 
	{
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		
		String val = f.getLabel();
		assertNull("Check we have no label qualifier", val);
		
		// Now add one
		f.addQualifierValues(new Qualifier("label","labelval"));
		val = f.getLabel();
		
		assertNotNull("Now we have added a label qualifier check that we can get its value", val);
		assertEquals("labelval", val);
	
	}
	
	@Test
	public void testGetGeneName() throws Exception 
	{
		final Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		
		String val = f.getGeneName();
		assertNull("Check we have no gene name qualifier", val);
		
		// Now add one
		f.addQualifierValues(new Qualifier("gene","geneName"));
		val = f.getGeneName();
		
		assertNotNull("Now we have added a gene name qualifier check that we can get its value", val);
		assertEquals("geneName", val);
	}
	
	@Test
	public void testGetIDString() throws Exception 
	{
		final Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		
		String val = f.getIDString();
		assertNotNull("Check we have an ID qualifier- name", val);
		assertEquals("PSOP24", val);
		
		// Now remove it
		f.removeQualifierByName("Name");
		val = f.getIDString();
		
		assertNotNull("Check we have an ID qualifier - should fall back to ID string now name has gone", val);
		assertEquals(GENE_FEATURE_ID, val);
		
	}
	
	@Test
	public void testGetSystematicName()
	{
		final Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		
		String val = f.getSystematicName();
		assertNotNull(val);
		assertEquals(GENE_FEATURE_ID, val);

	}
	
	@Test
	public void testGetNote() throws Exception
	{
		final Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		
		String val = f.getNote();
		assertNull("There is no note present on this feature", val);
		
		// Now add one
		f.addQualifierValues(new Qualifier(NOTE_QUALIFIER_NAME, "some notes"));
		val = f.getNote();
		assertNotNull("There should be a note present now", val);
		assertEquals("some notes", val);

	}
	
	@Test
	public void testGetProductString() throws Exception
	{
		Feature f = Utils.getFeatureByIdAndKey(POLYPEPTIDE_FEATURE_ID, POLYPEPTIDE_ID, features);
		
		String val = f.getProductString();
		assertNotNull("There should be a product present", val);
		assertEquals("stevor;", val);
	
		
		// Try a feature with no product
		f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		val = f.getProductString();
		assertNull("There should not be a product for this feature", val);

	}
	
	@Test
	public void testGetBaseCount() 
	{
		
		final Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		int count =  f.getBaseCount();
		
		assertEquals(3898, count);
		
	}
	
	@Test
	public void testGetBases() 
	{
		final Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		String bases = f.getBases();
		
		assertNotNull(bases);
		assertEquals(EXPECTED_GENE_BASES_1, bases);

	}
	
	@Test
	public void testGetTranslationBasesLength() 
	{
		Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		int len = f.getTranslationBasesLength();
		
		assertEquals(3897, len);
		
		f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		len = f.getTranslationBasesLength();
		
		assertEquals(5469, len);
		
		f = Utils.getFeatureByIdAndKey(POLYPEPTIDE_FEATURE_ID, POLYPEPTIDE_ID, features);
		len = f.getTranslationBasesLength();
		
		assertEquals(1038, len);

	}
	
	@Test
	public void testGetTranslationBases() 
	{
		Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		String translationBases = f.getTranslationBases();
		
		assertNotNull(translationBases);
		assertEquals(3897, translationBases.length());
	}
	
	@Test
	public void testGetTranslation() 
	{
		Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		AminoAcidSequence seq = f.getTranslation();
		
		// These are nice clean features...
		
		assertNotNull(seq);
		assertEquals( (3897 / 3), seq.length() );
		
		f = Utils.getFeatureByIdAndKey(POLYPEPTIDE_FEATURE_ID, POLYPEPTIDE_ID, features);
		seq = f.getTranslation();
		assertNotNull(seq);
		assertEquals( (1038 / 3), seq.length() );
		
		f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		seq = f.getTranslation();
		assertNotNull(seq);
		assertEquals( (5469 / 3), seq.length() );
	}
	
	@Test
	public void testGetAACount()
	{
		Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		int aaCount = f.getAACount();
		
		assertEquals((3897 / 3), aaCount);
		
		f = Utils.getFeatureByIdAndKey(POLYPEPTIDE_FEATURE_ID, POLYPEPTIDE_ID, features);
		aaCount = f.getAACount();
		
		assertEquals((1038 / 3), aaCount);
		
		f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		aaCount = f.getAACount();
		
		assertEquals((5469 / 3), aaCount);
	}
	
	@Test 
	public void testGetMolecularWeight() 
	{
		Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		float molWght = f.getMolecularWeight();
		
		assertEquals(151950.19f, molWght, 0.0001f);
		
		f = Utils.getFeatureByIdAndKey(POLYPEPTIDE_FEATURE_ID, POLYPEPTIDE_ID, features);
		molWght = f.getMolecularWeight();
		
		assertEquals(36548.04f, molWght, 0.0001f);
		
		f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		molWght = f.getMolecularWeight();
		
		assertEquals(219419.3f, molWght, 0.0001f);
	}
	
	@Test 
	public void testGetCodonCount()
	{
		Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		int numCodons = f.getCodonCount(0,0,1); // ttc
		
		assertEquals(5, numCodons);
		
		f = Utils.getFeatureByIdAndKey(POLYPEPTIDE_FEATURE_ID, POLYPEPTIDE_ID, features);
		numCodons = f.getCodonCount(0,0,1); //ttc
		
		assertEquals(9, numCodons);
	}
	
	@Test 
	public void  testGetResidueCount() throws Exception
	{
		final Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		
		for (int i = 0; i < TEST_GENE_AA_COUNTS.length; i++) 
		{
			int idx = f.getTranslation().getSymbolIndex(TEST_GENE_AA_COUNTS[i][0].charAt(0));
			int residueCount = f.getResidueCount(idx);
			
			assertEquals(
					"Check expected residue count for each amino acid in turn: " + TEST_GENE_AA_COUNTS[i][0], 
					Integer.parseInt(TEST_GENE_AA_COUNTS[i][1]), 
					residueCount);
		}

	}
	
	@Test
	public void testGetPositionalBaseCount()
	{
		final Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		
		for (Map.Entry<String, Integer> entry : TEST_GENE_BASE_POSITIONAL_COUNTS.entrySet())
		{
			char baseLetter  = entry.getKey().substring(0,1).charAt(0);
			Integer codonIdx = Integer.parseInt(entry.getKey().substring(1,2));
			
			int count = f.getPositionalBaseCount(codonIdx, Bases.getIndexOfBase(baseLetter));
				
			assertEquals(
					"Check each base/position in turn: " + entry.getKey(),
					entry.getValue().intValue(), 
					count);
		}
	}
	
	@Test
	public void testGetBaseCountByLetter()
	{
		final Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		
		for (int i = 0; i < TEST_GENE_BASE_COUNTS.length; i++) 
		{
			int count = f.getBaseCount(Bases.getIndexOfBase(TEST_GENE_BASE_COUNTS[i][0].charAt(0)));
			
			assertEquals(
					"Check expected count of each base letter in turn: " + TEST_GENE_BASE_COUNTS[i][0], 
					Integer.parseInt(TEST_GENE_BASE_COUNTS[i][1]), 
					count);
		}

	}
	
	@Test
	public void testGet12CorrelationScore() 
	{
		final Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		double correlationScore = f.get12CorrelationScore();
		
		assertEquals(55.2518090839107d, correlationScore, 0.001d);
		
	}
	
	@Test
	public void testGetPercentGC()
	{
		final Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		double gcPercent = f.getPercentGC();
		
		assertEquals(23.140071831708568d, gcPercent, 0.001d);
	}
	
	@Test
	public void testGetFirstBase()
	{
		final Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		int idx = f.getFirstBase();
		
		assertEquals(1231743, idx);
	}
	
	@Test
	public void testGetLastBase()
	{
		final Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		int idx = f.getLastBase();
		
		assertEquals(1235640, idx);
	}

	@Test 
	public void testLessThan() 
	{
		final Feature geneFeature = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		final Feature pepFeature = Utils.getFeatureByIdAndKey(POLYPEPTIDE_FEATURE_ID, POLYPEPTIDE_ID, features);
		final boolean result = geneFeature.lessThan(pepFeature);
		
		assertTrue(geneFeature.getFirstBase() + " should be less than " + pepFeature.getFirstBase(), result);
	}
	
	@Test 
	public void testGreaterThan() 
	{
		final Feature geneFeature = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		final Feature pepFeature = Utils.getFeatureByIdAndKey(POLYPEPTIDE_FEATURE_ID, POLYPEPTIDE_ID, features);
		boolean result = pepFeature.greaterThan(geneFeature);
		
		assertTrue(pepFeature.getFirstBase() + " should be greater than " + geneFeature.getFirstBase(), result);
	}
	
	@Test 
	public void testGetFirstBaseMarker() 
	{
		final Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		final Marker marker = f.getFirstBaseMarker();
		assertNotNull(marker);
		assertEquals(1231743, marker.getPosition());
	}
	
	@Test 
	public void testGetLastBaseMarker() 
	{
		final Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		final Marker marker = f.getLastBaseMarker();
		assertNotNull(marker);
		assertEquals(1235640, marker.getPosition());
	}
	
	@Test 
	public void testGetFirstCodingBaseMarker() 
	{
		final Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		final Marker marker = f.getFirstCodingBaseMarker();
		assertNotNull(marker);
		assertEquals(1231743, marker.getPosition());
	}
	
	@Test
	public void testGetPositionInSequence () throws Exception
	{
		Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		Marker marker = f.getPositionInSequence(10); // With respect to first base marker!!
		assertNotNull(marker);
		assertEquals(1231752, marker.getPosition());
		
		marker = f.getPositionInSequence(200); // With respect to first base marker!!
		assertNotNull(marker);
		assertEquals(1231942, marker.getPosition());
		
		try 
		{
			marker = f.getPositionInSequence(-1);
			org.junit.Assert.fail("Should have thrown an exception");
		}
		catch (OutOfRangeException e)
		{
			// Expected
		}
		
	}
	
	@Test
	public void testGetFeaturePositionFromMarker() throws Exception
	{
		Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		Marker marker = f.getStrand().makeMarker(1231752);
		int pos = f.getFeaturePositionFromMarker(marker);
		assertEquals(9, pos);
		
		marker = f.getStrand().makeMarker(1231942);
		pos = f.getFeaturePositionFromMarker(marker);
		assertEquals(199, pos);
		
		// out of range
		marker = f.getStrand().makeMarker(1231);
		pos = f.getFeaturePositionFromMarker(marker);
		assertEquals(-1, f.getFeaturePositionFromMarker(marker));
		
	}
	
	@Test
	public void testGetColour() throws Exception
	{
		// Gene feature with no colour qualifier
		Feature f = Utils.getFeatureByIdAndKey(GENE_FEATURE_ID, GENE_ID, features);
		Color colour = f.getColour();
		assertNull(colour);
		
		// CDS with colour qualifier
		f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		colour = f.getColour();
		assertNotNull(colour);
		assertEquals(
				"This feature has an attached colour qualifier", 
				Options.getOptions().getColorFromColourNumber(10).getRGB(), 
				colour.getRGB());
		
		// No default for this feature so must be taken from properties
		f = Utils.getFeatureByIdAndKey(REPEAT_REGION_FEATURE_ID, REPEAT_REGION_ID, features);
		colour = f.getColour();
		assertNotNull(colour);
		assertEquals(
				"This feature has no colour but should get the default from options", 
				Options.getOptions().getDefaultFeatureColour(new Key(REPEAT_REGION_ID)).getRGB(), 
				colour.getRGB());
		
		// Reset the value read from file to a garbage value. It should then fall back to default
		f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		f.removeQualifierByName(COLOUR_QUALIFIER_NAME);
		f.resetColour();
		f.addQualifierValues(new Qualifier(COLOUR_QUALIFIER_NAME,"ABC"));
		colour = f.getColour();
		assertNotNull(colour);
		assertEquals(
				"This feature has an invalid colour qualifier", 
				Options.getOptions().getColorFromColourNumber(5).getRGB(), 
				colour.getRGB());
	}
	
	@Test
	public void testGetValuesOfQualifier() throws Exception
	{
		Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);

		StringVector quals = f.getValuesOfQualifier(COLOUR_QUALIFIER_NAME);
		assertNotNull(quals);
		assertTrue(quals.size() == 1);
		assertNotNull(quals.get(0), "10");
	}
	
	@Test
	public void testGetQualifierByName() throws Exception
	{
		Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		Qualifier qual = f.getQualifierByName(COLOUR_QUALIFIER_NAME);
		
		assertNotNull(qual);
		assertEquals(qual.getName(), COLOUR_QUALIFIER_NAME);
		assertEquals(qual.getValues().size(), 1);
		assertEquals(qual.getValues().get(0), "10");
		
		f = Utils.getFeatureByIdAndKey(REPEAT_REGION_FEATURE_ID, REPEAT_REGION_ID, features);
		qual = f.getQualifierByName(COLOUR_QUALIFIER_NAME);
		assertNull(qual);
	}
	
	@Test
	public void testGetValueOfQualifier() throws Exception
	{
		Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		String colour = f.getValueOfQualifier(COLOUR_QUALIFIER_NAME);
		assertNotNull(colour);
		assertEquals(colour, "10");
		
		f = Utils.getFeatureByIdAndKey(REPEAT_REGION_FEATURE_ID, REPEAT_REGION_ID, features);
		colour = f.getValueOfQualifier(COLOUR_QUALIFIER_NAME);
		assertNull(colour);
	}
	
	@Test
	public void testRemoveFromEntry() throws Exception
	{
		Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		f.removeFromEntry();
		assertNull(f.getEntry());
	}
	
	@Test
	public void testGetSegments() 
	{
		Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		FeatureSegmentVector segments = f.getSegments();
		assertEquals(7, segments.size());
		
		f = Utils.getFeatureByIdAndKey(REPEAT_REGION_FEATURE_ID, REPEAT_REGION_ID, features);
		segments = f.getSegments();
		assertEquals(1, segments.size());
		assertEquals(555, segments.elementAt(0).getBases().length());
		
	}
	
	@Test
	public void testGetStrand()
	{
		Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		Strand strand = f.getStrand();
		
		assertNotNull(strand);
	}
	
	@Test
	public void testForcedDuplicate()  throws Exception
	{
		Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		int startingNumFeaturesInEntry = f.getEntry().getFeatureCount();
		Feature duplicate = f.duplicate(true);
		
		assertNotNull(duplicate);
		
		assertEquals("Check ID", f.getKey().getKeyString(), duplicate.getKey().getKeyString());
		assertEquals("Check Number of qualifiers", f.getQualifiers().size(), duplicate.getQualifiers().size());
		assertEquals("Check bases", f.getBases(), duplicate.getBases());
		assertEquals("Check amino acids", f.getTranslation().toString(), duplicate.getTranslation().toString());
		assertEquals("Check number of amino acids", f.getAACount(), duplicate.getAACount());
		assertEquals("Check number of bases", f.getBaseCount(), duplicate.getBaseCount());
		
		// Forced duplication adds a DUP prefix to the ID.
		assertEquals("Check ID string", "DUP1-"+f.getIDString(), duplicate.getIDString());
		assertEquals("Check molecular weight", f.getMolecularWeight(), duplicate.getMolecularWeight(), 0.001f);
		assertTrue("Check number of features has increased", (duplicate.getEntry().getAllFeatures().size() == startingNumFeaturesInEntry+1) );
	}
	
	@Test
	public void testUnforcedDuplicate()  throws Exception
	{
		Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		int startingNumFeaturesInEntry = f.getEntry().getFeatureCount();
		Feature duplicate = f.duplicate(false);
		
		assertNotNull(duplicate);
		
		assertEquals("Check ID", f.getKey().getKeyString(), duplicate.getKey().getKeyString());
		assertEquals("Check Number of qualifiers", f.getQualifiers().size(), duplicate.getQualifiers().size());
		assertEquals("Check bases", f.getBases(), duplicate.getBases());
		assertEquals("Check amino acids", f.getTranslation().toString(), duplicate.getTranslation().toString());
		assertEquals("Check number of amino acids", f.getAACount(), duplicate.getAACount());
		assertEquals("Check number of bases", f.getBaseCount(), duplicate.getBaseCount());
		
		// Forced duplication adds a DUP prefix to the ID.
		assertEquals("Check ID string", f.getIDString(), duplicate.getIDString());
		assertEquals("Check molecular weight", f.getMolecularWeight(), duplicate.getMolecularWeight(), 0.001f);
		assertTrue("Check number of features has increased", (duplicate.getEntry().getAllFeatures().size() == startingNumFeaturesInEntry+1) );
	}
	
	@Test
	public void testMoveTo()  throws Exception
	{
		
		// TODO Need to create more test data for this test.
		
		/*System.out.println("NUM ENTRIES: " + egrp.size());
		
		Entry entry = new Entry(null);
		Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		f.moveTo(entry,  true);
		
		assertTrue(entry.getFeatureCount()==1);*/
		
		
	}
	
	@Test
	public void testCopyTo()  throws Exception
	{
		// TODO Complete this test
	}
	
	@Test
	public void testIsReadOnly() 
	{
		Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		assertFalse(f.isReadOnly());
	}
	
	@Test 
	public void testSetLocation() throws Exception
	{
		Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		Location loc = f.getLocation();
		Range range = new Range(10,30);
		Location newLoc = new Location(range);
		f.setLocation(newLoc);
		
		assertTrue(f.getLocation().getFirstBase()==10);
		assertTrue(f.getLocation().getLastBase()==30);
	}
	
	@Test
	public void testAddSegment() throws Exception
	{
		Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		int numSegments = f.getSegments().size();
		f.addSegment(new Range(10,30));
        
		assertTrue("Number of segments should have increased", f.getSegments().size()==numSegments+1);

	}
	
	@Test
	public void testRemoveSegment() throws Exception
	{
		Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		int numSegments = f.getSegments().size();
		assertTrue(numSegments > 0);
		
		FeatureSegmentVector vec = f.getSegments();
		FeatureSegment segment = vec.elementAt(0);
		f.removeSegment(segment);
        
		assertTrue("Number of segments should have decreased", f.getSegments().size()==numSegments-1);
	}
	
	@Test
	public void testGetAllQualifierNames()
	{
		Feature f = Utils.getCDSFeatureByIdPrefix(TEST_FEATURE_ID, features);
		StringVector qualifiers = f.getAllQualifierNames(egrp.getAllFeatures());
		
		assertEquals("[Derives_from, ID, Name, Parent, codon_start, colour, label, note, product, synonym]", qualifiers.toString().trim());
	}
	
}
