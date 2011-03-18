
package uk.ac.sanger.artemis.components.variant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.InputStream;
import java.io.StringWriter;
import java.net.URL;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import javax.swing.JFrame;
import javax.swing.JPanel;


import junit.framework.Assert;

import org.junit.BeforeClass;
import org.junit.Before;
import org.junit.Test;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.components.EntryFileDialog;
import uk.ac.sanger.artemis.components.variant.VCFview;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.sequence.MarkerRange;
import uk.ac.sanger.artemis.sequence.Strand;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.DocumentFactory;

/**
 * Tests for writing out sequences based on VCF data.
 */
public class WriteVCFTest
{
  private static VCFview vcfView;
  private static VCFview bcfView;
  
  // TEST WRITING REGIONS
  
  /**
   * Write FASTA, from a range selection on the forward strand.
   * Variation is a SNP.
   * VCF line : test    120768  .       C       T
   */
  @Test
  public void testRegionWriteFwdFastaSNP()
  {
    StringWriter writer = getRegionWriter(true, 120768, 120777);
    
    StringBuffer buff = new StringBuffer("> test.embl.gz 120768:120777\n");
    buff.append("cttgtcaagg\n");
    buff.append(">test.vcf.gz test.embl.gz 120768:120777\n");
    buff.append("tttgtcaagg\n");
    assertEquals("Export FASTA range ", writer.toString(), buff.toString());
  }
  
  /**
   * Write FASTA, from a range selection on the reverse strand.
   * Variation is a SNP.
   * VCF line : test    120768  .       C       T
   */
  @Test
  public void testRegionWriteRevFastaSNP()
  {
    StringWriter writer = getRegionWriter(false, 120768, 120777);

    StringBuffer buff = new StringBuffer("> test.embl.gz 120768:120777 reverse\n");
    buff.append("ccttgacaag\n");
    buff.append(">test.vcf.gz test.embl.gz 120768:120777 reverse\n");
    buff.append("ccttgacaaa\n");
    assertEquals("Export FASTA range ", writer.toString(), buff.toString());
  }
  
  /**
   * Write FASTA, from a range selection on the forward strand.
   * Variation is a deletion.
   * test    396838  .       tt      t
   */
  @Test
  public void testRegionWriteFwdFastaDeletion()
  {
    StringWriter writer = getRegionWriter(true, 396835, 396845);

    StringBuffer buff = new StringBuffer("> test.embl.gz 396835:396845\n");
    buff.append("tttttaggtat\n");
    buff.append(">test.vcf.gz test.embl.gz 396835:396845\n");
    buff.append("tttt-aggtat\n");
    assertEquals("Export FASTA range ", writer.toString(), buff.toString());
  }
  
  /**
   * Write FASTA, from a range selection on the reverse strand.
   * Variation is a deletion.
   * test    396838  .       tt      t
   */
  @Test
  public void testRegionWriteRevFastaDeletion()
  {
    StringWriter writer = getRegionWriter(false, 396835, 396845);

    StringBuffer buff = new StringBuffer("> test.embl.gz 396835:396845 reverse\n");
    buff.append("atacctaaaaa\n");
    buff.append(">test.vcf.gz test.embl.gz 396835:396845 reverse\n");
    buff.append("ataccta-aaa\n");
    assertEquals("Export FASTA range ", writer.toString(), buff.toString());
  }
  
  /**
   * Write FASTA, from a range selection on the forward strand.
   * Variation is a insertion.
   * test    366787  .       t       tT
   */
  @Test
  public void testRegionWriteFwdFastaInsertion()
  {
    StringWriter writer = getRegionWriter(true, 366785, 366795);

    StringBuffer buff = new StringBuffer("> test.embl.gz 366785:366795\n");
    buff.append("tttcgcttttt\n");
    buff.append(">test.vcf.gz test.embl.gz 366785:366795\n");
    buff.append("tttTcgcttttt\n");
    assertEquals("Export FASTA range ", writer.toString(), buff.toString());
  }
  
  /**
   * Write FASTA, from a range selection on the reverse strand.
   * Variation is a insertion.
   * test    366787  .       t       tT
   */
  @Test
  public void testRegionWriteRevFastaInsertion()
  {
    StringWriter writer = getRegionWriter(false, 366785, 366795);

    StringBuffer buff = new StringBuffer("> test.embl.gz 366785:366795 reverse\n");
    buff.append("aaaaagcgaaa\n");
    buff.append(">test.vcf.gz test.embl.gz 366785:366795 reverse\n");
    buff.append("aaaaagcgaaaa\n");
    assertEquals("Export FASTA range ", writer.toString(), buff.toString());
  }
  
  
  /**
   * Write FASTA, from a range selection on the forward strand.
   * Variation is a multiple allele, (PL number indicates either G or T):
   * test    361978  .       G       T       77.76842        .       DP=11;AF1=0.95;CI95=0.5,1;DP4=1,0,10,0;MQ=58;PV4=1,9.6e-06,0.37,1       PL:DP:SP:GT:GQ  106,0,0:11:0:1/1:10
   */
  @Test
  public void testRegionWriteFwdFastaMultipleAllele()
  {
    StringWriter writer = getRegionWriter(true, 361975, 361985);

    StringBuffer buff = new StringBuffer("> test.embl.gz 361975:361985\n");
    buff.append("actgaaaaatt\n");
    buff.append(">test.vcf.gz test.embl.gz 361975:361985\n");
    buff.append("actkaaaaatt\n");
    assertEquals("Export FASTA range ", writer.toString(), buff.toString());
  }
  
  // TEST WRITING FEATURES
  /**
   * Write FASTA, from a feature selection on the forward strand.
   */
  @Test
  public void testFeatureWriteFastaSNP()
  {
    StringWriter writer = getFeatureWriter("SPN23F03630,SPN23F03800,SPN23F04290", vcfView);

    StringBuffer fastaBuff = new StringBuffer(">test.embl.gz\n");
    StringBuffer basesBuff = new StringBuffer();
    //SPN23F03630
    
    basesBuff.append("atgaagaaaactgtttataaaaaattgggtatttcaattattgcgagtactttattggct");
    basesBuff.append("agccagttatcgacagtatctgctttgagtgttatttctagtacaggtgaagaatatgag");
    basesBuff.append("gtaagtgagacactagaaaaaggtccagagtctaatgattcttcattatctgagatttca");
    basesBuff.append("ccaacgtatggttcatactaccaaaagcaatcagaagtattatcggtaatgatgatttga");
    
    //SPN23F03800
    basesBuff.append("atggcaagtatcacactcacaccaagcgaaaaggatattcaggcttttcttgaacactat");
    basesBuff.append("caaaccagtctggctcctagcaagaatccctatatccgctactttttgaaactacctcaa");
    basesBuff.append("gcaacggtttctatctatacttctggaaaaatcttgcttcagggtgaaggggctgaaaaa");
    basesBuff.append("tacgccagtttctttggctatcaagctgtagagcaaaccagcggacaaaatcttccttta");
    basesBuff.append("attgggacagatgaggtgggaaatggttcctactttggtgggcttgcagttgtggctgcc");
    basesBuff.append("tttgtcacacctgaccagcacgactttttacgaaaactcggtgtgggggattctaagact");
    basesBuff.append("ctgaccgaccaaaagatccgtcagattgctcctattctcaaggaaaaaattcagcaccag");
    basesBuff.append("gcactccttctctcacccagcaagtacaacgaggtcatcggagaccgctacaatgctgtt");
    basesBuff.append("tcggttaaggttgccctccataatcaggctatctatctcctccttcaaaaaggtgttcag");
    basesBuff.append("cctgagaaaattgtgattgatgcctttaccagtgctaaaaattatgacaagtacttggca");
    basesBuff.append("caagagaccaatcgtttcagcaatcctatcagcttagaagaaaaggctgagggcaaatac");
    basesBuff.append("ttggctgtcgcagtttcttctgtcattgcgcgtgatctctttctggaaaatcttgaaaat");
    basesBuff.append("ttgggacgagaactgggttatcagcttccaagtggagctggaacggcttctgacaaggtg");
    basesBuff.append("gctagccagattttgcaagcctatggtatgcagggactcaacttctgcgccaaattgcac");
    basesBuff.append("tttaaaaatactgaaaaagcgaaaaacgcttag");
    
    //SPN23F04290
    basesBuff.append("atgctttatgtgggcattgatatcgctaaaaataaacacgatgttacagccttgaatgtt");
    basesBuff.append("ccaggaaaaactgttcttaaaccactcactttttcaaataataaagctggttttgaactc");
    basesBuff.append("ttagatctgtctcttcgacagctcaaccaagactgtctcatcgctcttaaacttctttct");
    basesBuff.append("gatcccaatcgtgaacaatttcaacacgataatcggcaagtagacctaaaaatactggct");
    basesBuff.append("agacatattcatcgtctcaagaaaaaacagtctgattggaaagtacaatacactcgttgt");
    basesBuff.append("cttgatatcatctttcctgagttggataaaatcgttggaaagcattcagaatatacctac");
    basesBuff.append("caactcttgacgcgctaccctaatcctcagaaaaggattgaggcaggatttgataagctg");
    basesBuff.append("atagaaattaagcgattgaccgcttctaaaattcaggatatcctctcagtcgcacctcgt");
    basesBuff.append("tctatcgaaacaacatctcctgctcgtgaattcgaaatcatcgaaatcatcaaacattac");
    basesBuff.append("aagaggctcattgacaaggcggaaacatgtgtcaatgacttgatggctgagttcaactca");
    basesBuff.append("gtcatcacgacggttactgggattgggggtcgtttaggggcggtcattttagccgagatt");
    basesBuff.append("cgaaatattcatgcctttgataatcctgctcaattacaagctttcgctggactggattct");
    basesBuff.append("tctatttatcagtcaggtcagattgatttagctggaagaatgatcaaacggggttcccct");
    basesBuff.append("catctgcggtgggcactcatacaagctgccaaagcatgcgctcgcttttcacctgctttt");
    basesBuff.append("aaggcctatcttaagactaagttagaacaaggaaaacattacaatgtagccatcatccac");
    basesBuff.append("cttgcaaaaaaacttatccgaaccctgttttatatcctaaaaaagagctgccatttgacg");
    basesBuff.append("aacaaaaagtga");
    
    IOUtils.wrapString(basesBuff.toString(), fastaBuff);

    fastaBuff.append(">test.vcf.gz \n");
    
    basesBuff = new StringBuffer();
    //SPN23F03630
    basesBuff.append("atgaagaaaactatttataaaaaattgggtatttcaattattgcgagtactttattggct");
    basesBuff.append("agccagttatcgacagtatctgctttgagtgttatttctagtacaggtgaagaatatgag");
    basesBuff.append("gtaagtgagacaytagaaaaaggtccagagtctaatrattcttcattatctgagatttca");
    basesBuff.append("ccaacgtatggttcatactaccaaaagcaatcagaagtattatcggtaatgatgatttga");

    //SPN23F03800
    basesBuff.append("atggcaagtatcacactcacaccaagcgaaaaggatattcaggcttttcttgaacactat");
    basesBuff.append("caaaccagtctggctcctagcaagaatccctatatccgctactttttgaaactacctcaa");
    basesBuff.append("gcaacggtttctatctatacttctggaaaaatcttgcttcagggtgaaggggctgaaaaa");
    basesBuff.append("tacgccagtttctttggctatcaagctgtagagcaaaccagcggacaaaatcttccttta");
    basesBuff.append("attgggacagatgaggtgggaaatggttcctactttggtgggcttgcagttgtggctgcc");
    basesBuff.append("tttgtcacacctgaccagcacgactttttacgaaaactcggtgtgggggattctaagact");
    basesBuff.append("ctgaccgaccaaaagatccgtcagattgctcctattctcaaggaaaaaatccagcaccag");
    basesBuff.append("gcactccttctctcacccagcaagtacaacgaggtcatcggagaccgctacaatgctgtt");
    basesBuff.append("tcggttaaggttgccctccataatcaggctatctatctcctccttcaaaaaggtgttcag");
    basesBuff.append("cctgagaaaattgtgattgatgcctttaccagtgctaaaaattatgacaagtacttggca");
    basesBuff.append("caagaggccaatcgtttcagcaatcctatcagcttagaaaaaaaggctgagggcaaatac");
    basesBuff.append("ttggctgtcgcagtttcttctgtcattgcgcgtgatctctttctggaaaatcttgaaaat");
    basesBuff.append("ctgggacgagaactgggttatcagctcccaagtggagctggaacagcttctgacaaggtg");
    basesBuff.append("gctagccagattttgcaagcctatggtatgcagggactcagcttctgcgccaaattgcac");
    basesBuff.append("tttaaaaacactgaaaaagcgaaaaaaacgcttag");

    //SPN23F04290
    basesBuff.append("atgctttatgtgggcattgatatcgctaaaaataaacacgatgttacagccttgaatgtt");
    basesBuff.append("ccaggaaaaactgttcttaaaccactcactttttcaaataataaagctggttttgaactc");
    basesBuff.append("ttagatctgtctcttcgacagctcaaccaagactgtctcatcgctcttaaacttctttct");
    basesBuff.append("gaccccaatcgtgaacaatttcaacacgataatcggcaagtagaactaaaaatactggct");
    basesBuff.append("agacatattcatcgtctcaagaaaaaacagtctgattggaaagtacaatacactcgttgt");
    basesBuff.append("cttgatatcatctttcctgagttggataaaatcgttgaaaagcattcagaatatacctac");
    basesBuff.append("caactcttgacgcgctaccctaatcctcagaaaaggcttgaggcaggatttgataagctg");
    basesBuff.append("atagaaattaagcgattgaccgcttctaaaattcaggatatcctctcagttgcacctcgt");
    basesBuff.append("tctatcgmaacaacatctcctgctcgtgaattcgaaatcatcgaaatcatcaaacattac");
    basesBuff.append("aagaggctcattgacaaggcggaaacatgtgtcaatgacttgatggctgagttcaactcg");
    basesBuff.append("gtcatcacgacggtcactgggattgggaatcgtttagaggcggtcattttagccgagatt");
    basesBuff.append("cgaaatattcatgcctttgataatcctgctcaattacaagctttcgctggactggattct");
    basesBuff.append("tctatttatcagtcaggtcagattgatttagctggaagaatgatcaaacggggttcccct");
    basesBuff.append("catctgcggtgggcactcatacaagctgccaaagcatgccctcgcttttcacctgctttt");
    basesBuff.append("aaggcctatcttaagactaagttagaataaggaaaacattacaatgtagccatcatccac");
    basesBuff.append("cttgcaaaaaaacttatccgaaccctgttttatatcctaaaaaagagctgccatttgacg");
    basesBuff.append("aacaaaaagtga");
    
    IOUtils.wrapString(basesBuff.toString(), fastaBuff);

    assertEquals("Export FASTA feature ", writer.toString(), fastaBuff.toString());
  }
  
  
  /**
   * Write FASTA, from a feature selection on the forward strand.
   */
  @Test
  public void testBCFFeatureWriteFastaSNP()
  {
    StringWriter writer = getFeatureWriter("PFA0140c,PFA0475c", bcfView);

    StringBuffer fastaBuff = new StringBuffer(">MAL1.embl.gz\n");
    StringBuffer basesBuff = new StringBuffer();
    //PFA0140c
    basesBuff.append("atgagcaatataaatgataataatattcaaaacagcgatgtaaaggaaataaaaaatgat");
    basesBuff.append("gaatatatgggagaagggtattataatacacataatgaatgtgagggttttgaacaaaac");
    basesBuff.append("aattatatgaataataatgaaatcaatatggtagatggaaaacatgataatagtgatgat");
    basesBuff.append("aatataaataattccatgggagggataacacatatgggacatgtgaatgagtataaaaaa");
    basesBuff.append("ataaaatcattttgtaagttaagaaggaaatttgtatatataaggagaacgggttataat");
    basesBuff.append("tcatacatatataatagaaagaaaaaaaaaaataaaaataaaatagagactgaaaggaat");
    basesBuff.append("gtagatgaaatgtataatatgaataataacaataataataatcataacaacaacaacaat");
    basesBuff.append("aataataataataataataataataataacagatataataataatattaattgtttggat");
    basesBuff.append("aattttagcgatttatttttaaaatttataaaacaattcattgatcatataaatgtacat");
    basesBuff.append("ttgaataatgtaaaaaattattatcgaaattatgtgtgtaatgataatatgaaatatgta");
    basesBuff.append("ttttatttatttattatatatatgttgttaaatattatgattatctttttttcgatgttt");
    basesBuff.append("gtatatttctttgtatatttttattttatcccacaaaataagtatgtatttcctatagat");
    basesBuff.append("ttttcattggtaaaaaatcctatagaagattatttacaaaacgaaaaaagatataacatt");
    basesBuff.append("atgaataaaatgaataatatggacaatatggataatatggataatatgaatcatataaat");
    basesBuff.append("tatatgaataacatgaatcataataataataataataataataattgtaatataaagaat");
    basesBuff.append("aattgtaataatattaatatgatgaaaaattttgataagtctaatttcaatgagtcagaa");
    basesBuff.append("gatttctatataaaatatttcgattctttaaaagaaagcaggaaaagaaatgattattta");
    basesBuff.append("tataataataaattaaaaggtttttataaagattttcaaaataatattttaagaggtcat");
    basesBuff.append("atactatttcaaaataaattatattattcagaaaaagaagagttttataagaatttcaat");
    basesBuff.append("ccttttaattttatatttttttttaagaaaagaacaggaaaagataatttacttaatata");
    basesBuff.append("aagaaaggatataaaattgatatatttttaaatttttcttatatgaataatgaatataat");
    basesBuff.append("gacaaatttaattttattcaattagttactgaaatttttaataaaaataataatgtaatg");
    basesBuff.append("tttagaaatgaaaaattatatatgaataataataattatgaatttataaagaaattacat");
    basesBuff.append("ttatttttaaatgctcctttttatttttttaatatgtttaataataaaacaaaagaaatt");
    basesBuff.append("cgattagtacatgattataaatatacgacgaatttttataaaattcaaatatatatatac");
    basesBuff.append("ccaccgatacaaatacataaagcttctatcgttattcttgtctacgttaattttatatat");
    basesBuff.append("tattacatgtataattatccttttgtgttcttttatatttttgtttttttcttgagcttc");
    basesBuff.append("atcttgatatttttcaacactatccttttttttatattgagttgttattattattttaaa");
    basesBuff.append("aacggattataa");
    
    // PFA0475c
    basesBuff.append("atggaaggacagcataaagaaaaaaatacaaaaataaaaaaaagtaaaaaaccgttagtt");
    basesBuff.append("gttagtaatagaaaaccattcaacgtttttgaaaaaaaaaaaagtgcaaagccgctcgta");
    basesBuff.append("agagatcctcgtttttctgatttttcaggatcgttcaatgcgaacttttttagaaatgca");
    basesBuff.append("tacaagtttttatacgactcaagagaacaagaaaagaagattatagaaaaaaaattaaaa");
    basesBuff.append("agcaaaaatataacacaagaagaaaaggacgaattaaaaaaaaaatacaatgattataaa");
    basesBuff.append("agcactgatatattattaaaaaaaaaagaagaagaaaggaaattaaaagcagaattagta");
    basesBuff.append("aaacaggaaaaacaaaatattctaaccaaaaataaaaaaccttattattatagtgataga");
    basesBuff.append("aaaattaaaaaaatagttcaagaaaaattgtcgagttataagagtttaaagaaagttata");
    basesBuff.append("aaaaaagaaaaaaaaacattgcaaaaagagagaaagagaaacatcaaaccaacaaaaaaa");
    basesBuff.append("aaaatttttttaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa");
    basesBuff.append("aaaaaaaaaaaaaaaacataa");

    IOUtils.wrapString(basesBuff.toString(), fastaBuff);

    fastaBuff.append(">MAL1_8_16_24h.raw.bcf \n");
    
    basesBuff = new StringBuffer();
    //PFA0140c
    basesBuff.append("atgagcaatataaatgataataatattcaaaacagcgatgtaaaggaaataaaaaatgat");
    basesBuff.append("gaatatatgggagaagggtattataatacacataatgaatgtgagggttttgaacaaaac");
    basesBuff.append("aattatatgaataataatgaaatcaatatggtagatggaaaacatgataatagtgatgat");
    basesBuff.append("aatataaataattccatgggagggataacacatatgggacatgtgaatgagtataaaaaa");
    basesBuff.append("ataaaatcattttgtaagttaagaaggaaatttgtatatataaggagaacgggttataat");
    basesBuff.append("tcatacatatataatagaaagaaaaaaaaaaataaaaataaaatagagactgaaaggaat");
    basesBuff.append("gtagatgaaatgtataatatgaataataacaataataataatcataacaacaacaacaat");
    basesBuff.append("aataataataataataataataataataacagatataataataatattaattgtttggat");
    basesBuff.append("aattttagcgatttatttttaaaatttataaaacaattcattgatcatataaatgtacat");
    basesBuff.append("ttgaataatgtaaaaaattattatcgaaattatgtgtgtaatgataatatgaaatatgta");
    basesBuff.append("ttttatttatttattatatatatgttgttaaatattatgattatctttttttcgatgttt");
    basesBuff.append("gtatatttctttgtatatttttattttatcccacaaaataagtatgtatttcctatagat");
    basesBuff.append("ttttcattggtaaaaaatcctatagaagattatttacaaaacgaaaaaagatataacatt");
    basesBuff.append("atgaataaaatgaataatatggacaatatggataatatggataatatgaatcatataaat");
    basesBuff.append("tatatgaataacatgaatcataataataataataataataataattgtaatataaagaat");
    basesBuff.append("aattgtaataatattaatatgatgaaaaattttgataagtctaatttcaatgagtcagaa");
    basesBuff.append("gatttctatataaaatatttcgattctttaaaagaaagcaggaaaagaaatgattattta");
    basesBuff.append("tataataataaattaaaaggtttttataaagattttcaaaataatattttaagaggtcat");
    basesBuff.append("atactatttcaaaataaattatattattcagaaaaagaagagttttataagaatttcaat");
    basesBuff.append("ccttttaattttatatttttttttaagaaaagaacaggaaaagataatttacttaatata");
    basesBuff.append("aagaaaggatataaaattgatatatttttaaatttttcttatatgaataatgaatataat");
    basesBuff.append("gacaaatttaattttattcaattagttactgaaatttttaataaaaataataatgtaatg");
    basesBuff.append("tttagaaatgaaaaattatatatgaataataataattatgaatttataaagaaattacat");
    basesBuff.append("ttatttttaaatgctcctttttatttttttaatatgtttaataataacacaaaagaaatt");
    basesBuff.append("cgattagtacatgattataaatatacgacgaatttttataaaattcaaatatatatatac");
    basesBuff.append("ccaccgatacaaatacataaagcttctatcgttattcttgtctacgttaattttatatat");
    basesBuff.append("tattacatgtataattatccttttgtgttcttttatatttttgtttttttcttgagcttc");
    basesBuff.append("atcttgatatttttcaacactatccttttttttatattgagttgttattattattttaaa");
    basesBuff.append("aacggattataa");

    // PFA0475c
    basesBuff.append("atggaaggacagcataaagaaaaaaatacaaaaataaaaaaaagtaaaaaaccgttagtt");
    basesBuff.append("gttagtaatagaaaaccattcaacgtttttgaaaaaaaaaaaagtgcaaagccgctcgta");
    basesBuff.append("agagatcctcgtttttctgatttttcaggatcgttcaatgcgaacttttttagaaatgca");
    basesBuff.append("tacaagtttttatacgactcaagagaacaagaaaagaagattatagaaaaaaaattaaaa");
    basesBuff.append("agcaaaaatataacacaagaagaaaaggacgaattaaaaaaaaaatacaatgattataaa");
    basesBuff.append("agcactgatatattattaaaaaaaaaagaagaagaaaggaaattaaaagcagaattagta");
    basesBuff.append("aaacaggaaaaacaaaatattctaaccaaaaataaaaaaccttattattatagtgataga");
    basesBuff.append("aaaattaaaaaaatagttcaagaaaaattgtcgagtaataagagtttaaagaaagttata");
    basesBuff.append("aaaaaagaaaaaaaaacattgcaaaaagagagaaagagaaacatcataccaacaaaaaaa");
    basesBuff.append("aaaatttttttaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa");
    basesBuff.append("aaaaaaaaaaaaaaaacataa");

    IOUtils.wrapString(basesBuff.toString(), fastaBuff);

    assertEquals("Export FASTA feature ", writer.toString(), fastaBuff.toString());
  }

  /**
   * Get sequences for a selected region.
   * @param isFwd
   * @param start
   * @param end
   * @return
   */
  private StringWriter getRegionWriter(boolean isFwd, int start, int end)
  {
    StringWriter writer = new StringWriter();
    Selection selection = new Selection(null);
    
    FeatureVector features = vcfView.getEntryGroup().getAllFeatures();
    Feature feature = null;
    for(int i=0; i<features.size(); i++)
    {
      if(features.elementAt(i).isForwardFeature() && isFwd)
      {
        feature = features.elementAt(i);
        break;
      }
      else if(!features.elementAt(i).isForwardFeature() && !isFwd)
      {
        feature = features.elementAt(i);
        
        int length = vcfView.getEntryGroup().getSequenceEntry().getBases().getLength();
        int tmp = start;
        start = length - end + 1;
        end   = length - tmp + 1;

        break;
      }
    }

    try
    {
      selection.setMarkerRange(new MarkerRange(feature.getStrand(), start, end));
    }
    catch(uk.ac.sanger.artemis.util.OutOfRangeException e)
    {
      e.printStackTrace();
    }
    IOUtils.exportFastaByRange(vcfView, selection, false, writer);
    return writer;
  }
  
  /**
   * Get sequences for features.
   * @param id comma separated list of feature id's
   * @return
   */
  private StringWriter getFeatureWriter(String ids, final VCFview view)
  {
    StringWriter writer = new StringWriter();
    String id[] = ids.split(",");
    
    FeatureVector features = view.getEntryGroup().getAllFeatures();
    FeatureVector selectedFeatures = new FeatureVector();
    Feature feature = null;
    for(int i=0; i<features.size(); i++)
    {
      for(int j=0; j<id.length; j++)
      {
        if(features.elementAt(i).getIDString().equals(id[j]))
        {
          feature = features.elementAt(i);
          selectedFeatures.add(feature);
        }
      }
    }
    
    IOUtils.exportFasta(view, selectedFeatures, false, writer);
    return writer;
  }
  
  
  /**
   * Load in a VCF and BCF into separate VCFviews to be used 
   * separately in testing.
   */
  @BeforeClass
  public static void oneTimeSetUp()
  {
    // load VCF
    URL ref = WriteVCFTest.class.getResource("/data/test.embl.gz");
    URL vcf = WriteVCFTest.class.getResource("/data/test.vcf.gz");
  
    List<String> vcfFileList = new Vector<String>();
    vcfFileList.add(vcf.getFile());

    vcfView = new VCFview(null, new JPanel(), 
        vcfFileList, 
        5000, 100000000, null, 
        ref.getFile(), null);
    
    // load BCF
    ref = WriteVCFTest.class.getResource("/data/MAL1.embl.gz");
    vcf = WriteVCFTest.class.getResource("/data/MAL1_8_16_24h.raw.bcf");
    
    List<String> bcfFileList = new Vector<String>();
    bcfFileList.add(vcf.getFile());
    bcfView = new VCFview(null, new JPanel(), 
        bcfFileList, 
        5000, 100000000, null, 
        ref.getFile(), null);
  }

}