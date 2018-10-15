package uk.ac.sanger.artemis.components.variant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.net.URL;
import java.util.List;
import java.util.Vector;

import javax.swing.JPanel;

import org.junit.BeforeClass;
import org.junit.Test;

import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.components.variant.VCFview;
import uk.ac.sanger.artemis.sequence.Bases;


public class WriteVariantSitesTest
{ 
  /**
   * Test SNP's
   */
  @Test
  public void testSNPSites()
  {
    // load VCF
    final URL ref = WriteVCFTest.class.getResource("/data/test.embl.gz");
    final URL vcf = WriteVCFTest.class.getResource("/data/test1.vcf.gz");
    VCFview vcfView = getVcfView(ref, vcf);
    
    StringWriter sw = new StringWriter();
    PrintWriter pw  = new PrintWriter(sw);

    EntryGroup entryGroup = vcfView.getEntryGroup();
    
    IOUtils.exportVariantFasta(vcfView, pw, 
        entryGroup.getSequenceEntry().getBases().getLength(), 
        entryGroup.getAllFeatures(), 
        entryGroup.getSequenceEntry().getBases());

    StringBuilder fastaBuff = new StringBuilder(">1_sample\n");
    fastaBuff.append("TTTTGAACTTATACA\n");
    fastaBuff.append(">test.embl.gz\n");
    fastaBuff.append("CCCCATGTGGGCGTG\n");
    
    assertEquals("Export variants as FASTA ", sw.toString(), fastaBuff.toString());
  }
  
  /**
   * Test insertions
   */
  @Test
  public void testInsertionSites()
  {
    // load VCF
    final URL ref = WriteVCFTest.class.getResource("/data/test.embl.gz");
    final URL vcf = WriteVCFTest.class.getResource("/data/test2.vcf.gz");
    VCFview vcfView = getVcfView(ref, vcf);
    
    StringWriter sw = new StringWriter();
    PrintWriter pw  = new PrintWriter(sw);

    EntryGroup entryGroup = vcfView.getEntryGroup();
    
    IOUtils.exportVariantFasta(vcfView, pw, 
        entryGroup.getSequenceEntry().getBases().getLength(), 
        entryGroup.getAllFeatures(), 
        entryGroup.getSequenceEntry().getBases());

    StringBuilder fastaBuff = new StringBuilder(">1_sample\n");
    fastaBuff.append("tGgTttttttAaaTttaAGAaAGTAGa\n");
    fastaBuff.append(">test.embl.gz\n");
    fastaBuff.append("t-gttttt-ta-att-a---aa-----\n");
    
    assertEquals("Export variants as FASTA ", sw.toString(), fastaBuff.toString());
  }
  
  /**
   * Test Multi-allele
   */
  @Test
  public void testMultiAlleleSites()
  {
    // load VCF
    final URL ref = WriteVCFTest.class.getResource("/data/test.embl.gz");
    final URL vcf = WriteVCFTest.class.getResource("/data/test3.vcf.gz");
    VCFview vcfView = getVcfView(ref, vcf);
    
    StringWriter sw = new StringWriter();
    PrintWriter pw  = new PrintWriter(sw);

    EntryGroup entryGroup = vcfView.getEntryGroup();
    
    IOUtils.exportVariantFasta(vcfView, pw, 
        entryGroup.getSequenceEntry().getBases().getLength(), 
        entryGroup.getAllFeatures(), 
        entryGroup.getSequenceEntry().getBases());

    StringBuilder fastaBuff = new StringBuilder(">1_sample\n");
    fastaBuff.append("wwmrwr\n");
    fastaBuff.append(">test.embl.gz\n");
    fastaBuff.append("ggtgcc\n");
    
    assertEquals("Export variants as FASTA ", sw.toString(), fastaBuff.toString());
  }
  
  
  /**
   * Test deletion. Note: do not write out if just deletion.
   */
  @Test
  public void testDeletionSites()
  {
    // load VCF
    final URL ref = WriteVCFTest.class.getResource("/data/test.embl.gz");
    final URL vcf = WriteVCFTest.class.getResource("/data/test4.vcf.gz");
    VCFview vcfView = getVcfView(ref, vcf);
    
    StringWriter sw = new StringWriter();
    PrintWriter pw  = new PrintWriter(sw);

    EntryGroup entryGroup = vcfView.getEntryGroup();
    
    IOUtils.exportVariantFasta(vcfView, pw, 
        entryGroup.getSequenceEntry().getBases().getLength(), 
        entryGroup.getAllFeatures(), 
        entryGroup.getSequenceEntry().getBases());

    StringBuilder fastaBuff = new StringBuilder(">1_sample\n");
    fastaBuff.append(">test.embl.gz\n");
    
    assertEquals("Export variants as FASTA ", sw.toString(), fastaBuff.toString());
  }
  
  /**
   * Load in a VCF and BCF into separate VCFviews to be used 
   * separately in testing.
   */
  private VCFview getVcfView(final URL ref, final URL vcf)
  {
    List<String> vcfFileList = new Vector<String>();
    vcfFileList.add(vcf.getFile());

    return new VCFview(null, new JPanel(), 
        vcfFileList, 
        5000, 100000000, null, 
        ref.getFile(), null, null);
    
    // load BCF
/*    ref = WriteVCFTest.class.getResource("/data/MAL1.embl.gz");
    vcf = WriteVCFTest.class.getResource("/data/MAL1_8_16_24h.raw.bcf");
    
    List<String> bcfFileList = new Vector<String>();
    bcfFileList.add(vcf.getFile());
    bcfView = new VCFview(null, new JPanel(), 
        bcfFileList, 
        5000, 100000000, null, 
        ref.getFile(), null, null);*/
  }
}