
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
  private VCFview vcfView;
  
  /**
   * Write FASTA, from a range selection on the forward strand.
   * Variation is a SNP.
   * VCF line : test    120768  .       C       T
   */
  @Test
  public void testWriteFwdFastaSNP()
  {
    StringWriter writer = getWriter(true, 120768, 120777);
    
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
  public void testWriteRevFastaSNP()
  {
    StringWriter writer = getWriter(false, 120768, 120777);

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
  public void testWriteFwdFastaDeletion()
  {
    StringWriter writer = getWriter(true, 396835, 396845);

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
  public void testWriteRevFastaDeletion()
  {
    StringWriter writer = getWriter(false, 396835, 396845);

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
  public void testWriteFwdFastaInsertion()
  {
    StringWriter writer = getWriter(true, 366785, 366795);

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
  public void testWriteRevFastaInsertion()
  {
    StringWriter writer = getWriter(false, 366785, 366795);

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
  public void testWriteFwdFastaMultipleAllele()
  {
    StringWriter writer = getWriter(true, 361975, 361985);

    StringBuffer buff = new StringBuffer("> test.embl.gz 361975:361985\n");
    buff.append("actgaaaaatt\n");
    buff.append(">test.vcf.gz test.embl.gz 361975:361985\n");
    buff.append("actkaaaaatt\n");
    assertEquals("Export FASTA range ", writer.toString(), buff.toString());
  }
  
  
  /**
   * G
   * @param isFwd
   * @param start
   * @param end
   * @return
   */
  private StringWriter getWriter(boolean isFwd, int start, int end)
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
   * Open a flat file and create the components in the TransferAnnotaionTool
   * used to control the transfer of annotation.
   */
  @Before
  public void setup()
  {
    final URL ref = WriteVCFTest.class.getResource("/data/test.embl.gz");
    final URL vcf = WriteVCFTest.class.getResource("/data/test.vcf.gz");
  
    List<String> vcfFileList = new Vector<String>();
    vcfFileList.add(vcf.getFile());
    JFrame f = new JFrame();
    vcfView = new VCFview(f, (JPanel) f.getContentPane(), 
        vcfFileList, 
        5000, 100000000, null, 
        ref.getFile(), null);
    f.setVisible(false);
  }

}