package uk.ac.sanger.artemis.io;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.net.URL;
import java.util.LinkedHashMap;
import java.util.Map;

import org.apache.log4j.Level;
import org.junit.Test;

import junit.framework.Assert;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.components.EntryFileDialog;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.io.DocumentEntryFactory;
import uk.ac.sanger.artemis.io.Entry;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.FeatureVector;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.ValidateFeature;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.DocumentFactory;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.sequence.NoSequenceException;

public class ValidateFeatureTest
{
  @Test
  public void testGFF()
  {
    try 
    {
      final Entry entry = getEntry("/data/test.gff.gz");
      final EntryGroup egrp = new SimpleEntryGroup();
      egrp.add(new uk.ac.sanger.artemis.Entry(entry));
      ValidateFeature validate = new ValidateFeature(egrp);
      
      final FeatureVector features = entry.getAllFeatures();

      for(uk.ac.sanger.artemis.io.Feature f: features)
      {
        if(ValidateFeature.isGFF(f, null))
        {
          GFFStreamFeature gffFeature = (GFFStreamFeature)f;
          String id = GeneUtils.getUniqueName(gffFeature);
          
          assertTrue("Boundary "+id,  ValidateFeature.isBoundaryOK(gffFeature) == 0);
          assertTrue("Strand "+id,    ValidateFeature.isStrandOK(gffFeature));
          assertTrue("CDS phase "+id, ValidateFeature.isCDSPhaseOK(gffFeature));
          assertTrue("Attribute "+id, ValidateFeature.isAttributesOK(gffFeature).length() == 0);
          assertTrue("ID check "+id, ValidateFeature.isIdPrefixConsistent(gffFeature));
          assertTrue("Start_range check "+id, ValidateFeature.isPartialConsistent(gffFeature, "Start_range"));
          assertTrue("End_range check "+id, ValidateFeature.isPartialConsistent(gffFeature, "End_range"));

          if(ValidateFeature.isPartOfGene(gffFeature))
            assertTrue("Gene model "+id, ValidateFeature.isCompleteGeneModelOK(gffFeature) == 0);
        }

        assertTrue("Stop codon", validate.hasValidStop(f));
        assertTrue("Internal stop codon", !validate.isInternalStops(f));
      }
    }
    catch(OutOfRangeException e)
    {
      Assert.fail(e.getMessage());
    }
    catch(NoSequenceException e)
    {
      Assert.fail(e.getMessage());
    }
  }
  
  /**
   * Test the gene model boundary is consistent
   */
  @Test
  public void testGFFBoundary()
  {
    final Entry entry = getEntry("/data/test_boundary.gff.gz");
    final FeatureVector features = entry.getAllFeatures();

    for(uk.ac.sanger.artemis.io.Feature f: features)
    {
      if(ValidateFeature.isGFF(f, null))
      {
        GFFStreamFeature gffFeature = (GFFStreamFeature)f;
        String id = GeneUtils.getUniqueName(gffFeature);
        if(id.equals("PF3D7_0200100"))
          assertTrue("Boundary check: "+id, ValidateFeature.isBoundaryOK(gffFeature) != 0); // boundary not OK
        else
          assertTrue("Boundary check: "+id, ValidateFeature.isBoundaryOK(gffFeature) == 0); // boundary OK
      }
    }
  }
  
  /**
   * Test the gene model strand is consistent
   */
  @Test
  public void testGFFStrand()
  {
    final Entry entry = getEntry("/data/test_boundary.gff.gz");
    final FeatureVector features = entry.getAllFeatures();

    for(uk.ac.sanger.artemis.io.Feature f: features)
    {
      if(ValidateFeature.isGFF(f, null))
      {
        GFFStreamFeature gffFeature = (GFFStreamFeature)f;
        String id = GeneUtils.getUniqueName(gffFeature);
        if(id.equals("PF3D7_0200300"))
          assertTrue("Strand check: "+id, !ValidateFeature.isStrandOK(gffFeature)); // strand not ok
        else
          assertTrue("Strand check: "+id, ValidateFeature.isStrandOK(gffFeature));  // strand ok
      }
    }
  }
  
  /**
   * Test the CDS has a phase
   */
  @Test
  public void testGFFPhase()
  {
    final Entry entry = getEntry("/data/test_boundary.gff.gz");
    final FeatureVector features = entry.getAllFeatures();

    for(uk.ac.sanger.artemis.io.Feature f: features)
    {
      if(ValidateFeature.isGFF(f, null))
      {
        GFFStreamFeature gffFeature = (GFFStreamFeature)f;
        String id = GeneUtils.getUniqueName(gffFeature);
        if(id.equals("PF3D7_0200200.1:exon{2,1}"))
          assertTrue("CDS phase check: "+id, !ValidateFeature.isCDSPhaseOK(gffFeature)); // phase not ok
        else
          assertTrue("CDS phase check: "+id, ValidateFeature.isCDSPhaseOK(gffFeature));  // phase ok
      }
    }
  }
  
  /**
   * Test the gene models are complete
   */
  @Test
  public void testGFFCompleteGeneModel()
  {
    final Entry entry = getEntry("/data/test_boundary.gff.gz");
    final FeatureVector features = entry.getAllFeatures();

    for(uk.ac.sanger.artemis.io.Feature f: features)
    {
      if(ValidateFeature.isGFF(f, null))
      {
        GFFStreamFeature gffFeature = (GFFStreamFeature)f;
        if(ValidateFeature.isPartOfGene(gffFeature))
        {
          String id = GeneUtils.getUniqueName(gffFeature);
          if(id.startsWith("PF3D7_0200500")) 
            assertTrue("Complete gene model check: "+id,   // gene model missing mRNA
                ValidateFeature.isCompleteGeneModelOK(gffFeature) != 0);
          else
            assertTrue("Complete gene model check: "+id,   // gene model complete
                ValidateFeature.isCompleteGeneModelOK(gffFeature) == 0);
        }
      }
    }
  }
  
  
  /**
   * Test if the ID GFF3 attribute prefix is consistent within a gene model
   */
  @Test
  public void testGFFId()
  {
    final Entry entry = getEntry("/data/test_boundary.gff.gz");
    final FeatureVector features = entry.getAllFeatures();

    for(uk.ac.sanger.artemis.io.Feature f: features)
    {
      if(ValidateFeature.isGFF(f, null))
      {
        GFFStreamFeature gffFeature = (GFFStreamFeature)f;
        if(ValidateFeature.isPartOfGene(gffFeature))
        {
          String id = GeneUtils.getUniqueName(gffFeature);
          assertTrue("Complete gene model check: "+id,
              ValidateFeature.isIdPrefixConsistent(gffFeature));
        }
      }
    }
  }
  
  
  /**
   * Test if the Start_range and End_range are constistent within a gene model
   */
  @Test
  public void testGFFPartials()
  {
    final Entry entry = getEntry("/data/test_boundary.gff.gz");
    final FeatureVector features = entry.getAllFeatures();

    for(uk.ac.sanger.artemis.io.Feature f: features)
    {
      if(ValidateFeature.isGFF(f, null))
      {
        GFFStreamFeature gffFeature = (GFFStreamFeature)f;
        if(ValidateFeature.isPartOfGene(gffFeature))
        {
          String id = GeneUtils.getUniqueName(gffFeature);
          assertTrue("Start_range check: "+id,
              ValidateFeature.isPartialConsistent(gffFeature, "Start_range"));
          assertTrue("End_range check: "+id,
              ValidateFeature.isPartialConsistent(gffFeature, "End_range"));
        }
      }
    }
  }
  
  
  /**
   * Test if the Start_range and End_range are constistent within a gene model
   */
  @Test
  public void testGFFAttributes()
  {
    final Entry entry = getEntry("/data/test_boundary.gff.gz");
    final FeatureVector features = entry.getAllFeatures();

    for(uk.ac.sanger.artemis.io.Feature f: features)
    {
      if(ValidateFeature.isGFF(f, null))
      {
        GFFStreamFeature gffFeature = (GFFStreamFeature)f;
        if(ValidateFeature.isPartOfGene(gffFeature))
        {
          String id = GeneUtils.getUniqueName(gffFeature);
          assertTrue("Attributes : "+id,
              ValidateFeature.isAttributesOK(gffFeature).length() == 0);
        }
      }
    }
  }
  
  
  /**
   * Test stop codons for genes
   */
  @Test
  public void testGFFValidStop()
  {
    try
    {
      final Entry entry = getEntry("/data/test_boundary.gff.gz");
      final FeatureVector features = entry.getAllFeatures();
      final EntryGroup egrp = new SimpleEntryGroup();
      egrp.add(new uk.ac.sanger.artemis.Entry(entry));
      ValidateFeature validate = new ValidateFeature(egrp);

      for (uk.ac.sanger.artemis.io.Feature f : features)
      {
        GFFStreamFeature gffFeature = (GFFStreamFeature)f;
        String id = GeneUtils.getUniqueName(gffFeature);
        if( (id.startsWith("PF3D7_0200100") || id.startsWith("PF3D7_0200500.1:exon:2")) && 
            f.getKey().getKeyString().equals("CDS"))
          assertTrue("Stop codon "+id, !validate.hasValidStop(f)); // not valid
        else
          assertTrue("Stop codon "+id, validate.hasValidStop(f));
      }
    }
    catch (OutOfRangeException e)
    {
      Assert.fail(e.getMessage());
    }
    catch (NoSequenceException e)
    {
      Assert.fail(e.getMessage());
    }
  }

  /**
   * Test for internal stop codons
   */
  @Test
  public void testGFFInternalStop()
  {
    try
    {
      final Entry entry = getEntry("/data/test_boundary.gff.gz");
      final FeatureVector features = entry.getAllFeatures();
      final EntryGroup egrp = new SimpleEntryGroup();
      egrp.add(new uk.ac.sanger.artemis.Entry(entry));
      ValidateFeature validate = new ValidateFeature(egrp);

      for (uk.ac.sanger.artemis.io.Feature f : features)
      {
        GFFStreamFeature gffFeature = (GFFStreamFeature)f;
        String id = GeneUtils.getUniqueName(gffFeature);
        assertTrue("Internal stop codon "+id, !validate.isInternalStops(f));
      }
    }
    catch (OutOfRangeException e)
    {
      Assert.fail(e.getMessage());
    }
    catch (NoSequenceException e)
    {
      Assert.fail(e.getMessage());
    }
  }
  
  
  private Entry getEntry(final String gff)
  {
    try
    {
      URL gffFile = ValidateFeatureTest.class.getResource(gff);
      final Document doc = DocumentFactory.makeDocument(gffFile.getFile());
      return DocumentEntryFactory.makeDocumentEntry(
          Options.getArtemisEntryInformation(),doc,null);
    }
    catch(EntryInformationException e) 
    {
      Assert.fail(e.getMessage());
    }
    catch(IOException e) 
    {
      Assert.fail(e.getMessage());
    }
    return null;
  }
}