/* 
 * This file is part of Artemis
 *
 * Copyright (C) 2013  Genome Research Limited
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
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.net.URL;
import java.util.LinkedHashMap;
import java.util.Map;

import org.apache.log4j.Level;
import org.junit.Test;

import junit.framework.Assert;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
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
    testAll(Utils.getEntryGroup("/data/test.gff.gz"));
  }
  
  public static void testAll(final EntryGroup egrp)
  {
    ValidateFeature validate = new ValidateFeature(egrp);
    final uk.ac.sanger.artemis.FeatureVector features = egrp.getAllFeatures();

    for(int i=0; i<features.size(); i++)
    {
      Feature artFeature = features.elementAt(i);
      uk.ac.sanger.artemis.io.Feature f = artFeature.getEmblFeature();
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
          assertTrue("Gene model not complete "+id, ValidateFeature.isCompleteGeneModelOK(gffFeature) == 0);
      }

      assertTrue("Stop codon "+artFeature.getIDString(), validate.hasValidStop(f));
      assertTrue("Internal stop codon "+artFeature.getIDString(), !validate.isInternalStops(f));
    }
  }
  
  /**
   * Test the gene model boundary is consistent
   */
  @Test
  public void testGFFBoundary()
  {
    final Entry entry = Utils.getEntry("/data/test_boundary.gff.gz");
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
    final Entry entry = Utils.getEntry("/data/test_boundary.gff.gz");
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
    final Entry entry = Utils.getEntry("/data/test_boundary.gff.gz");
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
    final Entry entry = Utils.getEntry("/data/test_boundary.gff.gz");
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
    final Entry entry = Utils.getEntry("/data/test_boundary.gff.gz");
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
    final Entry entry = Utils.getEntry("/data/test_boundary.gff.gz");
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
    final Entry entry = Utils.getEntry("/data/test_boundary.gff.gz");
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
      final Entry entry = Utils.getEntry("/data/test_boundary.gff.gz");
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
      final Entry entry = Utils.getEntry("/data/test_boundary.gff.gz");
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
  
}