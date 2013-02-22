/* GFFTest
 *
 * created: 2013
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2013 Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or(at your option) any later version.
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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.List;
import java.util.regex.Pattern;

import javax.swing.JOptionPane;

import org.apache.log4j.Level;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureKeyPredicate;
import uk.ac.sanger.artemis.FeatureKeyQualifierPredicate;
import uk.ac.sanger.artemis.FeaturePredicate;
import uk.ac.sanger.artemis.FeaturePredicateConjunction;
import uk.ac.sanger.artemis.FilteredEntryGroup;
import uk.ac.sanger.artemis.GotoEventSource;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.components.BasePlotGroup;
import uk.ac.sanger.artemis.components.EntryFileDialog;
import uk.ac.sanger.artemis.components.FeatureListFrame;
import uk.ac.sanger.artemis.components.FileViewer;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.sequence.AminoAcidSequence;
import uk.ac.sanger.artemis.sequence.NoSequenceException;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.DocumentFactory;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.StringVector;

public class GFFTest
{
  private static String[] geneModelParts;
  
  //##sequence-region seqid start end
  private static Pattern HEADER_SEQ_REGION = Pattern.compile("##sequence-region \\S+ \\d+ \\d+");
  private static Pattern CAPITAL_START = Pattern.compile("^[A-Z]+\\S*");
  
  private static String[] RESERVED_TAGS =
    { "ID", "Name", "Alias", "Parent", 
      "Target", "Gap", "Derives_from", "Note", 
      "Dbxref", "Ontology_term", "Is_circular"};

  private static String[] OTHER_RESERVED_TAGS =
    { "GO", "EC_number", "EMBL_qualifier", "SignalP_prediction", "GPI_anchor_cleavage_site", "GPI_anchored" };
  
  private EntryGroup entryGrp;
  
  public GFFTest(EntryGroup entryGrp)
  {
    this.entryGrp = entryGrp;
  }

  public static void testHeader(final String headerTxt)
  {
    final BufferedReader reader = new BufferedReader(
        new StringReader(headerTxt)); 
    try 
    {
      String str;
      while ((str = reader.readLine()) != null) 
      {           
        if(str.startsWith("##sequence-region "))
        {
          if(!HEADER_SEQ_REGION.matcher(str).matches())
          {
            System.out.println("ERROR : HEADER "+str);
          }
        }
      }
    } 
    catch(IOException e) 
    {
      e.printStackTrace();
    }
  }
  
  /**
   * Check features
   * @param features
   */
  public void testFeatures(final FeatureVector features)
  {
    final FileViewer fv = new FileViewer("Validation Report");
    for(int i=0; i<features.size(); i++)
    {
      try
      {
        testFeature((GFFStreamFeature)features.elementAt(i), fv);
      }
      catch(ClassCastException e){}
    }
    
    fv.setVisible(true);
  }
  
  public void featureListErrors(final EntryGroup grp, 
                     final Selection sel,
                     final GotoEventSource gotoSrc,
                     final BasePlotGroup plotGrp)
  {
    showReport(ATTR_PREDICATE, "Attributes", grp, sel, gotoSrc, plotGrp);
    
    if(!GeneUtils.isDatabaseEntry(grp))
      showReport(CDS_PHASE_PREDICATE, "CDS phase", grp, sel, gotoSrc, plotGrp);
    showReport(STRAND_PREDICATE, "Gene Strand Errors", grp, sel, gotoSrc, plotGrp);
    showReport(BOUNDARY_PREDICATE, "Gene Boundary Errors", grp, sel, gotoSrc, plotGrp);
    showReport(COMPLETE_GENE_MODEL_PREDICATE, "Incomplete Gene Model", grp, sel, gotoSrc, plotGrp);
    showReport(INTERNAL_STOP, "Internal Stop Codons", grp, sel, gotoSrc, plotGrp);
  }
  
  private void showReport(final FeaturePredicate predicate,
                          final String title,
                          final EntryGroup grp, 
                          final Selection sel,
                          final GotoEventSource gotoSrc,
                          final BasePlotGroup plotGrp)
  {
    final FilteredEntryGroup fltrGrp = new FilteredEntryGroup (grp, predicate, title);
    if(fltrGrp.getAllFeaturesCount() < 1)
      return;
    
    final FeatureListFrame featureList =
        new FeatureListFrame (title, sel, gotoSrc, fltrGrp, plotGrp);
    featureList.setVisible(true);
  }
  
  
  /**
   * Check a single feature
   * @param gffFeature
   */
  public void testFeature(final GFFStreamFeature gffFeature, final FileViewer fv)
  {
    fv.appendString("\n"+GeneUtils.getUniqueName(gffFeature)+":\n", Level.INFO);
    if(isPartOfGene(gffFeature))
    {
      int compl = isCompleteGeneModelOK(gffFeature);
      
      switch (compl)
      {
        case 1: fv.appendString("No gene found\n", Level.FATAL); break;
        case 2: fv.appendString("Missing transcript\n", Level.FATAL); break;
        default: fv.appendString("Gene model is complete\n", Level.INFO); break;
      }
      
      int gb = isBoundaryOK(gffFeature);
      if(gb == 0)
        fv.appendString(getGeneBoundaryMsg(gb)+"\n", Level.INFO);
      else
        fv.appendString(getGeneBoundaryMsg(gb)+"\n", Level.FATAL);

      if(isStrandOK(gffFeature))
        fv.appendString("Gene features all on same strand\n", Level.INFO);
      else
        fv.appendString("Gene features found on different strand\n", Level.FATAL);
      
    }
    if(!isCDSPhaseOK(gffFeature))
      fv.appendString("CDS phase (codon_start) not set\n", Level.FATAL);
    
    if(isInternalStops(gffFeature))
      fv.appendString("Internal stop codon found\n", Level.FATAL);;

    final String attr = isAttributesOK(gffFeature);
    fv.appendString(attr, Level.FATAL);
  }
  
  /**
   * Check if the gene model is complete
   * @param gffFeature
   * @return 0 - complete gene model
   *         1 - no gene found
   *         2 - missing transcript
   */
  public static int isCompleteGeneModelOK(final GFFStreamFeature gffFeature)
  {
    final ChadoCanonicalGene gene = gffFeature.getChadoGene();
    if(gene == null)
    {
      //errorLog(gffFeature, "no gene found");
      return 1;
    }
    
    final List<uk.ac.sanger.artemis.io.Feature> transcripts = gene.getTranscripts();
    if(transcripts.size() < 1)
    {
      //errorLog(gffFeature, "gene feature with no transcript");
      return 2;
    }
    return 0;
  }
  
  public static int isBoundaryOK(final GFFStreamFeature gffFeature)
  {
    final ChadoCanonicalGene gene = gffFeature.getChadoGene();
    int gb = 0;
    if(gene != null && isGene(gffFeature) && (gb = GeneUtils.isBoundaryOK(gene)) > 0)
    {
      //errorLog(gffFeature, "gene boundary inconsistent : "+getGeneBoundaryMsg(gb));
    }
    return gb;
  }
  
  public static boolean isStrandOK(final GFFStreamFeature gffFeature)
  {
    final ChadoCanonicalGene gene = gffFeature.getChadoGene();
    if(gene != null && isGene(gffFeature) && !GeneUtils.isStrandOK(gene))
    {
      //errorLog(gffFeature, "gene model features on different strand");
      return false;
    }
    return true;
  }
  
  /**
   * The phase is REQUIRED for all CDS features.
   * @param gffFeature
   */
  public static boolean isCDSPhaseOK(final GFFStreamFeature gffFeature)
  {
    if(!gffFeature.getKey().getKeyString().equals("CDS"))
      return true;
    
    final Qualifier codonStart = gffFeature.getQualifierByName("codon_start");
    if(codonStart == null)
    {
      //errorLog(gffFeature, "cds phase (codon_start) absent");
      return false;
    }
    return true;
  }
 
  
  /**
   * Check attribute column
   * @param gffFeature
   */
  public static String isAttributesOK(final GFFStreamFeature gffFeature)
  {
    final StringBuilder str = new StringBuilder();
    final QualifierVector qualifiers = gffFeature.getQualifiers();
    for(Qualifier qualifier: qualifiers)
    {
      if(CAPITAL_START.matcher(qualifier.getName()).matches())
      {
        boolean found = false;
        for(String r: RESERVED_TAGS)
          if(qualifier.getName().equals(r))
            found = true;

        for(String r: OTHER_RESERVED_TAGS)
          if(qualifier.getName().equals(r))
            found = true;

        if(!found)
        {
          String msg = qualifier.getName()+" non-reserved attribute name begins with uppercase\n";
          //errorLog(gffFeature, msg);
          str.append(msg);
        }
      }
      
      // check format tag=value
      final StringVector values = qualifier.getValues();
      if(  values == null || values.size() < 1 ||
         ( values.size() == 1 && values.get(0).equals("")) )
      {
        String msg = qualifier.getName()+" atribute has no value\n";
        //errorLog(gffFeature, msg);
        str.append(msg);
      }
    }
    return str.toString();
  }
  
  public boolean isInternalStops(final uk.ac.sanger.artemis.io.Feature feature)
  {
    final FeaturePredicate cds_predicate;
    if(entryGrp != null && GeneUtils.isDatabaseEntry( entryGrp ))
      cds_predicate = new FeatureKeyPredicate(new Key(DatabaseDocument.EXONMODEL));
    else
      cds_predicate =
          new FeaturePredicateConjunction(
              new FeatureKeyQualifierPredicate(Key.CDS, "pseudo", false),
              new FeatureKeyQualifierPredicate(Key.CDS, "pseudogene", false),
              FeaturePredicateConjunction.AND);
    
    final Feature f = (Feature) feature.getUserData();
    if(!cds_predicate.testPredicate (f)) 
      return false;

    final AminoAcidSequence amino_acids = f.getTranslation ();
    if(amino_acids.containsStopCodon ())
      return true;
    else
      return false;
  }
  
/*  private  static void errorLog(final GFFStreamFeature gffFeature, final String msg)
  {
    System.err.println("ERROR : "+msg+" "+GeneUtils.getUniqueName(gffFeature));
  }*/
  
  private static String getGeneBoundaryMsg(int gb)
  {
    String msg = "Valid gene boundary";
    switch(gb)
    {
      case 1: msg = "Transcript start or end is outside gene range";
              break;
      case 2: msg = "Child feature of a transcript is outside the transcript range";
              break;
      case 3: msg = "Span of the children features does not match start and end of the transcript";
              break;
      case 4: msg = "Protein range does not match CDS";
              break;
      case 5: msg = "Gene range does not match the largest transcript range";
              break;
    }
    return msg;
  }
  
  private static boolean isPartOfGene(final GFFStreamFeature gffFeature)
  {
    final String keyStr = gffFeature.getKey().getKeyString();
    for(String part: getGeneModelParts())
      if( part.equals(keyStr) )
        return true;

    return false;
  }
  
  private static boolean isGene(final GFFStreamFeature gffFeature)
  {
    final String keyStr = gffFeature.getKey().getKeyString();
    return keyStr.equals("gene") || keyStr.equals("pseudogene");
  }
  
  private static String[] getGeneModelParts()
  {
    if(geneModelParts != null)
      return geneModelParts;

    final String[] geneModelKeys = { 
        "gene", "transcript", "mRNA", "CDS", "exon", "polypeptide", "three_prime_UTR", "five_prime_UTR", 
        "pseudogene", "pseudogenic_transcript", "pseudogenic_exon" };
    final String[] ncTranscripts = GeneUtils.getNonCodingTranscripts();
    geneModelParts = new String[
          geneModelKeys.length+ncTranscripts.length];
    
    for(int i=0; i<geneModelKeys.length; i++)
      geneModelParts[i] = geneModelKeys[i];
    
    for(int i=0; i<ncTranscripts.length; i++)
      geneModelParts[i+geneModelKeys.length] = ncTranscripts[i];
    
    return geneModelParts;
  }
  
  
  private static FeaturePredicate ATTR_PREDICATE = new FeaturePredicate() 
  {
    public boolean testPredicate(Feature feature)
    {
      return isAttributesOK((GFFStreamFeature) feature.getEmblFeature()).length() > 0;
    }
  };
  
  private static FeaturePredicate CDS_PHASE_PREDICATE = new FeaturePredicate() 
  {
    public boolean testPredicate(Feature feature)
    {
      return !isCDSPhaseOK((GFFStreamFeature) feature.getEmblFeature());
    }
  };

  private static FeaturePredicate STRAND_PREDICATE = new FeaturePredicate() 
  {
    public boolean testPredicate(Feature feature)
    {
      return !isStrandOK((GFFStreamFeature) feature.getEmblFeature());
    }
  };
  
  private static FeaturePredicate BOUNDARY_PREDICATE = new FeaturePredicate() 
  {
    public boolean testPredicate(Feature feature)
    {
      return isBoundaryOK((GFFStreamFeature) feature.getEmblFeature()) > 0;
    }
  };
  
  private static FeaturePredicate COMPLETE_GENE_MODEL_PREDICATE = new FeaturePredicate() 
  {
    public boolean testPredicate(Feature feature)
    {
      if(!isPartOfGene((GFFStreamFeature) feature.getEmblFeature()))
        return false;
      return isCompleteGeneModelOK((GFFStreamFeature) feature.getEmblFeature()) > 0;
    }
  };
  
  private FeaturePredicate INTERNAL_STOP = new FeaturePredicate() 
  {
    public boolean testPredicate(Feature feature)
    {
      return isInternalStops(feature.getEmblFeature());
    }
  };
  
 
  public static void main(final String args[])
  {
    System.setProperty("black_belt_mode","true");
    Options.getOptions();
    
    if (args.length > 0)
    {
      for(int i=0; i<args.length; i++)
      {
        System.out.println("\nTEST: "+args[i]);
        final Document entry_document = DocumentFactory.makeDocument(args[i]);

        final EntryInformation artemis_entry_information = 
            Options.getArtemisEntryInformation();
        final uk.ac.sanger.artemis.io.Entry entry = EntryFileDialog.getEntryFromFile(
            null, entry_document, artemis_entry_information, false);
        
        GFFTest gffTest = new GFFTest(null);
        gffTest.testFeatures(entry.getAllFeatures());
        System.out.println("Done");
      }
    }
    else
    {
      uk.ac.sanger.artemis.components.FileDialogEntrySource entrySource = 
          new uk.ac.sanger.artemis.components.FileDialogEntrySource(null, null);   
      final EntryGroup entryGrp = new SimpleEntryGroup();
      
      try
      {
        final Entry entry = entrySource.getEntry(true);
        entryGrp.add(entry);

        GFFTest gffTest = new GFFTest(entryGrp);
        gffTest.testFeatures(entry.getEMBLEntry().getAllFeatures());
        GFFTest.testHeader(entryGrp.getDefaultEntry().getHeaderText());
      }
      catch(OutOfRangeException e)
      {
        e.printStackTrace();
      }
      catch(NoSequenceException e)
      {
        JOptionPane.showMessageDialog(null, "No sequence found!", 
          "Sequence Missing", JOptionPane.WARNING_MESSAGE);
      }
    }
  }
  
}