/* ValidateFeature
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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;


import org.apache.log4j.Level;

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
import uk.ac.sanger.artemis.components.BasePlotGroup;
import uk.ac.sanger.artemis.components.EntryFileDialog;
import uk.ac.sanger.artemis.components.FeatureListFrame;
import uk.ac.sanger.artemis.components.FileViewer;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.sequence.AminoAcidSequence;
import uk.ac.sanger.artemis.sequence.Marker;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.DocumentFactory;
import uk.ac.sanger.artemis.util.StringVector;

public class ValidateFeature
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
    { "GO", "EC_number", "EMBL_qualifier", "SignalP_prediction", 
      "GPI_anchor_cleavage_site", "GPI_anchored", "PlasmoAP_score" };
  
  private EntryGroup entryGrp;
  private FeaturePredicate cds_predicate;

  public ValidateFeature(EntryGroup entryGrp)
  {
    this.entryGrp = entryGrp;
    this.cds_predicate = null;
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
  

  public void featureListErrors(final EntryGroup grp, 
                     final Selection sel,
                     final GotoEventSource gotoSrc,
                     final BasePlotGroup plotGrp)
  {
    if(GeneUtils.isGFFEntry( grp ))
    {
      showFeatureList(ATTR_PREDICATE, "Attributes", grp, sel, gotoSrc, plotGrp);
      if(!GeneUtils.isDatabaseEntry(grp))
        showFeatureList(CDS_PHASE_PREDICATE, "CDS phase", grp, sel, gotoSrc, plotGrp);

      showFeatureList(STRAND_PREDICATE, "Gene Strand Errors", grp, sel, gotoSrc, plotGrp);
      showFeatureList(BOUNDARY_PREDICATE, "Gene Boundary Errors", grp, sel, gotoSrc, plotGrp);
      showFeatureList(COMPLETE_GENE_MODEL_PREDICATE, "Incomplete Gene Model", grp, sel, gotoSrc, plotGrp);
    }
    
    showFeatureList(INTERNAL_STOP, "Internal Stop Codons", grp, sel, gotoSrc, plotGrp);
    showFeatureList(NO_VALID_STOP, "No Valid Stop Codons", grp, sel, gotoSrc, plotGrp);
  }
  
  private void showFeatureList(final FeaturePredicate predicate,
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
   * 
   * GFF - check complete gene model
   *     - check boundaries are valid
   *     - check all features are on the same strand
   *     - check CDS features have a phase
   *     - check attribute column 
   *          - qualifiers have a value (not empty)
   *          - only reserved tags start with uppercase
   * - CDS have no internal stop codon
   * - CDS have valid stop codon
   * 
   * @param f
   * @param fv
   * @param showOnlyFailures
   * @return true if passed
   */
  public boolean featureValidate(final uk.ac.sanger.artemis.io.Feature f, 
                                 final FileViewer fv, 
                                 final boolean showOnlyFailures)
  {
    boolean pass = true;
    
    final LinkedHashMap<String, Level> report = new LinkedHashMap<String, Level>();
    final String fTxt = featureTxt(f);
    final boolean isGFF = entryGrp == null || GeneUtils.isGFFEntry( entryGrp );
    if(isGFF)
    {
      final GFFStreamFeature gffFeature = (GFFStreamFeature)f;
      report.put("\n"+GeneUtils.getUniqueName(gffFeature)+" ("+fTxt+"):", Level.INFO);
      
      if(isPartOfGene(gffFeature))
      {
        int compl = isCompleteGeneModelOK(gffFeature);
        switch (compl)
        {
          case 1:  report.put("No gene found", Level.FATAL); pass = false; break;
          case 2:  report.put("Missing transcript", Level.FATAL); pass = false; break;
          default: report.put("Gene model is complete", Level.INFO); break;
        }
      
        final int gb = isBoundaryOK(gffFeature);
        if(gb == 0)
          report.put(getGeneBoundaryMsg(gb), Level.INFO);
        else
        {
          pass = false;
          report.put(getGeneBoundaryMsg(gb), Level.FATAL);
        }

        if(isStrandOK(gffFeature))
          report.put("Gene features all on same strand", Level.INFO);
        else
        {
          pass = false;
          report.put("Gene features found on different strand", Level.FATAL);
        }
      }
    
      if( (entryGrp == null || !GeneUtils.isDatabaseEntry(entryGrp)) && !isCDSPhaseOK(gffFeature))
      {
        pass = false;
        report.put("CDS phase (codon_start) not set", Level.FATAL);
      }

      String attr = isAttributesOK(gffFeature);
      if(!attr.equals(""))
      {
        pass = false;
        if(attr.endsWith("\n"))
          attr = attr.substring(0, attr.length()-1);
        report.put(attr, Level.FATAL);
      }
    }
    else
      report.put("\n"+((uk.ac.sanger.artemis.Feature)f.getUserData()).getIDString()+
          " ("+fTxt+"):", Level.INFO);

    if(isInternalStops(f))
    {
      pass = false;
      report.put("Internal stop codon found", Level.FATAL);
    }
    
    if(!hasValidStop(f))
    {
      pass = false;
      report.put("No valid stop codon found", Level.FATAL);
    }
    
    if(pass)
      report.put("PASS", Level.INFO);

    if(!pass || !showOnlyFailures)
    {
      for (Map.Entry<String, Level> entry : report.entrySet())
        fv.appendString(entry.getKey()+"\n", entry.getValue());
    }
    
    return pass;
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
      return 1;
    
    final List<uk.ac.sanger.artemis.io.Feature> transcripts = gene.getTranscripts();
    if(transcripts.size() < 1)
      return 2;

    return 0;
  }
  
  public static int isBoundaryOK(final GFFStreamFeature gffFeature)
  {
    final ChadoCanonicalGene gene = gffFeature.getChadoGene();
    int gb = 0;
    if(gene != null && isGene(gffFeature))
      gb = GeneUtils.isBoundaryOK(gene);
    return gb;
  }
  
  public static boolean isStrandOK(final GFFStreamFeature gffFeature)
  {
    final ChadoCanonicalGene gene = gffFeature.getChadoGene();
    if(gene != null && isGene(gffFeature) && !GeneUtils.isStrandOK(gene))
      return false;
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
      return false;
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
          str.append(msg);
        }
      }
      
      // check format tag=value
      final StringVector values = qualifier.getValues();
      if(  values == null || values.size() < 1 ||
         ( values.size() == 1 && values.get(0).equals("")) )
      {
        String msg = qualifier.getName()+" atribute has no value\n";
        str.append(msg);
      }
    }  
    return str.toString();
  }
  
  public boolean isInternalStops(final uk.ac.sanger.artemis.io.Feature feature)
  { 
    if(feature.getUserData() == null)
      return false;
    final uk.ac.sanger.artemis.Feature f = (uk.ac.sanger.artemis.Feature) feature.getUserData();
    final FeaturePredicate cds_predicate = getCodingFeaturePredicate();
    if(!cds_predicate.testPredicate (f)) 
      return false;

    final AminoAcidSequence aa = f.getTranslation ();
    return aa.containsStopCodon ();
  }
  
  public boolean hasValidStop(final uk.ac.sanger.artemis.io.Feature feature)
  {
    final FeaturePredicate cds_predicate = getCodingFeaturePredicate();
    if(feature.getUserData() == null)
      return true;
    final uk.ac.sanger.artemis.Feature f = (uk.ac.sanger.artemis.Feature) feature.getUserData();
    if(!cds_predicate.testPredicate (f)) 
      return true;

    return f.hasValidStopCodon (true);
  }
  
  private FeaturePredicate getCodingFeaturePredicate()
  {
    if(cds_predicate != null)
      return cds_predicate;

    if(entryGrp != null && GeneUtils.isDatabaseEntry( entryGrp ))
    {
      
      final FeaturePredicate codingPredicate = new FeaturePredicate(){
        private String nonCodingTranscripts[] = GeneUtils.getNonCodingTranscripts();
        public boolean testPredicate(Feature feature)
        {
          try
          {
            final GFFStreamFeature gffF = (GFFStreamFeature)feature.getEmblFeature();
            final ChadoCanonicalGene chadoGene = gffF.getChadoGene();
            if(chadoGene != null)
            {
              uk.ac.sanger.artemis.io.Feature transcriptF =
                chadoGene.getTranscriptFeatureFromName(GeneUtils.getUniqueName(gffF));
              final String transcriptKey = transcriptF.getKey().getKeyString();

              for(int i=0; i<nonCodingTranscripts.length; i++)
                if(nonCodingTranscripts[i].equals(transcriptKey))
                  return false;
            }
          }
          catch(Exception e){}
          return true;
        };
      };
      cds_predicate =
          new FeaturePredicateConjunction(
              new FeatureKeyPredicate(new Key(DatabaseDocument.EXONMODEL)),
              codingPredicate,
              FeaturePredicateConjunction.AND);
    }
    else
      cds_predicate =
          new FeaturePredicateConjunction(
              new FeatureKeyQualifierPredicate(Key.CDS, "pseudo", false),
              new FeatureKeyQualifierPredicate(Key.CDS, "pseudogene", false),
              FeaturePredicateConjunction.AND);
    
    return cds_predicate;
  }
  
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
  
  
  private String featureTxt(final uk.ac.sanger.artemis.io.Feature f)
  {
    final Feature feature = (Feature) (f.getUserData());
    if(feature == null)
      return f.getLocation().toStringShort();
    final Marker low_marker  = feature.getFirstBaseMarker();
    final Marker high_marker = feature.getLastBaseMarker();
    String low_pos, high_pos;
    if(low_marker == null || high_marker == null) 
    {
      low_pos  = "unknown";
      high_pos = "unknown";
    }
    else 
    {
      if(low_marker.getRawPosition() < high_marker.getRawPosition()) 
      {
        low_pos = String.valueOf(low_marker.getRawPosition());
        high_pos = String.valueOf(high_marker.getRawPosition());
      } 
      else
      {
        low_pos = String.valueOf(high_marker.getRawPosition());
        high_pos = String.valueOf(low_marker.getRawPosition());
      }
    }
    
    if(GeneUtils.isDatabaseEntry(entryGrp))
    {
      try
      {
        if(feature.getQualifierByName("isFminPartial") != null)
          low_pos = "<"+low_pos;
        if(feature.getQualifierByName("isFmaxPartial") != null)
          high_pos = ">"+high_pos;
      }
      catch (InvalidRelationException e){}
    }
    else
    {
      if(feature.getLocation().isPartial(true)) // 5prime
      {
        if(feature.isForwardFeature()) 
          low_pos = "<"+low_pos;
        else
          high_pos = ">"+high_pos;
      }
      if(feature.getLocation().isPartial(false)) // 3prime
      {
        if(feature.isForwardFeature())
          high_pos = ">"+high_pos;
        else
          low_pos = "<"+low_pos;
      }
    }
    
    final StringBuilder txt = new StringBuilder();
    txt.append(feature.getKey().getKeyString()).append(" ");
    txt.append(low_pos).append("..").append(high_pos);

    if(!feature.isForwardFeature()) 
      txt.append(" c");
    return txt.toString();
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
    public boolean testPredicate(uk.ac.sanger.artemis.Feature feature)
    {
      return isAttributesOK((GFFStreamFeature) feature.getEmblFeature()).length() > 0;
    }
  };
  
  private static FeaturePredicate CDS_PHASE_PREDICATE = new FeaturePredicate() 
  {
    public boolean testPredicate(uk.ac.sanger.artemis.Feature feature)
    {
      return !isCDSPhaseOK((GFFStreamFeature) feature.getEmblFeature());
    }
  };

  private static FeaturePredicate STRAND_PREDICATE = new FeaturePredicate() 
  {
    public boolean testPredicate(uk.ac.sanger.artemis.Feature feature)
    {
      return !isStrandOK((GFFStreamFeature) feature.getEmblFeature());
    }
  };
  
  private static FeaturePredicate BOUNDARY_PREDICATE = new FeaturePredicate() 
  {
    public boolean testPredicate(uk.ac.sanger.artemis.Feature feature)
    {
      return isBoundaryOK((GFFStreamFeature) feature.getEmblFeature()) > 0;
    }
  };
  
  private static FeaturePredicate COMPLETE_GENE_MODEL_PREDICATE = new FeaturePredicate() 
  {
    public boolean testPredicate(uk.ac.sanger.artemis.Feature feature)
    {
      if(!isPartOfGene((GFFStreamFeature) feature.getEmblFeature()))
        return false;
      return isCompleteGeneModelOK((GFFStreamFeature) feature.getEmblFeature()) > 0;
    }
  };
  
  private FeaturePredicate INTERNAL_STOP = new FeaturePredicate() 
  {
    public boolean testPredicate(uk.ac.sanger.artemis.Feature feature)
    {
      return isInternalStops(feature.getEmblFeature());
    }
  };
  
  
  private FeaturePredicate NO_VALID_STOP = new FeaturePredicate() 
  {
    public boolean testPredicate(uk.ac.sanger.artemis.Feature feature)
    {
      return !hasValidStop(feature.getEmblFeature());
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
        
        ValidateFeature gffTest = new ValidateFeature(null);
        FeatureVector features = entry.getAllFeatures();
        final FileViewer fv = new FileViewer("Validation Report :: "+features.size()+" feature(s)");
        for(int j=0; j<features.size();j++)
        {
          try
          {
            gffTest.featureValidate((uk.ac.sanger.artemis.io.Feature)features.elementAt(j), fv, true);
          }
          catch(ClassCastException e)
          {
            e.printStackTrace();
          }
        }
        
        fv.setVisible(true);

        System.out.println("Done");
      }
    }
  }
  
}