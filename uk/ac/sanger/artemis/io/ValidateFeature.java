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
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringReader;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.JOptionPane;

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
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.components.BasePlotGroup;
import uk.ac.sanger.artemis.components.EntryFileDialog;
import uk.ac.sanger.artemis.components.FeatureListFrame;
import uk.ac.sanger.artemis.components.FileViewer;
import uk.ac.sanger.artemis.components.database.DatabaseEntrySource;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.components.genebuilder.cv.GoBox;
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
  private static Pattern ID_PREFIX = Pattern.compile("^[^.:]+");
  
  private static String[] RESERVED_TAGS =
    { "ID", "Name", "Alias", "Parent", 
      "Target", "Gap", "Derives_from", "Note", 
      "Dbxref", "Ontology_term", "Is_circular",
      "Start_range", "End_range"};

  private static String[] OTHER_RESERVED_TAGS =
    { "GO", "EC_number", "EMBL_qualifier", "SignalP_prediction", 
      "GPI_anchor_cleavage_site", "GPI_anchored", "PlasmoAP_score" };
  
  private static String[] ALLOWED_TAGS_WITH_NO_VALUE =
    { "stop_codon_redefined_as_selenocysteine", "Name" };
  
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
      showFeatureList(PARTIAL_PREDICATE, "Check Partial Settings", grp, sel, gotoSrc, plotGrp);
      showFeatureList(ID_PREDICATE, "Check ID Settings", grp, sel, gotoSrc, plotGrp);
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
  
  protected static boolean isGFF(final uk.ac.sanger.artemis.io.Feature f, final EntryGroup entryGrp)
  {
    final boolean isGFF = entryGrp == null || GeneUtils.isGFFEntry( entryGrp );
    if(isGFF && f instanceof GFFStreamFeature)
      return true;
    return false;
  }

  /**
   * Check a single feature
   * GFF - check complete gene model
   *     - check boundaries are valid
   *     - check all features are on the same strand
   *     - check CDS features have a phase
   *     - check partial attributes consistent
   *     - check attribute column 
   *          - qualifiers have a value (not empty)
   *          - only reserved tags start with uppercase
   * - CDS have no internal stop codon
   * - CDS have valid stop codon
   * 
   * @param f
   * @param showOnlyFailures
   * @return report
   */
  private LinkedHashMap<String, Level> featureValidate(final uk.ac.sanger.artemis.io.Feature f, 
                                                      final boolean showOnlyFailures)
  {
    boolean pass = true;

    final LinkedHashMap<String, Level> report = new LinkedHashMap<String, Level>();
    final String fTxt = featureTxt(f);
    if(isGFF(f, entryGrp))
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

        if(!isPartialConsistent(gffFeature, "Start_range") || 
           !isPartialConsistent(gffFeature, "End_range"))
        {
          pass = false;
          report.put("Partial settings not consistent", Level.FATAL);
        }
        
        if(!isIdPrefixConsistent(gffFeature))
        {
          pass = false;
          report.put("Prefix of ID attribute not consistent within gene model", Level.FATAL);
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
    {
      if(f.getUserData() == null)
        return report;
      report.put("\n"+((uk.ac.sanger.artemis.Feature)f.getUserData()).getIDString()+
          " ("+fTxt+"):", Level.INFO);
    }

    if(isInternalStops(f))
    {
      pass = false;
      report.put("Internal stop codon found", Level.FATAL);
    }
    
    final EntryInformation eInfo;
    if(f.getEntry() == null)
      eInfo = ((Feature)f.getUserData()).getEntry().getEntryInformation();
    else 
      eInfo = f.getEntry().getEntryInformation();
    
    // test GO qualifier
    final String goErrMsg;
    if((goErrMsg = validateGO(f.getQualifiers(), eInfo)).length() > 0)
    {
      pass = false;
      report.put(goErrMsg, Level.FATAL);
    }

    if(!hasValidStop(f))
    {
      pass = false;
      report.put("No valid stop codon found", Level.FATAL);
    }
    
    if(pass)
    {
      if(showOnlyFailures)
        return null;
      report.put("PASS", Level.INFO);
    }   
    return report;
  }
  
  /**
   * Check a single feature and append validation report to the
   * viewer.
   * @return true if successfully validate
   */
  public boolean featureValidate(final uk.ac.sanger.artemis.io.Feature f,
                                 final FileViewer fv,
                                 final boolean showFailedFeaturesOnly)
  {
    final LinkedHashMap<String, Level> report = featureValidate(f, showFailedFeaturesOnly);
    boolean pass = false;
    if (report != null)
    {
      for (Map.Entry<String, Level> entry : report.entrySet())
      {
        if (!showFailedFeaturesOnly && entry.getKey().equals("PASS"))
          pass = true;
        fv.appendString(entry.getKey() + "\n", entry.getValue());
      }
    }
    else
      return true;
    return pass;
  }

  /**
   * Check GO annotation for:
   * 1) unexpected white space in with/from and dbxref columns
   * 2) with/from must be empty when using IDA, NAS, ND, TAS or EXP evidence code
   * 3) GO:0005515 can only have IPI evidence code
   * 4) IEP is not allowed for molecular_function and cellular_component terms
   * 5) with field must be filled when using ISS, ISA, ISO and ISM codes.
   * @param qualifiers
   * @param eInfo
   * @return errors found in GO annotation
   */
  public static String validateGO(QualifierVector qualifiers, EntryInformation eInfo)
  {
    final StringBuilder buff = new StringBuilder();
    final Qualifier q = qualifiers.getQualifierByName("GO");
    if (q != null)
    {
      try
      {
        final QualifierInfo qI = eInfo.getQualifierInfo("GO");
        final StringVector qualifierStrs = StreamQualifier.toStringVector(qI,q);
        for (String qualifierStr: qualifierStrs)
        {
          final String code = getEvidenceCodeAbbreviation(qualifierStr);
          final String goid = getField("GOid=", qualifierStr);
          final String with = getField("with=", qualifierStr);
          final String dbxref = getField("dbxref=", qualifierStr);

          // System.err.println(GeneUtils.getUniqueName(f)+" "+code+" "+qualifierStr);
          if (code != null)
          {
            if (code.equals("IDA") || code.equals("NAS") || code.equals("ND") || 
                code.equals("TAS") || code.equals("EXP"))
            {
              // with/from must be empty
              if (with.length() > 0)
                buff.append("GOid=" + goid
                  + ", the with/from must be empty when using " + code + "\n");
            }
            else if (code.equals("ISS") || code.equals("ISA") || 
                     code.equals("ISO") || code.equals("ISM"))
            {
              // with field must be filled
              if (with.length() == 0)
                buff.append("GOid=" + goid +
                  ", the with/from field must be filled when using " + code + "\n");
            }
          }

          // GO:0005515 can only have IPI evidence code
          if (goid.equals("GO:0005515") && (code == null || !code.equals("IPI")))
            buff.append("GOid=" + goid + ", can only have IPI evidence code\n");

          // IEP is not allowed for molecular_function and cellular_component terms
          if (code != null && code.equals("IEP")
              && !getField("aspect=", qualifierStr).equals("P"))
            buff.append("GOid=" + goid + ", IEP is restricted to Biological Process terms\n");

          if (with.indexOf(" ") > -1 || dbxref.indexOf(" ") > -1)
          {
            buff.append("GOid=" + goid + ", ");
            if (with.indexOf(" ") > -1)
              buff.append("with/from field (" + with + ") contains white space ");
            if (dbxref.indexOf(" ") > -1)
              buff.append("dbxref field (" + dbxref + ") contains white space");
            buff.append("\n");
          }
        }
      } catch(Exception e){}
    }

    return buff.toString();
  }
  
  /**
   * Extract the value of a field from a qualifier string
   * @param fieldName
   * @param qualifierString
   * @return
   */
  private static String getField(final String fieldName, final String qualifierStr)
  {
    String field = "";
    int ind1 = qualifierStr.toLowerCase().indexOf(fieldName.toLowerCase());
    int ind2 = qualifierStr.indexOf(";", ind1);

    int len = fieldName.length();
    if(ind2 > ind1 && ind1 > -1)
      field = qualifierStr.substring(ind1+len,ind2);
    else if(ind1 > -1)
      field = qualifierStr.substring(ind1+len);
    
    if(field.endsWith("\""))
      field = field.substring(0, field.length()-1);
    if(field.startsWith("\""))
      field = field.substring(1);
    return field;
  }
  
  private static String getEvidenceCodeAbbreviation(String goText)
  {
    final String evidence = getField("evidence=", goText); 
    for(int i=0; i<GoBox.evidenceCodes[2].length; i++)
      if(GoBox.evidenceCodes[2][i].equalsIgnoreCase(evidence))
        return GoBox.evidenceCodes[0][i];
    
    for(int i=0; i<GoBox.evidenceCodes[0].length; i++)
      if(GoBox.evidenceCodes[0][i].equalsIgnoreCase(evidence))
        return GoBox.evidenceCodes[0][i];
    return null;
  }
  
  /**
   * Check if the gene model is complete
   * @param gffFeature
   * @return 0 - complete gene model
   *         1 - no gene found
   *         2 - missing transcript
   */
  protected static int isCompleteGeneModelOK(final GFFStreamFeature gffFeature)
  {
    final ChadoCanonicalGene gene = gffFeature.getChadoGene();
    if(gene == null)
      return 1;
    
    final List<uk.ac.sanger.artemis.io.Feature> transcripts = gene.getTranscripts();
    if(transcripts.size() < 1)
      return 2;

    return 0;
  }
  
  protected static int isBoundaryOK(final GFFStreamFeature gffFeature)
  {
    final ChadoCanonicalGene gene = gffFeature.getChadoGene();
    int gb = 0;
    if(gene != null && isGene(gffFeature))
      gb = GeneUtils.isBoundaryOK(gene);
    return gb;
  }
  
  protected static boolean isStrandOK(final GFFStreamFeature gffFeature)
  {
    final ChadoCanonicalGene gene = gffFeature.getChadoGene();
    if(gene != null && isGene(gffFeature) && !GeneUtils.isStrandOK(gene))
      return false;
    return true;
  }
  
  /**
   * Test if the partial qualifiers are consistent within a gene model
   * @param gffFeature
   * @param partialKeyStr - either Start_range or End_range qualifier keys
   * @return
   */
  protected static boolean isPartialConsistent(final GFFStreamFeature gffFeature,
                                             final String partialKeyStr)
  {
    final ChadoCanonicalGene gene = gffFeature.getChadoGene();
    if(gene == null)
      return true;
    try
    {
      // is the gene marked as partial
      boolean partial = 
          (gene.getGene().getQualifierByName(partialKeyStr) == null ? false : true);
      
      if(partial)
      {
        if(!isPartialChildren(partial, gene, gene.getGene(), partialKeyStr))
          return false;
      }
      else
      {
        final List<uk.ac.sanger.artemis.io.Feature> transcripts = gene.getTranscripts();
        for(uk.ac.sanger.artemis.io.Feature f: transcripts)
        {
          // is the transcript marked as partial
          partial = (f.getQualifierByName(partialKeyStr) == null ? false : true);
          if(!isPartialChildren(partial, gene, f, partialKeyStr))
            return false;
        }
      }
    }
    catch (Exception e){ e.printStackTrace(); }
    return true;
  }
  
  /**
   * Return true if the child features agree with the parent 
   * feature for a given key
   * @param parentIsPartial
   * @param gene
   * @param f
   * @param partialKeyStr
   * @return
   * @throws InvalidRelationException
   */
  private static boolean isPartialChildren(
      final boolean parentIsPartial,
      final ChadoCanonicalGene gene,
      final uk.ac.sanger.artemis.io.Feature f,
      final String partialKeyStr) throws InvalidRelationException
  {
    final Set<uk.ac.sanger.artemis.io.Feature> children = gene.getChildren(f);
    for(uk.ac.sanger.artemis.io.Feature child: children)
    {
      String keyStr = child.getKey().getKeyString();
      boolean isComplement = f.getLocation().isComplement();
      if(keyStr.equals("three_prime_UTR"))
      {
        if(isComplement)
        {
          if(partialKeyStr.equals("End_range"))
            continue;
        }
        else
        {
          if(partialKeyStr.equals("Start_range"))
             continue;
        }
      }
      else if(keyStr.equals("five_prime_UTR"))
      {
        if(isComplement)
        {
          if(partialKeyStr.equals("Start_range"))
            continue;
        }
        else
        {
          if(partialKeyStr.equals("End_range"))
             continue;
        }
      }

      boolean partial = (child.getQualifierByName(partialKeyStr) == null ? false : true);
      if( (  partial && !parentIsPartial ) ||
          ( !partial &&  parentIsPartial ) )
        return false;
    }
    return true;
  }
  
  /**
   * Test if the ID GFF3 attribute prefix is consistent within a gene model
   * @param gffFeature
   * @return true if the prefix is the same within the gene model features
   */
  protected static boolean isIdPrefixConsistent(final GFFStreamFeature gffFeature)
  {
    final ChadoCanonicalGene gene = gffFeature.getChadoGene();
    if(gene == null)
      return true;

    try
    {
      if(gffFeature.getKey().getKeyString().endsWith("gene"))
        return (gene.getGene().getQualifierByName("ID") != null);

      if(gene.getGene().getQualifierByName("ID") == null)
        return true;
      if(gffFeature.getQualifierByName("ID") == null)
        return false;
      
      String id = gene.getGene().getQualifierByName("ID").getValues().elementAt(0);
      final Matcher m = ID_PREFIX.matcher(id);
      if(m.matches())
      {
        id = gffFeature.getQualifierByName("ID").getValues().elementAt(0);
        return id.startsWith( m.group() );
      }
    }
    catch (Exception e)
    {
      e.printStackTrace();
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
      return false;
    return true;
  }
 
  
  /**
   * Check attribute column
   * @param gffFeature
   */
  protected static String isAttributesOK(final GFFStreamFeature gffFeature)
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
        boolean allowed = false;
        
        if(GeneUtils.isDatabaseEntry(gffFeature))
        {
          for(String a: ALLOWED_TAGS_WITH_NO_VALUE)
            if(qualifier.getName().equals(a))
              allowed = true;
        }
        
        if(!allowed)
        {
          String msg = qualifier.getName()+" atribute has no value\n";
          str.append(msg);
        }
      }
    }  
    return str.toString();
  }
  
  protected boolean isInternalStops(final uk.ac.sanger.artemis.io.Feature feature)
  { 
    if(feature.getUserData() == null)
      return false;
    final uk.ac.sanger.artemis.Feature f = (uk.ac.sanger.artemis.Feature) feature.getUserData();
    final FeaturePredicate cds_predicate = getCodingFeaturePredicate();
    if(!cds_predicate.testPredicate (f)) 
      return false;

    if(feature instanceof GFFStreamFeature && GFFStreamFeature.isSelenocysteine(feature))
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

    if(entryGrp != null && GeneUtils.isGFFEntry( entryGrp ))
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
  
  protected static boolean isPartOfGene(final GFFStreamFeature gffFeature)
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
        if(feature.getQualifierByName("Start_range") != null)
          low_pos = "<"+low_pos;
        if(feature.getQualifierByName("End_range") != null)
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
  
  /**
   * Open a panel with the validation results
   * @param entry
   * @param seq
   */
  private void showReport(final Entry entry, final String seq)
  {
    final FeatureVector features = entry.getAllFeatures();
    final FileViewer fv = new FileViewer("Validation Report :: "+seq+" "+
                      features.size()+" feature(s)", false, false, true); 
    int nfail = 0;
    for(uk.ac.sanger.artemis.io.Feature f: features)
    {
      if(!featureValidate(f, fv, true))
        nfail++;
    }
    fv.setTitle(fv.getTitle()+" Fails:"+nfail);
    fv.setVisible(true);
  }
  
  /**
   * Write the validation results to a file
   * @param writer
   * @param entry
   * @param seq
   */
  private void writeReport(final BufferedWriter writer, final Entry entry, final String seq)
  {
    final FeatureVector features = entry.getAllFeatures();
    try
    {
      for(uk.ac.sanger.artemis.io.Feature f: features)
      {
        LinkedHashMap<String, Level> report = featureValidate(f, true);
        if (report != null)
        {
          for (Map.Entry<String, Level> ent : report.entrySet())
          {
            writer.append(ent.getKey());
            writer.newLine();
          }
        }
      }
    }
    catch(IOException e)
    {
      e.printStackTrace();
    }
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
  
  private static FeaturePredicate PARTIAL_PREDICATE = new FeaturePredicate() 
  {
    public boolean testPredicate(uk.ac.sanger.artemis.Feature feature)
    {
      if( isPartialConsistent((GFFStreamFeature) feature.getEmblFeature(), "Start_range") &&
          isPartialConsistent((GFFStreamFeature) feature.getEmblFeature(), "End_range") )
        return false;
      return true;
    }
  };
  
  private static FeaturePredicate ID_PREDICATE = new FeaturePredicate() 
  {
    public boolean testPredicate(uk.ac.sanger.artemis.Feature feature)
    {
      if( isIdPrefixConsistent((GFFStreamFeature) feature.getEmblFeature() ))
        return false;
      return true;
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
    if( (args != null && args.length == 1 && args[0].startsWith("-h")) ||
        (args == null || args.length < 1))
    {
      System.out.println("Artemis validation options:");
      System.out.println("-h\tshow help");
      System.out.println("-s\tspace separated list of sequences to read and write out");
      System.out.println("-c\tthe URL for your Chado database e.g. db.genedb.org:5432/snapshot?genedb_ro");
      System.out.println("-o\twrite report to specified file");
      System.exit(0);
    }

    BufferedWriter outfile = null;
    for(int i = 0; i < args.length; i++)
    {
      if (args[i].equalsIgnoreCase("-c"))
        System.setProperty("chado", args[i + 1]);
      else if (args[i].equalsIgnoreCase("-o"))
        try
        {
          outfile = new BufferedWriter(new FileWriter(args[i + 1]));
        }
        catch (IOException e)
        {
          JOptionPane.showMessageDialog(null, e.getMessage());
          e.printStackTrace();
          System.exit(1);
        }
    }

    final Vector<String> seqs = new Vector<String>();
    for(int i = 0; i < args.length; i++)
    {
      if(args[i].toLowerCase().equals("-s"))
      {
        for(int j = i + 1; j < args.length; j++)
        {
          if(args[j].startsWith("-"))
            break;
          seqs.add(args[j]);
        }
      }
      else if(args[i].startsWith("-"))
        i++;
      else
      {
        if(!seqs.contains(args[i]))
          seqs.add(args[i]);
      }
    }

    System.setProperty("black_belt_mode","true");
    Options.getOptions();
    
    if(System.getProperty("chado") != null)
    {
      DatabaseEntrySource entrySrc = null;
      try
      {
        for(String seq: seqs)
        {
          System.out.println("VALIDATING... "+seq);
          uk.ac.sanger.artemis.Entry entry = ReadAndWriteEntry.readEntryFromDatabase(seq, entrySrc);
          entrySrc = ReadAndWriteEntry.getEntrySource();
          final EntryGroup egrp = new SimpleEntryGroup();
          egrp.add(entry);

          ValidateFeature validate = new ValidateFeature(egrp);
          if(outfile == null)
            validate.showReport(entry.getEMBLEntry(), seq);
          else
            validate.writeReport(outfile, entry.getEMBLEntry(), seq);
        }
      }
      catch (Exception e)
      {
        e.printStackTrace();
        JOptionPane.showMessageDialog(null, e.getMessage());
      }
    }  
    else
    {
      for(String seq: seqs)
      {
        System.out.println("VALIDATING... "+seq);
        final Document doc = DocumentFactory.makeDocument(seq);
        final uk.ac.sanger.artemis.io.Entry entry = EntryFileDialog.getEntryFromFile(
            null, doc, Options.getArtemisEntryInformation(), false);
        
        ValidateFeature validate = new ValidateFeature(null);
        if(outfile == null)
          validate.showReport(entry, doc.getName());
        else
          validate.writeReport(outfile, entry, doc.getName());
      }
    }
    
    if(outfile != null)
      try
      {
        outfile.close();
      }
      catch (IOException e){}
  }

}