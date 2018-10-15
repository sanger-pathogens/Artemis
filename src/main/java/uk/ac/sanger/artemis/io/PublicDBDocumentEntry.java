/* PublicDBDocumentEntry.java
 *
 * created: Sat Sep 11 1999
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1999  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/PublicDBDocumentEntry.java,v 1.22 2008-12-15 14:07:57 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.FeatureSegmentVector;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.components.genebuilder.cv.GoBox;
import uk.ac.sanger.artemis.sequence.AminoAcidSequence;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.util.*;

import java.io.IOException;
import java.io.InputStream;
import java.util.Enumeration;
import java.util.List;
import java.util.Properties;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *  This class extends the Entry class with the data for the entry coming from
 *  a Document object.  The Document must contain an EMBL entry or a GENBANK
 *  entry.
 *
 *  @author Kim Rutherford
 *  @version $Id: PublicDBDocumentEntry.java,v 1.22 2008-12-15 14:07:57 tjc Exp $
 **/

public class PublicDBDocumentEntry extends SimpleDocumentEntry
    implements DocumentEntry 
{
  // database keys mapping to EMBL keys with any extra
  // qualifiers to add in.
  private static Object[][] DATABASE_MAP_KEYS;
  
  // database qualifiers mapping to EMBL qualifiers
  private static String[][] DATABASE_QUALIFIERS_TO_MAP;
  
  // database qualifiers to ignore when converting to embl
  private static Object[] DATABASE_QUALIFIERS_TO_REMOVE;
  
  public static boolean IGNORE_OBSOLETE_FEATURES = true;
  
  /**
   *  Create a new PublicDBDocumentEntry object associated with the given
   *  Document.
   *  @param entry_information The EntryInformation object of the new Entry.
   *  @param document This is the file that we will read from.  This is also
   *    used for saving the entry back to the file it came from and to give
   *    the new object a name.
   *  @param listener The object that will listen for ReadEvents.
   *  @exception IOException thrown if there is a problem reading the entry -
   *    most likely ReadFormatException.
   *  @exception EntryInformationException Thrown if force is false and if this
   *    Entry cannot contain the Key, Qualifier or Key/Qualifier combination of
   *    one of the features in the given Entry.
   **/
  PublicDBDocumentEntry(final EntryInformation entry_information,
                        final Document document, final ReadListener listener)
      throws IOException, EntryInformationException 
  {
    super(entry_information, document, listener);
  }

  /**
   *  Create a new PublicDBDocumentEntry that will be a copy of the given
   *  Entry and has no Document associated with it.  The new
   *  PublicDBDocumentEntry cannot be saved to a file with save () unless save
   *  (Document) has been called first.
   *  @param entry_information The EntryInformation object of the new Entry.
   *  @param force If true then invalid qualifiers and any features with
   *    invalid keys in the new Entry will be quietly thrown away.  "Invalid"
   *    means that the key/qualifier is not allowed to occur in an Entry of
   *    this type (probably determined by the EntryInformation object of this
   *    Entry).  If false an EntryInformationException will be thrown for
   *    invalid keys or qualifiers.
   *  @exception EntryInformationException Thrown if force is false and if this
   *    Entry cannot contain the Key, Qualifier or Key/Qualifier combination of
   *    one of the features in the given Entry.
   **/
  public PublicDBDocumentEntry(final EntryInformation entry_information,
                               final Entry new_entry, final boolean force)
      throws EntryInformationException 
  {
    super(entry_information, new_entry, force);
  }

  /**
   *  Create a new empty PublicDBDocumentEntry object that has no Document
   *  associated with it.  The new PublicDBDocumentEntry cannot be saved to a
   *  file with save () unless save (Document) has been called first.  The
   *  save (Document) method will assign a Document.
   *  @param entry_information The EntryInformation object of the Entry that
   *    will contain this Feature.
   **/
  public PublicDBDocumentEntry(final EntryInformation entry_information) 
  {
    super(entry_information);
  }

  /**
   *  If the given feature can be added directly to this Entry, then return
   *  it, otherwise create and return a new feature of the appropriate type.
   *  @param copy if true then always new a new copy of the Feature.
   **/
  protected Object makeNativeFeature(final Feature feature,
                                     final boolean copy) 
  {
    if (!copy && (feature instanceof EmblStreamFeature &&
                  this instanceof EmblDocumentEntry ||
                  feature instanceof GenbankStreamFeature &&
                  this instanceof GenbankDocumentEntry)) 
    {
      return (PublicDBStreamFeature) feature;
    } 
    else 
    {
      try
      {
        if(feature instanceof GFFStreamFeature)
          return mapGffToNativeFeature(feature);
        else if (this instanceof EmblDocumentEntry)
          return new EmblStreamFeature(feature);
        else
          return new GenbankStreamFeature(feature);
      }
      catch(NullPointerException npe)
      {
        System.err.println( 
            ((uk.ac.sanger.artemis.Feature)feature.getUserData()).getIDString() );
        throw npe;
      }
    }
  }
 

  /**
   * Map GFF features to EMBL/Genbank
   * @param feature
   * @return
   */
  private Object mapGffToNativeFeature(final Feature feature)
  {
    if(DATABASE_MAP_KEYS == null)
      initDatabaseMappings();
    
    Key key = feature.getKey();
    QualifierVector qualifiers = feature.getQualifiers().copy();
    
    // ignore if obsolete
    if(IGNORE_OBSOLETE_FEATURES)
    {
      Qualifier isObsoleteQualifier = qualifiers.getQualifierByName("isObsolete");
      if(isObsoleteQualifier != null)
      {
        String value = (String)isObsoleteQualifier.getValues().get(0);
        if(Boolean.parseBoolean(value))
          return null;
      }
    }
    
    key = map(key, qualifiers);
    if(getEntryInformation().isValidQualifier((String) DATABASE_QUALIFIERS_TO_REMOVE[0]))
    { 
      try
      {
        if(this instanceof EmblDocumentEntry)
          return new EmblStreamFeature (
              key, 
              feature.getLocation(), 
              qualifiers);
        else
          return new GenbankStreamFeature (
            key, 
            feature.getLocation(), 
            qualifiers);
      }
      catch(InvalidRelationException e)
      {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
    }
    
    Location location = joinUtrs(feature, key, qualifiers);
    if(location == null)
      return null;
    // flatten gene model - combining qualifiers
    if(key.getKeyString().equals(DatabaseDocument.EXONMODEL) &&
        ((GFFStreamFeature)feature).getChadoGene() != null)
    {
      ChadoCanonicalGene chadoGene = ((GFFStreamFeature)feature).getChadoGene();

      final String name = GeneUtils.getUniqueName(feature);
      final String transcriptName = chadoGene.getTranscriptFromName(name);
      
      StringVector sv = new StringVector();
      sv.add(transcriptName);
      final Feature transcript = chadoGene.containsTranscript(sv);
      
      if(transcript != null && GeneUtils.isNonCodingTranscripts(transcript.getKey()))
        return null;
      
      qualifiers.removeQualifierByName("ID");
      int ntranscripts = 0;
      // add transcript & protein qualifiers to CDS
      try
      {
        final Feature protein = chadoGene.getProteinOfTranscript(transcriptName);
        if(protein != null)
          combineQualifiers(qualifiers, protein.getQualifiers().copy(), false);

        if(transcript != null)
          ntranscripts = handleTranscripts(qualifiers, transcript, ntranscripts, chadoGene);
      }
      catch(NullPointerException npe){}
      
      // add gene qualifiers to CDS
      QualifierVector geneQualifiers = chadoGene.getGene().getQualifiers().copy();
      
      // multiple transcripts
      if(ntranscripts > 1 && geneQualifiers.getQualifierByName("ID") != null)
      {
      	Qualifier newIDQualifier = 
        	new Qualifier("shared_id", (String)geneQualifiers.getQualifierByName("ID").getValues().get(0));
    	addNewQualifier(qualifiers, newIDQualifier);
    	geneQualifiers.removeQualifierByName("ID");
      }
      else if(ntranscripts == 1)
      {
        if(qualifiers.getQualifierByName("Start_range") != null)
          addNewQualifier(qualifiers, qualifiers.getQualifierByName("Start_range"));
        if(qualifiers.getQualifierByName("End_range") != null)
          addNewQualifier(qualifiers, qualifiers.getQualifierByName("End_range"));
      }
      combineQualifiers(qualifiers, geneQualifiers, true);
    }
    else if(GeneUtils.isNonCodingTranscripts(key))
    {
      // use gene id for non-coding transcripts
      ChadoCanonicalGene chadoGene = ((GFFStreamFeature)feature).getChadoGene();
      if(chadoGene != null)
      {
        qualifiers.removeQualifierByName("ID");
        QualifierVector geneQualifiers = chadoGene.getGene().getQualifiers().copy();
        combineQualifiers(qualifiers, geneQualifiers, true);
        
        if(qualifiers.getQualifierByName("product") != null)
        {
          Qualifier newQualifier = new Qualifier("product", processProductValues(qualifiers.getQualifierByName("product")));
          qualifiers.setQualifier(newQualifier);
        }
      }
    }
    
    try
    {
      location = handlePartials(qualifiers, location);
      for(int i=0; i<DATABASE_QUALIFIERS_TO_MAP.length; i++)
      {
        if(!getEntryInformation().isValidQualifier(DATABASE_QUALIFIERS_TO_MAP[i][0]))
        {
          changeQualifierName(qualifiers, DATABASE_QUALIFIERS_TO_MAP[i][0], 
                                          DATABASE_QUALIFIERS_TO_MAP[i][1]);
        }
      }
      
      if(qualifiers.getQualifierByName("stop_codon_redefined_as_selenocysteine") != null)
      {
        handleSelenocysteine(qualifiers, feature);
      }
      
      for(int i=0; i<DATABASE_QUALIFIERS_TO_REMOVE.length; i++)
      {
        if(!getEntryInformation().isValidQualifier((String) DATABASE_QUALIFIERS_TO_REMOVE[i]))
          qualifiers.removeQualifierByName((String) DATABASE_QUALIFIERS_TO_REMOVE[i]);
      }
     
      if(key.getKeyString().equals("polypeptide"))
        return null;
      else if(key.getKeyString().equals("gene"))
        return null;
      else if(key.getKeyString().equals("centromere"))
        return null;
      else if(key.getKeyString().equals("transcript") || 
              key.getKeyString().equals("mRNA"))
        return null;
      
      
      if(this instanceof EmblDocumentEntry)
        return new EmblStreamFeature (
            key, location, 
            qualifiers);
      else
        return new GenbankStreamFeature (
            key, location, 
            qualifiers);
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
      if(feature instanceof DatabaseStreamFeature)
        return new EmblStreamFeature ();
      else
        return new GenbankStreamFeature ();
    }
  }
  
  /**
   * Handle UTR joins
   * @param feature
   * @param key
   * @param qualifiers
   * @return
   */
  private Location joinUtrs(Feature feature, Key key, QualifierVector qualifiers)
  {
    Location location = feature.getLocation();
    if(key.getKeyString().equals("5'UTR") ||
       key.getKeyString().equals("3'UTR"))
    {
      ChadoCanonicalGene gene = ((GFFStreamFeature)feature).getChadoGene();
      String utrName = GeneUtils.getUniqueName(feature);
      String transcriptName = gene.getTranscriptFromName(utrName);
      List<Feature> utrs;
      
      if(key.getKeyString().equals("5'UTR"))
        utrs = gene.get5UtrOfTranscript(transcriptName);
      else
        utrs = gene.get3UtrOfTranscript(transcriptName);
      
      if(utrs.size() > 1)
      {
        int start = Integer.MAX_VALUE;
        RangeVector ranges = new RangeVector();
        for(int i=0; i<utrs.size(); i++)
        {
          Feature utr = utrs.get(i);
          Range range = utr.getLocation().getTotalRange();
          if(start > range.getStart())
            start = range.getStart();
          ranges.add(range);
        }
        
        if(start != feature.getLocation().getTotalRange().getStart())
          return null;
        
        if(feature.getLocation().isComplement())
          ranges.reverse();
        location =
          new Location(ranges, feature.getLocation().isComplement());
      }

      int ntranscripts = gene.getTranscripts().size();
      if(ntranscripts == 1)
    	  transcriptName = gene.getGeneUniqueName();
      qualifiers.setQualifier(new Qualifier("locus_tag", transcriptName));
      qualifiers.removeQualifierByName("ID");
    }
    return location;
  }

  /**
   * Merge qualifiers
   * @param qualifiers
   * @param newQualifiers
   */
  private void combineQualifiers(final QualifierVector qualifiers,
                                 final QualifierVector newQualifiers,
                                 final boolean isGene)
  {
    for(int i=0; i<newQualifiers.size(); i++)
    {
      Qualifier newQualifier = (Qualifier) newQualifiers.get(i);
      
      if( newQualifier.getName().equals("ID") && !isGene )
      {
    	continue;
      }
      
      // convert GO evidence to codes (e.g. ND=No biological Data available)
      if(newQualifier.getName().equals("GO"))
      {
        final StringVector newValues = newQualifier.getValues();
        final StringVector tmpNewValues = new StringVector();
        for(int j=0; j<newValues.size(); j++)
        {
          String val = GoBox.getEvidenceCodeGoTextFromText((String)newValues.get(j));
          tmpNewValues.add(val);
        }
 
        newQualifier = new Qualifier("GO", tmpNewValues);
      }
      
      if(newQualifier.getName().equals("product"))
      {
        newQualifier = new Qualifier("product", processProductValues(newQualifier));
      }
      
      if(newQualifier.getName().equals("orthologous_to") ||
         newQualifier.getName().equals("paralogous_to"))
      {
        final StringVector newValues = newQualifier.getValues();
        final StringVector tmpNewValues = new StringVector();
        for(int j=0; j<newValues.size(); j++)
        {
          if(!newValues.get(j).equals(""))
            tmpNewValues.add(newValues.get(j));
        }
        if(tmpNewValues.size() == 0)
          continue;
        
        Pattern p = Pattern.compile("\\w+[ :]link=\\w+");
        for(int j=0; j<tmpNewValues.size(); j++)
        {
          String valueStr = (String)tmpNewValues.get(j);
          String newValueStr;
          int indexEnd = valueStr.indexOf(';');
          String endStr = "";
          if(indexEnd > -1)
            endStr = valueStr.substring(indexEnd);
          Matcher m = p.matcher(valueStr);
          while(m.find())
          {
            int index = valueStr.indexOf("link=", m.start());
            newValueStr = valueStr.substring(m.start(), index)+
                          valueStr.substring(index+5, m.end())+endStr;
            if(newQualifier.getName().equals("orthologous_to"))
              newQualifier = new Qualifier("orthologous_to", newValueStr);
            else
              newQualifier = new Qualifier("paralogous_to", newValueStr);
            qualifiers.addElement(newQualifier);
          }
        }
        continue;
      }
      
      addNewQualifier(qualifiers, newQualifier);
    }
  }
  
  /**
   * Process product qualifier
   * @param productQualifier
   * @return
   */
  private StringVector processProductValues(Qualifier productQualifier)
  {
    final StringVector values = productQualifier.getValues();
    final StringVector tmpNewValues = new StringVector();
    for(String val: values)
    {
      int ind = 0;
      if((ind=val.indexOf(";db_xref="))>-1 && this instanceof EmblDocumentEntry)
        val = val.substring(0,ind);

      if((ind=val.indexOf(";evidence="))>-1 && this instanceof EmblDocumentEntry)
        val = val.substring(0,ind);

      if((ind=val.indexOf(";with="))>-1 && this instanceof EmblDocumentEntry)
        val = val.substring(0,ind);

      if((ind=val.indexOf("rank="))>-1 && this instanceof EmblDocumentEntry)
      {
        int ind2 = val.indexOf(";", ind);
        if(ind == 0 && ind2 > -1)
          val = val.substring(ind2+1);
      }

      if(val.startsWith("term="))
        val = val.substring(5,  val.length());

      if(val.endsWith(";"))
        val = val.substring(0,  val.length()-1);

      tmpNewValues.add(val);
    }
    return tmpNewValues;
  }
  
  /**
   * Add a new qualifier to a list of qualifiers
   * @param qualifiers
   * @param newQualifier
   */
  private void addNewQualifier(QualifierVector qualifiers, Qualifier newQualifier)
  {
	Qualifier qualifier;
    if( (qualifier = qualifiers.getQualifierByName(newQualifier.getName())) != null)
    { 
      final StringVector newValues = newQualifier.getValues();
      final StringVector values = qualifier.getValues();
      
      if(newValues == null)
        return;
      for(int j=0; j<newValues.size(); j++)
      {
        String newValue = (String) newValues.get(j);
        if(!values.contains(newValue))
          qualifier.addValue(newValue);
      }
    }
    else
      qualifiers.addElement(newQualifier);
  }

  /**
   * Routine to combine transcript qualifiers and for multiple transcripts 
   * create links to the other transcripts (other_transcript) and 
   * to use the transcript ID.
   * @param qualifiers
   * @param transcript
   * @param ntranscripts
   * @param chadoGene
   */
  private int handleTranscripts(QualifierVector qualifiers, 
		                        Feature transcript,
		                        int ntranscripts,
		                        ChadoCanonicalGene chadoGene)
  {
	QualifierVector transcriptQualifiers = transcript.getQualifiers().copy();
    combineQualifiers(qualifiers, transcriptQualifiers, false);
    ntranscripts = chadoGene.getTranscripts().size();
    if(ntranscripts > 1)
    {
      addNewQualifier(qualifiers, transcriptQualifiers.getQualifierByName("ID"));
      List<Feature> transcripts = chadoGene.getTranscripts();
      for(int i=0;i<ntranscripts;i++)
      {
        Feature thisTranscript = (Feature)transcripts.get(i);
        String thisTranscriptName = GeneUtils.getUniqueName(thisTranscript);
        if(!thisTranscriptName.equals( GeneUtils.getUniqueName(transcript )))
        {
          Qualifier qualifier = new Qualifier("other_transcript", thisTranscriptName);
          addNewQualifier(qualifiers, qualifier);
        }
      }
    }
    return ntranscripts;
  }
  
  /**
   * Maps database (SO) keys to EMBl keys.
   * @param key
   * @return
   */
  protected static Key mapKeys(final Key key)
  {
    if(DATABASE_MAP_KEYS == null)
      initDatabaseMappings();
    for(int i=0; i<DATABASE_MAP_KEYS.length; i++)
    {
      if(key.getKeyString().equals(DATABASE_MAP_KEYS[i][0]))
        return new Key((String)DATABASE_MAP_KEYS[i][1]);
    }
    return key;
  }
  
  
  /**
   * Maps database (SO) keys to EMBl keys. It will add any extra qualifiers
   * found in the 3rd column of DATABASE_MAP_KEYS.
   * @param key
   * @return
   */
  private Key map(final Key key, final QualifierVector qualifiers)
  {
    if(DATABASE_MAP_KEYS == null)
      initDatabaseMappings();
    for(int i=0; i<DATABASE_MAP_KEYS.length; i++)
    {
      if(key.getKeyString().equals(DATABASE_MAP_KEYS[i][0]))
      {
        Key mappedKey = new Key((String)DATABASE_MAP_KEYS[i][1]);
        if(DATABASE_MAP_KEYS[i][2] != null)
        {
          
          Qualifier newQualifier = (Qualifier) DATABASE_MAP_KEYS[i][2];
          if(!getEntryInformation().isValidQualifier(mappedKey, newQualifier.getName()))
          {
            try
            {
              final int nvalues;
              if(newQualifier.getValues() == null)
                nvalues = 0;
              else
                nvalues = newQualifier.getValues().size();
              
              final int type;
              if(nvalues == 0)
                type = QualifierInfo.NO_VALUE;
              else
                type = QualifierInfo.QUOTED_TEXT;
              getEntryInformation().addQualifierInfo(
                  new QualifierInfo(newQualifier.getName(), 
                      type, null, null,
                      false));
            }
            catch(QualifierInfoException e){}
          }
          
          qualifiers.addQualifierValues(newQualifier);
        }
        
        return mappedKey;
      }
    }
    return key;
  }
  
  /**
   * Change the name of a qualifier
   * @param qualifiers
   * @param oldName
   * @param newName
   */
  private void changeQualifierName(QualifierVector qualifiers, 
                                   final String oldName, 
                                   final String newName)
  {
    QualifierVector tmpQualifiers = new QualifierVector();
    for(int i=0; i<qualifiers.size(); i++)
    {
      Qualifier qualifier = (Qualifier) qualifiers.elementAt(i);
      if(!qualifier.getName().equals(oldName))
      {
        tmpQualifiers.addElement(qualifier);
        continue;
      }

      Qualifier newQualifier = new Qualifier(newName, qualifier.getValues());
      tmpQualifiers.addQualifierValues(newQualifier);
    }
    qualifiers.removeAllElements();
    
    for(int i=0; i<tmpQualifiers.size(); i++)
      qualifiers.addElement(tmpQualifiers.elementAt(i));
  }
  
  /**
   *  If the given Sequence can be added directly to this Entry, then return a
   *  copy of it, otherwise create and return a new feature of the appropriate
   *  type for this Entry.
   **/
  protected StreamSequence makeNativeSequence (final Sequence sequence) 
  {
    if(this instanceof EmblDocumentEntry)
      return new EmblStreamSequence (sequence);
    else 
      return new GenbankStreamSequence (sequence);
  }
  
  public static Object[] getDatabaseQualifiersToRemove()  
  {
    initDatabaseMappings();
    return DATABASE_QUALIFIERS_TO_REMOVE;
  }
  
  /**
   *  Read key and qualifier mappings for CHADO to EMBL
   **/
  private static void initDatabaseMappings()
  {
    InputStream keyStream = 
      Options.class.getResourceAsStream("/key_mapping");
    if (keyStream == null)
      keyStream = 
        Options.class.getResourceAsStream("/etc/key_mapping");

    InputStream qualifierStream = 
      Options.class.getResourceAsStream("/qualifier_mapping");
    if (qualifierStream == null)
      qualifierStream = 
        Options.class.getResourceAsStream("/etc/qualifier_mapping");

    final Properties keyMapProperties = new Properties();
    final Properties qualifierMapProperties = new Properties();
    try
    {
      keyMapProperties.load(keyStream);
      qualifierMapProperties.load(qualifierStream);
      
      if(System.getProperty("nohistory") != null)
        qualifierMapProperties.setProperty("history", "");
    }
    catch(IOException e)
    {
      e.printStackTrace();
    }

    // parse the keyMapProperties
    DATABASE_MAP_KEYS = new Object[keyMapProperties.size()][3];
    final Enumeration keysenum = keyMapProperties.propertyNames();
    int n = 0;
    while(keysenum.hasMoreElements()) 
    {
      String current_map_name = (String) keysenum.nextElement();

      final StringVector property_values =
          Options.getPropertyValues(keyMapProperties, current_map_name);

      DATABASE_MAP_KEYS[n][0] = current_map_name;
      DATABASE_MAP_KEYS[n][1] = property_values.get(0);
      if(property_values.size() == 2)
      {
        String qualifierString[] = ((String)property_values.get(1)).split("=");
        final uk.ac.sanger.artemis.io.Qualifier qualifier;
        if(qualifierString.length == 2)
          qualifier = new uk.ac.sanger.artemis.io.Qualifier(qualifierString[0], qualifierString[1]);
        else
          qualifier = new uk.ac.sanger.artemis.io.Qualifier(qualifierString[0]);
        DATABASE_MAP_KEYS[n][2] = qualifier; 
      }
      else
        DATABASE_MAP_KEYS[n][2] = null;
      n++; 
    }
    
    // parse the qualifier mappings
    Enumeration qualifiersenum = qualifierMapProperties.propertyNames();
    n = 0;

    Vector<String> qualifiersToRemove = new Vector<String>();
    while(qualifiersenum.hasMoreElements()) 
    {
      String current_map_name = (String) qualifiersenum.nextElement();
      final StringVector property_values =
        Options.getPropertyValues(qualifierMapProperties, current_map_name);
      if(property_values == null || property_values.size() == 0)
        qualifiersToRemove.add(current_map_name);
      else
        n++;
    }
    
    DATABASE_QUALIFIERS_TO_MAP = new String[n][2];
    DATABASE_QUALIFIERS_TO_REMOVE = qualifiersToRemove.toArray();
    
    qualifiersenum = qualifierMapProperties.propertyNames();
    n = 0;

    while(qualifiersenum.hasMoreElements()) 
    {
      String current_map_name = (String) qualifiersenum.nextElement();
      final StringVector property_values =
        Options.getPropertyValues(qualifierMapProperties, current_map_name);
      if(property_values != null && property_values.size() > 0)
      {
        DATABASE_QUALIFIERS_TO_MAP[n][0] = current_map_name;
        DATABASE_QUALIFIERS_TO_MAP[n][1] = (String) property_values.get(0);
        n++;
      }
    }
  }
  
  /**
   * Use '<' and '>' signs in the location descriptors to
   * indicate that the sequence is partial.
   * @param qualifiers
   * @param location
   * @return
   */
  private Location handlePartials(QualifierVector qualifiers, Location location)
  {
    if(qualifiers.getQualifierByName("Start_range") != null)
    {
      try
      {
        location = new Location(location.toStringShort().replaceFirst("(\\d)", "<$1"));
      }
      catch (LocationParseException e)
      {
        e.printStackTrace();
      }
    }
    if(qualifiers.getQualifierByName("End_range") != null)
    {
      try
      {
        location = new Location(location.toStringShort().replaceAll("^(.*)(\\.)(.*)$","$1$2>$3"));
      }
      catch (LocationParseException e)
      {
        e.printStackTrace();
      }
    }
    return location;
  }
  
  /**
   * Change the stop_codon_redefined_as_selenocysteine SO qualifier
   * to the transl_except EMBL qualifier.
   * @param qualifiers
   * @param feature
   */
  private void handleSelenocysteine(QualifierVector qualifiers, Feature feature)
  {
    if(!feature.getKey().getKeyString().equals(DatabaseDocument.EXONMODEL))
      return;
    qualifiers.removeQualifierByName("stop_codon_redefined_as_selenocysteine");

    uk.ac.sanger.artemis.Feature f = ((uk.ac.sanger.artemis.Feature)feature.getUserData());
    
    int translatedBasePosion = 0;
    String aa = f.getTranslation().toString();
    for(int i=0; i<aa.length(); i++)
    {
      if(AminoAcidSequence.isStopCodon(aa.charAt(i)))
      {
        translatedBasePosion = i*3;
        break;
      }
    }

    FeatureSegmentVector segments = f.getSegments();
    int nbases = 0;
    int sequenceloc = 0;
    for(int i=0; i<segments.size(); i++)
    {
      int seglen = segments.elementAt(i).getBases().length();
      if(nbases+seglen > translatedBasePosion && sequenceloc == 0)
      {
        Bases bases = f.getStrand().getBases();
        sequenceloc = segments.elementAt(i).getStart().getPosition() +
                              (translatedBasePosion-nbases);
        
        if(!f.isForwardFeature())
          sequenceloc = bases.getComplementPosition(sequenceloc);
      }
      nbases += seglen;
    }
    
    String pos;
    if(f.isForwardFeature())
      pos = sequenceloc+".."+(sequenceloc+2);
    else
      pos = sequenceloc+".."+(sequenceloc-2);
    
    qualifiers.add(new Qualifier("transl_except","(pos:"+pos+",aa:Sec)"));
  }

}
