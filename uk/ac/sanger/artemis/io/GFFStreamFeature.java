/* GFFStreamFeature.java
 *
 * created: Tue Sep 14 1999
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1999,2000,2001  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/GFFStreamFeature.java,v 1.72 2009-08-28 10:33:12 tjc Exp $
 */

package uk.ac.sanger.artemis.io;


import java.util.Hashtable;
import java.util.Enumeration;
import java.util.List;
import java.util.StringTokenizer;
import java.util.Vector;
import java.io.IOException;
import java.io.Writer;
import java.sql.Timestamp;
import java.text.SimpleDateFormat;

import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.chado.ClusterLazyQualifierValue;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.components.genebuilder.ProteinMapPanel;
import uk.ac.sanger.artemis.components.genebuilder.ortholog.MatchPanel;
import uk.ac.sanger.artemis.util.LinePushBackReader;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.util.StringVector;


/**
 *  A StreamFeature that thinks it is a GFF feature.
 *
 *  @author Kim Rutherford
 *  @version $Id: GFFStreamFeature.java,v 1.72 2009-08-28 10:33:12 tjc Exp $
 **/

public class GFFStreamFeature extends SimpleDocumentFeature
                       implements DocumentFeature, StreamFeature, ComparableFeature 
{

  private static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(GFFStreamFeature.class);
  
  /**
   *  This is the line of GFF input that was read to get this
   *  GFFStreamFeature.  A GFFStreamFeature that was created from multiple GFF
   *  lines will have a gff_lines variable that contains multiple line.
   **/
  StringVector gff_lines = null;

  /** store for spliced features containing id and range of each segment */
  private Hashtable id_range_store;
  
  /** store a record of the new and old uniquenames that have been changed */
  private Hashtable newIdMapToOldId;

  /** store the Timestamp for the feature */
  private Timestamp timelastmodified;
  
  private ChadoCanonicalGene chadoGene;
  
  private boolean visible = true;
  
  /** combined feature_relationship.rank store for exons */
  private Hashtable feature_relationship_rank_store;
  
  /** first tabbed parameter  */
  private String gffSeqName;
  /** second tabbed parameter */
  private String gffSource;
  /** duplication count */
  private short duplicate = 0;
  private boolean lazyLoaded = false;
  private org.gmod.schema.sequence.Feature chadoLazyFeature;
  private boolean readOnlyFeature = false;
  

  private static String MAP_DECODE[][] = {
    { " ",  "%20" },  // white space
    { ",",  "%2C" },  // comma
    { ";",  "%3B" },  // semi-colon
    { "=",  "%3D" },  // equals
    { "\t", "%09" },  // tab
    { " ",  "+"   },  // white space
    { "+",  "%2B" },
    { "(",  "%28" },  // left bracket
    { ")",  "%29" },  // right bracket
    { "'", "\"" }
  };
  
  private static String MAP_ENCODE[][] = {
//  { " ",  "%20" },  // white space
    { ",",  "%2C" },  // comma 
    { ";",  "%3B" },  // semi-colon
    { "=",  "%3D" },  // equals
    { "\t", "%09" },  // tab
    { "+",  "%2B" },
    { " ",  "+"   },  // white space
    { "(",  "%28" },  // left bracket
    { ")",  "%29" },  // right bracket
    { "\n", "%5C" }   // new-line 
  };
  
  /**
   *  Create a new GFFStreamFeature object.  The feature should be added
   *  to an Entry (with Entry.add()).
   *  @param key The new feature key
   *  @param location The Location object for the new feature
   *  @param qualifiers The qualifiers for the new feature
   **/
  public GFFStreamFeature(final Key key, final Location location,
                          final QualifierVector qualifiers) 
  {
    super(null);
    
    try 
    {
      setKey(key);
      setLocation(location);
      setQualifiers(qualifiers);
      
      /*
      if(getQualifierByName("score") == null)
        setQualifier(new Qualifier("score", "."));
      
      if(getQualifierByName("gff_source") == null)
        setQualifier(new Qualifier("gff_source", "artemis"));
      
      if(getQualifierByName("gff_seqname") == null)
        setQualifier(new Qualifier("gff_seqname", "."));
      */
      if(getQualifierByName("ID") == null)
      {
        String idStr = null;
        StringVector v = Options.getOptions().getSystematicQualifierNames();
        for(int i=0; i<v.size(); i++)
        {
          final String sysName = (String)v.get(i);
          if(getQualifierByName(sysName) != null)
          {
            idStr = (String)getQualifierByName(sysName).getValues().get(0);
            break;
          }
        }
        // autogenerate ID
        if(idStr == null)
          idStr = key.getKeyString()+":"+location.toString();
        setQualifier(new Qualifier("ID", idStr));
      }
      
    } 
    catch(EntryInformationException e) 
    {
      // this should never happen because the feature will not be in an Entry
      throw new Error("internal error - unexpected exception: " + e);
    }
    catch(ReadOnlyException e) 
    {
      // this should never happen because the feature will not be in an Entry
      throw new Error("internal error - unexpected exception: " + e);
    } 
    catch(OutOfRangeException e) 
    {
      // this should never happen because the feature will not be in an Entry
      throw new Error("internal error - unexpected exception: " + e);
    }
  }

  public GFFStreamFeature(final Feature feature)
  {
    this(feature, false);
  }
  
  /**
   *  Create a new GFFStreamFeature with the same key, location and
   *  qualifiers as the given feature.  The feature should be added to an
   *  Entry (with Entry.add()).
   *  @param feature The feature to copy.
   **/
  public GFFStreamFeature(final Feature feature, final boolean isDuplicatedInChado) 
  {
    this(feature.getKey(), feature.getLocation(), feature.getQualifiers());
    
    if(feature instanceof GFFStreamFeature)
    {
      if(((GFFStreamFeature)feature).id_range_store != null)
        this.id_range_store = 
          (Hashtable)(((GFFStreamFeature)feature).id_range_store).clone();
      
      if(((GFFStreamFeature)feature).feature_relationship_rank_store != null)
        this.feature_relationship_rank_store = 
          (Hashtable)(((GFFStreamFeature)feature).feature_relationship_rank_store).clone();
      
      this.setGffSeqName(((GFFStreamFeature)feature).getGffSeqName());
      this.setGffSource(((GFFStreamFeature)feature).getGffSource());
      
      
      if(isDuplicatedInChado)
      {
        try
        {
          final String uniquename;
          final String duplicatePrefix;
          
          if(feature instanceof GFFStreamFeature)
          {
            ((GFFStreamFeature)feature).duplicate++;
            duplicatePrefix = "DUP"+Short.toString(((GFFStreamFeature)feature).duplicate)+"-";
          }
          else
            duplicatePrefix = "DUP";
          if(id_range_store != null)
          {
            final Hashtable new_id_range_store = new Hashtable(id_range_store.size());
            final Enumeration enumIdRangeStore = id_range_store.keys();
            while(enumIdRangeStore.hasMoreElements())
            {
              final String keyId = (String)enumIdRangeStore.nextElement();
              final Range range  = (Range)id_range_store.get(keyId);
              new_id_range_store.put(duplicatePrefix+keyId, range);
            }
            id_range_store.clear();
            this.id_range_store = (Hashtable) new_id_range_store.clone();
            
            uniquename = getSegmentID(getLocation().getRanges());
          }
          else
            uniquename = duplicatePrefix+ (String)getQualifierByName("ID").getValues().get(0);
          setQualifier(new Qualifier("ID", uniquename));
          
          if(getQualifierByName("Parent") != null)
          {
            final String parent =
              (String) getQualifierByName("Parent").getValues().get(0);
            setQualifier(new Qualifier("Parent", duplicatePrefix+parent));
          }
          
          if(getQualifierByName("Derives_from") != null)
          {
            final String derives_from =
              (String) getQualifierByName("Derives_from").getValues().get(0);
            setQualifier(new Qualifier("Derives_from", duplicatePrefix+derives_from));
          }
          
          // remove qualifiers that don't get transferred to duplicate
          final String removeQualifierNames[] = 
            {  "feature_id", 
          		 "timelastmodified", 
          		 "feature_relationship_rank", 
          		 ProteinMapPanel.POLYPEPTIDE_DOMAIN,
          		 ProteinMapPanel.TMHMM[0],
          		 ProteinMapPanel.TMHMM[1],
          		 ProteinMapPanel.TMHMM[2],
          		 ProteinMapPanel.TMHMM[3],
          		 MatchPanel.ORTHOLOG,
          		 MatchPanel.ORTHOLOG
          	};
          
          for(int i=0;i<removeQualifierNames.length; i++)
            removeQualifierByName(removeQualifierNames[i]);
        }
        catch(ReadOnlyException e){}
        catch(EntryInformationException e){}
      }
      else
      {
        chadoGene = ((GFFStreamFeature)feature).chadoGene;
      }
    }
  }

  /**
   *  Create a new GFFStreamFeature from the given line.  The String should be
   *  in gene finder format.
   **/
  public GFFStreamFeature(final String line)
      throws ReadFormatException 
  {
    super(null);

    final StringVector line_bits = StringVector.getStrings(line, "\t", true);

    if(line_bits.size() < 8) 
      throw new ReadFormatException("invalid GFF line: 8 fields needed " +
                                    "(got " + line_bits.size () +
                                    " fields) from: " + line);

    final String start_base_string = ((String)line_bits.elementAt(3)).trim();
    final String end_base_string   = ((String)line_bits.elementAt(4)).trim();

    final int start_base;
    final int end_base;

    try 
    {
      start_base = Integer.parseInt(start_base_string);
      end_base   = Integer.parseInt(end_base_string);
    } 
    catch(NumberFormatException e)
    {
      throw new ReadFormatException("Could not understand the start or end base " +
                                    "of a GFF feature: " + start_base_string + 
                                    " " + end_base_string);
    }

    // start of qualifier parsing and setting
    try 
    {
      final boolean complement_flag;

      if(((String)line_bits.elementAt(6)).equals("+")) 
        complement_flag = false;
      else if(((String)line_bits.elementAt(6)).equals("-"))
        complement_flag = true;
      else 
      {
        // must be unstranded
        complement_flag = false;

        // best we can do
        //final String note_string = "this feature is unstranded";
        //setQualifier(new Qualifier("note", note_string));
      }

      if(line_bits.size() == 9) 
      {
        final String rest_of_line = (String)line_bits.elementAt(8); 

        // parse the rest of the line as ACeDB format attributes
        final Hashtable attributes = parseAttributes(rest_of_line);
//      final String type = (String)line_bits.elementAt(2);

        for(final java.util.Enumeration attribute_enum = attributes.keys();
            attribute_enum.hasMoreElements();)
        {
          String name = (String)attribute_enum.nextElement();

          final StringVector values = (StringVector)attributes.get(name);

          if(MatchPanel.isClusterTag(name))
          {
            List lazyValues = new Vector();
            for(int i=0; i<values.size(); i++)
              lazyValues.add(
                  new ClusterLazyQualifierValue( (String)values.get(i), name,
                                         this ));
            
            setQualifier(new QualifierLazyLoading(name, lazyValues));
          }
          else
          {
            if(values.size() == 0)
              setQualifier(new Qualifier(name));
            else
              setQualifier(new Qualifier(name, values));
          }
        }
      }

      /*if( !((String)line_bits.elementAt(0)).equals("null") )
      {
        final Qualifier gff_seqname =
          new Qualifier("gff_seqname", decode((String)line_bits.elementAt(0)));

        setQualifier(gff_seqname);
      }*/
      if( !((String)line_bits.elementAt(0)).equals("null") )
        setGffSeqName( decode((String)line_bits.elementAt(0)) );
      
      final Key key = new Key((String)line_bits.elementAt(2));
      setKey(key);

      /*final Qualifier source_qualifier =
        new Qualifier("gff_source", (String)line_bits.elementAt(1));
      setQualifier(source_qualifier);*/
      this.setGffSource((String)line_bits.elementAt(1));
      
      if( !((String)line_bits.elementAt(5)).equals(".") )
      {
        final Qualifier score_qualifier =
          new Qualifier("score", (String)line_bits.elementAt(5));
        setQualifier(score_qualifier);
      }
      
      String frame = (String)line_bits.elementAt(7);

      if(frame.equals ("0"))
        frame = "1";
      else if(frame.equals("1"))
        frame = "2";
      else if(frame.equals("2")) 
        frame = "3";
      else
        frame = ".";

      if(!frame.equals(".")) 
      {
        final Qualifier codon_start_qualifier =
          new Qualifier("codon_start", frame);

        setQualifier(codon_start_qualifier);
      }

      if(start_base > end_base) 
        throw new ReadFormatException("start position is greater than end " +
                                      "position: " + start_base + " > " +
                                      end_base+"\n"+line);

      if(start_base < 0)
        throw new ReadFormatException("start position must be positive: " +
                                      start_base); 
      
      final Range location_range = new Range(start_base, end_base);
      final RangeVector location_ranges = new RangeVector(location_range);
      setLocation(new Location(location_ranges, complement_flag));
    }
    catch(ReadOnlyException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    } 
    catch(EntryInformationException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    } 
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }

    this.gff_lines = new StringVector(line);
  }

  /**
  *
  * Store for spliced regions of segments ID's and ranges.
  *
  */
  public void setSegmentRangeStore(Hashtable id_range_store)
  {
    this.id_range_store = id_range_store;
  }

  public Hashtable getSegmentRangeStore()
  {
    if(id_range_store == null)
    {
      id_range_store = new Hashtable();
      id_range_store.put((String)this.getQualifierByName("ID").getValues().get(0), 
                         this.getLocation().getTotalRange());
    }
    return id_range_store;
  }
  
  public Hashtable getNewIdMapToOldId()
  {
    return newIdMapToOldId;
  }
  
  /**
   * Used when changing spliced feature uniquenames
   * @param newIdMapToOldId
   */
  public void setNewIdMapToOldId(Hashtable newIdMapToOldId)
  {
    this.newIdMapToOldId = newIdMapToOldId;
  }
  
  /**
   * Store for ID's and CHADO feature_relationship.rank
   * @param feature_relationship_rank_store
   */
  public void setFeature_relationship_rank_store(
      Hashtable feature_relationship_rank_store)
  {
    this.feature_relationship_rank_store = feature_relationship_rank_store;
  }
  
  /**
   * Store for ID's and CHADO feature_relationship.rank
   * @return
   */
  public Hashtable getFeature_relationship_rank_store()
  {
    return feature_relationship_rank_store;
  }
  
 
  /**
   * Get the chado uniquename 
   * @param r
   * @return
   */
  public String getSegmentID(final Range r)
  {
    if(id_range_store != null)
    {
      Enumeration enum_ranges = id_range_store.keys();
      //Iterator it = id_range_store.values().iterator();
      while(enum_ranges.hasMoreElements())
      //while(it.hasNext())
      {
        String key  = (String)enum_ranges.nextElement();
        Range range = (Range)id_range_store.get(key);
        if(range.getStart() == r.getStart() &&
           range.getEnd()   == r.getEnd())
          return key;
      }
    }
    else if (getQualifierByName("ID") != null)
    {
      return (String)getQualifierByName("ID").getValues().get(0);
    }

    logger4j.warn("RANGE NOT FOUND "+r.toString());

    return null;
  }

  /**
   * Get the feature ID based on the segments chado 
   * uniquename's.
   * @param rv
   * @return
   */
  public String getSegmentID(RangeVector rv)
  {
    String id = "";
    if(id_range_store != null)
    {
      String id_new;
      Range range;
      int index;
      for(int i=0; i<rv.size(); i++)
      {
        range  = (Range)rv.get(i);
        id_new = getSegmentID(range);
        
        String prefix[] = getPrefix(id_new, ':');
        if(prefix[0] != null)
        {
          index = id.indexOf(prefix[0]);
          if(id.equals("") || index < 0)
          {
            if(!id.equals(""))
              id = id +",";
            id = id+prefix[0] + "{" + prefix[1] + "}";
            continue;
          }
          
          index = id.indexOf('}', index);
          id = id.substring(0,index) + "," + 
               prefix[1] + id.substring(index);
        }
        else if(id_new != null)
        {
          if(!id.equals(""))
            id = id +",";
          id = id+id_new;
        }
      }
    }
    
    return id;
  }
  
  /**
   * Get the ID prefix, e.g. for SPAC1556.06.1:exon:2
   * returns SPAC1556.06.1:exon as the prefix and 2 as the
   * index.
   * @param id
   * @return
   */
  public String[] getPrefix(final String id,
                            final char separator)
  {
    String prefix[] = new String[2];
    int index = id.lastIndexOf(separator);

    if(index > -1)
    {
      prefix[0] = id.substring(0,index);
      prefix[1] = id.substring(index+1);
    }
    return prefix;
  }
 
  /**
   * Used to automatically generate
   * @param prefix
   * @return
   */
  public int getAutoNumber(final String prefix,
                           final char separator)
  {
    int auto   = 1;
    String val = prefix + separator + auto;
    while(id_range_store.containsKey(val))
    {
      auto++;
      val = prefix + separator + auto;
    }
    return auto;
  }
  

  
  /**
  * For gff-version 3:
  * http://song.sourceforge.net/gff3-jan04.shtml
  *
  * Remove URL escaping rule (e.g. space="%20" or "+")
  */
  public static String decode(String s)
  {
    int ind;
    String enc;
    String dec;

    for(int i=0; i<MAP_DECODE.length; i++)
    {
      enc = MAP_DECODE[i][1];
      dec = MAP_DECODE[i][0];
      while( (ind = s.indexOf(enc)) > -1)
        s = s.substring(0,ind) + dec + s.substring(ind+enc.length());
    }

    return s;
  }

  
  /**
  * For gff-version 3:
  * http://song.sourceforge.net/gff3-jan04.shtml
  *
  * Add URL escaping rule (e.g. space="%20" or "+")
  */
  public static String encode(String s)
  {
    int ind;
    String enc;
    String dec;

    for(int i=0; i<MAP_ENCODE.length; i++)
    {
      enc = MAP_ENCODE[i][1];
      dec = MAP_ENCODE[i][0];
      while( (ind = s.indexOf(dec)) > -1 )
        s = s.substring(0,ind) + enc + s.substring(ind+1);
    }

    return s;
  }

   
  /**
   *  Return the reference of a new copy of this Feature.
   **/
  public Feature copy() 
  {
    final Feature return_value = new GFFStreamFeature(this);
    return return_value;
  }

  /**
   *  Read and return a GFFStreamFeature from a stream.  A feature must be the
   *  next thing in the stream.
   *  @param stream the Feature is read from this stream
   *  @exception IOException thrown if there is a problem reading the Feature -
   *    most likely ReadFormatException.
   *  @exception InvalidRelationException Thrown if this Feature cannot contain
   *    the given Qualifier.
   *  @return null if in_stream is at the end of file when the method is
   *    called
   */
  protected static GFFStreamFeature readFromStream(LinePushBackReader stream)
      throws IOException, InvalidRelationException 
  {
    String line = stream.readLine();
    if(line == null) 
      return null;

    try
    {
      final GFFStreamFeature new_feature = new GFFStreamFeature(line);
      return new_feature;
    } 
    catch(ReadFormatException exception) 
    {
      // re-throw the exception with the line number added
      final String new_error_string = exception.getMessage();

      throw new ReadFormatException(new_error_string,
                                    stream.getLineNumber());
    }
  }

  /**
   *  Read the details of a feature from an EMBL stream into the current
   *  object.
   *  @param entry_information The EntryInformation object of the Entry that
   *    will contain the Feature.
   *  @param in_stream the Feature is read from this stream
   *  @exception IOException thrown if there is a problem reading the Feature -
   *    most likely ReadFormatException if the stream does not contain GFF
   *    feature.
   **/
  public void setFromStream(final EntryInformation entry_information,
                            final LinePushBackReader in_stream)
      throws IOException, InvalidRelationException, ReadOnlyException 
  {
    throw new ReadOnlyException();
  }

  protected static Hashtable contig_ranges;

  /**
   *  Write this Feature to the given stream.
   *  @param writer The stream to write to.
   *  @exception IOException thrown if there is an io problem while writing
   *    the Feature.
   **/
  public void writeToStream(final Writer writer)
      throws IOException 
  {
    final RangeVector ranges = getLocation().getRanges();
    final int ranges_size = ranges.size();

//  final Hashtable contig_ranges = SimpleDocumentEntry.getContigRanges();
    for(int i = 0; i < ranges_size; ++i) 
    {
      Range this_range = (Range)ranges.elementAt(i);

      String seqname = getGffSeqName();
      String source  = getGffSource();
      Qualifier score   = getQualifierByName("score");
      Qualifier group   = getQualifierByName("group");

      // source becomes a Dbxref in chado
      String source_str = null;
      if(getQualifierByName("Dbxref") != null)
      {
        source_str = getDbxrefGFFSource(getQualifierByName("Dbxref"));
      }
      
      if(seqname == null && ((GFFDocumentEntry)getEntry()).getDocument() != null) 
        seqname = ((GFFDocumentEntry)getEntry()).getDocument().getName();
      if(seqname == null)
      {
        try
        {
          seqname = ((GFFStreamFeature)(getEntry().getAllFeatures().elementAt(0))).getGffSeqName();
        }
        catch(Exception e) {}
        if(seqname == null)
          seqname = "gff_seqname";
      }
      
      if(source == null) 
        source = "artemis";

      if(score == null) 
        score = new Qualifier("score", ".");

      int start = this_range.getStart();
      int end   = this_range.getEnd();
      if(seqname != null && contig_ranges != null &&
         contig_ranges.containsKey(seqname))
      {
        Range offset_range = (Range)contig_ranges.get(seqname);
        start = start-offset_range.getStart()+1;
        end   = end-offset_range.getStart()+1;
      }

      if(group == null || group.getValues() == null ||
         group.getValues().elementAt(0).equals(""))
      {
        final Qualifier gene = getQualifierByName("gene");

        if(gene == null) 
          group = new Qualifier("group", "");
        else 
          group = gene;
      }

      String frame = ".";

      final Qualifier codon_start = getQualifierByName("codon_start");

      if(codon_start != null && i == 0) 
      {
        frame = (String)(codon_start.getValues()).elementAt(0);

        if(frame.equals ("1")) 
          frame = "0";
        else if(frame.equals("2"))
          frame = "1";
        else if(frame.equals("3"))
          frame = "2";
        else
          frame = ".";
      }

      final String myId = getSegmentID(this_range);
      String attribute_string = unParseAttributes(myId);

      if(source_str == null && source != null)
       source_str = source;

      String key = getKey().getKeyString();

      String translation = getTranslation();
      if(translation != null)
        attribute_string = attribute_string + ";" + translation;
      writer.write(seqname + "\t" +
                   source_str + "\t" +
                   key + "\t" +
                   start + "\t" +
                   end + "\t" +
                   score.getValues() .elementAt(0)+ "\t" +
                   (getLocation().isComplement() ? "-\t" : "+\t") +
                   frame + "\t" +
                   attribute_string + "\n");
    }

  }

  /**
   * Get the GFF_source value of a Dbxref qualifier.
   * @param qualifier
   * @return  the gff_source value or NULL
   */
  /*
  private String getDbxrefGFFSource(final Qualifier qualifier)
  {
    StringVector qualifier_strings =
      StreamQualifier.toStringVector(null, qualifier);
    
    for(int i=0; i<qualifier_strings.size(); i++)
    {
      String qualifier_string = (String)qualifier_strings.elementAt(i);
      
      if(qualifier_string.indexOf("GFF_source:") >-1)
      {
        int index = qualifier_string.indexOf(":")+1;
        int len = qualifier_string.length();
        if(qualifier_string.endsWith("\""))
          len--;
        return qualifier_string.substring(index, len);
      }
    }
    return null;
  }
  */
  
  /**
   *  Return a String containing the qualifiers of this feature in a form
   *  suitable for using as the last field of a GFF line.  The codon_start
   *  attribute is not included since GFF has a frame field.  gff_seqname,
   *  gff_source and score aren't included since they have corresponding
   *  fields.
   **/
  private String unParseAttributes(final String myId) 
  {
    final StringBuffer buffer = new StringBuffer();
    final QualifierVector qualifiers = getQualifiers();

    final String names[] = { "ID", "Name", "Alias", "Parent",
                             "Derives_from",
                             "Target", "Gap", "Note", 
                             "Dbxref", "Ontology_term" };
    int count = 0;
    Qualifier this_qualifier;
    final int names_length = names.length;

    if(myId != null)
    {
      buffer.append("ID=");
      buffer.append(encode(myId));
      count++;
    }
    
    for(int i=1; i<names_length; i++)
    {
      this_qualifier = (Qualifier)qualifiers.getQualifierByName(names[i]);
      
      if(this_qualifier == null) 
        continue;
      
      // GSV :: see new getQualifierString signature
      // this qualifier is one of the reserved qualifiers 
      final String this_qualifier_str = getQualifierString(this_qualifier, true);
      if(this_qualifier_str == null)
        continue;

      if(count != 0)
        buffer.append(";");
      buffer.append(this_qualifier_str);
      count++;
    }
    
    boolean lname;
    final int qualifiers_size = qualifiers.size();
    for(int i = 0; i < qualifiers_size; i++) 
    {
      this_qualifier = (Qualifier)qualifiers.elementAt(i);

      lname = false;
      for(int j=0; j<names_length; j++)
        if(this_qualifier.getName().equals(names[j]))
          lname = true;

      if(lname)
        continue;
      
      if( (this_qualifier.getName().equals("private") && System.getProperty("noprivate") != null) ||
          (this_qualifier.getName().equals("history") && System.getProperty("nohistory") != null) )
        continue;

      // GSV :: see new getQualifierString signature
      // this qualifier is NOT one of the reserved qualifiers 
      String this_qualifier_str = getQualifierString(this_qualifier, false);
      
      if(this_qualifier_str == null)
        continue;
      
      if(count != 0)
        buffer.append(";");
      buffer.append(this_qualifier_str);
    }

    return buffer.toString();
  }

  
  /**
   * Get the translation qualifier string for polypeptide features.
   */
  private String getTranslation() 
  {
    if (! getKey().getKeyString().equals("polypeptide"))
        return null;
    if (chadoGene != null)
    {
      if(getUserData() == null)
      {
        uk.ac.sanger.artemis.Feature f = new uk.ac.sanger.artemis.Feature(this);
      }
      // the above line constructs the appropriate userData within this current GFFStreamFeature object, 
      // which is required by the following GeneUtils.deriveResidues()     
      String residues = GeneUtils.deriveResidues(this);
      if (residues != null)
        return "translation="+residues;
    } 
    return null;
  }
  
  /**
   * Used to write out the GFF attributes.
   * @param q the qualifier to represent as a <code>String</code>
   * @param reserved indicate if this is one of the reserved tags or not
   * @return  the <code>String</code> representation
   * 
   * GSV: modified the signature to force the caller to declare if this 
   * qualifier is one of the reserved ones.
   */
  private String getQualifierString(Qualifier q, boolean reserved )
  {
    StringBuffer buffer = new StringBuffer();
    final String name = q.getName();

    if(name.equals("codon_start") || name.equals("gff_source") ||
       name.equals("gff_seqname") || name.equals("score"))
      return null;

    final StringVector values = q.getValues();
    
    /* 
     * GSV :
     * 
     * The Bio::FeatureIO perl module falls over if there are Uppercased 
     * attribute names for tags which aren't part of the standard reserved 
     * set. So we lowercase these, since in the specification it says :
     * 
     * "All attributes that begin with an uppercase letter are reserved for
     * later use.  Attributes that begin with a lowercase letter can be used
     * freely by applications."
     * 
     * see http://www.sequenceontology.org/gff3.shtml
     */
    String nameToBuffer = encode(name);
    if (! reserved)
    	nameToBuffer = Character.toLowerCase(nameToBuffer.charAt(0)) + nameToBuffer.substring(1);
    buffer.append(nameToBuffer);
    
    if(values != null)
    {
      buffer.append('=');
      for(int value_index = 0; value_index < values.size();
          ++value_index)
      {
        final String this_value;
        if(name.equals("class"))
        {
          int index = ((String)values.elementAt(value_index)).indexOf("::");
          if(index > -1)
            this_value = encode(((String)values.elementAt(value_index)).substring(0,index));
          else
            this_value = encode((String)values.elementAt(value_index));
        }
        else
          this_value = encode((String)values.elementAt(value_index));
        
        if(value_index>0)
          buffer.append("%2C");
        
        if(name.equals("Parent"))
          buffer.append(this_value);
        else
        {
          try
          {
            buffer.append(Integer.valueOf(this_value));
          }
          catch(NumberFormatException _)
          {
            // not an integer
            try
            {
              buffer.append(Double.valueOf(this_value));
            }
            catch (NumberFormatException __)
            {
              // not a double or integer so quote it
              buffer.append(this_value);
            }
          }
        }
      }
    }
    return buffer.toString();
  }

  /**
   *  Parse the given String as ACeDB format attributes.
   *  Adapted from code by Matthew Pocock for the BioJava project.
   *
   *  Modified for gff-version 3.
   *
   *  @return Return a Hashtable.  Each key is an attribute name and each value
   *    of the Hashtable is a StringVector containing the attribute values.
   *    If the attribute has no value then the Hashtable value will be a zero
   *    length vector.
   **/
  private Hashtable parseAttributes(final String att_val_list) 
  {
    Hashtable attributes = new Hashtable();

//  StringTokenizer tokeniser = new StringTokenizer(att_val_list, ";", false);
//  while(tokeniser.hasMoreTokens()) 
//  {
//    final String this_token = tokeniser.nextToken().trim();

    int ind_start = 0;
    int ind_end;
    while( (ind_end = att_val_list.indexOf(";",ind_start)) > -1 || 
           ind_start < att_val_list.length() )
    {
      if(ind_end < 0)
        ind_end = att_val_list.length();

      final String this_token = decode(att_val_list.substring(ind_start, ind_end).trim());
      ind_start = ind_end+1;
      
      /*if(this_token.startsWith("feature_relationship_rank="))
      {
        setFeature_relationship_rank( 
            Integer.parseInt(this_token.substring(26)) );
        continue;
      }*/

      int index_of_first_space = this_token.indexOf(" ");
       
      String att_name;
      StringVector att_values = new StringVector();

      if( this_token.indexOf("=") > -1 &&
         (this_token.indexOf("=") < index_of_first_space ||
          index_of_first_space == -1) )
      {
        index_of_first_space = this_token.indexOf("=");
        att_name = this_token.substring(0, index_of_first_space);
        att_values.add(this_token.substring(index_of_first_space+1).trim());
      }
      else if(index_of_first_space == -1) 
        att_name = this_token;
      else 
      {
        att_name = this_token.substring(0, index_of_first_space);

        String rest_of_token =
          this_token.substring(index_of_first_space+1).trim();

        while(rest_of_token.length() > 0) 
        {
          if(rest_of_token.startsWith("\""))
          {
            int quote_index = 0;
            do 
            {
              quote_index++;
              quote_index = rest_of_token.indexOf("\"", quote_index);
            } while(quote_index > -1 &&
                    rest_of_token.charAt(quote_index - 1) == '\\');

            if(quote_index < 0) 
            {
              // no closing quote - panic
              final Hashtable panic_attributes = new Hashtable();
              final StringVector notes = new StringVector();
              notes.add(att_val_list);
              panic_attributes.put("note", notes);

              return panic_attributes;
            }

            final String next_bit = rest_of_token.substring(1, quote_index);
            att_values.add(next_bit);
            rest_of_token = rest_of_token.substring(quote_index + 1).trim();
          } 
          else
          {
            final int index_of_next_space = rest_of_token.indexOf(" ");

            if(index_of_next_space == -1) 
            {
              att_values.add(rest_of_token);
              rest_of_token = "";
            } 
            else 
            {
              final String next_bit =
                rest_of_token.substring(0, index_of_next_space);

              att_values.add(next_bit);
              rest_of_token =
                rest_of_token.substring(index_of_next_space).trim();
            }
          }
        }

        if(!rest_of_token.equals(""))
          att_values.add(rest_of_token);
      }

      if(att_name.equals("Dbxref") || att_name.equals("Alias")) // convert to multi-line
      {
        StringTokenizer stok = 
            new StringTokenizer((String)att_values.get(0), ",");
        StringVector str_values = new StringVector();
        while(stok.hasMoreTokens())
          str_values.add(stok.nextElement());

        att_values = str_values;
      }
      
      if(att_name.equals("timelastmodified"))
      {
        try
        {
          this.timelastmodified = 
                  new Timestamp( Long.parseLong((String)att_values.get(0)) );
          SimpleDateFormat date_format = 
                  new SimpleDateFormat("dd.MM.yyyy hh:mm:ss z");
          att_values.set(0,date_format.format(timelastmodified));
        }
        catch(NumberFormatException e)
        {
          att_values.set(0,(String)att_values.get(0));
        }
      }
      
      if(attributes.get(att_name) != null) 
        ((StringVector)attributes.get(att_name)).add(att_values);
      else 
        attributes.put(att_name, att_values);
    }

    return attributes;
  }

  /**
   * Get the feature time last modified timestamp.
   * @return
   */
  public Timestamp getLastModified()
  {
    return timelastmodified;
  }
  
  /**
   * Get the GFF_source value of a Dbxref qualifier.
   * @param qualifier
   * @return  the gff_source value or NULL
   */
  private String getDbxrefGFFSource(final Qualifier qualifier)
  {
    StringVector qualifier_strings =
      StreamQualifier.toStringVector(null, qualifier);
    
    for(int i=0; i<qualifier_strings.size(); i++)
    {
      String qualifier_string = (String)qualifier_strings.elementAt(i);
      
      if(qualifier_string.indexOf("GFF_source:") >-1)
      {
        int index = qualifier_string.indexOf(":")+1;
        int len = qualifier_string.length();
        if(qualifier_string.endsWith("\""))
          len--;
        return qualifier_string.substring(index, len);
      }
    }
    return null;
  }
  
  /**
   * Set the feature time last modified timestamp.
   * @param timelastmodified
   */
  public void setLastModified(final Timestamp timelastmodified)
  {
    this.timelastmodified = timelastmodified;
    
    // now update the qualifier value itself
    QualifierVector qualifiers = getQualifiers();
    Qualifier qualifier = qualifiers.getQualifierByName("timelastmodified");
    SimpleDateFormat date_format = 
      new SimpleDateFormat("dd.MM.yyyy hh:mm:ss z");
    
    if(qualifier != null)
      qualifier.removeValue((String)qualifier.getValues().get(0));
    else
    {
      try
      {
        qualifier = new Qualifier("timelastmodified", 
               date_format.format(timelastmodified));
        setQualifier(qualifier);
        return;
      }
      catch(EntryInformationException eie)
      {}
      catch(ReadOnlyException roe)
      {}
    }
    
    qualifier.addValue(date_format.format(timelastmodified));
  }
  
  /**
   *  Returns true if and only if this Feature can't be changed or can't be
   *  removed from it's entry.
   **/
  public boolean isReadOnly () 
  {
    if(readOnlyFeature)
      return true;
    return super.isReadOnly();
  }

  public void setReadOnlyFeature(boolean readOnlyFeature)
  {
    this.readOnlyFeature = readOnlyFeature;
  }
  
  public ChadoCanonicalGene getChadoGene()
  {
    return chadoGene;
  }

  public void setChadoGene(ChadoCanonicalGene chadoGene)
  {
    this.chadoGene = chadoGene;
  }
  
  public boolean isVisible()
  {
    return visible;
  }

  public void setVisible(boolean visible)
  {
    this.visible = visible;
  }

  public String getGffSeqName()
  {
    return gffSeqName;
  }

  public void setGffSeqName(String gffSeqName)
  {
    this.gffSeqName = gffSeqName;
  }

  public String getGffSource()
  {
    return gffSource;
  }

  public void setGffSource(String gffSource)
  {
    this.gffSource = gffSource;
  }

  public boolean isLazyLoaded()
  {
    return lazyLoaded;
  }

  public void setLazyLoaded(boolean lazyLoaded)
  {
    this.lazyLoaded = lazyLoaded;
  }

  public org.gmod.schema.sequence.Feature getChadoLazyFeature()
  {
    return chadoLazyFeature;
  }

  public void setChadoLazyFeature(
      org.gmod.schema.sequence.Feature chadoLazyFeature)
  {
    this.chadoLazyFeature = chadoLazyFeature;
  }
  
  protected static boolean isGTF(Feature feature)
  {
    if(!(feature instanceof GFFStreamFeature))
      return false;
    
    final String names[] = { "ID", "Name", "Alias", "Parent",
        "Derives_from",
        "Target", "Gap", "Note", 
        "Dbxref", "Ontology_term" };
    
    for(String name: names)
    {
      if(feature.getQualifiers().getQualifierByName(name) != null)
        return false;
    }
    
    if(feature.getQualifiers().getQualifierByName("gene_id") != null && 
       feature.getQualifiers().getQualifierByName("transcript_id") != null)
    {
      logger4j.debug(feature.getEntry().getName()+" is in GTF format");
      return true;
    }
    return false;
  }
  
  public static void main(String args[])
  {
    Key key = new Key("region");
    try
    {
      final EntryInformation entry_information =
        SimpleEntryInformation.getDefaultEntryInformation ();
      GFFDocumentEntry entry = new GFFDocumentEntry(entry_information);

      Location location = new Location("1003..1222");
      QualifierVector qualifiers = new QualifierVector();
      GFFStreamFeature f = new GFFStreamFeature(key, location, qualifiers);
      entry.add(f);
      
      java.io.File aFile = new java.io.File("x");
      java.io.FileWriter writer = new java.io.FileWriter(aFile);
      f.writeToStream(writer);
      writer.close();
    }
    catch (LocationParseException e)
    {
      e.printStackTrace();
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
    catch (EntryInformationException e)
    {
      e.printStackTrace();
    }
  }
}
