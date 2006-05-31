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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/GFFStreamFeature.java,v 1.36 2006-05-31 10:38:48 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;

import java.io.*;
import java.util.Hashtable;
import java.util.Enumeration;
import java.util.StringTokenizer;
import java.sql.Timestamp;
import java.text.SimpleDateFormat;


/**
 *  A StreamFeature that thinks it is a GFF feature.
 *
 *  @author Kim Rutherford
 *  @version $Id: GFFStreamFeature.java,v 1.36 2006-05-31 10:38:48 tjc Exp $
 **/

public class GFFStreamFeature extends SimpleDocumentFeature
                       implements DocumentFeature, StreamFeature, ComparableFeature 
{

  /**
   *  This is the line of GFF input that was read to get this
   *  GFFStreamFeature.  A GFFStreamFeature that was created from multiple GFF
   *  lines will have a gff_lines variable that contains multiple line.
   **/
  StringVector gff_lines = null;

  /** store for spliced features containing id and range of each segment */
  private Hashtable id_range_store;

  /** store the Timestamp for the feature */
  private Timestamp timelastmodified;
  
  private ChadoCanonicalGene chadoGene;
  
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
      if(getQualifierByName("score") == null)
        setQualifier(new Qualifier("score", "."));
      
      if(getQualifierByName("gff_source") == null)
        setQualifier(new Qualifier("gff_source", "artemis"));
      
      if(getQualifierByName("gff_seqname") == null)
        setQualifier(new Qualifier("gff_seqname", "."));

      if(getQualifierByName("ID") == null)
        setQualifier(new Qualifier("ID", "to_be_set"));
      
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

  /**
   *  Create a new GFFStreamFeature with the same key, location and
   *  qualifiers as the given feature.  The feature should be added to an
   *  Entry (with Entry.add()).
   *  @param feature The feature to copy.
   **/
  public GFFStreamFeature(final Feature feature) 
  {
    this(feature.getKey(), feature.getLocation(), feature.getQualifiers());
    
  //  if(feature instanceof GFFStreamFeature)
  //  {
  //    this.id_range_store = ((GFFStreamFeature)feature).id_range_store;
  //  }
  }

  /**
   *  Create a new GFFStreamFeature from the given line.  The String should be
   *  in gene finder format.
   **/
  private GFFStreamFeature(final String line)
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

          if(values.size() == 0)
            setQualifier(new Qualifier(name));
          else
            setQualifier(new Qualifier(name, values));
        }
      }

      final Qualifier gff_seqname =
        new Qualifier("gff_seqname", decode((String)line_bits.elementAt(0)));

      setQualifier(gff_seqname);

      final Key key = new Key((String)line_bits.elementAt(2));
      setKey(key);

      final Qualifier source_qualifier =
        new Qualifier("gff_source", (String)line_bits.elementAt(1));
      setQualifier(source_qualifier);

      final Qualifier score_qualifier =
        new Qualifier("score", (String)line_bits.elementAt(5));
      setQualifier(score_qualifier);

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
                                      end_base);

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


  public String getSegmentID(Range r)
  {
    if(id_range_store != null)
    {
      Enumeration enum_ranges = id_range_store.keys();
      while(enum_ranges.hasMoreElements())
      {
        Range range = (Range)enum_ranges.nextElement();
        if(range.getStart() == r.getStart() ||
           range.getEnd()   == r.getEnd())
        return (String)id_range_store.get(range);
      }
    }

    return null;
  }

  public String getSegmentID(RangeVector rv)
  {
    String id = "";
    if(id_range_store != null)
    {
      String id_new;
      Range range;
      for(int i=0; i<rv.size(); i++)
      {
        range  = (Range)rv.get(i);
        id_new = getSegmentID(range);
        if(id_new != null)
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
  *
  * For gff-version 3:
  * http://song.sourceforge.net/gff3-jan04.shtml
  *
  * Remove URL escaping rule (e.g. space="%20" or "+")
  *
  */
  public static String decode(String s)
  {
    final String map[][] = {
                             { " ",  "%20" },  // white space
                             { ",",  "%2C" },  // comma
                             { ";",  "%3B" },  // semi-colon
                             { "=",  "%3D" },  // equals
                             { "\t", "%09" },  // tab
                             { " ",  "+"   },  // white space
                             { "(",  "%28" },  // left bracket
                             { ")",  "%29" }   // right bracket )
                           };

    int ind;
    String enc;
    String dec;

    for(int i=0; i<map.length; i++)
    {
      enc = map[i][1];
      dec = map[i][0];
      while( (ind = s.indexOf(enc)) > -1)
        s = s.substring(0,ind) + dec + s.substring(ind+enc.length());
    }

    return s;
  }

  /**
  *
  * For gff-version 3:
  * http://song.sourceforge.net/gff3-jan04.shtml
  *
  * Add URL escaping rule (e.g. space="%20" or "+")
  *
  */
  public static String encode(String s)
  {
    final String map[][] = {
//                           { " ",  "%20" },  // white space
                             { ",",  "%2C" },  // comma 
                             { ";",  "%3B" },  // semi-colon
                             { "=",  "%3D" },  // equals
                             { "\t", "%09" },  // tab
                             { " ",  "+"   },  // white space
                             { "(",  "%28" },  // left bracket
                             { ")",  "%29" }   // right bracket )
                           };

    int ind;
    String enc;
    String dec;

    for(int i=0; i<map.length; i++)
    {
      enc = map[i][1];
      dec = map[i][0];
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

      Qualifier seqname = getQualifierByName("gff_seqname");
      Qualifier source  = getQualifierByName("gff_source");
      Qualifier score   = getQualifierByName("score");
      Qualifier group   = getQualifierByName("group");

      // source becomes a Dbxref in chado
     // String source_str = null;
     // if(getQualifierByName("Dbxref") != null)
     //   source_str = getDbxrefGFFSource(getQualifierByName("Dbxref"));
      
      if(seqname == null) 
        seqname = new Qualifier("gff_seqname", ".");

      if(source == null) 
        source = new Qualifier("source", "artemis");

      if(score == null) 
        score = new Qualifier("score", ".");

      int start = this_range.getStart();
      int end   = this_range.getEnd();
      if(seqname != null && contig_ranges != null &&
         contig_ranges.containsKey(seqname.getValues().elementAt(0)))
      {
        Range offset_range = (Range)contig_ranges.get(seqname.getValues().elementAt(0));
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

      final String attribute_string = unParseAttributes();

      //if(source_str == null)
      String source_str = (String)source.getValues().elementAt(0);
      
      writer.write(seqname.getValues().elementAt(0) + "\t" +
                   source_str + "\t" +
                   getKey() + "\t" +
                   start + "\t" +
                   end + "\t" +
                   score.getValues() .elementAt(0)+ "\t" +
                   (getLocation().isComplement() ? "-\t" : "+\t") +
                   frame + "\t" +
                   attribute_string + "\n");
    }

 // for(int i = 0 ; i < gff_lines.size() ; ++i) 
 //   writer.write(gff_lines.elementAt(i) + "\n");
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
  private String unParseAttributes() 
  {
    final StringBuffer buffer = new StringBuffer();
    final QualifierVector qualifiers = getQualifiers();

    final String names[] = { "ID", "Name", "Alias", "Parent",
                             "Target", "Gap", "Note", 
                             "Dbxref", "Ontology_term" };
    int count = 0;
    Qualifier this_qualifier;
    final int names_length = names.length;

    for(int i=0; i<names_length; i++)
    {
      this_qualifier = (Qualifier)qualifiers.getQualifierByName(names[i]);
 
      if(this_qualifier == null)
        continue;

      final String this_qualifier_str = getQualifierString(this_qualifier);
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
      
      String this_qualifier_str = getQualifierString(this_qualifier);
      
      if(this_qualifier_str == null)
        continue;
      
      if(count != 0)
        buffer.append(";");
      buffer.append(this_qualifier_str);
    }

    return buffer.toString();
  }

  /**
   * Used to write out the GFF attributes.
   * @param q the qualifier to represent as a <code>String</code>
   * @return  the <code>String</code> representation
   */
  private String getQualifierString(Qualifier q)
  {
    StringBuffer buffer = new StringBuffer();
    final String name = q.getName();

    if(name.equals("codon_start") || name.equals("gff_source") ||
       name.equals("gff_seqname") || name.equals("score"))
      return null;

    final StringVector values = q.getValues();
    buffer.append(encode(name));

    if(values != null)
    {
      buffer.append('=');
      for(int value_index = 0; value_index < values.size();
          ++value_index)
      {
        final String this_value = encode((String)values.elementAt(value_index));
        if(value_index>0)
          buffer.append("%2C");
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
        this.timelastmodified = 
                  new Timestamp( Long.parseLong((String)att_values.get(0)) );
        SimpleDateFormat date_format = 
                  new SimpleDateFormat("dd.MM.yyyy hh:mm:ss z");
        att_values.set(0,date_format.format(timelastmodified));
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

  public ChadoCanonicalGene getChadoGene()
  {
    return chadoGene;
  }

  public void setChadoGene(ChadoCanonicalGene chadoGene)
  {
    this.chadoGene = chadoGene;
  }
}
