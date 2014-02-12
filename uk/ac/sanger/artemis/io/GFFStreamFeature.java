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
import java.util.HashSet;
import java.util.Enumeration;
import java.util.List;
import java.util.Set;
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
import uk.ac.sanger.artemis.io.GFF3Encoder;

/**
 * A StreamFeature that thinks it is a GFF feature.
 * 
 * @author Kim Rutherford
 **/
public class GFFStreamFeature extends SimpleDocumentFeature implements
    DocumentFeature, StreamFeature, ComparableFeature {

  private static org.apache.log4j.Logger    logger4j        = org.apache.log4j.Logger
                                                                .getLogger(GFFStreamFeature.class);

  /** store for spliced features containing id and range of each segment */
  private Hashtable<String, Range>          id_range_store;

  /** store a record of the new and old uniquenames that have been changed */
  private Hashtable<String, String>         newIdMapToOldId;

  /** store the Timestamp for the feature */
  private Timestamp                         timelastmodified;

  private ChadoCanonicalGene                chadoGene;

  private boolean                           visible         = true;

  /** combined feature_relationship.rank store for exons */
  private Hashtable<String, Integer>        feature_relationship_rank_store;

  /** first tabbed parameter */
  private String                            gffSeqName;
  /** second tabbed parameter */
  private String                            gffSource;
  /** duplication count */
  private short                             duplicate       = 0;

  protected static Hashtable<String, Range> contig_ranges;
  private boolean                           lazyLoaded      = false;
  private org.gmod.schema.sequence.Feature  chadoLazyFeature;
  private boolean                           readOnlyFeature = false;

  private static Set<String>                attrs_to_filter = new HashSet<String>();

  /**
   * Registers an attribute not to be included in the GFF3 output for
   * GFFStreamFeatures
   * 
   * @param attr
   *          The GFF3 attribute to remove
   **/
  public static void removeAttribute(String attr) {
    attrs_to_filter.add(attr);
  }

  /**
   * Registers an attribute to be included in the GFF3 output for
   * GFFStreamFeatures
   * 
   * @param attr
   *          The GFF3 attribute to include
   **/
  public static void includeAttribute(String attr) {
    attrs_to_filter.remove(attr);
  }

  /**
   * Create a new GFFStreamFeature object. The feature should be added to an
   * Entry (with Entry.add()).
   * 
   * @param key
   *          The new feature key
   * @param location
   *          The Location object for the new feature
   * @param qualifiers
   *          The qualifiers for the new feature
   **/
  public GFFStreamFeature(final Key key, final Location location,
      final QualifierVector qualifiers) {
    super(null);

    try {
      setKey(key);
      setLocation(location);
      setQualifiers(qualifiers);

      if (getQualifierByName("ID") == null) {
        String idStr = null;
        StringVector v = Options.getOptions().getSystematicQualifierNames();
        for (int i = 0; i < v.size(); i++) {
          final String sysName = (String) v.get(i);
          if (getQualifierByName(sysName) != null) {
            idStr = (String) getQualifierByName(sysName).getValues().get(0);
            break;
          }
        }
        // autogenerate ID
        if (idStr == null)
          idStr = key.getKeyString() + ":" + location.toString();
        setQualifier(new Qualifier("ID", idStr));
      }

    } catch (EntryInformationException e) {
      // this should never happen because the feature will not be in an Entry
      throw new Error("internal error - unexpected exception: " + e);
    } catch (ReadOnlyException e) {
      // this should never happen because the feature will not be in an Entry
      throw new Error("internal error - unexpected exception: " + e);
    } catch (OutOfRangeException e) {
      // this should never happen because the feature will not be in an Entry
      throw new Error("internal error - unexpected exception: " + e);
    }
  }

  public GFFStreamFeature(final Feature feature) {
    this(feature, false);
  }

  /**
   * Create a new GFFStreamFeature with the same key, location and qualifiers as
   * the given feature. The feature should be added to an Entry (with
   * Entry.add()).
   * 
   * @param feature
   *          The feature to copy.
   **/
  @SuppressWarnings("unchecked")
  public GFFStreamFeature(final Feature feature,
      final boolean isDuplicatedInChado) {
    this(feature.getKey(), feature.getLocation(), feature.getQualifiers());

    if (feature instanceof GFFStreamFeature) {
      if (((GFFStreamFeature) feature).id_range_store != null)
        this.id_range_store = (Hashtable<String, Range>) (((GFFStreamFeature) feature).id_range_store)
            .clone();

      if (((GFFStreamFeature) feature).feature_relationship_rank_store != null)
        this.feature_relationship_rank_store = (Hashtable<String, Integer>) (((GFFStreamFeature) feature).feature_relationship_rank_store)
            .clone();

      this.setGffSeqName(((GFFStreamFeature) feature).getGffSeqName());
      this.setGffSource(((GFFStreamFeature) feature).getGffSource());

      if (isDuplicatedInChado) {
        try {
          final String uniquename;
          final String duplicatePrefix;

          if (feature instanceof GFFStreamFeature) {
            ((GFFStreamFeature) feature).duplicate++;
            duplicatePrefix = "DUP"
                + Short.toString(((GFFStreamFeature) feature).duplicate) + "-";
          } else
            duplicatePrefix = "DUP";
          if (id_range_store != null) {
            final Hashtable<String, Range> new_id_range_store = new Hashtable<String, Range>(
                id_range_store.size());
            final Enumeration<String> enumIdRangeStore = id_range_store.keys();
            while (enumIdRangeStore.hasMoreElements()) {
              final String keyId = enumIdRangeStore.nextElement();
              final Range range = id_range_store.get(keyId);
              new_id_range_store.put(duplicatePrefix + keyId, range);
            }
            id_range_store.clear();
            this.id_range_store = (Hashtable<String, Range>) new_id_range_store
                .clone();

            if (getLocation().getRanges().size() > 1)
              uniquename = getSegmentID(getLocation().getRanges());
            else {
              if (((String) getQualifierByName("ID").getValues().get(0))
                  .endsWith("}"))
                uniquename = id_range_store.keys().nextElement();
              else
                uniquename = duplicatePrefix
                    + (String) getQualifierByName("ID").getValues().get(0);
            }
          } else
            uniquename = duplicatePrefix
                + (String) getQualifierByName("ID").getValues().get(0);
          setQualifier(new Qualifier("ID", uniquename));

          if (getQualifierByName("Parent") != null) {
            final String parent = (String) getQualifierByName("Parent")
                .getValues().get(0);
            setQualifier(new Qualifier("Parent", duplicatePrefix + parent));
          }

          if (getQualifierByName("Derives_from") != null) {
            final String derives_from = (String) getQualifierByName(
                "Derives_from").getValues().get(0);
            setQualifier(new Qualifier("Derives_from", duplicatePrefix
                + derives_from));
          }

          // remove qualifiers that don't get transferred to duplicate
          final String removeQualifierNames[] = { "feature_id",
              "timelastmodified", "feature_relationship_rank",
              ProteinMapPanel.POLYPEPTIDE_DOMAIN, ProteinMapPanel.TMHMM[0],
              ProteinMapPanel.TMHMM[1], ProteinMapPanel.TMHMM[2],
              ProteinMapPanel.TMHMM[3], MatchPanel.ORTHOLOG,
              MatchPanel.ORTHOLOG };

          for (int i = 0; i < removeQualifierNames.length; i++)
            removeQualifierByName(removeQualifierNames[i]);
        } catch (ReadOnlyException e) {
        } catch (EntryInformationException e) {
        }
      } else {
        chadoGene = ((GFFStreamFeature) feature).chadoGene;
      }
    }
  }

  /**
   * Create a new GFFStreamFeature from the given line. The String should be in
   * gene finder format.
   **/
  public GFFStreamFeature(final String line) throws ReadFormatException {
    super(null);

    final StringVector line_bits = StringVector.getStrings(line, "\t", true);
    if (line_bits.size() < 8)
      throw new ReadFormatException("invalid GFF line: 8 fields needed "
          + "(got " + line_bits.size() + " fields) from: " + line);

    final String start_base_str = line_bits.elementAt(3).trim();
    final String end_base_str = line_bits.elementAt(4).trim();

    final int start_base;
    final int end_base;
    try {
      start_base = Integer.parseInt(start_base_str);
      end_base = Integer.parseInt(end_base_str);
    } catch (NumberFormatException e) {
      throw new ReadFormatException(
          "Could not understand the start or end base " + "of a GFF feature: "
              + start_base_str + " " + end_base_str);
    }

    // start of qualifier parsing and setting
    try {
      final boolean complement_flag;
      if (line_bits.elementAt(6).equals("+"))
        complement_flag = false;
      else if (line_bits.elementAt(6).equals("-"))
        complement_flag = true;
      else {
        // must be unstranded
        complement_flag = false;
      }

      if (line_bits.size() == 9) {
        final String rest_of_line = line_bits.elementAt(8);
        final Hashtable<String, StringVector> attributes = parseAttributes(rest_of_line);
        for (final Enumeration<String> attribute_enum = attributes.keys(); attribute_enum
            .hasMoreElements();) {
          String name = attribute_enum.nextElement();
          final StringVector values = attributes.get(name);

          if (MatchPanel.isClusterTag(name)) {
            List<ClusterLazyQualifierValue> lazyValues = new Vector<ClusterLazyQualifierValue>();
            for (int i = 0; i < values.size(); i++)
              lazyValues.add(new ClusterLazyQualifierValue((String) values
                  .get(i), name, this));
            setQualifier(new QualifierLazyLoading(name, lazyValues));
          } else {
            if (values.size() == 0)
              setQualifier(new Qualifier(name));
            else
              setQualifier(new Qualifier(name, values));
          }
        }
      }

      if (!line_bits.elementAt(0).equals("null"))
        setGffSeqName(GFF3Encoder.decode(line_bits.elementAt(0)));

      setKey(new Key(line_bits.elementAt(2)));
      setGffSource(line_bits.elementAt(1));

      if (!line_bits.elementAt(5).equals(".")) {
        final Qualifier score_qualifier = new Qualifier("score",
            line_bits.elementAt(5));
        setQualifier(score_qualifier);
      }

      String frame = line_bits.elementAt(7);

      if (frame.equals("0"))
        frame = "1";
      else if (frame.equals("1"))
        frame = "2";
      else if (frame.equals("2"))
        frame = "3";
      else
        frame = ".";

      if (!frame.equals(".")) {
        final Qualifier codon_start_qualifier = new Qualifier("codon_start",
            frame);

        setQualifier(codon_start_qualifier);
      }

      if (start_base > end_base)
        throw new ReadFormatException("start position is greater than end "
            + "position: " + start_base + " > " + end_base + "\n" + line);

      if (start_base < 0)
        throw new ReadFormatException("start position must be positive: "
            + start_base);

      final Range location_range = new Range(start_base, end_base);
      final RangeVector location_ranges = new RangeVector(location_range);
      setLocation(new Location(location_ranges, complement_flag));
    } catch (ReadOnlyException e) {
      throw new Error("internal error - unexpected exception: " + e);
    } catch (EntryInformationException e) {
      throw new Error("internal error - unexpected exception: " + e);
    } catch (OutOfRangeException e) {
      throw new Error("internal error - unexpected exception: " + e);
    }

    // this.gff_lines = new StringVector(line);
  }

  /**
   * 
   * Store for spliced regions of segments ID's and ranges.
   * 
   */
  public void setSegmentRangeStore(Hashtable<String, Range> id_range_store) {
    this.id_range_store = id_range_store;
  }

  public Hashtable<String, Range> getSegmentRangeStore() {
    if (id_range_store == null) {
      id_range_store = new Hashtable<String, Range>();
      id_range_store.put((String) this.getQualifierByName("ID").getValues()
          .get(0), this.getLocation().getTotalRange());
    }
    return id_range_store;
  }

  public Hashtable<String, String> getNewIdMapToOldId() {
    return newIdMapToOldId;
  }

  /**
   * Used when changing spliced feature uniquenames
   * 
   * @param newIdMapToOldId
   */
  public void setNewIdMapToOldId(Hashtable<String, String> newIdMapToOldId) {
    this.newIdMapToOldId = newIdMapToOldId;
  }

  /**
   * Store for ID's and CHADO feature_relationship.rank
   * 
   * @param feature_relationship_rank_store
   */
  public void setFeature_relationship_rank_store(
      Hashtable<String, Integer> feature_relationship_rank_store) {
    this.feature_relationship_rank_store = feature_relationship_rank_store;
  }

  /**
   * Store for ID's and CHADO feature_relationship.rank
   * 
   * @return
   */
  public Hashtable<String, Integer> getFeature_relationship_rank_store() {
    return feature_relationship_rank_store;
  }

  /**
   * Get the chado uniquename
   * 
   * @param r
   * @return
   */
  public String getSegmentID(final Range r) {
    if (id_range_store != null) {
      int offset = 0;
      if (getGffSeqName() != null && contig_ranges != null
          && contig_ranges.containsKey(getGffSeqName())) {
        // adjust for coordinates in multi-sequence GFF
        Range offset_range = contig_ranges.get(getGffSeqName());
        offset = offset_range.getStart() - 1;
      }

      Enumeration<String> enum_ranges = id_range_store.keys();
      while (enum_ranges.hasMoreElements()) {
        String key = enum_ranges.nextElement();
        Range range = id_range_store.get(key);
        if (range.getStart() == r.getStart() - offset
            && range.getEnd() == r.getEnd() - offset)
          return key;
      }
    } else if (getQualifierByName("ID") != null) {
      return (String) getQualifierByName("ID").getValues().get(0);
    }

    logger4j.warn("RANGE NOT FOUND " + r.toString());

    return null;
  }

  /**
   * Get the feature ID based on the segments chado uniquename's.
   * 
   * @param rv
   * @return
   */
  public String getSegmentID(final RangeVector rv) {
    String id = "";
    if (id_range_store != null) {
      String id_new;
      Range range;
      int index;
      for (int i = 0; i < rv.size(); i++) {
        range = (Range) rv.get(i);
        id_new = getSegmentID(range);

        String prefix[] = getPrefix(id_new, ':');
        if (prefix[0] != null) {
          index = id.indexOf(prefix[0]);
          if (id.equals("") || index < 0) {
            if (!id.equals(""))
              id = id + ",";
            id = id + prefix[0] + "{" + prefix[1] + "}";
            continue;
          }

          index = id.indexOf('}', index);
          id = id.substring(0, index) + "," + prefix[1] + id.substring(index);
        } else if (id_new != null) {
          if (!id.equals(""))
            id = id + ",";
          id = id + id_new;
        }
      }
    }

    return id;
  }

  /**
   * Get the ID prefix, e.g. for SPAC1556.06.1:exon:2 returns SPAC1556.06.1:exon
   * as the prefix and 2 as the index.
   * 
   * @param id
   * @return
   */
  public String[] getPrefix(final String id, final char separator) {
    String prefix[] = new String[2];
    int index = id.lastIndexOf(separator);

    if (index > -1) {
      prefix[0] = id.substring(0, index);
      prefix[1] = id.substring(index + 1);
    }
    return prefix;
  }

  /**
   * Used to automatically generate
   * 
   * @param prefix
   * @return
   */
  public int getAutoNumber(final String prefix, final char separator) {
    int auto = 1;
    String val = prefix + separator + auto;
    while (id_range_store.containsKey(val)) {
      auto++;
      val = prefix + separator + auto;
    }
    return auto;
  }

  /**
   * Return the reference of a new copy of this Feature.
   **/
  public Feature copy() {
    final Feature return_value = new GFFStreamFeature(this);
    return return_value;
  }

  /**
   * Read and return a GFFStreamFeature from a stream. A feature must be the
   * next thing in the stream.
   * 
   * @param stream
   *          the Feature is read from this stream
   * @exception IOException
   *              thrown if there is a problem reading the Feature - most likely
   *              ReadFormatException.
   * @exception InvalidRelationException
   *              Thrown if this Feature cannot contain the given Qualifier.
   * @return null if in_stream is at the end of file when the method is called
   */
  protected static GFFStreamFeature readFromStream(LinePushBackReader stream)
      throws IOException, InvalidRelationException {
    final String line = stream.readLine();
    if (line == null)
      return null;

    try {
      return new GFFStreamFeature(line);
    } catch (ReadFormatException exception) {
      // re-throw the exception with the line number added
      final String new_error_string = exception.getMessage();

      throw new ReadFormatException(new_error_string, stream.getLineNumber());
    }
  }

  /**
   * Read the details of a feature from an EMBL stream into the current object.
   * 
   * @param entry_information
   *          The EntryInformation object of the Entry that will contain the
   *          Feature.
   * @param in_stream
   *          the Feature is read from this stream
   * @exception IOException
   *              thrown if there is a problem reading the Feature - most likely
   *              ReadFormatException if the stream does not contain GFF
   *              feature.
   **/
  public void setFromStream(final EntryInformation entry_information,
      final LinePushBackReader in_stream) throws IOException,
      InvalidRelationException, ReadOnlyException {
    throw new ReadOnlyException();
  }

  /**
   * Write this Feature to the given stream.
   * 
   * @param writer
   *          The stream to write to.
   * @exception IOException
   *              thrown if there is an io problem while writing the Feature.
   **/
  public void writeToStream(final Writer writer) throws IOException {
    final RangeVector ranges = getLocation().getRanges();
    final int ranges_size = ranges.size();

    // final Hashtable contig_ranges = SimpleDocumentEntry.getContigRanges();
    for (int i = 0; i < ranges_size; ++i) {
      Range this_range = (Range) ranges.elementAt(i);

      String seqname = getGffSeqName();
      String source = getGffSource();
      Qualifier score = getQualifierByName("score");
      Qualifier group = getQualifierByName("group");

      // source becomes a Dbxref in chado
      String source_str = null;
      if (getQualifierByName("Dbxref") != null) {
        source_str = getDbxrefGFFSource(getQualifierByName("Dbxref"));
      }

      int start = this_range.getStart();
      int end = this_range.getEnd();

      if (seqname == null
          && ((GFFDocumentEntry) getEntry()).getDocument() != null)
        seqname = ((GFFDocumentEntry) getEntry()).getDocument().getName();
      if (seqname == null)
        seqname = deriveSeqName(start);

      if (source == null)
        source = "artemis";

      if (score == null)
        score = new Qualifier("score", ".");

      if (seqname != null && contig_ranges != null
          && contig_ranges.containsKey(seqname)) {
        Range offset_range = contig_ranges.get(seqname);
        start = start - offset_range.getStart() + 1;
        end = end - offset_range.getStart() + 1;
      }

      if (group == null || group.getValues() == null
          || group.getValues().elementAt(0).equals("")) {
        final Qualifier gene = getQualifierByName("gene");

        if (gene == null)
          group = new Qualifier("group", "");
        else
          group = gene;
      }

      String frame = ".";
      final Qualifier codon_start = getQualifierByName("codon_start");

      if (codon_start != null) {
        frame = (String) (codon_start.getValues()).elementAt(0);

        if (frame.equals("1"))
          frame = "0";
        else if (frame.equals("2"))
          frame = "1";
        else if (frame.equals("3"))
          frame = "2";
        else
          frame = ".";
      }

      // phase is REQUIRED for all CDS features
      if (getKey().equals("CDS") && frame.equals("."))
        frame = "0";

      final String myId = getSegmentID(this_range);
      String attribute_string = unParseAttributes(myId);

      if (source_str == null && source != null)
        source_str = source;

      final String translation = getTranslation();
      if (translation != null)
        attribute_string = attribute_string + ";" + translation;
      writer.write(seqname + "\t" + source_str + "\t" + getKey().getKeyString()
          + "\t" + start + "\t" + end + "\t" + score.getValues().elementAt(0)
          + "\t" + (getLocation().isComplement() ? "-\t" : "+\t") + frame
          + "\t" + attribute_string + "\n");
    }
  }

  /**
   * If the seqname is not set for this feature try to derive the
   * contig/chromosome it is located on
   * 
   * @param start
   * @return
   */
  private String deriveSeqName(int start) {
    String seqname = null;
    if (contig_ranges != null) {
      final Enumeration<String> contigEnum = contig_ranges.keys();
      while (contigEnum.hasMoreElements()) {
        final String key = contigEnum.nextElement();
        final Range r = contig_ranges.get(key);
        if (r.getStart() > start)
          continue;
        if (r.getEnd() > start)
          return key;
      }
    } else {
      try {
        seqname = ((GFFStreamFeature) (getEntry().getAllFeatures().elementAt(0)))
            .getGffSeqName();
      } catch (Exception e) {
      }
    }

    if (seqname == null)
      seqname = "gff_seqname";
    return seqname;
  }

  /**
   * Return a String containing the qualifiers of this feature in a form
   * suitable for using as the last field of a GFF line.
   **/
  private String unParseAttributes(final String myId) {
    final QualifierVector qualifiers = getQualifiers();
    GFF3AttributeBuilder abuf = new GFF3AttributeBuilder();
    prepareProcessors(abuf);
    
    for (String attr : attrs_to_filter) {
      abuf.ignore(attr);
    }

    final int names_length = abuf.reserved_a.length;

    // add ID attribute
    if (myId != null) {
     abuf.add("ID", myId);
    }
    
    // build reserved attributes
    for (int i = 1; i < names_length; i++) {
      Qualifier this_qualifier = qualifiers.getQualifierByName(abuf.reserved_a[i]);
      
      if (this_qualifier == null)
        continue;

      abuf.add(this_qualifier.getName(), this_qualifier.getValues());
    }

    // build remaining attributes
    boolean lname;
    for (Qualifier this_qualifier : qualifiers) {
      lname = false;

      // skip reserved names
      for (int j = 0; j < names_length; j++)
        if (this_qualifier.getName().equals(abuf.reserved_a[j]))
          lname = true;
      if (lname)
        continue;

      // skip internal qualifiers
      if ((this_qualifier.getName().equals("private") && System
          .getProperty("noprivate") != null)
          || (this_qualifier.getName().equals("history") && System
              .getProperty("nohistory") != null))
        continue;

      abuf.add(this_qualifier.getName(), this_qualifier.getValues());
    }

    return abuf.toString();
  }

  void prepareProcessors(GFF3AttributeBuilder abuf) {
    GFF3AttributeAggregator productProc = new GFF3AttributeAggregator() {
      @Override
      public String process(StringVector values) {
        StringBuilder buffer = new StringBuilder();
        if (values != null && values.size() > 0) {
          for (int value_index = 0; value_index < values.size(); ++value_index) {
            final String this_value;
            int index = values.elementAt(value_index).indexOf("term=");
            // strip off the 'term=' etc
            if (index > -1)
              this_value = GFF3Encoder.encode(values.elementAt(value_index)
                  .substring(index + 5,
                      values.elementAt(value_index).length() - 1));
            else
              this_value = GFF3Encoder.encode(values.elementAt(value_index));
            if (value_index > 0 && value_index < (values.size())) {
              buffer.append(",");
            }
            buffer.append(this_value);
          }
        }
        return buffer.toString();
      }
    };

    GFF3AttributeAggregator ecProc = new GFF3AttributeAggregator() {
      @Override
      public String process(StringVector values) {
        StringBuilder buffer = new StringBuilder();
        if (values != null && values.size() > 0) {
          for (int value_index = 0; value_index < values.size(); ++value_index) {
            final String this_value = "EC:"
                + GFF3Encoder.encode(values.elementAt(value_index));
            if (value_index > 0 && value_index < (values.size())) {
              buffer.append(",");
            }
            buffer.append(this_value);
          }
        }
        return buffer.toString();
      }
    };

    GFF3AttributeAggregator psysIDProc = new GFF3AttributeAggregator() {
      @Override
      public String process(StringVector values) {
        StringBuilder buffer = new StringBuilder();
        if (values != null && values.size() > 0) {
          for (int value_index = 0; value_index < values.size(); ++value_index) {
            final String this_value;
            int index = values.elementAt(value_index).indexOf(";current=");
            if (index > -1)
              this_value = GFF3Encoder.encode(values.elementAt(value_index)
                  .substring(0, index));
            else
              this_value = GFF3Encoder.encode(values.elementAt(value_index));
            if (value_index > 0 && value_index < (values.size())) {
              buffer.append(",");
            }
            buffer.append(this_value);
          }
        }
        return buffer.toString();
      }
    };

    GFF3AttributeAggregator classProc = new GFF3AttributeAggregator() {
      @Override
      public String process(StringVector values) {
        StringBuilder buffer = new StringBuilder();
        if (values != null && values.size() > 0) {
          for (int value_index = 0; value_index < values.size(); ++value_index) {
            final String this_value;
            int index = values.elementAt(value_index).indexOf("::");
            if (index > -1)
              this_value = GFF3Encoder.encode(values.elementAt(value_index)
                  .substring(0, index));
            else
              this_value = GFF3Encoder.encode(values.elementAt(value_index));
            if (value_index > 0 && value_index < (values.size())) {
              buffer.append(",");
            }
            buffer.append(this_value);
          }
        }
        return buffer.toString();
      }
    };

    GFF3AttributeAggregator startEndRangeProc = new GFF3AttributeAggregator() {
      @Override
      public String process(StringVector values) {
        StringBuilder buffer = new StringBuilder();
         if (values != null && values.size() > 0) {
          for (int value_index = 0; value_index < values.size(); ++value_index) {
            if (value_index > 0 && value_index < (values.size())) {
              buffer.append(",");
            }
            buffer.append(values.elementAt(value_index));
          }
        }
        return buffer.toString();
      }
    };

    GFF3AttributeAggregator goProc = new GFF3AttributeAggregator() {
      @Override
      public String process(StringVector values) {
        StringBuilder buffer = new StringBuilder();
        if (values != null && values.size() > 0) {
          for (int value_index = 0; value_index < values.size(); ++value_index) {
            int goindex = values.elementAt(value_index).indexOf("GOid=");
            int termindex = values.elementAt(value_index).indexOf(";term=");
            if (goindex > -1 && termindex > -1) {
              buffer.append(GFF3Encoder.encode(values.elementAt(value_index)
                  .substring(goindex + 5, termindex)));
              if (value_index < (values.size()) - 1)
                buffer.append(",");
            }
          }
        }
        return buffer.toString();
      }
    };

    GFF3AttributeAggregator curcomProc = new GFF3AttributeAggregator() {
      @Override
      public String process(StringVector values) {
        StringBuilder buffer = new StringBuilder();
        if (values != null && values.size() > 0) {
          for (int value_index = 0; value_index < values.size(); ++value_index) {
            buffer.append(GFF3Encoder.encode(values.elementAt(value_index)));
            if (value_index < (values.size()) - 1)
              buffer.append(" ");
          }
        }
        return buffer.toString();
      }
    };

    // map GO -> full_GO
    abuf.setMapping("GO", "full_GO");
    abuf.setGlue("full_GO", ",");

    // merge curation and comment
    abuf.setMapping("curation", "comment");
    abuf.setGlue("comment", " ");
    abuf.setAggregator("comment", curcomProc);
    
    // also put GOs in Ontology_term
    abuf.setClone("full_GO", "Ontology_term");
    abuf.setAggregator("Ontology_term", goProc);
    abuf.setGlue("Ontology_term", ",");

    // class
    abuf.setAggregator("class", classProc);

    // EC numbers go into Dbxref
    abuf.setMapping("EC_number", "Dbxref");
    abuf.setAggregator("EC_number", ecProc);
    abuf.setGlue("Dbxref", ",");

    // start/end ranges
    abuf.setAggregator("Start_range", startEndRangeProc);
    abuf.setAggregator("End_range", startEndRangeProc);

    // previous_systematic_id
    abuf.setAggregator("previous_systematic_id", psysIDProc);
    
    // product
    abuf.setAggregator("product", productProc);
  }

  /**
   * Get the translation qualifier string for polypeptide features.
   */
  private String getTranslation() {
    if (!getKey().getKeyString().equals("polypeptide"))
      return null;
    if (chadoGene != null) {
      if (getUserData() == null)
        new uk.ac.sanger.artemis.Feature(this);
      // the above line constructs the appropriate userData within this current
      // GFFStreamFeature object,
      // which is required by the following GeneUtils.deriveResidues()
      String residues = GeneUtils.deriveResidues(this);
      if (residues != null)
        return "translation=" + residues;
    }
    return null;
  }

  /**
   * Parse the given String as ACeDB format attributes. Adapted from code by
   * Matthew Pocock for the BioJava project.
   * 
   * Modified for gff-version 3.
   * 
   * @return Return a Hashtable. Each key is an attribute name and each value of
   *         the Hashtable is a StringVector containing the attribute values. If
   *         the attribute has no value then the Hashtable value will be a zero
   *         length vector.
   **/
  private Hashtable<String, StringVector> parseAttributes(
      final String att_val_list) {
    final Hashtable<String, StringVector> attr = new Hashtable<String, StringVector>();

    int ind_start = 0;
    int ind_end;
    while ((ind_end = att_val_list.indexOf(";", ind_start)) > -1
        || ind_start < att_val_list.length()) {
      if (ind_end < 0)
        ind_end = att_val_list.length();

      final String this_token = GFF3Encoder.decode(att_val_list.substring(ind_start,
          ind_end).trim());
      ind_start = ind_end + 1;

      int index_of_first_space = this_token.indexOf(" ");

      final String att_name;
      StringVector att_values = new StringVector();

      if (this_token.indexOf("=") > -1
          && (this_token.indexOf("=") < index_of_first_space || index_of_first_space == -1)) {
        index_of_first_space = this_token.indexOf("=");
        att_name = this_token.substring(0, index_of_first_space);
        att_values.add(this_token.substring(index_of_first_space + 1).trim());
      } else if (index_of_first_space == -1)
        att_name = this_token;
      else {
        att_name = this_token.substring(0, index_of_first_space);

        String rest_of_token = this_token.substring(index_of_first_space + 1)
            .trim();

        while (rest_of_token.length() > 0) {
          if (rest_of_token.startsWith("\"")) {
            int quote_index = 0;
            do {
              quote_index++;
              quote_index = rest_of_token.indexOf("\"", quote_index);
            } while (quote_index > -1
                && rest_of_token.charAt(quote_index - 1) == '\\');

            if (quote_index < 0) {
              // no closing quote - panic
              final Hashtable<String, StringVector> panic_attributes = new Hashtable<String, StringVector>();
              final StringVector notes = new StringVector();
              notes.add(att_val_list);
              panic_attributes.put("note", notes);

              return panic_attributes;
            }

            final String next_bit = rest_of_token.substring(1, quote_index);
            att_values.add(next_bit);
            rest_of_token = rest_of_token.substring(quote_index + 1).trim();
          } else {
            final int index_of_next_space = rest_of_token.indexOf(" ");

            if (index_of_next_space == -1) {
              att_values.add(rest_of_token);
              rest_of_token = "";
            } else {
              final String next_bit = rest_of_token.substring(0,
                  index_of_next_space);

              att_values.add(next_bit);
              rest_of_token = rest_of_token.substring(index_of_next_space)
                  .trim();
            }
          }
        }

        if (!rest_of_token.equals(""))
          att_values.add(rest_of_token);
      }

      if (att_name.equals("Dbxref") || att_name.equals("Alias")) // convert to
                                                                 // multi-line
      {
        StringTokenizer stok = new StringTokenizer((String) att_values.get(0),
            ",");
        StringVector str_values = new StringVector();
        while (stok.hasMoreTokens())
          str_values.add(stok.nextToken());

        att_values = str_values;
      }

      if (att_name.equals("timelastmodified")) {
        try {
          this.timelastmodified = new Timestamp(
              Long.parseLong((String) att_values.get(0)));
          SimpleDateFormat date_format = new SimpleDateFormat(
              "dd.MM.yyyy hh:mm:ss z");
          att_values.set(0, date_format.format(timelastmodified));
        } catch (NumberFormatException e) {
          att_values.set(0, (String) att_values.get(0));
        }
      }

      if (attr.get(att_name) != null)
        attr.get(att_name).add(att_values);
      else
        attr.put(att_name, att_values);
    }

    return attr;
  }

  /**
   * Get the feature time last modified timestamp.
   * 
   * @return
   */
  public Timestamp getLastModified() {
    return timelastmodified;
  }

  /**
   * Get the GFF_source value of a Dbxref qualifier.
   * 
   * @param qualifier
   * @return the gff_source value or NULL
   */
  private String getDbxrefGFFSource(final Qualifier qualifier) {
    StringVector qualifier_strings = StreamQualifier.toStringVector(null,
        qualifier);

    for (int i = 0; i < qualifier_strings.size(); i++) {
      String qualifier_string = (String) qualifier_strings.elementAt(i);

      if (qualifier_string.indexOf("GFF_source:") > -1) {
        int index = qualifier_string.indexOf(":") + 1;
        int len = qualifier_string.length();
        if (qualifier_string.endsWith("\""))
          len--;
        return qualifier_string.substring(index, len);
      }
    }
    return null;
  }

  /**
   * Set the feature time last modified timestamp.
   * 
   * @param timelastmodified
   */
  public void setLastModified(final Timestamp timelastmodified) {
    this.timelastmodified = timelastmodified;

    // now update the qualifier value itself
    QualifierVector qualifiers = getQualifiers();
    Qualifier qualifier = qualifiers.getQualifierByName("timelastmodified");
    SimpleDateFormat date_format = new SimpleDateFormat("dd.MM.yyyy hh:mm:ss z");

    if (qualifier != null)
      qualifier.removeValue((String) qualifier.getValues().get(0));
    else {
      try {
        qualifier = new Qualifier("timelastmodified",
            date_format.format(timelastmodified));
        setQualifier(qualifier);
        return;
      } catch (EntryInformationException eie) {
      } catch (ReadOnlyException roe) {
      }
    }

    qualifier.addValue(date_format.format(timelastmodified));
  }

  /**
   * Returns true if and only if this Feature can't be changed or can't be
   * removed from it's entry.
   **/
  public boolean isReadOnly() {
    if (readOnlyFeature)
      return true;
    return super.isReadOnly();
  }

  public void setReadOnlyFeature(boolean readOnlyFeature) {
    this.readOnlyFeature = readOnlyFeature;
  }

  public ChadoCanonicalGene getChadoGene() {
    return chadoGene;
  }

  public void setChadoGene(ChadoCanonicalGene chadoGene) {
    this.chadoGene = chadoGene;
  }

  public boolean isVisible() {
    return visible;
  }

  public void setVisible(boolean visible) {
    this.visible = visible;
  }

  public String getGffSeqName() {
    return gffSeqName;
  }

  public void setGffSeqName(String gffSeqName) {
    this.gffSeqName = gffSeqName;
  }

  public String getGffSource() {
    return gffSource;
  }

  public void setGffSource(String gffSource) {
    this.gffSource = gffSource;
  }

  public boolean isLazyLoaded() {
    return lazyLoaded;
  }

  public void setLazyLoaded(boolean lazyLoaded) {
    this.lazyLoaded = lazyLoaded;
  }

  public org.gmod.schema.sequence.Feature getChadoLazyFeature() {
    return chadoLazyFeature;
  }

  public void setChadoLazyFeature(
      org.gmod.schema.sequence.Feature chadoLazyFeature) {
    this.chadoLazyFeature = chadoLazyFeature;
  }

  protected static boolean isGTF(Feature feature) {
    if (!(feature instanceof GFFStreamFeature))
      return false;

    final String names[] = { "ID", "Name", "Alias", "Parent", "Derives_from",
        "Target", "Gap", "Note", "Dbxref", "Ontology_term" };

    for (String name : names) {
      if (feature.getQualifiers().getQualifierByName(name) != null)
        return false;
    }

    if (feature.getQualifiers().getQualifierByName("gene_id") != null
        && feature.getQualifiers().getQualifierByName("transcript_id") != null) {
      if (feature.getEntry() != null)
        logger4j.debug(feature.getEntry().getName() + " is in GTF format");
      return true;
    }
    return false;
  }
}
