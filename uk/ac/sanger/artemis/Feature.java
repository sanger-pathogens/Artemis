/* Feature.java
 *
 * created: Sun Oct 11 1998
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1998,1999,2000  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/Feature.java,v 1.35 2009-02-03 11:36:39 tjc Exp $
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.io.OutOfDateException;
import uk.ac.sanger.artemis.io.LocationParseException;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.PartialSequence;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.RangeVector;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.DateStampFeature;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.EmblStreamFeature;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.FastaStreamSequence;
import uk.ac.sanger.artemis.io.StreamFeature;
import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.plot.*;

import java.awt.Color;
import java.util.Vector;
import java.io.*;
import java.util.Date;
import java.util.regex.Pattern;

/**
 *  This class extends an embl.Feature with the other information needed by
 *  Diana.  It also able to send events to other objects that are interested
 *  in changes to this object.  (see FeatureChangeEvent details of the
 *  possible change events.)  To make changes to the feature it calls methods
 *  in the Entry class.  Changes to this object will update the underlying
 *  embl.Feature and embl.Entry objects.
 *
 *  @author Kim Rutherford
 *  @version $Id: Feature.java,v 1.35 2009-02-03 11:36:39 tjc Exp $
 **/

public class Feature
    implements EntryChangeListener, Selectable, SequenceChangeListener,
               MarkerChangeListener, OptionChangeListener
{

  /**
   *  The Entry that controls/contains this feature.  This is the Entry object
   *  that was passed to the constructor.
   **/
  private Entry entry;

  /**
   *  This is a reference to the low level Feature object that this class is
   *  providing a wrapper for.  This is the embl.Feature object that was
   *  passed to constructor.
   **/
  private uk.ac.sanger.artemis.io.Feature embl_feature;

  /**
   *  The segment of this Feature.  All features have at least one segment.
   **/
  private FeatureSegmentVector segments = null;

  /**
   *  A vector of those objects listening for feature change events.
   **/
  private final Vector feature_listener_list = new Vector();

  /**
   *  The translation of the bases of this feature.
   **/
  private AminoAcidSequence amino_acids = null;

  /**
   *  The bases of this feature.
   **/
  private String bases = null;

  /**
   *  The count of the number of amino acids in this feature. (set by
   *  getAACount() and resetCache()).
   **/
  private int aa_count = -1;

  /**
   *  The count of the bases in this feature (set by getBaseCount() and
   *  resetCache()).
   **/
  private int base_count = -1;

  /**
   *  This array contains counts of the occurrences of each three bases.  The
   *  first index is the first base, the second index is the second base, etc.
   *  The indices are the same as those for Bases.letter_index.
   *  (set by resetCache()).
   **/
  private int [][][] codon_counts = null;

  /**
   *  This array contains counts of the occurrences of each amino acids in the
   *  translation of the bases of the current feature.
   *  (set by resetCache()).
   **/
  private int [] residue_counts = null;

  /**
   *  This array contains counts of number of each base that appear in each
   *  codon position.  The first index is the codon position and the second is
   *  the base (indexed in the same way as Bases.letter_index).
   *  (set by resetCache()).
   **/
  private int [][] positional_base_counts = null;

  /**
   *  This array contains the counts of the total number of each base in the
   *  feature.
   *  The indices are the same as those for Bases.letter_index.
   *  (set by resetCache()).
   **/
  private int [] base_counts = null;

  /**
   *  The current Location reference is saved each time setLocation() is
   *  called so that if the reference changes resetCache() can
   *  be called.  This is needed to handle RWCorbaFeature objects which can
   *  change their Location at arbitrary times.
   **/
  private Location old_location = null;

  /**
   *  Incremented when startListening() is called - decremented when
   *  stopListening() is called.
   **/
  private int listen_count = 0;
 
  private Color colour = null;
 
  /**
   *  Create a new Feature object.
   *  @param entry The uk.ac.sanger.artemis.Entry object that contains this Feature.
   *  @param embl_feature The embl.Feature object that this class is
   *    providing a wrapper for.
   **/
  public Feature(uk.ac.sanger.artemis.io.Feature embl_feature) 
  {
    this.embl_feature = embl_feature;
    embl_feature.setUserData(this);
    old_location = embl_feature.getLocation();
  }

  /**
   *  Implementation of the EntryChangeListener interface.  We listen to
   *  EntryChange events so that if the feature is deleted we can destroy any
   *  objects that use it.
   *
   *  NOTE: not active currently
   **/
  public void entryChanged(EntryChangeEvent event) 
  {
  }

  /**
   *  This method fixes up the location of this Feature when a Marker changes.
   **/
  public void markerChanged(final MarkerChangeEvent event) 
  {
    try 
    {
      final Location old_location = getLocation();
      updateEMBLFeatureLocation();
      locationChanged(old_location);
    }
    catch(ReadOnlyException e) {}
  }


  /**
   *  This method fixes up the location of this Feature when a sequence
   *  changes.
   **/
  public void sequenceChanged(final SequenceChangeEvent event) 
  {
    try 
    {
      // we don't send a FeatureChangeEvent because the logical location
      // hasn't changed

      // all the Markers for the ranges of this feature will have changed
      // before this method is called because the markers are added as
      // SequenceChangeListener with a higher priority than this Feature

      if(event.getType() == SequenceChangeEvent.REVERSE_COMPLEMENT) 
        reverseComplement(getEntry().getBases().getLength());
      else if(event.getType() == SequenceChangeEvent.CONTIG_REVERSE_COMPLEMENT)
      {
         final Location old_location = getLocation();

        // if the event is contained within this feature then the feature
        // sequence may have changed
        final Range this_feature_range = getMaxRawRange();

        boolean feature_changed = false;
        Range eventRange = event.getRange();
        
        // reverse complement feature if within contig region
        if(eventRange.getStart() <= this_feature_range.getStart() &&
           eventRange.getEnd() >= this_feature_range.getEnd())
        {
          try
          {
            final Location new_location =
               getLocation().reverseComplement(event.getLength(), eventRange.getStart());
     
            setLocationInternal(new_location);
          }
          catch(OutOfRangeException e)
          {
            throw new Error("internal error - inconsistent location: " + e);
          }
        } 
      }
      else if(event.getType() == SequenceChangeEvent.CONTIG_REORDER)
      {
        final Location old_location = getLocation();

        final int new_base_pos = event.getPosition();
        final int range_start  = event.getRange().getStart(); 
        final int range_end = event.getRange().getEnd();
        final Range this_feature_range = getMaxRawRange();
        
        // check if feature is effected
        if( (this_feature_range.getStart() >= new_base_pos ||
             this_feature_range.getStart() >= range_start) &&
            (this_feature_range.getStart() < new_base_pos ||
             this_feature_range.getStart() < range_end) )
        {
          try
          {
            final int diff;
            
            if(range_start <= this_feature_range.getStart() &&
               range_end   >= this_feature_range.getEnd())
            {
              if(new_base_pos < range_start)
                diff = new_base_pos-range_start;
              else
                diff = new_base_pos-range_end-1;
            }
            else
            {
              if(this_feature_range.getStart() < range_start)
                diff = range_end-range_start+1;
              else
                diff = range_start-range_end-1;
            }

            Location new_location = moveSegments(diff);
            setLocationInternal(new_location);
          }
          catch(OutOfRangeException e)
          {
            throw new Error("internal error - inconsistent location: " + e);
          }
        }
      }
      else 
      {
        final Location old_location = getLocation();

        // if the event is contained within this feature then the feature
        // sequence may have changed
        final Range this_feature_range = getMaxRawRange();

        boolean feature_changed = false;
        int eventPosition = event.getPosition();

        if(eventPosition >= this_feature_range.getStart() &&
           eventPosition <= this_feature_range.getEnd() + 1) 
        {
          // check each segment
          final FeatureSegmentVector segments = getSegments();
          int seg_size = segments.size();

          for(int i = 0; i < seg_size; ++i) 
          {
            final FeatureSegment this_segment = segments.elementAt(i);

            final Range this_segment_range = this_segment.getRawRange();

            if(eventPosition >= this_segment_range.getStart() &&
               eventPosition <= this_segment_range.getEnd() + 1) 
              feature_changed = true;
          }
        }

        updateEMBLFeatureLocation();

        if(feature_changed) 
        {
          resetCache();
          locationChanged(old_location);
        }
      }
    }
    catch(ReadOnlyException e)
    {
      throw new Error("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Invoked when an Option is changed.
   **/
  public void optionChanged(OptionChangeEvent event) 
  {
    // if the eukaryotic mode option changes the sequence may change (ie. the
    // start codon may need to be translated differently) - see
    // getTranslation()
    locationChanged(getLocation());
  }

  /**
   *  Returns the embl feature that was passed to the constructor.
   **/
  public uk.ac.sanger.artemis.io.Feature getEmblFeature() 
  {
    return embl_feature;
  }

  /**
   *  Write this Feature to the given stream in the native format of this
   *  Feature.  The output will be in Genbank format if this Feature is a
   *  Genbank feature and EMBL format otherwise.
   *  @param writer The record is written to this Writer.
   **/
  public void writeNative(final Writer writer)
      throws IOException 
  {
    if(getEmblFeature() instanceof StreamFeature) 
      ((StreamFeature)getEmblFeature()).writeToStream(writer);
    else 
    {
      final EntryInformation entry_info = getEntry().getEntryInformation();

      // this is a hack to make the correct EntryInformation object available
      // to the feature writing code.
      final uk.ac.sanger.artemis.io.EmblDocumentEntry document_entry =
        new uk.ac.sanger.artemis.io.EmblDocumentEntry(entry_info);

      final uk.ac.sanger.artemis.io.Feature returned_feature =
        document_entry.forcedAdd(new EmblStreamFeature(getEmblFeature()));

      ((EmblStreamFeature)returned_feature).writeToStream(writer);
    }
  }

  /**
   *  Write a PIR database record of this feature to the given Writer.
   *  @param writer The record is written to this Writer.
   **/
  public void writePIROfFeature(final Writer writer) 
  {
    final String gene_name = getGeneName();

    final String pir_name;

    if(gene_name == null) 
    {
      if(getLabel() == null) 
        pir_name = getKey().toString();
      else 
        pir_name = getLabel();
    }
    else 
      pir_name = gene_name;

    final String header_line =
      ">BL;" + pir_name + ", " +
      getEntry().getName() + " " +
      getWriteRange() +
      " MW:" + (int) getMolecularWeight();

    final PrintWriter print_writer = new PrintWriter(writer);

    print_writer.println(header_line);

    final String translation_string =
      getTranslation().toString().toUpperCase();

    wrapAndWrite(print_writer, translation_string, 80);
    print_writer.println("*");
    print_writer.flush();
  }

  /**
   *  Write the bases of this feature to the given Writer.
   *  @param writer The bases are written to this Writer.
   **/
  public void writeBasesOfFeature(final Writer writer)
      throws IOException 
  {
    final FastaStreamSequence stream_sequence =
      new FastaStreamSequence(getBases(),
                              getIDString() + ", " +
                              getWriteRange());

    stream_sequence.writeToStream(writer);
  }

  /**
   *  Write the amino acid symbols of this feature to the given Writer.
   *  @param writer The amino acids are written to this Writer.
   **/
  public void writeAminoAcidsOfFeature(final Writer writer)
      throws IOException 
  {
    final StringBuffer header_buffer = new StringBuffer(">");

    header_buffer.append(getSystematicName());
    header_buffer.append(" ");
    header_buffer.append(getIDString());
    header_buffer.append(" ");

    final String product = getProductString();

    if(product == null)
      header_buffer.append("undefined product");
    else
      header_buffer.append(product);

    header_buffer.append(" ").append(getWriteRange()).append(" MW:");
    header_buffer.append((int) getMolecularWeight());

    final PrintWriter print_writer = new PrintWriter(writer);

    print_writer.println(header_buffer);

    final String translation_string =
      getTranslation().toString().toUpperCase();

    wrapAndWrite(print_writer, translation_string, 60);

    print_writer.flush();
  }

  /**
   *  Return a String containing the bases upstream of this feature.
   *  @param count The number of (immediately) upstream bases to return.
   **/
  public String getUpstreamBases(final int count) 
  {
    final int feature_start_base = getFirstBase();
    final int start_base;
    final int end_base;

    if(feature_start_base == 1) 
    {
      // there are no bases before this feature
      return "";
    }
    else
    {
      end_base = feature_start_base - 1;

      if(feature_start_base > count) 
        start_base = feature_start_base - count;
      else 
        start_base = 1;
    }

    final String bases_string;

    try 
    {
      bases_string =
        getStrand().getSubSequence(new Range(start_base, end_base));
    } 
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }

    return bases_string;
  }

  /**
   *  Return a String containing the bases downstream of this feature.
   *  @param count The number of (immediately) downstream bases to return.
   **/
  public String getDownstreamBases(final int count) 
  {
    final int feature_end_base = getLastBase();
    final int start_base;
    final int end_base;
    final int sequenceLength = getSequenceLength();
    
    if(feature_end_base == sequenceLength)
    {
      // there are no bases after this feature
      return "";
    }
    else 
    {
      start_base = feature_end_base + 1;

      if(sequenceLength - feature_end_base > count) 
        end_base = feature_end_base + count;
      else 
        end_base = sequenceLength;
    }

    final String bases_string;

    try 
    {
      bases_string =
        getStrand().getSubSequence(new Range(start_base, end_base));
    } 
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }

    return bases_string;
  }

  /**
   *  Helper method for writePIROfFeature().
   *  @return A string of the form "100:200 reverse" or "1:2222 forward".
   **/
  public String getWriteRange()
  {
    String partial = " ";
    if(isPartial(true))
      partial += "partial 5' ";
    if(isPartial(false))
      partial += "partial 3' ";

    return (isForwardFeature() ?
       getFirstCodingBaseMarker().getRawPosition() + ":" +
       getLastBaseMarker().getRawPosition() + partial + "forward" :
       getLastBaseMarker().getRawPosition() + ":" +
       getFirstCodingBaseMarker().getRawPosition() + partial + "reverse");
  }
  
  /**
   * If lookAt5prime is set to true then only return true if the 5' end is 
   * partial otherwise only return true if the 3' end is partial.
   * @param lookAt5prime
   * @return
   */
  private boolean isPartial(final boolean lookAt5prime)
  {
    try
    {
      boolean isDatabaseFeature = GeneUtils.isDatabaseEntry(getEmblFeature());
      if(isDatabaseFeature)
      {
        if(lookAt5prime)
        {
          if(isForwardFeature())
          {
            if(getQualifierByName("Start_range") != null)
              return true;
          }
          else if(getQualifierByName("End_range") != null)
              return true;
        }
        else
        {
          if(isForwardFeature())
          {
            if(getQualifierByName("End_range") != null)
              return true;
          }
          else if(getQualifierByName("Start_range") != null)
              return true;
        }
      }
    } catch (Exception e) {}
    return getLocation().isPartial(lookAt5prime);
  }

  /**
   *  Write a String object to a PrintWriter object, wrapping it a the given
   *  coloumn.
   *  @param writer The String will be written to this object.
   *  @param string The String to write.
   *  @param wrap_column The lines of output will be no longer than this many
   *    characters.
   **/
  private void wrapAndWrite(final PrintWriter writer,
                            final String string,
                            final int wrap_column) 
  {
    String remaining_string = string;

    while(remaining_string.length() > 0) 
    {
      int last_index = wrap_column;

      if(wrap_column > remaining_string.length()) 
        last_index = remaining_string.length();

      final String write_string = remaining_string.substring(0, last_index);
      writer.println(write_string);
      remaining_string = remaining_string.substring(last_index);
    }
  }

  /**
   *  Return this Feature as a EMBL, Genbank or GFF formatted String
   *  (depending on the type of this Feature).
   **/
   public String toString() 
   {
     final StringWriter string_writer = new StringWriter();

     try 
     {
       writeNative(string_writer);
     } 
     catch(IOException e) 
     {
       throw new Error("internal error - unexpected exception: " + e);
     }

     return string_writer.toString() ;
   }

  /**
   *  Return a Reader object that gives an EMBL format version of this
   *  Feature when read.
   **/
  public Reader toReader()
  {
    return new StringReader(toString());
  }

  /**
   *  Return true if and only if this Feature is on the forward strand.
   **/
  public boolean isForwardFeature() 
  {
    if(getLocation().isComplement()) 
      return false;
    else
      return true;
  }

  /**
   *  Return the key of this Feature.
   **/
  public Key getKey()
  {
    return getEmblFeature().getKey();
  }

  /**
   *  Return true if and only if the key of this feature is protein feature -
   *  one that should be displayed on the translation lines of the
   *  FeatureDisplay component rather than on the forward or backward strand.
   **/
  public boolean isProteinFeature() 
  {
    if(getKey().toString().startsWith("CDS") ||
       getKey().equals(DatabaseDocument.EXONMODEL) ||
       getKey().equals("exon") ||
       getKey().equals("BLASTCDS") || 
       getKey().equals("polypeptide"))
      return true;
    else
      return false;
  }

  /**
   *  Return true if and only if the key of this feature is CDS feature.
   **/
  public boolean isCDS() 
  {
    if(getKey().equals("CDS")) 
      return true;
    else
      return false;
  }

  /**
   *  Return true if and only if the key of this feature is CDS feature and
   *  the feature has a /partial qualifier.
   **/
  private boolean isPartialCDS() 
  {
    try 
    {
      if(getKey().equals("CDS") && getQualifierByName("partial") != null) 
        return true;
      else
        return false;
    }
    catch(InvalidRelationException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Return the Location of this Feature.
   **/
  public Location getLocation()
  {
    final Location current_location = getEmblFeature().getLocation();
    return current_location;
  }

  /**
   *  Return a Vector containing the qualifiers of this Feature.
   *  XXX - FIXME - should return a copy or be private
   **/
  public QualifierVector getQualifiers()
  {
    // return the embl.QualifierVector from the underlying embl.Feature object
    return getEmblFeature().getQualifiers();
  }

  /**
   *  Return the Entry that owns this object.
   **/
  public Entry getEntry() 
  {
    return entry;
  }

  /**
   *  Return the value of the first /codon_start qualifier as an integer or
   *  returns 1 if codon_start doesn't exist or if it makes no sense.
   **/
  public int getCodonStart()
  {
    try 
    {
      final String codon_start_string = getValueOfQualifier("codon_start");

      if(codon_start_string == null) 
      {
        // default value
        return 1;
      }

      if(codon_start_string.equals("2")) 
        return 2;
      else 
      {
        if(codon_start_string.equals("3")) 
          return 3;
        else 
          return 1;  // default value
      }
    } 
    catch(InvalidRelationException e) 
    {
      return 1;
    }
  }

  /**
   *  Return the value of the first /score qualifier of this feature as an
   *  integer.
   *  @return The score or -1 if the feature has no /score qualifier.
   **/
  public int getScore() 
  {
    try
    {
      final String score_string = getValueOfQualifier("score");

      if(score_string == null) 
      {
        // default value
        return -1;
      }

      try 
      {
        final int score_int = Float.valueOf(score_string).intValue();

        if(score_int > 100) 
          return 100;

        if(score_int < 0)
          return 0;

        return score_int;
      }  
      catch(NumberFormatException e) 
      {
        // assume there is no /score

        return -1;
      }
    } 
    catch(InvalidRelationException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Set the Entry that owns this object.
   **/
  public void setEntry(final Entry entry) 
  {
    if(this.entry != null)
      stopListening();

    final Entry old_entry = this.entry;

    this.entry = entry;

    if(old_entry == entry) 
      return;
    else 
    {
      if(old_entry != null)
        removeFeatureChangeListener(old_entry);

      if(entry != null)
      {
        // the Entry object acts as a proxy for FeatureChange events, other
        // objects can can addFeatureChangeListener() once on the Entry object
        // instead of calling it for every Feature.
        addFeatureChangeListener(getEntry());
      }
    }

    if(this.entry != null)
      startListening();
  }

  /**
   *  Change this Feature and it's underlying embl Feature object to the given
   *  key, location and qualifiers.
   *  @exception InvalidRelationException Thrown if this Feature cannot contain
   *    the given Qualifier.
   *  @exception OutOfRangeException Thrown if the location is out of
   *    range for this Entry.
   *  @exception ReadOnlyException Thrown if this Feature cannot be changed.
   **/
  public void set(Key new_key,
                  Location new_location,
                  QualifierVector new_qualifiers)
      throws EntryInformationException, OutOfRangeException,
             ReadOnlyException
  {
    try
    {
      set((Date)null, new_key, new_location, new_qualifiers);
    } 
    catch(OutOfDateException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Change this Feature and it's underlying embl Feature object to the given
   *  key, location and qualifiers.
   *  @exception InvalidRelationException Thrown if this Feature cannot contain
   *    the given Qualifier.
   *  @exception OutOfRangeException Thrown if the location is out of
   *    range for this Entry.
   *  @exception ReadOnlyException Thrown if this Feature cannot be changed.
   *  @exception OutOfDateException If the key has changed in the server since
   *    the time given by datestamp.  If datestamp argument is null then this
   *    exception will never be thrown.
   **/
  public void set(final Date datestamp,
                  final Key new_key,
                  final Location new_location,
                  final QualifierVector new_qualifiers)
      throws EntryInformationException, OutOfRangeException,
             ReadOnlyException, OutOfDateException 
  {
    final Key old_key = getKey();
    old_location = getLocation();
    final QualifierVector old_qualifiers = getQualifiers().copy();

    final int sequence_length;
    
    if(!(getEntry().getBases().getSequence() instanceof PartialSequence))
    {
      sequence_length = getEntry().getBases().getLength();

      if(new_location != null)
      {
        final Range span = new_location.getTotalRange();

        if(span.getEnd() > sequence_length || span.getStart() < 1)
        {
          throw new OutOfRangeException(new_location.toString());
        }
      }
    }
  
    if(datestamp == null ||
       !(getEmblFeature() instanceof DateStampFeature)) 
    {
      getEmblFeature().set(new_key, new_location, new_qualifiers);
    } 
    else
    {
      ((DateStampFeature)getEmblFeature()).set(datestamp,
                                               new_key,
                                               new_location,
                                               new_qualifiers);
    }

    resetCache();

    if(new_location != old_location && segments != null) 
      reexamineSegments();

    // now inform the listeners that a change has occured
    final FeatureChangeEvent event =
      new FeatureChangeEvent(this,
                             this,
                             old_key,
                             old_location,
                             old_qualifiers,
                             FeatureChangeEvent.ALL_CHANGED);

    fireAction(feature_listener_list, event);
  }

  /**
   *  This method will send a FeatureChangeEvent with type LOCATION_CHANGED to
   *  all the FeatureChangeEvent listeners when a marker changes.  It also
   *  calls resetCache(), because changing the location will change the
   *  translation and bases of the feature.
   **/
  private void locationChanged(final Location old_location,
                               final QualifierVector qualifiers,
                               int type) 
  {
    resetCache();

    // now inform the listeners that a change has occured
    final FeatureChangeEvent feature_change_event =
      new FeatureChangeEvent(this,
                             this,
                             null,
                             old_location,
                             qualifiers,
                             type);

    this.old_location = getLocation();

    fireAction(feature_listener_list, feature_change_event);
  }

  private void locationChanged(final Location old_location) 
  {
    locationChanged(old_location, null,
                    FeatureChangeEvent.LOCATION_CHANGED);  
  }
  
  /**
   *  Add the values from the given qualifier to the Qualifier object with the
   *  same name in this Feature or if there is no Qualifier with that name
   *  just add a copy of the argument.
   *  @param qualifier This object contians name and values to add.
   *  @return The Qualifier that was changed or created.
   **/
  public Qualifier addQualifierValues(Qualifier qualifier)
      throws EntryInformationException, ReadOnlyException 
  {

    final QualifierVector old_qualifiers = getQualifiers().copy();
    final Qualifier return_qualifier;

    final Qualifier current_qualifier =
      getEmblFeature().getQualifierByName(qualifier.getName());

    if(current_qualifier == null)
      return_qualifier = qualifier.copy();
    else
    {
      return_qualifier = current_qualifier.copy();
      return_qualifier.addValues(qualifier.getValues());
    }

    setQualifier(return_qualifier);

    // now inform the listeners that a change has occured
    // THIS IS ALREADY DONE IN setQualifier()
/*    final FeatureChangeEvent event =
      new FeatureChangeEvent(qualifier,
                             this,
                             null,
                             null,
                             old_qualifiers,
                             FeatureChangeEvent.QUALIFIER_CHANGED);

    fireAction(feature_listener_list, event);
*/
    return return_qualifier;
  }

  /**
   *  Add the given Qualifier to this Feature, replacing any exisiting
   *  Qualifier that have the same name.
   *  @param qualifier The Qualifier to add.
   **/
  public void setQualifier(final Qualifier qualifier)
      throws EntryInformationException, ReadOnlyException 
  {
    final QualifierVector old_qualifiers = getQualifiers().copy();

    getEmblFeature().setQualifier(qualifier);

    // now inform the listeners that a change has occured
    final FeatureChangeEvent event =
      new FeatureChangeEvent(qualifier,
                             this,
                             null,
                             null,
                             old_qualifiers,
                             FeatureChangeEvent.QUALIFIER_CHANGED);

    if(qualifier.getName().equals("transl_except"))
    {
      // discard cache
      amino_acids = null;
    }

    fireAction(feature_listener_list, event);
  }

  /**
   *  Remove the Qualifier with the given name.  If there is no Qualifier with
   *  that name then return immediately.
   *  @param name The qualifier name to look for.
   **/
  public void removeQualifierByName(final String name)
      throws EntryInformationException, ReadOnlyException, OutOfDateException 
  {
    if(getEmblFeature().getQualifierByName(name) == null)
    {
      // nothing to remove
      return;
    }

    final QualifierVector old_qualifiers = getQualifiers().copy();

    getEmblFeature().removeQualifierByName(name);

    // now inform the listeners that a change has occured
    final FeatureChangeEvent event =
      new FeatureChangeEvent(this,
                             this,
                             null,
                             null,
                             old_qualifiers,
                             FeatureChangeEvent.QUALIFIER_CHANGED);

    if(name.equals("transl_except"))
    {
      // discard cache
      amino_acids = null;
    }

    fireAction(feature_listener_list, event);
  }

  /**
   *  Return the time when this feature last changed or null if this Feature
   *  doesn't support datestamps.
   **/
  public Date getDatestamp() 
  {
    if(getEmblFeature() instanceof DateStampFeature) 
      return ((DateStampFeature)getEmblFeature()).getDatestamp();
    else 
      return null;
  }

  
  /**
   *  Return true if and only if any qualifier in this feature contains the
   *  given text string.
   *  @param search_text The text to search for.
   *  @param fold_case If true then the text comparisons will ignore case.
   *  @param match_substring If true then matches to substrings are allowed.
   *  @param qualifier_names If null search all qualifiers, otherwise just
   *    search these names,
   **/
  public boolean containsText(final String search_text,
                              final boolean fold_case,
                              final boolean match_substring,
                              final StringVector qualifier_names) 
  {
    return findOrReplaceText(search_text, fold_case, match_substring, false,
                             qualifier_names, null);
  }

  
  /**
   *  Return true if and only if any qualifier in this feature contains the
   *  given text string.
   *  @param search_text The text to search for.
   *  @param fold_case If true then the text comparisons will ignore case.
   *  @param match_substring If true then matches to substrings are allowed.
   *  @param qualifier_names If null search all qualifiers, otherwise just
   *    search these names,
   *  @param replaceText text to replace all qualifier value matches. If null
   *    then returns true if text found.
   **/
  public boolean findOrReplaceText(final String search_text,
                              final boolean fold_case,
                              final boolean match_substring,
                              final boolean deleteQualifier,
                              final StringVector qualifier_names,
                              final String replaceText) 
  {
    final String real_search_text;

    if(fold_case) 
      real_search_text = search_text.toLowerCase();
    else 
      real_search_text = search_text;

    final QualifierVector qualifiers = getQualifiers();
    QualifierVector newQualifiers = null;
    Vector<String> qualifiersToDelete = null;
    int qual_size =  qualifiers.size();

    boolean hasReplacedOrDeletedText = false;
    
    for(int i = 0; i  < qual_size; ++i) 
    {
      final Qualifier this_qualifier = (Qualifier)qualifiers.elementAt(i);
   
      if(qualifier_names != null &&
         !qualifier_names.contains(this_qualifier.getName())) 
        continue;

      final StringVector values = this_qualifier.getValues();

      if(values != null)
      {
        StringVector newValues = null;
        StringVector deleteValues = null;
        int val_size = values.size();
        
        for(int values_index = 0; values_index < val_size; 
             ++values_index) 
        {
          String this_value_string = (String)values.elementAt(values_index);

          if(this_value_string == null) 
            continue;

          if(fold_case) 
            this_value_string = this_value_string.toLowerCase();

          if(! match_substring &&
             this_value_string.equals(real_search_text) ||
             match_substring &&
             this_value_string.indexOf(real_search_text) != -1)
          {
            if(deleteQualifier)
            {
              if(deleteValues == null)
                deleteValues = new StringVector();
              deleteValues.add(this_value_string);
              continue;
            }
            // found the match & return true if replace function off
            if(replaceText == null)
              return true;            
            
            String new_text = replaceText;
            if(match_substring)
            {
              final String value_string = (String)values.elementAt(values_index);
              new_text =
                  Pattern.compile(real_search_text, Pattern.CASE_INSENSITIVE)
                  .matcher(value_string).replaceAll(replaceText);
            }
            if(newValues == null)
              newValues = this_qualifier.getValues();
            
            newValues.setElementAt(new_text, values_index);
          }
        }
        
        if(newValues != null)
        {
          if(newQualifiers == null)
            newQualifiers = new QualifierVector();
          newQualifiers.setQualifier(
              new Qualifier(this_qualifier.getName(),newValues));
          hasReplacedOrDeletedText = true;
        }
        
        if(deleteQualifier && deleteValues != null)
        {
          newValues = this_qualifier.getValues();
          newValues.removeAll(deleteValues);
          
          if(newValues.size() < 1)
          {
            if(qualifiersToDelete == null)
              qualifiersToDelete = new Vector<String>();
            qualifiersToDelete.add(this_qualifier.getName());
          }
          else
          {
            newQualifiers = new QualifierVector();
            newQualifiers.setQualifier(
              new Qualifier(this_qualifier.getName(),newValues));
          }
          hasReplacedOrDeletedText = true;
        }
      }
    }
    
    // delete qualifiers
    if(qualifiersToDelete != null)
    {
      try
      {
        for (int i = 0; i < qualifiersToDelete.size(); i++)
          removeQualifierByName(qualifiersToDelete.get(i));
      }
      catch (Exception e)
      {
        e.printStackTrace();
      }
    }
    
    if(newQualifiers != null)
    {
      for(int i=0; i<newQualifiers.size(); i++)
      {
        try
        {
          setQualifier((Qualifier) newQualifiers.elementAt(i));
        }
        catch(ReadOnlyException e)
        {
          e.printStackTrace();
        }
        catch(EntryInformationException e)
        {
          e.printStackTrace();
        }
      }
    }
    
    return hasReplacedOrDeletedText;
  }

  public boolean hasValidStartCodon() 
  {
    return hasValidStartCodon(false);
  }
  
  /**
   *  Return true if and only if this feature ends in a valid stop codon.
   **/
  public boolean hasValidStartCodon(final boolean force) 
  {
    if(!isCDS() && !force)
      return true;

    try 
    {
      if(getQualifierByName("codon_start") != null &&
         isPartialCDS() &&
         getFirstBase() == 1) 
        return true;
    }
    catch(InvalidRelationException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }

    final AminoAcidSequence translation = getTranslation();

    if(translation.length() < 1) 
      return false;

    final String first_codon = getTranslationBases().substring(0, 3);

    final StringVector start_codons = Options.getOptions().getStartCodons();

//  if(Options.getOptions().isEukaryoticMode()) 
//    start_codons = Options.getOptions().getEukaryoticStartCodons();
//  else
//    start_codons = Options.getOptions().getProkaryoticStartCodons();

    if(start_codons.contains(first_codon)) 
      return true;
    else 
      return false;
  }
  
  public boolean hasValidStopCodon() 
  {
    return hasValidStopCodon(false);
  }

  /**
   *  Return true if and only if this feature ends in a valid stop codon.
   **/
  public boolean hasValidStopCodon(final boolean force) 
  {
    if(!isCDS() && !force)
      return true;

    if(isPartialCDS() &&
       getLastBase() == getEntry().getBases().getLength()) 
      return true;

    final int bases_length = getBaseCount() - getCodonStart() + 1;

    if(bases_length % 3 != 0) 
      return false;

    if(bases_length < 3) 
      return false;

    final String codon_string = getBases().substring(getBaseCount() - 3);

    final char last_codon_translation =
      AminoAcidSequence.getCodonTranslation(codon_string);

    if(AminoAcidSequence.isStopCodon(last_codon_translation))
      return true;
    else
      return false;
  }

  /**
   *  Return true if and only if this feature has a valid EMBL key (rather
   *  than an Artemis extension).
   **/
  public boolean hasValidEMBLKey() 
  {
    if(getEntry().getEntryInformation().isValidKey(getKey()))
      return true;
    else
      return false;
  }

  /**
   *  Return true if and only if this feature has all the required EMBL
   *  qualifiers (for submission to EMBL).
   **/
  public boolean hasRequiredQualifiers()
  {
    final StringVector required_qualifiers =
      getEntry().getEntryInformation().getRequiredQualifiers(getKey());

    if(required_qualifiers == null)
      return true;

    try
    {
      int reqd_size = required_qualifiers.size();
      for(int i = 0; i < reqd_size; ++i) 
      {
        if(getQualifierByName((String)required_qualifiers.elementAt(i)) == null) 
          return false;
      }
    } 
    catch(InvalidRelationException e)
    {
      throw new Error("internal error - unexpected exception: " + e);
    }

    return true;
  }

  /**
   *  This method will fix the end location of the feature so that the feature
   *  ends with a stop codon.
   *  @return true if and only if the feature has a stop codon or if a stop
   *    codon was successfully added (by moving the end by three bases).
   **/
  public boolean fixStopCodon() 
      throws ReadOnlyException 
  {
    final Location old_location = getLocation();
    final int codon_start = getCodonStart();

    // there is no possible way to fix this feature
    if((getBaseCount() - codon_start + 1) % 3 != 0)
      return false;

    final FeatureSegment last_segment = getSegments().lastElement();
    final Marker last_base_marker = last_segment.getEnd();
    final Marker last_codon_marker;

    try
    {
      last_codon_marker = last_base_marker.moveBy(-2);
    } 
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }

    final Range last_codon_range;

    try 
    {
      last_codon_range = new Range(last_codon_marker.getPosition(),
                                   last_codon_marker.getPosition() + 2);
    } 
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }

    final String last_codon_string =
      last_codon_marker.getStrand().getSubSequence(last_codon_range);

    final char last_amino_acid_char =
      AminoAcidSequence.getCodonTranslation(last_codon_string);

    if(AminoAcidSequence.isStopCodon(last_amino_acid_char))
      return true;

    final Range next_codon_range;

    try 
    {
      next_codon_range =
        new Range(last_codon_marker.getPosition() + 3,
                  last_codon_marker.getPosition() + 5);
    } 
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }

    final String next_codon_string =
      last_codon_marker.getStrand().getSubSequence(next_codon_range);

    final char next_amino_acid_char =
      AminoAcidSequence.getCodonTranslation(next_codon_string);

    if(AminoAcidSequence.isStopCodon(next_amino_acid_char)) 
    {
      try 
      {
        // move the end marker to the next codon
        final Marker new_end_marker = last_base_marker.moveBy(3);

        last_segment.setEndPosition(new_end_marker.getPosition());
      } 
      catch(OutOfRangeException e) 
      {
        throw new Error("internal error - unexpected exception: " + e);
      }

      updateEMBLFeatureLocation();
      locationChanged(old_location);

      return true;
    }

    return false;
  }

  /**
   *  This method will move the start of the first segment of this feature so
   *  that the first codon is a start codon.  If the segment already starts at
   *  a start codon it will return true immediately.
   *  @param trim_to_any If true then the features will be trimmed to the next
   *    codon with the base pattern: ATG, GTG or TTG, otherwise the features
   *    will be trimmed to ATG (Met) only.
   *  @param trim_to_next If true then the features will be trimmed to the
   *    next start codon (dependent on trim_to_any) regardless of whether the
   *    feature currently start on a start codon.  If false then the feature
   *    will only be trimmed if the feature doesn't start on a start codon.
   *  @return true if and only if the trim worked.  It will fail if the first
   *    segment is less than 6 bases long or if there is no start codon in the
   *    first segment or in the first 30% of the feature.
   **/
  public boolean trimStart(final boolean trim_to_any,
                           final boolean trim_to_next)
      throws ReadOnlyException 
  {
    final Location old_location = getLocation();
    final FeatureSegment first_segment = getSegments().elementAt(0);

    // too short to trim
    if(first_segment.getBaseCount() < 6) 
      return false;

    final BasePattern search_pattern;

    try 
    {
      if(trim_to_any) 
        search_pattern = new BasePattern("dtg");
      else 
        search_pattern = new BasePattern("atg");
    } 
    catch(BasePatternFormatException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }

    final String first_codon_bases = getTranslationBases().substring(0, 3);

    Marker current_marker = first_segment.getStart();

    try 
    {
      current_marker = current_marker.moveBy(getCodonStart() - 1);
    } 
    catch(OutOfRangeException e) 
    {
      return false;
    }

    if(search_pattern.matches(first_codon_bases)) 
    {
      // the segment already starts on a start codon
      if(trim_to_next) 
      {
        // move passed the current start codon
        try 
        {
          current_marker = current_marker.moveBy(3);
        } 
        catch(OutOfRangeException e) 
        {
          return false;
        }
      } 
      else 
        return true;
    }

    final int end_of_segment_position =
      first_segment.getEnd().getPosition();

    final int start_of_feature_position =
      first_segment.getStart().getPosition();

    while(true) 
    {
      try 
      {
        final String current_codon = Strand.getCodonAtMarker(current_marker);

        if(search_pattern.matches(current_codon)) 
          break;

        current_marker = current_marker.moveBy(3);

        final int current_marker_position = current_marker.getPosition();

        if(current_marker_position > end_of_segment_position - 2 ||
           (!trim_to_next &&
            1.0 * (current_marker_position - start_of_feature_position) /
            getTranslationBasesLength() > 0.3)) 
        {
          // the current_marker is past the end of the first segment or is
          // more than 30% into the feature
          return false;
        }
      }
      catch(OutOfRangeException e) 
      {
        return false;
      }
    }

    try
    {
      first_segment.setStartPosition(current_marker.getPosition());
    } 
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }

    try 
    {
      // remove /codon_start (if it exists) since the new location will be in
      // the correct frame
      removeQualifierByName("codon_start");
    } 
    catch(OutOfDateException _) 
    {
      // do nothing
    } 
    catch(EntryInformationException _) 
    {
      // do nothing
    }

    updateEMBLFeatureLocation();
    locationChanged(old_location);

    return true;
  }

  /**
   *  Return the label (the value of the first /label qualifier) of this
   *  Feature, but return null if there is no label or if the label is "*".
   **/
  public String getLabel() 
  {
    try
    {
      return getValueOfQualifier("label");
    }
    catch(InvalidRelationException e)
    {
      return null;
    }
  }

  /**
   *  Return the gene name (the value of the first /gene qualifier) of this
   *  Feature or null if there it has no gene name.
   **/
  public String getGeneName() 
  {
    try 
    {
      return getValueOfQualifier("gene");
    }
    catch(InvalidRelationException e) 
    {
      return null;
    }
  }

  /**
   *  Return the name of this Feature to use in the display.  The search order
   *  is /primary_name, /synonym, /systematic_id, /temporary_systematic_id,
   *  /gene, /locus_tag, /label.  If none of the qualifiers exists return
   *  the feature key of this feature.
   **/
  public String getIDString() 
  {
    final String picked_name =
      pickName(Options.getOptions().getDisplayQualifierNames());

    if(picked_name == null) 
      return getKey().toString();
    else 
      return picked_name;
  }

  /**
   *  Return the systematic name of this feature.
   **/
  public String getSystematicName() 
  {
    final String picked_name =
      pickName(Options.getOptions().getSystematicQualifierNames());

    if(picked_name == null) 
      return getIDString();
    else 
      return picked_name;
  }


  /**
   *  Look at the qualifier_names one-by-one and return the first value of the
   *  first qualifier found.
   **/
  private String pickName(final StringVector qualifier_names) 
  {
    int qn_size = qualifier_names.size();
    for(int i = 0; i < qn_size; ++i)
    {
      try
      {
        final Qualifier qualifier =
          getQualifierByName((String)qualifier_names.elementAt(i));

        if(qualifier != null)
        {
          final StringVector values = qualifier.getValues();

          if(values != null && values.size() > 0)
          {
            for(int j=0; j<values.size(); j++)
            {
              final String value = (String)values.elementAt(j);
              if(value != null && !value.endsWith("current=false") && !value.equals(""))
                return value;
            }
          }
        }
      } 
      catch(InvalidRelationException e){}
    }

    return null;
  }

  /**
   *  Return the note (the value of the first /note qualifier) of this
   *  Feature or null if there is no /note.
   **/
  public String getNote()
  {
    try
    {
      return getValueOfQualifier("note");
    } 
    catch(InvalidRelationException e)
    {
      return null;
    }
  }

  /**
   *  Return the product (the value of the first /product qualifier) of this
   *  Feature or null if there is no /product.
   **/
  public String getProductString()
  {
    try 
    {
      final String product = getValueOfQualifier("product");
      if( product != null && 
          getEmblFeature() instanceof GFFStreamFeature &&
          product.startsWith("term=") )
        return product.substring(5);
      
      return product;
    }
    catch(InvalidRelationException e) 
    {
      return null;
    }
  }

  /**
   *  Return the number of bases in this feature not inclding introns (total of
   *  all segments).
   **/
  public int getBaseCount() 
  {
    if(base_count == -1) 
    {
      int new_base_count = 0;
      int seg_size = getSegments().size(); 

      for(int i = 0; i < seg_size; ++i) 
        new_base_count += getSegments().elementAt(i).getBaseCount();

      base_count = new_base_count;
    }

    return base_count;
  }

  /**
   *  Return the bases in this feature (all bases in all segments).
   **/
  public String getBases() 
  {
    if(bases == null) 
    {
      final StringBuffer buffer = new StringBuffer();
      int seg_size = getSegments().size();

      for(int i = 0; i < seg_size; ++i) 
        buffer.append(getSegments().elementAt(i).getBases());

      bases = buffer.toString();
    }

    return bases;
  }

  /**
   *  Return the number of coding bases in this feature, that is, the number
   *  of codons divided by three.
   **/
  public int getTranslationBasesLength()
  {
    final int codon_start = getCodonStart();
    final int start_index = codon_start - 1;

    int end_index = getBases().length();

    final int mod_value = (end_index - start_index) % 3;

    // add 1 or 2 if necessary to make the range a multiple of 3
    end_index -= mod_value;

    final int new_length = end_index - start_index;

    if(new_length >= 3) 
    {
      // remove the stop codon (if present)

      final String last_codon = getBases().substring(end_index - 3);

      final char amino_acid_char =
        AminoAcidSequence.getCodonTranslation(last_codon);

      if(AminoAcidSequence.isStopCodon(amino_acid_char)) 
        return new_length - 3;
    }

    return new_length;
  }

  /**
   *  Return the bases in this feature (all bases in all segments), taking
   *  the /start_codon qualifier into consideration.
   **/
  public String getTranslationBases() 
  {
    final int codon_start = getCodonStart();
    final int start_index = codon_start - 1;
    final int end_index = getTranslationBasesLength() + start_index;

    return getBases().substring(start_index, end_index);
  }

  /**
   *  Return an AminoAcidSequence containing the translated sequence of the
   *  feature with changes from /transl_except qualifiers added.  Return null
   *  if there are no valid /transl_except qualifiers (ie. parsable and with
   *  coordinates in range).
   **/
  private AminoAcidSequence fixTranslationExceptions() 
  {
    final String amino_acids_string = amino_acids.toString();
    String new_amino_acids_string = amino_acids_string;

    try
    {
      final Qualifier except_qualifier =
        getQualifierByName("transl_except");

      if(except_qualifier != null) 
      {
        final StringVector values = except_qualifier.getValues();

        if(values != null) 
        {
          for(int i = 0; i < values.size(); ++i) 
          {
            final String value = (String)values.elementAt(i);

            final String START_STRING = "(pos:";
            final String COMMA_STRING = ",aa:";

            if(value.startsWith(START_STRING) && value.endsWith(")")) 
            {
              final int comma_pos = value.lastIndexOf(COMMA_STRING);

              if(comma_pos >= 0) 
              {
                final String location_part =
                  value.substring(START_STRING.length(), comma_pos);

                final String aa_part =
                  value.substring(comma_pos + COMMA_STRING.length(),
                                  value.length() - 1);

                char aa_part_one_letter_code =
                  AminoAcidSequence.getOneLetterCode(aa_part);

                // aa_part is probably "OTHER"
                if(aa_part_one_letter_code == (char)-1) 
                  aa_part_one_letter_code = '.';

                final Location location = new Location(location_part);

                int start_base_in_feature;

                if(isForwardFeature())
                {
                  start_base_in_feature =
                    location.getFirstBase() - getRawFirstBase() -
                    (getCodonStart() - 1);
                }
                else 
                {
                  start_base_in_feature =
                    getRawLastBase() - location.getLastBase() -
                    (getCodonStart() - 1);
                }

                if(start_base_in_feature >= 0 &&
                   start_base_in_feature <
                   getBaseCount() - getCodonStart() + 1) 
                {
                  final int start_aa_in_feature = start_base_in_feature / 3;

                  new_amino_acids_string =
                    new_amino_acids_string.substring(0, start_aa_in_feature) +
                    aa_part_one_letter_code +
                    new_amino_acids_string.substring(start_aa_in_feature + 1);
                }
              }
            }
          }
        }
      }
    }
    catch(InvalidRelationException e) {}
    catch(LocationParseException e) {}

    if(new_amino_acids_string == amino_acids_string) 
      return null;
    else
      return new AminoAcidSequence(new_amino_acids_string);
  }

  /**
   *  Return the translation of the bases in this feature.
   **/
  public AminoAcidSequence getTranslation() 
  {
    if(amino_acids == null)
    {
      amino_acids =
        AminoAcidSequence.getTranslation(getTranslationBases(), true);

      // a very short feature
      if(amino_acids.length() == 0) 
        return amino_acids;

      final AminoAcidSequence fixed_amino_acids =
        fixTranslationExceptions();

      if(fixed_amino_acids != null) 
        amino_acids = fixed_amino_acids;

      if(isCDS() && !isPartialCDS() && hasValidStartCodon())
      {
        if(amino_acids.elementAt(0) != 'm') 
        {
          // translation should always start with M
          final String amino_acids_string = amino_acids.toString();

          final String new_amino_acids_string =
            'M' + amino_acids_string.substring(1);

          amino_acids = new AminoAcidSequence(new_amino_acids_string);
        }
      }
    }
    return amino_acids;
  }

  /**
   *  Return the number of amino acids in this feature (total of all segments).
   **/
  public int getAACount() 
  {
    if(aa_count == -1) 
      aa_count = getTranslationBasesLength() / 3;
    
    return aa_count;
  }

  /**
   *  Return the total molecular weight of the amino acid translation of the
   *  bases of this feature.
   **/
  public float getMolecularWeight() 
  {
    return getTranslation().getMolecularWeight();
  }

  /**
   *  Return the number of occurrences of the given codon in the frame
   *  starting at the correct frame for this feature (see the
   *  getTranslationBases() method).  The first index is the first base, the
   *  second index is the second base, etc.  The indices are the same as those
   *  for Bases.letter_index.  An example: calling getCodonCount(0, 0, 1)
   *  will return the count of the codon with the sequence t,t,c.
   **/
  public int getCodonCount(final int first,
                           final int second,
                           final int third) 
  {
    if(codon_counts == null) 
      setArrays();
    
    return codon_counts[first][second][third];
  }

  /**
   *  Return a count of the occurrences of each amino acids in the translation
   *  of the bases of the current feature.  The index should be the same as
   *  the index that is returned by AminoAcidSequence.getSymbolIndex().
   **/
  public int getResidueCount(final int amino_acid_index) 
  {
    if(residue_counts == null) 
      setArrays();
    
    return residue_counts[amino_acid_index];
  }

  /**
   *  Return the count of each base that appears in each codon position (1, 2
   *  or 3).  The first index is the codon position and the second is the base
   *  (indexed in the same way as Bases.letter_index).  (set by resetCache ()).
   **/
  public int getPositionalBaseCount(final int codon_base_position,
                                    final int base_index) 
  {
    if(positional_base_counts == null)
      setArrays();
   
    return positional_base_counts[codon_base_position][base_index];
  }

  /**
   *  Return the count each base in the feature.  The indices are the same as
   *  those for Bases.letter_index. (set by resetCache()).
   **/
  public int getBaseCount(final int base_index) 
  {
    if(base_counts == null)
      setArrays();
    
    return base_counts[base_index];
  }

  /**
   *  Return the codon position 1 and 2 correlation score for the bases of
   *  this feature.
   **/
  public double get12CorrelationScore() 
  {
    final int t1_count =
      getPositionalBaseCount(0, Bases.getIndexOfBase('t'));
    final int c1_count =
      getPositionalBaseCount(0, Bases.getIndexOfBase('c'));
    final int a1_count =
      getPositionalBaseCount(0, Bases.getIndexOfBase('a'));
    final int g1_count =
      getPositionalBaseCount(0, Bases.getIndexOfBase('g'));

    final int t2_count =
      getPositionalBaseCount(1, Bases.getIndexOfBase('t'));
    final int c2_count =
      getPositionalBaseCount(1, Bases.getIndexOfBase('c'));
    final int a2_count =
      getPositionalBaseCount(1, Bases.getIndexOfBase('a'));
    final int g2_count =
      getPositionalBaseCount(1, Bases.getIndexOfBase('g'));

    final int c3_count =
      getPositionalBaseCount(2, Bases.getIndexOfBase('c'));
    final int g3_count =
      getPositionalBaseCount(2, Bases.getIndexOfBase('g'));

    final int base_total = getTranslationBases().length();

    final double cor1_2_score =
      3.0 * t1_count/base_total *
      Codon12CorrelationAlgorithm.correlation_score_factors_1[0] +
      3.0 * c1_count/base_total *
      Codon12CorrelationAlgorithm.correlation_score_factors_1[1] +
      3.0 * a1_count/base_total *
      Codon12CorrelationAlgorithm.correlation_score_factors_1[2] +
      3.0 * g1_count/base_total *
      Codon12CorrelationAlgorithm.correlation_score_factors_1[3] +
      3.0 * t2_count/base_total *
      Codon12CorrelationAlgorithm.correlation_score_factors_2[0] +
      3.0 * c2_count/base_total *
      Codon12CorrelationAlgorithm.correlation_score_factors_2[1] +
      3.0 * a2_count/base_total *
      Codon12CorrelationAlgorithm.correlation_score_factors_2[2] +
      3.0 * g2_count/base_total *
      Codon12CorrelationAlgorithm.correlation_score_factors_2[3] +
      0.5;                 // add 0.5 because that is what the old uk.ac.sanger.artemis did

    return cor1_2_score;
  }

  /**
   *  Return the percentage GC content of the bases of this feature.
   **/
  public double getPercentGC() 
  {
    final String bases = getBases();

    if(bases.length() > 0)
    {
      int gc_count = 0;

      final char[] sequence_chars = new char[bases.length()];

      bases.getChars(0, bases.length(), sequence_chars, 0);

      for(int i = 0 ; i < bases.length() ; ++i) 
      {
        final char this_char = sequence_chars[i];
        if(this_char == 'g' || this_char == 'c') 
          ++gc_count;
      }

      return 100.0 * gc_count / bases.length();
    } 
    else
      return 0.0;
  }

  /**
   *  Return the first (lowest) base number of this Feature.  The lowest base
   *  of the Feature is the minimum (with respect to the Feature's Strand) of
   *  the lowest bases of all the segments of this Feature.
   *  @return The first base of this feature or 0 if this feature doesn't know
   *    where it is.
   **/
  public int getFirstBase()
  {
    final Marker first_base_marker = getFirstBaseMarker();

    if(first_base_marker == null) 
      return 0;
    else 
      return first_base_marker.getPosition ();
  }

  /**
   *  Return the last (highest) base number of this Feature.  The highest base
   *  of the Feature is the maximum (with respect to the Feature's Strand) of
   *  the highest bases of all the segments of this Feature.
   *  @return The last base of this feature or 0 if this feature doesn't know
   *    where it is.
   **/
  public int getLastBase()
  {
    final Marker last_base_marker = getLastBaseMarker();

    if(last_base_marker == null)
      return 0;
    else
      return last_base_marker.getPosition ();
  }

  /**
   *  Return true if and only if the first base of this Feature is less than
   *  the argument Feature with respect to the Feature's Strand.
   **/
  public boolean lessThan(Feature other_feature) 
  {
    if(getFirstBase() < other_feature.getFirstBase())
      return true;
    else
      return false;
  }

  /**
   *  Return true if and only if the first base of this Feature is greater
   *  than the argument Feature with respect to the Feature's Strand.
   **/
  public boolean greaterThan(Feature other_feature) 
  {
    if(getFirstBase() > other_feature.getFirstBase())
    {
      return true;
    } 
    else 
    {
      return false;
    }
  }

  /**
   *  Return the first (lowest) base number of this Feature.  The lowest base
   *  of the Feature is the minimum (with respect to the raw sequence) of the
   *  lowest bases of all the segments of this Feature.
   *  @return The first base of this feature or 0 if this feature doesn't know
   *    where it is.
   **/
  public int getRawFirstBase()
  {
    final int A_BIG_NUMBER = 0x7fffffff; // largest integer
    int minimum = A_BIG_NUMBER;

    for(int i = 0; i < getSegments().size(); ++i) 
    {
      final int current_minimum;
      if(isForwardFeature()) 
      {
        current_minimum =
          getSegments().elementAt(i).getStart().getRawPosition();
      } 
      else 
      {
        current_minimum =
          getSegments().elementAt(i).getEnd().getRawPosition();
      }

      if(current_minimum < minimum)
        minimum = current_minimum;
    }

    if(minimum == A_BIG_NUMBER) 
      return 0;
    else
      return minimum;
  }

  /**
   *  Return the last (highest) base number of this Feature.  The highest base
   *  of the Feature is the maximum (with respect to the raw sequence) of the
   *  highest bases of all the segments of this Feature.
   *  @return The last base of this feature or 0 if this feature doesn't know
   *    where it is.
   **/
  public int getRawLastBase()
  {
    final int A_SMALL_NUMBER = -1;
    int maximum = A_SMALL_NUMBER;
    int seg_size = getSegments().size();

    for(int i = 0; i < seg_size; ++i)
    {
      final int current_maximum;
      if(isForwardFeature())
      {
        current_maximum =
          getSegments().elementAt(i).getEnd().getRawPosition();
      } 
      else
      {
        current_maximum =
          getSegments().elementAt(i).getStart().getRawPosition();
      }

      if(current_maximum > maximum) 
        maximum = current_maximum;
    }

    if(maximum == A_SMALL_NUMBER) 
      return 0;
    else
      return maximum;
  }

  /**
   *  Return the maximum extent of this feature.  The range returned starts at
   *  getRawFirstBase() and ends at getRawLastBase().
   **/
  public Range getMaxRawRange()
  {
    try 
    {
      return new Range(getRawFirstBase(), getRawLastBase());
    } 
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }
  }


  /**
   *  Return true if and only if the first base of this Feature is less than
   *  the argument Feature with to the raw sequence.
   **/
  public boolean rawLessThan(Feature other_feature) 
  {
    if(getRawFirstBase() < other_feature.getRawFirstBase()) 
      return true;
    else
      return false;
  }

  /**
   *  Return true if and only if the first base of this Feature is greater
   *  than the argument Feature with to the raw sequence.
   **/
  public boolean rawGreaterThan(Feature other_feature) 
  {
    if(getRawFirstBase() > other_feature.getRawFirstBase())
      return true;
    else
      return false;
  }

  /**
   *  Return the Marker of the first (lowest) base number of this Feature.
   *  The lowest base of the Feature is the minimum of the lowest bases of all
   *  the segments of this Feature.
   *  @return The Marker of the first base of this feature or null if this
   *    feature doesn't know where it is.
   **/
  public Marker getFirstBaseMarker() 
  {
    final int A_BIG_NUMBER = 0x7fffffff; // largest integer
    int minimum = A_BIG_NUMBER;
    Marker minimum_marker = null;

    for(int i = 0 ; i < getSegments().size() ; ++i)
    {
      final Marker current_marker = getSegments().elementAt(i).getStart();
      final int current_minimum = current_marker.getPosition();

      if(current_minimum < minimum)  
      {
        minimum = current_minimum;
        minimum_marker = current_marker;
      }
    }

    return minimum_marker;
  }

  /**
   *  Return the Marker of the last (highest) base number of this Feature.
   *  The highest base of the Feature is the maximum of the highest bases of
   *  all the segments of this Feature.
   *  @return The Marker of the last base of this feature or null if this
   *    feature doesn't know where it is.
   **/
  public Marker getLastBaseMarker() 
  {
    final long A_SMALL_NUMBER = -1;
    long maximum = A_SMALL_NUMBER;
    Marker maximum_marker = null;

    for(int i = 0; i < getSegments().size(); ++i) 
    {
      final Marker current_marker = getSegments().elementAt(i).getEnd();
      final int current_maximum = current_marker.getPosition();

      if(current_maximum > maximum) 
      {
        maximum = current_maximum;
        maximum_marker = current_marker;
      }
    }

    return maximum_marker;
  }

  /**
   *  Return the Marker of the first base of the first (coding) codon of this
   *  feature.  This is the minimum of the lowest bases of all the segments of
   *  this Feature plus codon_start minus one.
   *  @return The Marker of the first coding base of this feature or null if
   *    this feature doesn't know where it is.
   **/
  public Marker getFirstCodingBaseMarker() 
  {
    final Marker first_base_marker = getFirstBaseMarker();

    try 
    {
      return first_base_marker.moveBy(getCodonStart() - 1);
    } 
    catch(OutOfRangeException e) 
    {
      // if we are at the end of the sequence then just return the first base
      // marker
      return first_base_marker;
    }
  }

  /**
   *  Given a base position in this feature, return the Marker of the
   *  corresponding base on the Strand object.
   *  @param position The position within the feature - not counting introns.
   *  @exception OutOfRangeException Thrown if the new position is less than 1
   *    or greater than the length of the sequence.
   **/
  public Marker getPositionInSequence (final int position)
      throws OutOfRangeException 
  {
    // the position is outside of the feature
    if(position < 1) 
      throw new OutOfRangeException("position: " + position);

    // subtract one because positions are numbered from 1 not 0
    int bases_remaining = position - 1;

    for(int i = 0; i < segments.size(); ++i) 
    {
      final FeatureSegment this_segment = segments.elementAt(i);

      if(bases_remaining < this_segment.getBaseCount()) 
        return this_segment.getStart().moveBy(bases_remaining);
      else
        bases_remaining -= this_segment.getBaseCount();
    }

    // the position is outside of the feature
    throw new OutOfRangeException("position: " + position);
  }

  /**
   *  Return the (base) position in this feature of the given Marker or -1 if
   *  the given Marker does not correspond to any bases in this Feature.
   **/
  public int getFeaturePositionFromMarker(final Marker base_marker) 
  {
    // counts the bases in the previous segments as we loop over them
    int bases_so_far = 0;
    final int marker_position = base_marker.getPosition();
    int seg_size = segments.size();

    for(int i = 0; i < seg_size; ++i) 
    {
      final FeatureSegment this_segment = segments.elementAt(i);
      final int this_segment_start = this_segment.getStart().getPosition();

      final int marker_start_diff = marker_position - this_segment_start;

      if(marker_start_diff >= 0 &&
         marker_start_diff < this_segment.getBaseCount())
      {
        return bases_so_far + marker_start_diff;
      }

      bases_so_far += this_segment.getBaseCount();
    }

    return -1;
  }

  public void resetColour()
  {
    colour = null;
  }

  /**
   *  Return the colour (the value of the /colour qualifier) of this Feature
   *  or null if there is no colour qualifier and no default colour.
   **/
  public Color getColour()
  {
    if(colour != null)
      return colour;

    String colour_qualifier;
    try 
    {
      colour_qualifier = getValueOfQualifier("colour");

      // it's international "be nice to Americans day":
      if(colour_qualifier == null) 
        colour_qualifier = getValueOfQualifier("color");
    }
    catch(InvalidRelationException e) 
    {
      colour_qualifier = null;
    }

    // use default colour for this type of feature
    if(colour_qualifier == null) 
    {
      colour = Options.getOptions().getDefaultFeatureColour(getKey());
      return colour;
    }

    final StringVector colours = StringVector.getStrings(colour_qualifier);

    if(colours.size() < 1) 
    {
      colour = Options.getOptions().getDefaultFeatureColour(getKey());
      return colour;
    }

    try 
    {
      if(colours.size() == 3)
      {
        int red   = Integer.parseInt((String)colours.elementAt(0));
        int green = Integer.parseInt((String)colours.elementAt(1));
        int blue  = Integer.parseInt((String)colours.elementAt(2));

        if(red < 0) 
          red = 0;

        if(red > 255) 
          red = 255;

        if(green < 0)
          green = 0;

        if(green > 255) 
          green = 255;

        if(blue < 0)
          blue = 0;

        if(blue > 255)
          blue = 255;

        colour = new Color(red, green, blue);
        return colour;
      } 
      else
      {
        final String colour_string = (String)colours.elementAt(0);

        final int colour_number;

        colour_number = Integer.parseInt(colour_string);
        colour = Options.getOptions().getColorFromColourNumber(colour_number);
        return colour;
      }
    } 
    catch(NumberFormatException e) 
    {
      // use default colour for this type of feature
      colour = Options.getOptions().getDefaultFeatureColour(getKey());
      return colour;
    }
  }

  /**
   *  Return a vector of String objects containing all the values of the
   *  qualifier with the given name.  If the Feature has two qualifiers:
   *    /db_xref="PID:e1256466" /db_xref="SPTREMBL:O43016"
   *  then this method will return a vector containing "PID:e1256466" and
   *  "SPTREMBL:O43016".
   *  Returns null if the is no Qualifier with the given name.
   **/
  public StringVector getValuesOfQualifier(String qualifier_name)
      throws InvalidRelationException 
  {
    final Qualifier this_qualifier = getQualifierByName(qualifier_name);

    if(this_qualifier == null) 
      return null;
    else
      return this_qualifier.getValues();
  }

  /**
   *  Return the Qualifier in this Feature with the given name or null if
   *  there no such Qualifier.
   *  @param name The qualifier name to look for.
   **/
  public Qualifier getQualifierByName(final String name)
      throws InvalidRelationException 
  {
    return getEmblFeature().getQualifierByName(name);
  }

  /**
   *  Return the value of the first occurrence of the given qualifier in this
   *  Feature or null if this feature doesn't have this qualifier.  If this
   *  feature has the given qualifier but the qualifier doesn't have a value
   *  then it will return "" - a zero length String.
   *  @param qualifier_name The name of the qualifier to search for.
   **/
  public String getValueOfQualifier(String qualifier_name)
      throws InvalidRelationException 
  {
    final StringVector values = getValuesOfQualifier(qualifier_name);

    if(values == null) 
      return null;
    else
    {
      if(values.size() == 0) 
        return "";
      else
      {
        final String first_element = (String)values.elementAt(0);

        if(first_element == null)
          return "";
        else
          return first_element;
      }
    }
  }

  /**
   *  Remove this Feature from the Entry object that contains it.  It will
   *  also be removed from the underlying embl.Entry object.
   **/
  public void removeFromEntry()
      throws ReadOnlyException 
  {
    getEntry().remove(this);
  }

  /**
   *  Return the vector containing the references of the FeatureSegment
   *  objects of this Feature.
   **/
  public FeatureSegmentVector getSegments() 
  {
    final Location current_location = getEmblFeature().getLocation();

    if(segments == null) 
      createSegments();
    else 
    {
      // see comment on old_location
      if(current_location != old_location) 
      {
        reexamineSegments();
        if(!old_location.equals(current_location)) 
        {
          locationChanged(old_location);
          old_location = current_location;
        }
      }
    }

    return segments;
  }

  /**
   *  Examine each Range in the Location of this Feature.  Remove
   *  FeatureSegments that don't correspond to Ranges in the Location.  Add
   *  a FeatureSegment for each extra Range that has appeared.  This is
   *  necessary because the feature location may change at any time in a
   *  client/server environment.
   **/
  private void reexamineSegments()
  {
    stopSegmentsListening();
    
    final Location current_location = getEmblFeature().getLocation();
    final RangeVector ranges = current_location.getRanges();
    final Vector new_segments = new Vector();

    int ranges_size = ranges.size();
    new_segments.setSize(ranges_size);
  
    FeatureSegmentVector old_segments =
      (FeatureSegmentVector)segments.clone();

    // if the feature has changed strand, give up and recreate all
    // segments
    if(current_location.isComplement() != old_location.isComplement())
      old_segments.removeAllElements();

    // first find exact matches
EXACT_SEGMENTS:
    for(int old_segment_index = old_segments.size() - 1; old_segment_index >= 0;
        --old_segment_index) 
    {
      final FeatureSegment old_segment =
        old_segments.elementAt(old_segment_index);

      for(int range_index = 0; range_index < ranges_size;
          ++range_index)
      {
        final Range new_range = (Range)ranges.elementAt(range_index);

        // there is a Range in the new Location that exactly matches a Range
        // in an old FeatureSegment so reuse the old FeatureSegment
        if(old_segment.getRawRange().equals(new_range)) 
        {
          old_segment.setRange(new_range);
          new_segments.setElementAt(old_segment, range_index);
          old_segments.removeElementAt(old_segment_index);
          continue EXACT_SEGMENTS;
        }
      }
    }

    // find matches where the start or end has changed
CHANGED_END:
    for(int old_segment_index = old_segments.size() - 1; old_segment_index >= 0;
        --old_segment_index) 
    {
      final FeatureSegment old_segment =
        old_segments.elementAt(old_segment_index);

      for(int range_index = 0; range_index < ranges_size;
          ++range_index) 
      {
        final Range new_range = (Range)ranges.elementAt(range_index);

        if(old_segment.getRawRange().getStart() ==
           new_range.getStart() ||
           old_segment.getRawRange().getEnd() ==
           new_range.getEnd()) 
        {
          // there is a Range in the new Location that partly matches a Range
          // in an old FeatureSegment so reuse the old FeatureSegment
          old_segment.setRange(new_range);
          new_segments.setElementAt(old_segment, range_index);
          old_segments.removeElementAt(old_segment_index);
          continue CHANGED_END;
        }
      }
    }

    // create a segment for each range that we don't have a segment for
    for(int new_segment_index = 0; new_segment_index < ranges_size;
        ++new_segment_index) 
    {
      final Range missing_range = (Range)ranges.elementAt(new_segment_index);
      
      if(new_segments.elementAt(new_segment_index) == null) 
      {
        // new Range -> create a segment
        new_segments.setElementAt(makeSegment(missing_range),
                                   new_segment_index);
      }
    }

    segments = new FeatureSegmentVector();

    for(int i = 0; i < ranges.size(); ++i)
      segments.addElementAtEnd((FeatureSegment) new_segments.elementAt(i));
    
    startSegmentsListening();
  }

  /**
   *  Return the reference of the Strand that this Feature and it's
   *  FeatureSegment objects is associated with.
   **/
  public Strand getStrand() 
  {
    // get the strand of the first entry in the entry group
    if(isForwardFeature())
      return entry.getBases().getForwardStrand();
    else
      return entry.getBases().getReverseStrand();
  }
  
  /**
   *  Return the reference of a new copy of this Feature.  This method will
   *  update the underlying embl.Entry and the new Feature will be installed
   *  in the same Entry object as this one.
   *  @return The reference of the new feature.
   **/
  public Feature duplicate()
                 throws ReadOnlyException
  {
    return duplicate(false);
  }

  /**
   * Return the reference of a new copy of this Feature.  This method will
   *  update the underlying embl.Entry and the new Feature will be installed
   *  in the same Entry object as this one.
   * @param isDuplicatedInChado  if true then create new ID and update in chado
   * @return
   * @throws ReadOnlyException
   */
  public Feature duplicate(final boolean isDuplicatedInChado)
      throws ReadOnlyException
  {
    uk.ac.sanger.artemis.io.Feature new_embl_feature;

    if(getEmblFeature() instanceof GFFStreamFeature)
      new_embl_feature = new GFFStreamFeature(getEmblFeature(), isDuplicatedInChado);
    else
      new_embl_feature = new EmblStreamFeature(getEmblFeature());

    final Feature return_feature = new Feature(new_embl_feature);

    try 
    {
      getEntry().add(return_feature, !isDuplicatedInChado, false);
    } 
    catch(EntryInformationException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    } 
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }

    return return_feature;
  }

  /**
   *  Move this Feature to the given Entry.  The move is achieved by removing
   *  the feature from it's current Entry and adding to the the
   *  destination_entry.
   *  @param destination_entry The Feature will be moved to this Entry.
   *  @param force If true then invalid qualifiers will be quietly thrown away
   *    when saving.  Features with invalid keys won't be moved.  "Invalid"
   *    means that the key/qualifier is not allowed to occur in the
   *    destination type (probably determined by the default EntryInformation
   *    object of the destination type).  If false an
   *    EntryInformationException will be thrown for invalid keys or
   *    qualifiers.
   *  @exception EntryInformationException Thrown if force is false and if the
   *    destination type cannot contain the Key, Qualifier or Key/Qualifier
   *    combination of the given feature.
   **/
  public void moveTo(final Entry destination_entry, final boolean force)
      throws EntryInformationException, OutOfRangeException,
      ReadOnlyException 
  {
    if(destination_entry.isReadOnly()) 
      throw new ReadOnlyException();
    else 
    {
      // this (hack) causes the calls to stopListening() and startListening()
      // in setEntry() to return immediately which improves the speed a lot
      // when many features are loaded
      startListening();

      try 
      {
        final Entry old_entry = getEntry();

        getEntry().remove(this);
        
        try 
        {
          destination_entry.add(this, force);
        } 
        catch(EntryInformationException e) 
        {
          // put the feature back where it was
          old_entry.add(this, true);
          
          // re-throw
          throw e;
        }
      } 
      finally 
      {
        // see note above on startListening()
        stopListening();
      }
    }
  }

  /**
   *  Make a copy of the this Feature and add it to the given Entry.
   *  @return The reference of the new feature.
   *  @exception EntryInformationException Thrown if the destination Entry
   *    cannot contain a Feature with this Key, Qualifier or Key/Qualifier
   *    combination.
   **/
  public Feature copyTo(final Entry destination_entry)
      throws EntryInformationException, OutOfRangeException,
      ReadOnlyException 
  {
    if(destination_entry.isReadOnly())
      throw new ReadOnlyException();

    uk.ac.sanger.artemis.io.Feature new_embl_feature =
                         new EmblStreamFeature(getEmblFeature());

    final Feature return_feature = new Feature(new_embl_feature);

    destination_entry.add(return_feature, false);

    return return_feature;
  }

  /**
   *  Adds the specified event listener to receive feature change events from
   *  this object.
   *  @param l the change event listener.
   **/
  public void addFeatureChangeListener(FeatureChangeListener l) 
  {
    feature_listener_list.addElement(l);
  }

  /**
   *  Removes the specified event listener so that it no longer receives
   *  feature change events from this object.
   *  @param l the change event listener.
   **/
  public void removeFeatureChangeListener(FeatureChangeListener l) 
  {
    feature_listener_list.removeElement(l);
  }

  /**
   *  Returns true if and only if this feature is read only or is in a read
   *  only entry.
   **/
  public boolean isReadOnly()
  {
    return getEmblFeature().isReadOnly();
  }

  /**
   *  This method examines each FeatureSegment and uses that information to
   *  set the location of the underlying embl.Feature object.  This must be
   *  called any time a segment position changes, or when there is a change to
   *  the sequence (the sequence as represented by the Strand and Bases
   *  objects).  This method does not fire off any FeatureChange events.
   **/
  private void updateEMBLFeatureLocation()
      throws ReadOnlyException 
  {
    final boolean complement = getLocation().isComplement();
    final RangeVector ranges = new RangeVector();

    for(int i = 0; i < segments.size(); ++i) 
      ranges.addElement(segments.elementAt (i).getRawRange ());

    try 
    {
      final Location new_location = new Location(ranges, complement);

      getEmblFeature().setLocation(new_location);
      old_location = new_location;
    }
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - inconsistent location information: " +
                       e.getMessage());
    }
  }


  /**
   *  Send an event to those object listening for it.
   *  @param listeners A Vector of the objects that the event should be sent
   *    to.
   *  @param event The event to send
   **/
  private void fireAction(Vector listeners, ChangeEvent event) 
  {
    final Vector targets;
    // copied from a book - synchronising the whole method might cause a
    // deadlock
    synchronized(this)
    {
      targets = (Vector)listeners.clone();
    }

    for(int i = 0 ; i < targets.size() ; ++i) 
    {
      ChangeListener target = (ChangeListener)targets.elementAt(i);

      final FeatureChangeListener feature_change_listener =
                                           (FeatureChangeListener)target;
      feature_change_listener.featureChanged((FeatureChangeEvent) event);
    }
  }

  /**
   *  Reset these arrays: segments, amino_acids, bases, codon_counts,
   *  residue_counts, positional_base_counts and base_counts.
   **/
  private void resetCache()
  {
    amino_acids = null;
    bases = null;
    codon_counts = null;
    residue_counts = null;
    positional_base_counts = null;
    base_counts = null;
    aa_count = -1;
    base_count = -1;
  }

  /**
   *  Update the values stored in residue_counts, codon_counts, base_counts
   *  and positional_base_counts.
   **/
  private void setArrays() 
  {
    final String translation_bases = getTranslationBases();
    final AminoAcidSequence translation = getTranslation();

    final String translation_string = translation.toString();

    codon_counts = new int[4][4][4];
    residue_counts = new int[AminoAcidSequence.symbol_count];
    base_counts = new int[4];
    positional_base_counts = new int[3][4];

    // zero the residue_counts array
    for(int i = 0; i < residue_counts.length; ++i) 
      residue_counts[i] = 0;

    int trans_len = translation_string.length();
    for(int i = 0; i < trans_len; ++i) 
    {
      final int symbol_index =
        AminoAcidSequence.getSymbolIndex(translation_string.charAt(i));
      ++residue_counts[symbol_index];
    }

    // zero the codon_counts array
    for(int first = 0 ; first < 4 ; ++first) 
    {
      for(int second = 0 ; second < 4 ; ++second) 
      {
        for(int third = 0 ; third < 4 ; ++third)
          codon_counts[first][second][third] = 0;
      }
    }

    int base_len = translation_bases.length();
    for(int base_index = 0; base_index < base_len;
        ++base_index)
    {
      // this is 0, 1, 2 or 3
      final int index_of_base =
        Bases.getIndexOfBase(translation_bases.charAt(base_index));

      if(index_of_base < 4) 
        ++base_counts[index_of_base];
    }

    for(int i = 0; i < translation_bases.length() / 3; ++i) 
    {
      final int first_base_index =
        Bases.getIndexOfBase(translation_bases.charAt(i * 3));
      final int second_base_index =
        Bases.getIndexOfBase(translation_bases.charAt(i * 3 + 1));
      final int third_base_index =
        Bases.getIndexOfBase(translation_bases.charAt(i * 3 + 2));

      if(first_base_index < 4) 
        ++positional_base_counts[0][first_base_index];
      
      if(second_base_index < 4) 
        ++positional_base_counts[1][second_base_index];
      
      if(third_base_index < 4)
        ++positional_base_counts[2][third_base_index];
      
      if(first_base_index < 4 && second_base_index < 4 && third_base_index < 4) 
        ++codon_counts[first_base_index][second_base_index][third_base_index];
    }
  }

  /**
   *  Remove all FeatureSegments from all listener lists.
   **/
  private void stopSegmentsListening() 
  {
    if(segments != null) 
    {
      for(int i = 0 ; i < segments.size(); ++i) 
      {
        segments.elementAt(i).stopListening();
        segments.elementAt(i).removeMarkerChangeListener(this);
      }
    }
  }

  /**
   *  Remove this Feature and its FeatureSegments from all listener lists.
   **/
  private void stopListening()
  {
    if(listen_count == 1) 
    {
      if(getEntry() != null && getEntry().getBases() != null)
      {
        final Bases bases = getEntry().getBases();
        bases.removeSequenceChangeListener(this);
      }

      stopSegmentsListening();

      Options.getOptions().removeOptionChangeListener(this);
    }

    --listen_count;

    if(listen_count < 0) 
      throw new Error("Feature.listen_count < 0");
  }

  /**
   *  Add all FeatureSegments to all appropriate listener lists.
   **/
  private void startSegmentsListening() 
  {
    if(segments != null)
    {
      int seg_size = segments.size();
      for (int i = 0; i < seg_size; ++i) 
      {
        segments.elementAt(i).startListening();
        segments.elementAt(i).addMarkerChangeListener(this);
      }
    }
  }

  /**
   *  Add this Feature and its FeatureSegments to the all appropriate listener
   *  lists.
   **/
  private void startListening() 
  {
    if(listen_count == 0) 
    {
      if(getEntry() != null && getEntry().getBases() != null) 
      {
        final Bases bases = getEntry().getBases();
        startSegmentsListening();

        // it doesn't matter what the priority is as long as it is lower than
        // the FeatureSegment priority, so that all FeatureSegment objects get
        // updated before we try to update the location of this feature
        final int PRIORITY = Marker.LISTENER_PRIORITY - 2;
        bases.addSequenceChangeListener(this, PRIORITY);
      }

      Options.getOptions().addOptionChangeListener(this);
    }

    ++listen_count;
  }

  /**
   *  Create one FeatureSegment object for each segment/exon of this Feature.
   **/
  private void createSegments() 
  {
    // this call will remove each segment from the MarkerChangeListener list
    // of the start and end Marker objects of the segment.  If we don't do
    // this the segments will not be garbage collected.
    stopSegmentsListening();

    segments = new FeatureSegmentVector();
    final Location location = getLocation();
    final RangeVector ranges = location.getRanges();

    int ranges_size = ranges.size();
    for(int i = 0; i < ranges_size; ++i) 
    {
      final Range this_range = (Range)ranges.elementAt(i);
      segments.add(makeSegment(this_range));
    }

    startSegmentsListening();
  }

  /**
   *  Make a single FeatureSegment.
   *  @param range Provides the start and end positions of the segment.
   **/
  private FeatureSegment makeSegment(final Range range) 
  {
    return new FeatureSegment(this, range);
  }

  /**
   *  Set the Location of this feature and reset old_location and
   *  segments.  Don't send any FeatureChangeEvents.
   **/
  private void setLocationInternal(final Location new_location)
      throws ReadOnlyException, OutOfRangeException 
  {
    getEmblFeature().setLocation(new_location, getEntry().getEMBLEntry());

    if(new_location != old_location) 
      reexamineSegments();

    old_location = new_location;
    resetCache();
  }

  /**
   *  Set the Location of this feature and notify any FeatureChangeListeners.
   **/
  public void setLocation(final Location new_location)
      throws ReadOnlyException, OutOfRangeException 
  {
    final Location old_location = getLocation();
    setLocationInternal(new_location);
    locationChanged(old_location);
  }
  
  public void addSegment(final Range range)
         throws ReadOnlyException 
  {
    final QualifierVector old_qualifiers = getQualifiers().copy();
    addSegment(range, old_qualifiers);
  }

  /**
   *  Add a FeatureSegment to this Feature.
   *  @param range Provides the start and end positions of the segment.
   **/
  public void addSegment(final Range range,
                         final QualifierVector old_qualifiers)
         throws ReadOnlyException 
  {
    final Location old_location = getLocation();
    final Location new_location = old_location.addRange(range);
    
    try
    {
      setLocationInternal(new_location);
    } 
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - inconsistent location " +
                      "information: " + e);
    }

    reexamineSegments();
    locationChanged(old_location, old_qualifiers, FeatureChangeEvent.SEGMENT_CHANGED);
  }

  private Location moveSegments(final int shift)
      throws OutOfRangeException, ReadOnlyException
  {
    Location location = getLocation();
    RangeVector ranges = old_location.getRanges();
    for(int i=0; i<ranges.size(); i++)
    {
      old_location = location;
      final Range seg_range = (Range)ranges.elementAt(i);
      final Range seg_new_range = new Range(seg_range.getStart() + shift,
                                            seg_range.getEnd() + shift);
      location = location.changeRange(seg_range, seg_new_range);
   
//    System.out.println("DIFF "+ shift+" "+seg_new_range.getStart()+".."+
//                                         seg_new_range.getEnd()+"  was "+
//                                          seg_range.getStart()+".."+
//                                        seg_range.getEnd());
    }

    return location;
  }

  /**
   *  Delete a FeatureSegment from this object and update the location
   *  underlying embl.Feature object.
   *  @param segment The FeatureSegment to delete.
   *  @exception LastSegmentException thrown if an attempt is made to remove
   *    the last segment.
   **/
  public void removeSegment(FeatureSegment segment)
      throws ReadOnlyException, LastSegmentException 
  {
    if(getSegments().size() > 1) 
    {
      QualifierVector qualifiers = getQualifiers().copy();
      final Location old_location = getLocation();

      final Location new_location =
        old_location.removeRange(segment.getRawRange());

      try 
      {
        setLocationInternal(new_location);
      } 
      catch(OutOfRangeException e) 
      {
        throw new Error("internal error - inconsistent location " +
                        "information: " + e);
      }

      reexamineSegments();

      locationChanged(old_location, qualifiers, FeatureChangeEvent.SEGMENT_CHANGED);
    } 
    else 
      throw new LastSegmentException();
  }

  /**
   *  Return a list of all the names of the qualifiers in the given features.
   **/
  public static StringVector getAllQualifierNames(final FeatureVector features)
  {
    final StringVector qualifier_names = new StringVector();
    int feat_size = features.size();

    for(int i = 0; i < feat_size; ++i)
    {
      final Feature feature = features.elementAt(i);
      final QualifierVector qualifiers = feature.getQualifiers();

      int qualifiers_size = qualifiers.size();
      for(int qualifier_index = 0; qualifier_index < qualifiers_size;
          ++qualifier_index)
      {
        final Qualifier this_qualifier = (Qualifier)qualifiers.elementAt(qualifier_index);
        final String name = this_qualifier.getName();

        if(!qualifier_names.contains(name)) 
          qualifier_names.add(name);
      }
    }

    qualifier_names.sort();

    return qualifier_names;
  }

  /**
   *  Reverse and complement the location of this feature.
   *  @param sequence_length The length of the sequence that this Feature is
   *    associated with.
   **/
  private void reverseComplement(final int sequence_length)
      throws ReadOnlyException 
  {
    try 
    {
      final Location new_location =
        getLocation().reverseComplement(sequence_length);

      setLocationInternal(new_location);
    } 
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - inconsistent location: " + e);
    }
  }

  /**
   *  Return the sequence length of the Strand that this feature is attached
   *  to.
   **/
  private int getSequenceLength()
  {
    FeatureSegment first_seg = getSegments().elementAt(0);
    return first_seg.getStart().getStrand().getSequenceLength();
  }

  /**
   *  Set the embl.Feature reference of this Feature to be the given reference.
   **/
  void setEmblFeature(final uk.ac.sanger.artemis.io.Feature new_embl_feature)
  {
    embl_feature = new_embl_feature;
  }

}
