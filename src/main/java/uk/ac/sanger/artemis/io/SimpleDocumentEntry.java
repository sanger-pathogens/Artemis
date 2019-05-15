/* SimpleDocumentEntry.java
 *
 * created: Tue Feb 15 2000
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2019  Genome Research Limited
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

import uk.ac.sanger.artemis.util.*;

import java.io.*;
import java.util.Date;
import java.util.HashSet;
import java.util.Set;
import java.util.Vector;
import java.util.Hashtable;
import java.util.List;

/**
 *  This class contains the methods common to all DocumentEntry objects.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 **/

abstract public class SimpleDocumentEntry
                      implements DocumentEntry 
{

  /**
   *  A vector of ReadOnlyEmblStreamFeature objects - one for each fasta
   *  record in the DocumentEntry.  These won't be written when the
   *  DocumentEntry is written.
   **/
  private FeatureVector fake_fasta_features = new FeatureVector();
 
  /** EntryInformation object that was passed to constructor */
  final private EntryInformation entry_information;
                                                                                                     
  /** collection to send ReadEvents to */
  private Vector<ReadListener> listeners = new Vector<ReadListener>();

  /** 
   *  The Document object that was passed to the constructor.  This should be
   *  the document that this SimpleDocumentEntry was read from.
   **/
  private Document document = null;
                                                                                                               
  /**
   *  This contains all the lines(stored as LineGroup objects) from the entry
   *  stream that was passed to the constructor.
   **/
  protected LineGroupCache line_groups = new LineGroupCache();
                                                                                                               
  /**
   *  The DocumentEntryAutosaveThread that is started when the first call is
   *  made to setDirtyFlag().
   **/
  private Thread autosave_thread = null;
                                                                                                               
  /**
   *  The Date when this Entry last changed or null if this Entry
   *  hasn't changed since the last save.  Set to null by save().
   **/
  private java.util.Date last_change_time = null;
                                                                                                               
  /**
   *  Set to true in the constructor while features are added.  setDirtyFlag()
   *  will do nothing while this is true.
   **/
  protected boolean in_constructor = false;
  
  protected Hashtable<String, Range> contig_ranges;

  /**
   *  Create a new SimpleDocumentEntry from the given Document.
   *  @param entry_information The EntryInformation object of the new Entry.
   *  @param document This is the file that we will read from.  This is also
   *    used for saving the entry back to the file it came from and to give
   *    the new object a name.
   *  @exception IOException thrown if there is a problem reading the entry -
   *    most likely ReadFormatException.
   *  @exception EntryInformationException Thrown if force is false and if this
   *    Entry cannot contain the Key, Qualifier or Key/Qualifier combination of
   *    one of the features in the given Entry.
   **/
  public SimpleDocumentEntry(final EntryInformation entry_information,
                              final Document document,
                              final ReadListener read_listener)
      throws IOException, EntryInformationException 
  {
    this.document = document;
    this.entry_information = new SimpleEntryInformation(entry_information);
    this.in_constructor    = true;  // flag used by setDirtyFlag()

    if(read_listener != null)
      addReadListener(read_listener);

    final LinePushBackReader pushback_reader =
              getDocument().getLinePushBackReader();

    LineGroup new_line_group;

    final int MAX_LOOP = 9999;
    
    while((new_line_group =
            LineGroup.readNextLineGroup(pushback_reader, this)) != null) 
    {       
      if(new_line_group instanceof SimpleDocumentFeature)
      {
        final SimpleDocumentFeature new_feature =
                   (SimpleDocumentFeature)new_line_group;
        
        // try several times because adding the Feature may cause more than
        // one exception
        int i;
        EntryInformationException saved_error = null;

        for(i = 0; i<MAX_LOOP; ++i) 
        {
          try 
          {
            addInternal(new_feature, true);
            break;
          } 
          catch(EntryInformationException e) 
          {
            getEntryInformation().fixException(e);
            saved_error = e;
          }
        }

        if(i == MAX_LOOP) 
          throw new Error("internal error - too many exceptions: " +
                           saved_error.getMessage());
      }
      else 
        addLineGroup(new_line_group);
      
      if(new_line_group instanceof IndexFastaStream)
        break;
    }

    pushback_reader.close();
    
    // we added some features above hence:
    last_change_time = null;

    final Sequence sequence = getSequence();

    if(sequence != null && sequence instanceof FastaStreamSequence) 
    {
      // add a feature for each FASTA record if there are more 
      // than one record in the FASTA sequence
      final FastaStreamSequence fasta_sequence =
                                (FastaStreamSequence)sequence;

      final String[] header_strings = fasta_sequence.getFastaHeaderStrings();

      if(header_strings.length > 1) 
      {
        final int[] header_positions =
                      fasta_sequence.getFastaHeaderPositions();

        getFeatureTable();

        for(int i = 0 ; i < header_strings.length ; ++i) 
        {
          try
          {
            final Range new_range;

            if(i == header_strings.length - 1) 
            {
              if(header_positions[i] == fasta_sequence.length()) 
                throw new ReadFormatException("empty FASTA record: >" +
                                               header_strings[i]);

              new_range = new Range(header_positions[i] + 1,
                                    fasta_sequence.length());
            }  
            else
            {
              if(header_positions[i] == header_positions[i+1]) 
                throw new ReadFormatException("empty FASTA record: >" +
                                               header_strings[i]);

              new_range = new Range(header_positions[i] + 1,
                                    header_positions[i+1]);
            }

            String thisHeader[] = header_strings[i].split("\\s");
            final QualifierVector qualifiers = new QualifierVector();

            qualifiers.setQualifier(new Qualifier("note",
                                                    header_strings[i]));
            qualifiers.setQualifier(new Qualifier("label", thisHeader[0]));
            if(i % 2 == 0) 
              qualifiers.setQualifier(new Qualifier("colour", "10"));
            else 
              qualifiers.setQualifier(new Qualifier("colour", "11"));

            //ReadOnlyEmblStreamFeature
            final EmblStreamFeature new_feature =
              new EmblStreamFeature(new Key("fasta_record"),
                                             new Location(new_range),
                                             qualifiers);

            fake_fasta_features.add(new_feature);

            // record coordinates to adjust feature coordinates
            if(contig_ranges == null)
              contig_ranges = new Hashtable<String, Range>();           

            // find the sequence id from the header
            contig_ranges.put(thisHeader[0], new_range);
          }
          catch(InvalidRelationException e) 
          {
            throw new Error("internal error - unexpected exception: " + e);
          }
          catch(OutOfRangeException e)
          {
            throw new Error("internal error - unexpected exception: " + e);
          }
        }
        addFakeFeatures();
      }
    }

    this.in_constructor = false;
  }

  /**
   *  Create a new SimpleDocumentEntry with no Document associated with it.
   *  @param entry_information The EntryInformation object of the new Entry.
   **/
  public SimpleDocumentEntry(final EntryInformation entry_information)
  {
    this.entry_information = new SimpleEntryInformation(entry_information);
  }

  /**
   *  Create a new SimpleDocumentEntry that will be a copy of the given Entry
   *  and has no Document associated with it.  The new SimpleDocumentEntry
   *  cannot be saved to a file with save() unless save(Document) has been
   *  called first.  Some qualifier and location information will be lost if
   *  the argument Entry is a different type to this class.
   *  @param entry_information The EntryInformation object of the Entry that
   *    will contain this Feature.
   *  @param force If true then invalid qualifiers and any features with
   *    invalid keys in the new Entry will be quietly thrown away.  "Invalid"
   *    means that the key/qualifier is not allowed to occur in an Entry of
   *    this type(probably determined by the EntryInformation object of this
   *    Entry).  If false an EntryInformationException will be thrown for
   *    invalid keys or qualifiers.
   **/
  public SimpleDocumentEntry(final EntryInformation entry_information,
                              final Entry new_entry, final boolean force)
      throws EntryInformationException
  {
    this.entry_information = new SimpleEntryInformation(entry_information);

    if(new_entry.getClass().equals(this.getClass()) ||
        (this instanceof GFFDocumentEntry && new_entry instanceof DatabaseDocumentEntry)) 
    {
      try
      {
        setHeaderText(new_entry.getHeaderText());
      } 
      catch(IOException e) 
      {
        System.err.println(e);
        // if it doesn't work just ignore it
      }
    }

    final FeatureEnumeration feature_enum = new_entry.features();

    Set<String> failed = null;
    while(feature_enum.hasMoreFeatures()) 
    {
      final Feature new_feature = feature_enum.nextFeature();

      try
      {
        if(force) 
        {
          Feature f = (SimpleDocumentFeature)makeNativeFeature(new_feature, true);
          
          if(f != null)
          {
            if(forcedAdd(f) == null)
            {
              if(failed == null)
                failed = new HashSet<String>();
              failed.add(new_feature.getKey().getKeyString());
            }
          }
        }
        else
        {
          final Object docFeature = makeNativeFeature(new_feature, true);
          if(docFeature instanceof SimpleDocumentFeature[])
          {
            SimpleDocumentFeature[] docFeatures = (SimpleDocumentFeature[])docFeature;
            for(int i=0; i<docFeatures.length; i++)
              add((SimpleDocumentFeature)docFeatures[i]);
          }
          else if(docFeature != null)
            add((SimpleDocumentFeature)docFeature);
        }
      } 
      catch(ReadOnlyException e) 
      {
        throw new Error("internal error - unexpected exception: " + e);
      }
    }

    if(failed != null)
      UI.warn("Failed to use the following keys\n"+failed.toString(), "Warning - unknown keys");
    
    final Sequence new_sequence = new_entry.getSequence();

    if(new_sequence != null) 
      setSequence(makeNativeSequence(new_sequence));
  }

  /**
   *  Return the text of the header of this Entry or null if there is no
   *  header.
   **/
  public String getHeaderText() 
  {
	return line_groups.getHeadersAsText();
  }

  /**
   *  Set the header of this Entry to be the given text.
   *  @return true if and only if the header was successfully set.  Not all
   *    Entry objects can change their header, so it is up to the calling
   *    function to check the return value.  If null the current header will be
   *    removed.
   *  @exception IOException thrown if there is a problem reading the header
   *    from the String - most likely ReadFormatException.
   **/
  public boolean setHeaderText(final String new_header)
      throws IOException 
  {

    final Vector<LineGroup> new_line_groups = new Vector<LineGroup>();

    if(new_header != null) 
    {
      final Reader reader = new StringReader(new_header);

      final LinePushBackReader pushback_reader =
                          new LinePushBackReader(reader);

      LineGroup new_line_group;

      try 
      {
        while((new_line_group =
                LineGroup.readNextLineGroup(pushback_reader, this)) != null)
        {
          if(new_line_group instanceof MiscLineGroup) 
            new_line_groups.addElement(new_line_group);
          else
            throw new ReadFormatException("the header must contain only " +
                                           "header lines");
        }
      } 
      catch(InvalidRelationException e) 
      {
        throw new ReadFormatException("the header must contain only " +
                                       "header lines");
      }
    }

    // now remove all the EmblMisc and GenbankMisc LineGroup objects from
    // this Entry
    line_groups.clearMiscLineGroups();

    if(new_header != null) 
    {
      // then add the new LineGroup objects
      for(int i = 0 ; i < new_line_groups.size() ; ++i) 
    	  line_groups.add(new_line_groups.elementAt(i));
    }

    setDirtyFlag();

    return true;
  }

  /**
   *  Write this Entry to the given stream.
   *  @param writer The stream to write to.
   *  @exception IOException thrown if there is a problem writing the entry.
   **/
  public void writeToStream(final Writer writer)
      throws IOException 
  {
    removeFakeFeatures();

    try
    {
      line_groups.writeToStream(this, writer);
    } 
    finally
    {
      addFakeFeatures();
    }
  }

  /**
   *  Return a count of the number of Feature objects in this Entry.
   **/
  public int getFeatureCount() 
  {
    final FeatureTable feature_table = line_groups.getFeatureTable();

    if(feature_table == null) 
      return 0;
    else 
      return feature_table.getFeatureCount();
  }

  /**
   *  The method is identical to add() except that it doesn't check the
   *  read_only flag before adding and it doesn't throw ReadOnlyExceptions.
   *  This is used by the constructor to avoid exceptions.
   *
   *  Add the given Feature to this Entry.  If the Feature is already in an
   *  Entry then Entry.remove() should be called on that Entry before calling
   *  Entry.add().  An Error will be thrown otherwise.
   *  @param throw_entry_info_exceptions if true throw
   *    EntryInformationExceptions, otherwise just send a ReadEvent.
   *  @exception EntryInformationException Thrown if this Entry
   *    cannot contain the Key, Qualifier or Key/Qualifier combination of the
   *    given Feature or if a required qualifier is missing.
   *  @return A reference that was passed to add(), if that Feature can be
   *    stored directly in this Entry, otherwise returns a reference to a new
   *    Feature, that is a copy of the argument.  The argument reference
   *    should not be used after the call to add(), unless the return
   *    reference happens to be the same as the argument.
   **/
  private Feature addInternal(final Feature feature,
                               final boolean throw_entry_info_exceptions)
      throws EntryInformationException 
  {
    if(feature.getEntry() != null) 
      throw new Error("internal error - a feature must have one owner");

    final EntryInformation entry_information = getEntryInformation();

    if(!entry_information.isValidKey(feature.getKey()))  
    {
      final String message = feature.getKey() + " is not a valid key";

      fireEvent(new ReadEvent(this, message));

      throw new InvalidKeyException(feature.getKey() + " is not a valid " +
                                     "key for this entry", feature.getKey());
    }

    final QualifierVector new_qualifiers = feature.getQualifiers();

    final Key new_key = feature.getKey();

    // check the qualifiers
    int numQual = new_qualifiers.size();
    for(int i = 0 ; i < numQual ; ++i) 
    {
      final Qualifier this_qualifier = (Qualifier)new_qualifiers.elementAt(i);
      final String this_qualifier_name = this_qualifier.getName();

      if(!entry_information.isValidQualifier(new_key, this_qualifier_name)) 
      {
        final String message = new_key + " can't have " + this_qualifier_name 
                                       + " as a qualifier";

        fireEvent(new ReadEvent(this, message));

        if(throw_entry_info_exceptions) 
          throw new InvalidRelationException(message, new_key,
                                             this_qualifier);
      }
    }

    final SimpleDocumentFeature native_feature =
      (SimpleDocumentFeature)makeNativeFeature(feature, false);

    final FeatureTable feature_table = getFeatureTable();
    feature_table.add(native_feature);

    try 
    {
      native_feature.setDocumentEntry(this);
    }
    catch(ReadOnlyException e) 
    {
      // makeNativeFeature() should never return a read only feature
      throw new Error("internal error - unexpected exception: " + e);
    }

    setDirtyFlag();

    return native_feature;
  }

  /**
   *  Add the given Feature to this Entry.  If the Feature is already in an
   *  Entry then Entry.remove() should be called on that Entry before calling
   *  Entry.add().  An Error will be thrown otherwise.
   *  @exception ReadOnlyException If this entry is read only.
   *  @exception EntryInformationException Thrown if this Entry
   *    cannot contain the Key, Qualifier or Key/Qualifier combination of the
   *    given Feature or if a required qualifier is missing.
   *  @return A reference that was passed to add(), if that Feature can be
   *    stored directly in this Entry, otherwise returns a reference to a new
   *    Feature, that is a copy of the argument.  The argument reference
   *    should not be used after the call to add(), unless the return
   *    reference happens to be the same as the argument.
   **/
  public Feature add(final Feature feature)
      throws EntryInformationException, ReadOnlyException 
  {
    if(isReadOnly()) 
      throw new ReadOnlyException();

    return addInternal(feature, true);
  }

  /**
   *  Add the given Feature to this Entry.  If the Feature is already in an
   *  Entry then Entry.remove() should be called on that Entry before calling
   *  Entry.forcedAdd()(An Error will be thrown otherwise).  Invalid
   *  qualifiers will be quietly thrown away.  Features with invalid keys will
   *  not be added(and null will be returned).  "Invalid" means that the
   *  key/qualifier is non allowed to occur in an Entry of this type(probably
   *  determined by the EntryInformation object of this Entry).  Any
   *  qualifiers that are required for this Entry will be quietly added(with
   *  a zero-length string as the value).
   *  @exception ReadOnlyException If this entry is read only.
   *  @return A reference that was passed to add(), if that Feature can be
   *    stored directly in this Entry, otherwise returns a reference to a new
   *    Feature, that is a copy of the argument.  The argument reference
   *    should not be used after the call to add(), unless the return
   *    reference happens to be the same as the argument.  Returns null if and
   *    only if the new Feature has a key that is invalid for this Entry.
   **/
  public Feature forcedAdd(final Feature feature)
      throws ReadOnlyException 
  {
    if(isReadOnly()) 
      throw new ReadOnlyException();

    if(feature.getEntry() != null) 
      throw new Error("internal error - a feature must have one owner");

    final EntryInformation entry_information = getEntryInformation();

    if(!entry_information.isValidKey(feature.getKey())) 
      return null;

    final QualifierVector feature_qualifiers = feature.getQualifiers();
    final QualifierVector fixed_qualifiers = new QualifierVector();
    final Key new_key = feature.getKey();

    // set to true if there is an invalid qualifier
    boolean qualifiers_fixed = false;

    // check the qualifiers
    for(int i = 0 ; i < feature_qualifiers.size() ; ++i)
    {
      final Qualifier this_qualifier = (Qualifier)feature_qualifiers.elementAt(i);

      final String this_qualifier_name = this_qualifier.getName();

      if(entry_information.isValidQualifier(new_key, this_qualifier_name)) 
        fixed_qualifiers.setQualifier(this_qualifier);
      else 
        qualifiers_fixed = true;
    }

    final SimpleDocumentFeature native_feature =
      (SimpleDocumentFeature)makeNativeFeature(feature, false);

    if(qualifiers_fixed) 
    {
      try
      {
        native_feature.setQualifiers(fixed_qualifiers);
      }
      catch(EntryInformationException e) 
      {
        throw new Error("internal error - unexpected exception: " + e);
      }
    } 
    else
    {
      // otherwise use the feature as is
    }

    final FeatureTable feature_table = getFeatureTable();
    feature_table.add(native_feature);
    native_feature.setDocumentEntry(this);
    setDirtyFlag();

    return native_feature;
  }

  /**
   *  The method is identical to remove() except that it doesn't check the
   *  read_only flag before removing and it doesn't throw ReadOnlyExceptions.
   *
   *  Remove the given Feature from this Entry.
   *  @return true if and only if the Feature was in this Entry.
   **/
  protected boolean removeInternal(Feature feature) 
  {
    final FeatureTable feature_table = line_groups.getFeatureTable();

    if(feature_table == null) 
      return false;
    else 
    {
      final SimpleDocumentFeature feature_from_table =
       (SimpleDocumentFeature) feature_table.remove(feature);

      if(feature_from_table == null) 
        return false;
      else  
      {
        try 
        {
          feature_from_table.setDocumentEntry(null);
        } 
        catch(ReadOnlyException e) 
        {
          throw new Error("internal error - unexpected exception: " + e);
        }

        // get rid of the feature table
        if(feature_table.getFeatureCount() == 0) 
          removeLineGroup(feature_table);

        setDirtyFlag();

        return true;
      }
    }
  }

  /**
   *  Remove the given Feature from this Entry.
   *  @return true if and only if the Feature was in this Entry.
   **/
  public boolean remove(Feature feature)
      throws ReadOnlyException 
  {
    if(isReadOnly() || feature.isReadOnly()) 
      throw new ReadOnlyException();

    return removeInternal(feature);
  }

  /**
   *  Return the ith Feature from this Entry.  This Features are returned in a
   *  consistent order, sorted by the first base of each Feature.
   **/
  public Feature getFeatureAtIndex(int i) 
  {
    final FeatureTable feature_table = line_groups.getFeatureTable();

    if(feature_table == null) 
      return null;
    else 
      return feature_table.getFeatureAtIndex(i);
  }

  /**
   *  Return the index of the given Feature.  This does the reverse of
   *  getFeatureAtIndex().
   **/
  public int indexOf(final Feature feature) 
  {
    final FeatureTable feature_table = line_groups.getFeatureTable();

    if(feature_table == null) 
      return -1;
    else 
      return feature_table.indexOf(feature);
  }

  /**
   *  Returns true if and only if this Entry contains the given feature.
   **/
  public boolean contains(final Feature feature) 
  {
    final FeatureTable feature_table = line_groups.getFeatureTable();

    if(feature_table == null) 
      return false;
    else
      return feature_table.contains(feature);
  }

  /**
   *  Create a new Feature object of an appropriate type in this Entry.
   *  @param key The new feature key
   *  @param location The Location object for the new feature
   *  @param qualifiers The qualifiers for the new feature(can be null if
   *    there are no qualifiers).
   **/
  public Feature createFeature(final Key key,
                                final Location location,
                                final QualifierVector qualifiers)
      throws EntryInformationException, ReadOnlyException,
             OutOfRangeException 
  {
    if(isReadOnly())
      throw new ReadOnlyException();

    final Feature new_feature;

    if(this instanceof EmblDocumentEntry) 
      new_feature = new EmblStreamFeature(key, location, qualifiers);
    else if(this instanceof DatabaseDocumentEntry)
      new_feature = new DatabaseStreamFeature(key, location, qualifiers);
    else if(this instanceof GFFDocumentEntry)
      new_feature = new GFFStreamFeature(key, location, qualifiers);
    else
      new_feature = new GenbankStreamFeature(key, location, qualifiers);

    add(new_feature);
    setDirtyFlag();

    return new_feature;
  }

  /**
   *  Return a vector containing the references of the Feature objects within
   *  the given range.
   *  @param range Return features that overlap this range - ie the start of
   *    the feature is less than or equal to the end of the range and the end
   *    of the feature is greater than or equal to the start of the range.
   *  @return The features of this feature table the are within
   *    the given range.  The returned object is a copy - changes will not
   *    effect the FeatureTable object itself.
   **/
  public FeatureVector getFeaturesInRange(Range range) 
  {
    final FeatureTable feature_table = line_groups.getFeatureTable();

    if(feature_table == null) 
      return new FeatureVector();
    else 
      return feature_table.getFeaturesInRange(range);
  }

  /**
   *  Return a vector containing the references of all the Feature objects in
   *  this Entry.
   *  @return The features of this Entry.  The returned object
   *    is a copy - changes will not effect the Entry object itself.
   **/
  public FeatureVector getAllFeatures() 
  {
    final FeatureTable feature_table = line_groups.getFeatureTable();

    if(feature_table == null) 
      return new FeatureVector();
    else 
      return feature_table.getAllFeatures();
  }

  /**
   *  Returns an enumeration of the Feature objects in this
   *  SimpleDocumentEntry.  The returned Enumeration object will generate
   *  all features in this object in turn. The first item generated is the
   *  item at index 0, then the item at index 1, and so on.
   **/
  public FeatureEnumeration features() 
  {
    final FeatureTable feature_table = line_groups.getFeatureTable();

    if(feature_table == null)
    {
      return new FeatureEnumeration() {
        public boolean hasMoreFeatures() {
          return false;
        }

        public Feature nextFeature() {
          return null;
        }
      };
    }
    else 
      return feature_table.features();
  }

  /**
   *  Return the Sequence object from this entry or null if it does not
   *  contain one.
   *  @return a Sequence object for this Entry.  the returned object is
   *    not a copy - changes to it will change the Entry object itself
   **/
  public Sequence getSequence() 
  {
    return (Sequence)line_groups.getSequence();
  }

  /**
   *  Add a new LineGroup object to this Entry.
   *  @param new_line_group A new LineGroup to add.
   **/
  protected void addLineGroup(final LineGroup new_line_group) 
  {
    if(new_line_group instanceof GFFMisc)
    {
      String line = ((GFFMisc)new_line_group).getLine();
      if(line.indexOf("FASTA") > -1)  // ignore
        return;
    }

    line_groups.add(new_line_group);
  }

  /**
   *  Return the name of this Entry or null if it has no name.
   **/
  public String getName() 
  {
    if(getDocument() == null || getDocument().getName() == null) 
      return null;
    else 
      return getDocument().getName();
  }

  /**
   *  Set the name of this Entry - if possible(the return value will let the
   *  caller know).
   *  @return true if and only if the name was successfully set.  Not all
   *    Entry objects can change there name, so it is up to the calling
   *    function to check the return value.
   **/
  public boolean setName(final String name) 
  {
    if(name == null || name.length() == 0) 
      return false;

    setDocument(new FileDocument(new File(name)));

    return true;
  }

  /**
   *  Return the EntryInformation object for this Entry.
   **/
  public EntryInformation getEntryInformation() 
  {
    return entry_information;
  }

  /**
   *  Return the FeatureTable object from this entry(which will be created if
   *  this object does not already contain one).  This method has package
   *  scope because other package should add, remove and query Features in the
   *  FeatureTable indirectly(through the Entry class).
   *  @return a FeatureTable object for this Entry.  the returned object is
   *    not a copy - changes to it will change the Entry object itself
   **/
  protected FeatureTable getFeatureTable() 
  {
    final FeatureTable found_feature_table = line_groups.getFeatureTable();

    if(found_feature_table == null) 
      return line_groups.createFeatureTable();
    else
      return found_feature_table;
  }

  /**
   *  Remove the given LineGroup from this SimpleDocumentEntry.
   **/
  protected void removeLineGroup(final LineGroup line_group) 
  {
	LineGroup group = line_groups.remove(line_group);
	if (group != null)
	{
	  setDirtyFlag();	
	}
  }

  /**
   *  Return the File reference that was passed to the constructor or null if
   *  none was passed.
   **/
  public Document getDocument() 
  {
    return document;
  }


  /**
   *  Set the document to use when saving this DocumentEntry.
   **/
  public void setDocument(final Document document) 
  {
    this.document = document;
  }


  /**
   *  Write this Entry to the Document it came from.  This method uses the
   *  current Document reference(as given by getDocument()) and will throw a
   *  NullPointerException if the is no current Document ie. if getDocument()
   *  returns null.  Use save(Document) to save the current Entry
   *  to a different Document.
   *  @exception IOException thrown if there is an IO problem saving the entry.
   **/
  public void save()
      throws IOException
  {
    save(getDocument());
  }


  /**
   *  Write this Entry to the given Document.
   *  @param document This is the file that we will write to.
   *  @exception IOException thrown if there is an IO problem saving the entry.
   **/
  public void save(final Document document)
      throws IOException 
  {
    final Writer out_file;
    try
    {
      out_file = document.getWriter();
    }
    catch(NullPointerException npe)
    {
      return;
    }

    writeToStream(out_file);
    out_file.close();

    last_change_time = null;
  }


  /**
   *  Returns true if and only if this entry is read only.
   **/
  public boolean isReadOnly() 
  {
    return false;
  }


  /**
   *  Returns true if and only if there have been some changes to this Entry
   *  since the last save.
   **/
  public boolean hasUnsavedChanges() 
  {
    return last_change_time != null;
  }


  /**
   *  Set last_change_time so that hasUnsavedChanges() will return true.
   *  last_change_time can be read with getLastChangeTime().
   **/
  public void setDirtyFlag() 
  {
    if(in_constructor) 
    {
      // features are being added, but the entry hasn't really changed
    }
    else
    {
       if(autosave_thread == null && getName() != null) 
       {
         // this is the first change so start autosaving
         autosave_thread = new DocumentEntryAutosaveThread(this);
         autosave_thread.start();
       }

      final java.util.Calendar calendar = java.util.Calendar.getInstance();
      last_change_time = calendar.getTime();
    }
  }


  /**
   *  Return the Date when this Entry last changed.
   **/
  public Date getLastChangeTime()
  {
    return last_change_time;
  }


  /**
   *  Add the given Sequence object to this Entry.
   **/
  public void setSequence(final StreamSequence sequence) 
  {
    final Sequence old_sequence = getSequence();

    if(old_sequence != null) 
      removeLineGroup((StreamSequence) old_sequence);

    addLineGroup(sequence);
  }


  /**
   *  If the given feature can be added directly to this Entry, then return
   *  it, otherwise create and return a new feature of the appropriate type
   * (eg. return an EmblStreamFeature if this Entry is a EmblDocumentEntry).
   *  The feature should not be in an Entry when this method is called.
   *  @param copy if true then always make a new copy of the Feature.
   **/
  protected abstract Object
    makeNativeFeature(final Feature feature,
                       final boolean copy);

  /**
   *  If the given Sequence can be added directly to this Entry, then return a
   *  copy of it, otherwise create and return a new feature of the appropriate
   *  type for this Entry.
   **/
  protected abstract StreamSequence
    makeNativeSequence(final Sequence sequence);

  /**
   *  Add the elements of fake_fasta_features to the feature table.
   **/
  protected void addFakeFeatures() 
  {
    final FeatureTable feature_table = getFeatureTable();

    for(int i = 0 ; i < fake_fasta_features.size() ; ++i) 
      feature_table.add(fake_fasta_features.featureAt(i));
  }

  /**
   *  Remove the elements of fake_fasta_features from the feature table.
   **/
  protected void removeFakeFeatures() 
  {
    final FeatureTable feature_table = getFeatureTable();

    for(int i = 0 ; i < fake_fasta_features.size() ; ++i) 
      feature_table.remove(fake_fasta_features.featureAt(i));
  }

  /**
   *  Add an object that will receive ReadEvents.
   **/
  private void addReadListener(final ReadListener listener) 
  {
    listeners.addElement(listener);
  }

  /**
   *  Send the given event to all the listeners.
   **/
  private void fireEvent(final ReadEvent event)
  {
    for(int i = 0 ; i < listeners.size() ; ++i) 
     ((ReadListener)listeners.elementAt(i)).notify(event);
  }

  public void dispose()
  {
    line_groups.clear();
  }
}
