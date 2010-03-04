/* Entry.java
 *
 * created: Sun Oct 11 1998
 *
 * This file is part of Artemis
 *
 * Copyright(C) 1998,1999,2000  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/Entry.java,v 1.11 2008-06-20 10:00:25 tjc Exp $
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.sequence.*;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.EmblStreamFeature;
import uk.ac.sanger.artemis.io.DocumentEntry;
import uk.ac.sanger.artemis.io.EmblDocumentEntry;
import uk.ac.sanger.artemis.io.GFFDocumentEntry;
import uk.ac.sanger.artemis.io.PartialSequence;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.RangeVector;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.LocationParseException;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.DocumentEntryFactory;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.EntryInformationException;

import java.util.NoSuchElementException;
import java.util.Vector;
import java.io.IOException;
import java.io.File;


/**
 *  This class is a wrapper for the io.Entry class which contains the state
 *  and handles the events needed for editing one embl entry in Diana.  The
 *  state includes an EMBL.Entry object and a list of EntryChange listeners
 *  and FeatureChange listeners.  Other objects can register as listeners for
 *  changes to the entry.  For that reason, all changes to the embl.Entry
 *  objects should go through this class. (see ChangeEvent for details of the
 *  possible events.)
 *
 *  @author Kim Rutherford
 *  @version $Id: Entry.java,v 1.11 2008-06-20 10:00:25 tjc Exp $
 **/

public class Entry implements FeatureChangeListener, Selectable 
{

  /**
   *  The embl.Entry object that was passed to the constructor
   **/
  private uk.ac.sanger.artemis.io.Entry embl_entry;

  /**
   *  A vector of those objects listening for entry change events.
   **/
  final private Vector entry_listener_list = new Vector();

  /**
   *  A vector of those objects listening for feature change events.
   **/
  final private Vector feature_listener_list = new Vector();

  /**
   *  This is the Bases reference that was passed to the constructor.
   **/
  /*final*/ private Bases bases;

  /**
   *  Create a new Entry object.
   *  @param bases The Bases object which contains the Strand objects that will
   *    be used by the features of this Entry.
   *  @param embl_entry a reference to an embl.Entry object containing the
   *    underlying data for new object
   *  @exception OutOfRangeException Thrown if one of the features in
   *    embl_entry is out of range of the Bases object.
   **/
  public Entry(final Bases bases,
               final uk.ac.sanger.artemis.io.Entry embl_entry)
      throws OutOfRangeException 
  {
    this.embl_entry = embl_entry;
    this.bases = bases;

    checkLocations();
    createDianaFeatures();
  }

  /**
   *  Create a new Entry object from a uk.ac.sanger.artemis.io.Entry
   *  object.  A new Bases object will be created for the new Entry.
   *  @param embl_entry a reference to an embl.Entry object containing the
   *    underlying data for new object
   *  @exception OutOfRangeException Thrown if one of the features in
   *    embl_entry is out of range of the Bases object.
   *  @exception NoSequenceException Thrown if embl_entry contains no sequence.
   **/
  public Entry(final uk.ac.sanger.artemis.io.Entry embl_entry)
      throws OutOfRangeException, NoSequenceException 
  {
    this.embl_entry = embl_entry;

    if(embl_entry.getSequence() == null ||
       embl_entry.getSequence().length() == 0) 
      throw new NoSequenceException();
    else 
      this.bases = new Bases(embl_entry.getSequence());

    checkLocations();
    createDianaFeatures();
  }
  
  
  /**
   *  Returns true if and only if this entry is read only.
   **/
  public boolean isReadOnly() 
  {
    return getEMBLEntry().isReadOnly();
  }

  /**
   *  Save the changes to this Entry back to where the Entry came from.  If
   *  this object is read only a ReadOnlyException will be thrown.
   *  @param destination_type Should be one of EMBL_FORMAT, GENBANK_FORMAT or
   *    ANY_FORMAT.  If ANY_FORMAT then the Entry will be saved in the
   *    same format it was created, otherwise it will be saved in the given
   *    format.
   *  @exception EntryInformationException Thrown if the destination type
   *    cannot contain the Key, Qualifier or Key/Qualifier combinations of one
   *    of the features in this Entry.
   *  @exception IOException Thrown if an IO error occurs while saving.  One
   *    possibility is ReadOnlyException.
   **/
  public void save(final int destination_type)
      throws IOException, EntryInformationException 
  {
    if(destination_type == DocumentEntryFactory.ANY_FORMAT) 
      getEMBLEntry().save();
    else 
    {
      if(getEMBLEntry() instanceof DocumentEntry) 
        getEMBLEntryAsNewType(getEntryInformation(),
                              destination_type, false).save();
      else 
        throw new ReadOnlyException("operation cannot be " +
                                    "applied to this entry");
    }
  }

  /**
   *  Save the changes to this Entry back to where the Entry came from(like
   *  save()) but without any uk.ac.sanger.artemis extensions.  If this object is read only
   *  or doesn't support this a ReadOnlyException will be thrown.
   *  @param destination_type Should be one of EMBL_FORMAT, GENBANK_FORMAT or
   *    ANY_FORMAT.  If ANY_FORMAT then the Entry will be saved in the
   *    same format it was created, otherwise it will be saved in the given
   *    format.
   *  @exception EntryInformationException Thrown if the destination type
   *    cannot contain the Key, Qualifier or Key/Qualifier combinations of one
   *    of the features in this Entry.
   **/
  public void saveStandardOnly(final int destination_type)
      throws IOException,EntryInformationException 
  {
    if(destination_type == DocumentEntryFactory.ANY_FORMAT) 
      getEMBLEntry().save();
    else 
    {
      if(getEMBLEntry() instanceof DocumentEntry) 
      {
        final EntryInformation entry_information =
          Options.getDBEntryInformation();

        getEMBLEntryAsNewType(entry_information,
                               destination_type, false).save();
      } 
      else 
        throw new ReadOnlyException("operation cannot be " +
                                    "applied to this entry");
    }
  }

  /**
   *  Save the changes to this Entry to the file.  If this object is read only
   *  a ReadOnlyException will be thrown.  This method only works when the
   *  underlying embl.Entry object is a DocumentEntry.
   *  @param destination_type Should be one of EMBL_FORMAT, GENBANK_FORMAT,
   *    GFF_FORMAT or ANY_FORMAT.  If ANY_FORMAT then the Entry will
   *    be saved in the same format it was created, otherwise it will be saved
   *    in the given format.
   *  @param force If true then invalid qualifiers and any features with
   *    invalid keys will be quietly thrown away when saving.  "Invalid" means
   *    that the key/qualifier is not allowed to occur in the destination type
   *   (probably determined by the default EntryInformation object of the
   *    destination type).  If false an EntryInformationException will be
   *    thrown for invalid keys or qualifiers.
   *  @exception EntryInformationException Thrown if force is false and if the
   *    destination type cannot contain the Key, Qualifier or Key/Qualifier
   *    combination of one of the features in this Entry.
   **/
  public void save(final File file, final int destination_type,
                   final boolean force)
      throws IOException, EntryInformationException 
  {
    final EntryInformation artemis_entry_information;
    
    if((getEMBLEntry() instanceof DatabaseDocumentEntry ||
        getEMBLEntry() instanceof GFFDocumentEntry) &&
       destination_type == DocumentEntryFactory.EMBL_FORMAT)
      artemis_entry_information = Options.getArtemisEntryInformation();
    else
      artemis_entry_information = getEntryInformation();
    
    save(file, destination_type, force, artemis_entry_information);
  }
  
  
  /**
   *  Save the changes to this Entry to the file.  If this object is read only
   *  a ReadOnlyException will be thrown.  This method only works when the
   *  underlying embl.Entry object is a DocumentEntry.
   *  @param destination_type Should be one of EMBL_FORMAT, GENBANK_FORMAT,
   *    GFF_FORMAT or ANY_FORMAT.  If ANY_FORMAT then the Entry will
   *    be saved in the same format it was created, otherwise it will be saved
   *    in the given format.
   *  @param force If true then invalid qualifiers and any features with
   *    invalid keys will be quietly thrown away when saving.  "Invalid" means
   *    that the key/qualifier is not allowed to occur in the destination type
   *   (probably determined by the default EntryInformation object of the
   *    destination type).  If false an EntryInformationException will be
   *    thrown for invalid keys or qualifiers.
   *  @exception EntryInformationException Thrown if force is false and if the
   *    destination type cannot contain the Key, Qualifier or Key/Qualifier
   *    combination of one of the features in this Entry.
   **/
  public void save(final File file, final int destination_type,
                   final boolean force, final EntryInformation artemis_entry_information)
      throws IOException, EntryInformationException 
  {
    final FileDocument file_document = new FileDocument(file);
    final DocumentEntry document_entry =
      (DocumentEntry)getEMBLEntryAsNewType(artemis_entry_information,
                                           destination_type, force);
    document_entry.save(file_document);
  }

  /**
   *  Save the changes to this Entry to the file(like save(Document)) but
   *  without any uk.ac.sanger.artemis extensions.  If this object is read only a
   *  ReadOnlyException will be thrown.  This method only works when the
   *  underlying embl.Entry object is a DocumentEntry.
   *  @param destination_type Should be one of EMBL_FORMAT, GENBANK_FORMAT or
   *    ANY_FORMAT.  If ANY_FORMAT then the Entry will be saved in the
   *    same format it was created, otherwise it will be saved in the given
   *    format.
   *  @param force If true then invalid qualifiers and any features with
   *    invalid keys will be quietly thrown away when saving.  "Invalid" means
   *    that the key/qualifier is not allowed to occur in the destination type
   *   (probably determined by the default EntryInformation object of the
   *    destination type).  If false an EntryInformationException will be
   *    thrown for invalid keys or qualifiers.
   *  @exception EntryInformationException Thrown if force is false and if the
   *    destination type cannot contain the Key, Qualifier or Key/Qualifier
   *    combination of one of the features in this Entry.
   **/
  public void saveStandardOnly(final File file,
                               final int destination_type,
                               final boolean force)
      throws IOException, EntryInformationException 
  {
    final FileDocument file_document = new FileDocument(file);

    final DocumentEntry document_entry =
     (DocumentEntry) getEMBLEntryAsNewType(Options.getDBEntryInformation(),
                                           destination_type, force);
    document_entry.save(file_document);
  }

  /**
   *  Returns true if and only if there have been some changes to this Entry
   *  since the last save.
   **/
  public boolean hasUnsavedChanges() 
  {
    return getEMBLEntry().hasUnsavedChanges();
  }

  /**
   *  Create a new Entry in the given EntryGroup.
   *  @return A reference to the new Entry.
   **/
  public static Entry newEntry(final Bases bases) 
  {
    final uk.ac.sanger.artemis.io.Entry entry =
      new EmblDocumentEntry(Options.getArtemisEntryInformation());

    try 
    {
      return new Entry(bases, entry);
    } 
    catch(OutOfRangeException e) 
    {
      // an empty Entry cannot have any out of range features
      throw new Error("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Return the name of this Entry or null if it has no name.
   **/
  public String getName() 
  {
    return getEMBLEntry().getName();
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
    final String orig_entry_name = getName();

    if(getEMBLEntry().setName(name)) 
    {
      if(orig_entry_name == null || ! name.equals(orig_entry_name)) 
      {
        // now inform the listeners that the name has changed
        final EntryChangeEvent event =
          new EntryChangeEvent(name, this, EntryChangeEvent.NAME_CHANGED);

        fireAction(entry_listener_list, event);
      }
      return true;
    } 
    else 
      return false;
  }

  /**
   *  Return the text of the EMBL header of this Entry or null if there is no
   *  header.
   **/
  public String getHeaderText() 
  {
    return getEMBLEntry().getHeaderText();
  }

  /**
   *  Set the header of this Entry to be the given text.
   *  @return true if and only if the header was successfully set.  Not all
   *    Entry objects can change their header, so it is up to the calling
   *    function to check the return value.
   *  @exception IOException thrown if there is a problem reading the header
   *    from the String - most likely ReadFormatException.
   **/
  public boolean setHeaderText(final String new_header)
      throws IOException 
  {
    if(getEMBLEntry().setHeaderText(new_header)) 
    {
      final EntryChangeEvent event =
        new EntryChangeEvent(new_header,
                              this,
                              EntryChangeEvent.HEADER_CHANGED);

      fireAction(entry_listener_list, event);

      return true;
    }
    else 
      return false;
  }

  /**
   *  Return the path of the directory that this entry is in.
   **/
  public Document getRootDocument() 
  {
    if(getEMBLEntry() instanceof DocumentEntry &&  // XXX FIXME:
       ((DocumentEntry) getEMBLEntry()).getDocument() != null) 
    {
      return((DocumentEntry)getEMBLEntry()).getDocument().getParent();
    } 
    else 
      return null;
  }

  /**
   *  Return a vector containing the references of the Feature objects within
   *  the given range.
   *  @param range Return features that overlap this range - ie the start of
   *    the feature is less than or equal to the end of the range and the end
   *    of the feature is greater than or equal to the start of the range.
   *  @return The non-source key features of this feature table the are within
   *    the given range.  The returned object is a copy - changes will not
   *    effect the FeatureTable object itself.
   **/
  public FeatureVector getFeaturesInRange(Range range)
      throws OutOfRangeException 
  {
    final FeatureVector return_features = new FeatureVector();
    final uk.ac.sanger.artemis.io.FeatureVector embl_features =
      getEMBLEntry().getFeaturesInRange(range);

//    System.err.println("starting getFeaturesInRange()");
    Feature diana_feature;
    for(int i = 0 ; i < embl_features.size() ; ++i)
    {
      diana_feature = getFeatureOf(embl_features.featureAt(i));

//       System.err.println("getFeaturesInRange() - diana_feature:" +
//                           diana_feature + " " +
//                           embl_features.elementAt(i));

      return_features.add(diana_feature);
    }

//    System.err.println("ending getFeaturesInRange()"); 

    return return_features;
  }

  /**
   *  Return a vector containing the references of the Feature objects in this
   *  Entry.
   *  @return The non-source key features of this Entry.  The
   *    returned object is a copy - changes will not effect the Entry
   *    object itself.
   **/
  public FeatureVector getAllFeatures() 
  {
    final FeatureVector return_features = new FeatureVector();
    final uk.ac.sanger.artemis.io.FeatureVector embl_features =
      getEMBLEntry().getAllFeatures();

    Feature diana_feature;
    for(int i = 0; i < embl_features.size(); ++i) 
    {
      diana_feature = getFeatureOf(embl_features.featureAt(i));
      return_features.add(diana_feature);
    }

    return return_features;
  }
  
  public void remove(final Feature feature)
      throws ReadOnlyException 
  {
    remove(feature, false);
  }

  /**
   *  Delete a Feature from this object and the underlying embl.Entry object.
   *  @param feature The Feature to delete.
   **/
  public void remove(final Feature feature,
                     final boolean duplicate)
      throws ReadOnlyException 
  {
    synchronized(getBases()) 
    {
      // get the embl.Feature object out of the uk.ac.sanger.artemis.Feature object then
      // remove it
      if(!getEMBLEntry().remove(feature.getEmblFeature())) 
        throw new Error("internal error - remove failed");

      // now inform the listeners that a deletion has occured
      final EntryChangeEvent event =
        new EntryChangeEvent(this, feature, duplicate, EntryChangeEvent.FEATURE_DELETED);

      fireAction(entry_listener_list, event);
      feature.setEntry(null);
    }
  }

  /**
   *  Remove all the features in this Entry.
   **/
  public void removeAllFeatures()
      throws ReadOnlyException 
  {
    // remove the first feature
    while(getFeatureCount() > 0) 
      remove(getFeature(0));
  }

  /**
   *  Add a Feature to the Entry.  The new feature will be inserted in order
   *  in the feature vector.  The features in the vector are ordered by the
   *  first base of each feature.
   *  @param new_feature The new feature to add.  It should not be a member
   *    of the Entry object.
   *  @param force If true then invalid qualifiers will be quietly thrown away
   *    and a feature with invalid keys will not be added.  "Invalid" means
   *    that the key/qualifier is non allowed to occur in an Entry of this type
   *   (probably determined by the EntryInformation object of this Entry).
   *    If false an EntryInformationException will be thrown for invalid keys
   *    or qualifiers.
   *  @exception EntryInformationException Thrown force is false if this Entry
   *    cannot contain the Key, Qualifier or Key/Qualifier combination of the
   *    given Feature
   **/
  public void add(final Feature new_feature,
                  final boolean force)
        throws EntryInformationException, OutOfRangeException,
        ReadOnlyException 
  {
    add(new_feature, false, force);
  }
  
  /**
   *  Add a Feature to the Entry.  The new feature will be inserted in order
   *  in the feature vector.  The features in the vector are ordered by the
   *  first base of each feature.
   *  @param new_feature The new feature to add.  It should not be a member
   *    of the Entry object.
   *  @param force If true then invalid qualifiers will be quietly thrown away
   *    and a feature with invalid keys will not be added.  "Invalid" means
   *    that the key/qualifier is non allowed to occur in an Entry of this type
   *   (probably determined by the EntryInformation object of this Entry).
   *    If false an EntryInformationException will be thrown for invalid keys
   *    or qualifiers.
   *  @exception EntryInformationException Thrown force is false if this Entry
   *    cannot contain the Key, Qualifier or Key/Qualifier combination of the
   *    given Feature
   **/
  public void add(final Feature new_feature, final boolean duplicate,
                  final boolean force)
      throws EntryInformationException, OutOfRangeException,
             ReadOnlyException 
  {
    try
    {
      if(new_feature.getEntry() != null)
        throw new Error("internal error - Feature has a parent");

      new_feature.setEntry(this);

      final uk.ac.sanger.artemis.io.Feature old_embl_reference =
        new_feature.getEmblFeature();

      // this call will add the new_feature in the correct place in the
      // embl.Feature vector of the FeatureTable object
      final uk.ac.sanger.artemis.io.Feature new_embl_reference;

      if(force) 
      {
        new_embl_reference =
          getEMBLEntry().forcedAdd(new_feature.getEmblFeature());
      } 
      else 
      {
        new_embl_reference =
          getEMBLEntry().add(new_feature.getEmblFeature());
      }

      if(new_embl_reference != old_embl_reference) 
      {
        new_feature.getEmblFeature().setUserData(null);
        new_feature.setEmblFeature(new_embl_reference);
        new_embl_reference.setUserData(new_feature);
      }

      // now inform the listeners that a addition has occured
      final EntryChangeEvent event =
        new EntryChangeEvent(this, new_feature, duplicate,
                             EntryChangeEvent.FEATURE_ADDED);

      fireAction(entry_listener_list, event);

    } 
    catch(EntryInformationException e) 
    {
      new_feature.setEntry(null);
      throw e;
    }
  }

  /**
   *  Forget all the features in this entry and inform all the EntryChange
   *  listeners that all the features have gone away.  The garbage collector
   *  will really dispose of the features later.
   **/
  public void dispose() 
  {
    final FeatureEnumeration feature_enum = features();

    while(feature_enum.hasMoreFeatures()) 
    {
      final Feature current_feature = feature_enum.nextFeature();

      // now inform the listeners that a deletion has occured
      final EntryChangeEvent event =
        new EntryChangeEvent(this, current_feature,
                              EntryChangeEvent.FEATURE_DELETED);

      fireAction(entry_listener_list, event);

      current_feature.setEntry(null);
    }
  }

  /**
   *  Returns the Feature at the given index in the vector of Features.
   **/
  public Feature getFeature(int i) 
  {
    return getFeatureOf(getEMBLEntry().getFeatureAtIndex(i));
  }

  /**
   *  Create a new Feature in this Entry with no qualifiers, a default key and
   *  a location of "1..max_base"
   **/
  public Feature createFeature()
      throws ReadOnlyException 
  {
    final Key new_key;
    final Location new_location;
    final QualifierVector new_qualifiers = new QualifierVector();

    try 
    {
      new_key = getEntryInformation().getDefaultKey();
      new_location = new Location("1" + ".." + getBases().getLength());

      return createFeature(new_key, new_location, new_qualifiers);
    } 
    catch(LocationParseException e) 
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
  }

  /**
   *  Create and return a new Feature in this Entry with the given key,
   *  location and qualifiers.  This method creates a new embl.Feature object
   *  and then wraps it in a uk.ac.sanger.artemis.Feature object and 
   *  then inserts it in thisEntry object.
   *  @param new_key The Key to use for the new Feature.
   *  @param new_location The Location to use for the new Feature.
   *  @param new_qualifiers The a vector containing the Qualifier objects to
   *    use for the new Feature.
   *  @exception EntryInformationException Thrown force is false if this Entry
   *    cannot contain a feature with the given Key, Qualifier or
   *    Key/Qualifier combination.
   **/
  public Feature createFeature(Key new_key,
                                Location new_location,
                                QualifierVector new_qualifiers)
      throws EntryInformationException, ReadOnlyException,
      OutOfRangeException
  {
    final uk.ac.sanger.artemis.io.Feature new_embl_feature =
      getEMBLEntry().createFeature(new_key, new_location, new_qualifiers);

    final Feature new_feature = getFeatureOf(new_embl_feature);

    // now inform the listeners that a addition has occured
    final EntryChangeEvent event =
      new EntryChangeEvent(this, new_feature, EntryChangeEvent.FEATURE_ADDED);

    fireAction(entry_listener_list, event);
    return new_feature;
  }

  /**
   *  Create and return a new Feature in this Entry with the given key,
   *  location and no qualifiers.  This method creates a new embl.Feature
   *  object and then wraps it in a uk.ac.sanger.artemis.Feature object and then inserts it
   *  in this Entry object.
   *  @param new_key The Key to use for the new Feature.
   *  @param new_location The Location to use for the new Feature.
   **/
  public Feature createFeature(Key new_key,
                                Location new_location)
      throws EntryInformationException, ReadOnlyException,
      OutOfRangeException 
  {
    final QualifierVector new_qualifiers = new QualifierVector();
    return createFeature(new_key, new_location, new_qualifiers);
  }

  /**
   *  Return true if this Entry contains the given Feature.
   **/
  public boolean contains(Feature feature) 
  {
    return getEMBLEntry().contains(feature.getEmblFeature());
  }

  /**
   *  Return the index of the given Feature in the feature table of this Entry
   *  or -1 if the feature isn't in this Entry.
   **/
  public int indexOf(Feature feature) 
  {
    return getEMBLEntry().indexOf(feature.getEmblFeature());
  }

  /**
   *  Return the index of a feature in the embl feature table of this entry.
   *  The indices start at zero.
   *  @param feature The Feature to lookup.
   *  @return The index of the feature or -1 if the feature isn't in this
   *    Entry.
   **/
  public int getIndexOfFeature(Feature feature) 
  {
    return getEMBLEntry().indexOf(feature.getEmblFeature());
  }

  /**
   *  Returns an enumeration of the Feature objects in this Entry. The
   *  returned FeatureEnumeration object will generate all features in this
   *  object in turn. The first item generated is the item at index 0, then
   *  the item at index 1, and so on.
   **/
  public FeatureEnumeration features() 
  {
    return new FeatureEnumerator();
  }

  /**
   *  An Enumeration of Feature objects.
   **/
  public class FeatureEnumerator implements FeatureEnumeration 
  {
    /**
     *  Create a new FeatureEnumeration that will enumerate the enclosing
     *  Entry object.  The Entry object must not be changed while the
     *  enumeration is active.
     **/
    public FeatureEnumerator()
    {
      feature_enumerator = getEMBLEntry().features();
    }

    /**
     *  See the FeatureEnumeration interface for details.
     **/
    public boolean hasMoreFeatures() 
    {
      return feature_enumerator.hasMoreFeatures();
    }

    /**
     *  See the FeatureEnumeration interface for details.
     **/
    public Feature nextFeature()
        throws NoSuchElementException 
    {

      final uk.ac.sanger.artemis.io.Feature this_feature =
        feature_enumerator.nextFeature();

      return getFeatureOf(this_feature);
    }

    /**
     *  The Enumeration for the current entry
     **/
    private uk.ac.sanger.artemis.io.FeatureEnumeration feature_enumerator;
  }

  /**
   *  Returns the number of features in this entry.
   **/
  public int getFeatureCount() 
  {
    return getEMBLEntry().getFeatureCount();
  }

  /**
   *  Return the base length of the sequence of this entry.
   **/
  private int getSequenceLength() 
  {
    return getBases().getLength();
  }

  /**
   *  Return the object representing the forward sequence of bases for this
   *  Entry.
   **/
  private Strand getForwardStrand() 
  {
    return getBases().getForwardStrand();
  }

  /**
   *  Return the object representing the backward(reverse complemented)
   *  sequence of bases for this Entry.
   **/
  private Strand getReverseStrand() 
  {
    return getBases().getReverseStrand();
  }

  /**
   *  Return the Bases object that was passed to the constructor.
   **/
  public Bases getBases() 
  {
    return bases;
  }

  /**
   *  This method translates the start and end of every each Range in every
   *  Location into another coordinate system.  The Ranges will be truncated
   *  if necessary.
   *  @param constraint This contains the start and end base of the new
   *    coordinate system.  The position given by constraint.getStart() will
   *    be at postion/base 1 in the new coordinate system.
   *  @return a copy of the Entry which has been translated into the new
   *    coordinate system.
   **/
  public Entry truncate(final Bases new_bases, final Range constraint) 
  {
    final uk.ac.sanger.artemis.io.Entry new_embl_entry =
      new EmblDocumentEntry(getEntryInformation());

    final Entry new_entry;

    try 
    {
      new_entry = new Entry(new_bases, new_embl_entry);
    } 
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }

    final FeatureEnumeration feature_enum = features();

    while(feature_enum.hasMoreFeatures()) 
    {
      final Feature current_feature = feature_enum.nextFeature();
      final Location current_feature_location = current_feature.getLocation();

      final Location new_location =
        current_feature_location.truncate(constraint);

      if(new_location != null) 
      {
        final Key current_key = current_feature.getKey();
        final QualifierVector current_qualifiers =
          current_feature.getQualifiers();

        try 
        {
          final uk.ac.sanger.artemis.io.Feature new_embl_feature =
            new EmblStreamFeature(current_key,
                                   new_location,
                                   current_qualifiers);

          final Feature new_feature = new Feature(new_embl_feature);

          new_entry.add(new_feature, false);
        }
        catch(EntryInformationException e) 
        {
          throw new Error("internal error - unexpected exception: " + e);
        } 
        catch(OutOfRangeException e) 
        {
          throw new Error("internal error - unexpected exception: " + e);
        } 
        catch(ReadOnlyException e) 
        {
          throw new Error("internal error - unexpected exception: " + e);
        }
      }
    }

    return new_entry;
  }

  /**
   *  Return a Vector containing those features that don't have valid EMBL
   *  key.
   **/
  public FeatureVector checkForNonEMBLKeys() 
  {
    final FeatureVector non_embl_features = new FeatureVector();
    final FeatureEnumeration feature_enum = features();

    while(feature_enum.hasMoreFeatures()) 
    {
      final Feature current_feature = feature_enum.nextFeature();

      if(!current_feature.hasValidEMBLKey()) 
        non_embl_features.add(current_feature);
    }

    return non_embl_features;
  }

  /**
   *  Returns a vector containing those CDS features that have no /pseudo
   *  qualifier and do not have a valid start codon.  The returned features
   *  may be illegal in EMBL submissions.
   **/
  public FeatureVector checkFeatureStartCodons() 
  {
    // get all the CDS features that do not have a /pseudo qualifier
    final FeatureKeyQualifierPredicate predicate =
      new FeatureKeyQualifierPredicate(Key.CDS, "pseudo", false);

    final FeatureVector non_embl_features = new FeatureVector();

    final FeatureEnumeration feature_enum = features();

    while(feature_enum.hasMoreFeatures()) 
    {
      final Feature current_feature = feature_enum.nextFeature();

      if(predicate.testPredicate(current_feature)) 
      {
        if(!current_feature.hasValidStartCodon()) 
          non_embl_features.add(current_feature);
      }
    }

    return non_embl_features;
  }

  /**
   *  Returns a vector containing those CDS features that have no /pseudo
   *  qualifier and do not have a valid stop codon.  The returned features may
   *  be illegal in EMBL submissions.
   **/
  public FeatureVector checkFeatureStopCodons() 
  {
    // get all the CDS features that do not have a /pseudo qualifier
    final FeatureKeyQualifierPredicate predicate =
      new FeatureKeyQualifierPredicate(Key.CDS, "pseudo", false);

    final FeatureVector non_embl_features = new FeatureVector();

    final FeatureEnumeration feature_enum = features();

    while(feature_enum.hasMoreFeatures()) 
    {
      final Feature current_feature = feature_enum.nextFeature();

      if(predicate.testPredicate(current_feature)) 
      {
        if(!current_feature.hasValidStopCodon()) 
          non_embl_features.add(current_feature);
      }
    }

    return non_embl_features;
  }

  /**
   *  Returns a vector containing those features that have the same key and
   *  location as one or more other features.  These features are not allowed
   *  in EMBL submissions.
   **/
  public FeatureVector checkForEMBLDuplicates() 
  {
    final FeatureVector non_embl_features = new FeatureVector();
    final FeatureEnumeration feature_enum = features();
    Feature last_feature = null;

    while(feature_enum.hasMoreFeatures()) 
    {
      final Feature current_feature = feature_enum.nextFeature();
      final Key current_key = current_feature.getKey();
      final Location current_location = current_feature.getLocation();

      if(last_feature != null &&
          last_feature.getKey().equals(current_key) &&
          last_feature.getLocation().equals(current_location)) 
      {
        if(! non_embl_features.contains(last_feature)) 
          non_embl_features.add(last_feature);
        
        non_embl_features.add(current_feature);
        continue;
      }
      last_feature = current_feature;
    }

    return non_embl_features;
  }

  /**
   *  Returns a vector containing those features that have the same key and
   *  location as one or more other features.  These features are not allowed
   *  in EMBL submissions.
   **/
  public FeatureVector checkForOverlappingCDSs() 
  {
    final FeatureVector cds_features = new FeatureVector();
    final FeatureEnumeration feature_enum = features();

    while(feature_enum.hasMoreFeatures()) 
    {
      final Feature current_feature = feature_enum.nextFeature();

      if(current_feature.isCDS()) 
        cds_features.add(current_feature);
    }

    final FeatureVector return_features = new FeatureVector();

    for(int i = 0 ; i + 1 < cds_features.size() ; ++i) 
    {
      final Feature this_feature = cds_features.elementAt(i);
      final Feature next_feature = cds_features.elementAt(i + 1);

      final int this_feature_end = this_feature.getRawLastBase();
      final int next_feature_start = next_feature.getRawFirstBase();

      if(this_feature_end >= next_feature_start)
        return_features.add(this_feature);
    }

    return return_features;
  }

  /**
   *  Returns a vector containing those features that have a feature that is
   *  missing a required qualifier.  These features are not allowed in EMBL
   *  submissions.
   **/
  public FeatureVector checkForMissingQualifiers() 
  {
    final FeatureVector non_embl_features = new FeatureVector();

    final FeatureEnumeration feature_enum = features();

    while(feature_enum.hasMoreFeatures()) 
    {
      final Feature current_feature = feature_enum.nextFeature();

      if(!current_feature.hasRequiredQualifiers()) 
        non_embl_features.add(current_feature);
    }

    return non_embl_features;
  }

  /**
   *  Adds the specified event listener to receive entry change events from
   *  this object.
   *  @param l the event change listener.
   **/
  public void addEntryChangeListener(EntryChangeListener l) 
  {
    entry_listener_list.addElement(l);
  }

  /**
   *  Removes the specified event listener so that it no longer receives
   *  entry change events from this object.
   *  @param l the event change listener.
   **/
  public void removeEntryChangeListener(EntryChangeListener l) 
  {
    entry_listener_list.removeElement(l);
  }

  /**
   *  Adds the specified event listener to receive feature change events from
   *  this object.
   *  @param l the event change listener.
   **/
  public void addFeatureChangeListener(FeatureChangeListener l)
  {
    feature_listener_list.addElement(l);
  }

  /**
   *  Removes the specified event listener so that it no longer receives
   *  feature change events from this object.
   *  @param l the event change listener.
   **/
  public void removeFeatureChangeListener(FeatureChangeListener l) 
  {
    feature_listener_list.removeElement(l);
  }

  /**
   *  Implementation of the FeatureChangeListener interface.  We need to
   *  listen to feature change events from the Features in this object.
   *  @param event The change event.
   **/
  public void featureChanged(FeatureChangeEvent event) 
  {
    // pass the action straight through
    fireAction(feature_listener_list, event);
  }

  /**
   *  Return the EntryInformation object of this Entry.
   **/
  public EntryInformation getEntryInformation() 
  {
    return getEMBLEntry().getEntryInformation();
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

      if(event instanceof EntryChangeEvent) 
      {
        final EntryChangeListener entry_change_listener =
         (EntryChangeListener) target;
        entry_change_listener.entryChanged((EntryChangeEvent) event);
      } 
      else
      {
        final FeatureChangeListener feature_change_listener =
         (FeatureChangeListener) target;
        feature_change_listener.featureChanged((FeatureChangeEvent) event);
      }
    }
  }

  /**
   *  Return the uk.ac.sanger.artemis.Feature object of the given embl.Feature object.  This
   *  method will create an appropriate uk.ac.sanger.artemis.Feature if none exists.
   **/
  private Feature
    getFeatureOf(uk.ac.sanger.artemis.io.Feature embl_feature) 
  {
    final Feature test_feature =(Feature) embl_feature.getUserData();
 
    if(test_feature == null) 
    {
      final Feature new_feature = new Feature(embl_feature);
      new_feature.setEntry(this);
      // hack to create a FeatureSegment for each exon, which has to be done
      // aftre setting the Entry of the Feature
      new_feature.getSegments();
      return new_feature;
    } 
    else 
      return test_feature;
  }

  /**
   *  Check that all features in the embl.Entry object are in range for the
   *  Bases object that was passed to the constructor.
   **/
  private void checkLocations() throws OutOfRangeException
  {
    final uk.ac.sanger.artemis.io.FeatureEnumeration feature_enumerator =
      getEMBLEntry().features();

    while(feature_enumerator.hasMoreFeatures()) 
    {
      final uk.ac.sanger.artemis.io.Feature feature =
        feature_enumerator.nextFeature();

      final Location location = feature.getLocation();
      final RangeVector ranges = location.getRanges();

      for(int i = 0 ; i < ranges.size() ; ++i) 
      {
        if(((Range)ranges.elementAt(i)).getEnd() > getBases().getLength())
        {
          if(!(getBases().getSequence() instanceof PartialSequence))
            throw new OutOfRangeException(location.toString());
        }
      }
    }
  }

  /**
   *  Create a uk.ac.sanger.artemis.Feature object for each embl.Feature in embl_entry.
   **/
  private void createDianaFeatures() 
  {
    // enumerating the features will call getFeatureOf() for each
    // embl.Feature, which will in turn create a uk.ac.sanger.artemis.Feature

    final FeatureEnumeration feature_enumerator = features();

    while(feature_enumerator.hasMoreFeatures())  
    {
      final Feature this_feature = feature_enumerator.nextFeature();
    }
  }

  /**
   *  Return the embl.Entry object that was passed to the constructor.
   **/
  public uk.ac.sanger.artemis.io.Entry getEMBLEntry() 
  {
    return embl_entry;
  }

  /**
   *  Return the embl.Entry object that was passed to the constructor.
   *  @param entry_information The EntryInformation object for the new Entry.
   *  @param new_type Should be one of EMBL_FORMAT, GENBANK_FORMAT or
   *    ANY_FORMAT.  If ANY_FORMAT then the Entry will be returned as is.
   *  @param force If true then invalid qualifiers and any features with
   *    invalid keys in the new Entry will be quietly thrown away.  "Invalid"
   *    means that the key/qualifier is not allowed to occur in an Entry of
   *    this type(probably determined by the EntryInformation object of this
   *    Entry).  If false an EntryInformationException will be thrown for
   *    invalid keys or qualifiers.
   *  @exception EntryInformationException Thrown if an Entry using the given
   *    EntryInformation object cannot contain the Key, Qualifier or
   *    Key/Qualifier combination of one of the features in the Document.
   **/
  private uk.ac.sanger.artemis.io.Entry
    getEMBLEntryAsNewType(final EntryInformation entry_information,
                           final int new_type, final boolean force)
      throws EntryInformationException 
  {
    return DocumentEntryFactory.makeDocumentEntry(entry_information,
                                                  embl_entry, new_type,
                                                  force);
  }

}
