/* CorbaFeature.java
 *
 * created: Sun May 30 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/RWCorbaFeature.java,v 1.4 2005-11-28 16:46:38 tjc Exp $
 **/

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;

import nsdb.NucFeatureWriter;
import nsdb.NucFeature;
import nsdb.Datestamp;
import type.NoResult;
import nsdb.NucFeaturePackage.QualifierValue_u;

import java.util.Date;

/**
 *  This class implements the Feature interface by reading and writing to
 *  CORBA using a nsdb.NucFeatureWriter object.
 *
 *  @author Kim Rutherford
 *  @version $Id: RWCorbaFeature.java,v 1.4 2005-11-28 16:46:38 tjc Exp $
 **/

public class RWCorbaFeature extends EMBLObject implements DateStampFeature {
  /**
   *  Create a new CorbaFeature object.
   *  @param feature_handle The corba handle of the feature
   *  @param entry_information The EntryInformation object of the new Feature.
   *  @exception EntryInformationException Thrown if this
   *    Entry cannot contain the Key, Qualifier or Key/Qualifier combination of
   *    one of the features in the given Entry.
   **/
  public RWCorbaFeature (final EntryInformation entry_information,
                         final NucFeatureWriter feature_handle)
      throws NoResult, EntryInformationException {
    this.feature_handle = feature_handle;
    this.entry_information = entry_information;
  }

  /**
   *  Set the value of this object.
   *  @param key The new feature key
   *  @param location The Location object for the new feature
   *  @param qualifiers The qualifiers for the new feature
   *  @exception InvalidRelationException Thrown if this Feature cannot contain
   *    any of the given Qualifier objects.
   *  @exception ReadOnlyException Thrown if this Feature cannot be changed.
   *  @exception OutOfRangeException Thrown if any part of the location is out
   *    of range for this sequence.
   **/
  public void set (Key key,
                   Location location,
                   QualifierVector qualifiers)
      throws InvalidRelationException, OutOfRangeException, ReadOnlyException {
    setKey (key);
    setLocation (location);
    setQualifiers (qualifiers);
  }

  /**
   *  Set the value of this object.
   *  @param key The new feature key
   *  @param location The Location object for the new feature
   *  @param qualifiers The qualifiers for the new feature
   *  @exception InvalidRelationException Thrown if this Feature cannot contain
   *    any of the given Qualifier objects.
   *  @exception OutOfDate If the key has changed in the server since the time
   *    given by datestamp.
   *  @exception ReadOnlyException Thrown if this Feature cannot be changed.
   *  @exception OutOfRangeException Thrown if any part of the location is out
   *    of range for the sequence.
   **/
  public void set (final Date datestamp,
                   final Key key,
                   final Location location,
                   final QualifierVector qualifiers)
      throws InvalidRelationException, OutOfRangeException, ReadOnlyException,
             OutOfDateException {
    setKey (datestamp, key);
    
    // pass null because if the setKey() works the server datestamp will
    // change
    setLocation (null, location);
    setQualifiers (null, qualifiers);
  }

  
  /**
   *  Set the owning Entry of this Feature.  Other objects should call
   *  setCorbaEntry () to change the owner.
   *  @param entry The Entry that now owns this Feature.
   **/
  private void setEntry (RWCorbaEntry entry) {
    this.entry = entry;
  }

  /**
   *  Set the owning CorbaEntry of this Feature.
   *  @param entry The CorbaEntry that now owns this Feature.
   **/
  void setRWCorbaEntry (final RWCorbaEntry entry) {
    setEntry (entry);
  }

  /**
   *  Set the key field of this object.
   *  @param key The new feature key
   *  @exception InvalidRelationException Throw if this Feature cannot have
   *    the given Key and one of it's current qualifiers.
   *  @exception ReadOnlyException Thrown if this Feature cannot be changed.
   **/
  public void setKey (Key key)
      throws InvalidRelationException, ReadOnlyException {
    try {
      setKey (null, key);
    } catch (OutOfDateException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Set the key field of this object.
   *  @param key The new feature key
   *  @exception ReadOnlyException Thrown if this Feature cannot be changed.
   *  @exception OutOfDate If the key has changed in the server since the time
   *    given by datestamp.
   **/
  public void setKey (final Date datestamp, final Key key)
      throws InvalidRelationException, ReadOnlyException,
             OutOfDateException {
    if (key == null) {
      return;
    }

    try {
      feature_handle.setKey (makeServerStamp (datestamp), key.toString ());
      this.key = key;
    } catch (nsdb.ReadOnlyException e) {
      throw new ReadOnlyException ();
    } catch (nsdb.OutOfDate e) {
      throw new OutOfDateException ();
    } catch (nsdb.InvalidKey e) {
      throw new Error ("internal error - unexpected exception: " + e);
    } catch (type.InvalidRelation e) {
      throw new InvalidRelationException (e.reason, key);
    }
  }

  public void setLocation (Location location, Entry entry)
      throws OutOfRangeException, ReadOnlyException {
    setLocation(location);
  }

  /**
   *  Set the location of this object.
   *  @param location The Location object for the new feature
   *  @exception ReadOnlyException Thrown if this Feature cannot be changed.
   *  @exception OutOfRangeException Thrown if any part of the location is out
   *    of range for this sequence.
   **/
  public void setLocation (Location location)
      throws OutOfRangeException, ReadOnlyException {
    try {
      setLocation (null, location);
    } catch (OutOfDateException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Set the location of this object.
   *  @param location The Location object for the new feature
   *  @exception ReadOnlyException Thrown if this Feature cannot be changed.
   *  @exception OutOfDate If the key has changed in the server since the time
   *    given by datestamp.
   *  @exception OutOfRangeException Thrown if any part of the location is out
   *    of range for this sequence.
   **/
  public void setLocation (final Date datestamp,
                           final Location location)
      throws ReadOnlyException, OutOfDateException, OutOfRangeException {
    if (location == null) {
      return;
    }

    try {
      feature_handle.setLocation (makeServerStamp (datestamp),
                                  location.toString ());

      this.location = location;
      old_location  = location;
    } catch (nsdb.ReadOnlyException e) {
      throw new ReadOnlyException ();
    } catch (nsdb.OutOfDate e) {
      throw new OutOfDateException ();
    } catch (nsdb.LocationParse e) {
      throw new Error ("internal error - unexpected exception: " + e);
    } catch (type.IndexOutOfRange e) {
      throw new OutOfRangeException (location.toString ());
    }
  }

  /**
   *  Set the qualifiers of this object discarding all the current qualifiers.
   *  @param qualifiers The qualifiers for the new feature
   *  @exception InvalidRelationException Throw if this Feature cannot contain
   *    any of the given Qualifier objects.
   *  @exception ReadOnlyException Thrown if this Feature cannot be changed.
   **/
  public void setQualifiers (QualifierVector qualifiers)
      throws InvalidRelationException, ReadOnlyException {
    try {
      setQualifiers (null, qualifiers);
    } catch (OutOfDateException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Set the qualifiers of this object discarding all the current qualifiers.
   *  @param qualifiers The qualifiers for the new feature
   *  @exception InvalidRelationException Throw if this Feature cannot contain
   *    any of the given Qualifier objects.
   *  @exception ReadOnlyException Thrown if this Feature cannot be changed.
   *  @exception OutOfDateException If the key has changed in the server since
   *    the time given by datestamp.
   **/
  public void setQualifiers (final Date datestamp,
                             final QualifierVector qualifiers)
      throws InvalidRelationException, ReadOnlyException, OutOfDateException {
    if (qualifiers == null) {
      return;
    }

    try {
      final nsdb.NucFeaturePackage.Qualifier [] qualifier_list =
        new nsdb.NucFeaturePackage.Qualifier [qualifiers.size ()];

      for (int i = 0 ; i < qualifier_list.length ; ++i) {
        qualifier_list [i] = getStructFromQualifier((Qualifier)qualifiers.elementAt (i));
      }

      feature_handle.setQualifiers (makeServerStamp (datestamp),
                                    qualifier_list);

      this.qualifiers = qualifiers;
    } catch (nsdb.ReadOnlyException e) {
      throw new ReadOnlyException ();
    } catch (nsdb.OutOfDate e) {
      throw new OutOfDateException ();
    } catch (type.InvalidRelation e) {
      throw new InvalidRelationException (e.reason, getKey ());
    } catch (nsdb.InvalidQualifier e) {
      throw new Error ("internal error - unexpected exception: " + e);
    } catch (nsdb.QualifierParse e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Add the given Qualifier to this Feature.  If this Feature contains a
   *  Qualifier with the same name as the new Qualifier it will be replaced.
   *  @param qualifier The new qualifier to add.
   *  @exception InvalidRelationException Throw if this Feature cannot contain
   *    the given Qualifier.
   *  @exception ReadOnlyException Thrown if this Feature cannot be changed.
   *  @exception OutOfDateException If the key has changed in the server since
   *    the time given by datestamp.
   **/
  public void setQualifier (Qualifier qualifier)
      throws InvalidRelationException, ReadOnlyException {
    try {
      final nsdb.NucFeaturePackage.Qualifier new_qualifier=
        getStructFromQualifier (qualifier);

      feature_handle.setQualifier (new Datestamp (0), new_qualifier);
      qualifiers = null;
    } catch (nsdb.ReadOnlyException e) {
      throw new ReadOnlyException ();
    } catch (nsdb.OutOfDate e) {
      System.err.println ("setQualifier (): OutOfDate");
    } catch (type.InvalidRelation e) {
      throw new InvalidRelationException (e.reason, getKey(), qualifier);
    } catch (nsdb.InvalidQualifier e) {
      throw new Error ("internal error - unexpected exception: " + e);
    } catch (nsdb.QualifierParse e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Remove the Qualifier with the given name.  If there is no Qualifier with
   *  that name then return immediately.
   *  @param name The qualifier name to look for.
   *  @exception ReadOnlyException Thrown if this Feature cannot be changed.
   **/
  public void removeQualifierByName (final String name)
      throws ReadOnlyException {

  }

  /**
   *  Add the values from the given qualifier to the Qualifier object with the
   *  same name in this Feature or otherwise add a copy of the argument.
   *  @param qualifier This object contians name and values to add.
   *  @param values The values to add to the Qualifier.
   *  @return The Qualifier that was changed or created.
   *  @exception InvalidRelationException Thrown if this Feature cannot contain
   *    the given Qualifier.
   *  @exception ReadOnlyException Thrown if this Feature cannot be changed.
   **/
  public Qualifier addQualifierValues (Qualifier qualifier)
      throws InvalidRelationException, ReadOnlyException {
    // throw new ReadOnlyException ();
    return null;
  }

  /**
   *  Return the unique identifier of this StreamFeature.
   **/
  long getNumericID () {
    return id;
  }

  /**
   *  Return the key of this feature, as passed to the constructor.
   **/
  public Key getKey () {
    checkStamp ();

    if (key == null) {
      key = new Key (feature_handle.getKey ());
    }

    return key;
  }

  /**
   *  Return the Location of this Feature, as passed to the constructor.
   **/
  public Location getLocation () {
    checkStamp ();
    if (location == null) {
      return grabLocation ();
    } else {
      return location;
    }
  }

  /**
   *  Return a QualifierVector object containing the qualifiers for this
   *  feature.  This method does not return a copy of the qualifier vector so
   *  changing the vector will change the feature.
   **/
  public QualifierVector getQualifiers () {
    checkStamp ();

    if (qualifiers == null) {
      qualifiers = new QualifierVector ();

      try {
        final nsdb.NucFeaturePackage.Qualifier [] corba_qualifier_handles =
          feature_handle.getQualifiers ();

        for (int i = 0 ; i < corba_qualifier_handles.length ; ++i) {
          final nsdb.NucFeaturePackage.Qualifier this_corba_qualifier =
            corba_qualifier_handles[i];

          final Qualifier new_qualifier =
            getQualifierFromStruct (this_corba_qualifier);

          qualifiers.setQualifier (new_qualifier);
        }
      } catch (NoResult e) {

      }
    }

    return qualifiers;
  }

  /**
   *  Return the Qualifier in this Feature with the given name or null if
   *  there no such Qualifier.
   **/
  public Qualifier getQualifierByName (final String name)
      throws InvalidRelationException {
    checkStamp ();

    return getQualifiers ().getQualifierByName (name);
  }

  /**
   *  Given a Qualifier reference return a nsdb.NucFeaturePackage.Qualifier
   *  object.
   **/
  private nsdb.NucFeaturePackage.Qualifier
    getStructFromQualifier (final Qualifier qualifier) {

    final StringVector qualifier_values = qualifier.getValues ();

    if (qualifier_values == null) {
      return new nsdb.NucFeaturePackage.Qualifier (qualifier.getName (),
                                                   new QualifierValue_u [0]);
    } else {
      final QualifierValue_u [] return_values =
        new QualifierValue_u [qualifier_values.size ()];

      for (int i = 0 ; i < return_values.length ; ++i) {
        return_values[i] = new QualifierValue_u ();
        return_values[i].text ((String)qualifier_values.elementAt (i));
      }

      return new nsdb.NucFeaturePackage.Qualifier (qualifier.getName (),
                                                   return_values);
    }
  }

  /**
   *  Given a nsdb.NucFeaturePackage.Qualifier reference return a
   *  embl.Qualifier object.
   **/
  private Qualifier
    getQualifierFromStruct (final nsdb.NucFeaturePackage.Qualifier handle) {
    try {
      final StringVector values = getQualifierStringValues (handle);

      EntryInformation entry_info = getEntryInformation ();
      
      if (!entry_info.isValidQualifier (getKey (), handle.name)) {
        final String message =
          getKey () + " cannot have " + handle.name +
          " as a qualifier";
        throw new InvalidRelationException (message, getKey (),
                                            new Qualifier (handle.name));
      }

      if (values.size () == 0) {
        return new Qualifier (handle.name);
      } else {
        return new Qualifier (handle.name, values);
      }
    } catch (InvalidRelationException e) {
      System.out.println ("exception reading from corba: " + e);
    }

    return null;
  }

  /**
   *  Grab the location from the server and set the location field.
   **/
  Location grabLocation () {
    try {
      final String location_string =
        feature_handle.getLocation ().getLocationString ();

      if (location == null ||
          ! location.toString ().equals (location_string)) {
        if (old_location != null &&
            location_string.equals (old_location.toString ())) {
          location = old_location;
        } else {
          location = new Location (location_string);
        }
      } else {
        // location hasn't changed
      }
    } catch (LocationParseException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    } catch (NoResult e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }

    return location;
  }

  /**
   *  Return the first base of this feature.
   **/
  public int getFirstBase () {
    return getLocation ().getFirstBase ();
  }

  /**
   *  Return the last base of this feature.
   **/
  public int getLastBase () {
    return getLocation ().getLastBase ();
  }

  /**
   *  Return the Entry object that contains this Feature.
   **/
  public Entry getEntry () {
    return entry;
  }

  /**
   *  Return the reference of a new copy of this Feature.  The new Feature
   *  will not be in any Entry, so it should be explicitly added with a call
   *  to Entry.add ().
   **/
  public Feature copy () {
    final Feature new_feature = new EmblStreamFeature (this);

    return new_feature;
  }

  /**
   *  Return the datestamp of this feature - that is, the time it was last
   *  changed.
   **/
  public Date getDatestamp () {
    checkServerStamp ();
    return datestamp;
  }

  /**
   *  Return  the values of the given qualifier as String objects.
   **/
  private StringVector
    getQualifierStringValues (nsdb.NucFeaturePackage.Qualifier qualifier) {

    final QualifierValue_u [] corba_values =
      qualifier.values;

    final StringVector return_array = new StringVector ();

    for (int value_index = 0 ;
         value_index < corba_values.length ;
         ++value_index) {

      final QualifierValue_u corba_value =
        corba_values[value_index];

      /*final*/ String this_value_string;

      switch (corba_value.discriminator ()) {
      case NucFeature.string_qtc:
        return_array.add (corba_value.text ());
        break;
      case NucFeature.integer_qtc:
        return_array.add (String.valueOf (corba_value.integer ()));
        break;
      case NucFeature.real_qtc:
        return_array.add (String.valueOf (corba_value.real ()));
        break;
      case NucFeature.TranslationException_qtc:
        return_array.add (String.valueOf (corba_value.translation_exception ()));
        break;
      case NucFeature.CodonTranslation_qtc:
        return_array.add (String.valueOf (corba_value.codon_translation ()));
        break;
      case NucFeature.Anticodon_qtc:
        return_array.add (String.valueOf (corba_value.anti_codon ()));
        break;
      case NucFeature.SpliceConsensus_qtc:
        return_array.add (String.valueOf (corba_value.splice_consensus ()));
        break;
      case NucFeature.RepeatUnit_qtc:
        return_array.add (String.valueOf (corba_value.repeat_unit ()));
        break;
      case NucFeature.DbXref_qtc:
        return_array.add (String.valueOf (corba_value.db_xref ()));
        break;
      }
    }

    return return_array;
  }

  /**
   *  Return the NucFeatureWriter reference that was passed to the constructor.
   **/
  NucFeatureWriter getNucFeatureWriter () {
    return feature_handle;
  }

  /**
   *  Used by checkStamp () - the server won't be queried if the last call to
   *  checkStamp () (and server query) happened less then 10 seconds ago.
   **/
  private Date last_query_time = null;

  /**
   *  Get the Datestamp from the corba object of this feature, if the
   *  datestamp hasn't been retrieved in the last n seconds.  If it is newer
   *  than the current datestamp (stored in the datestamp field of this
   *  Feature), then all cached information will be forgotten
   **/
  private void checkStamp () {
    final java.util.Calendar calendar = java.util.Calendar.getInstance ();

    final java.util.Date current_time = calendar.getTime ();

    if (last_query_time == null ||
        current_time.getTime () - last_query_time.getTime () > 10000) {

      checkServerStamp ();

      last_query_time = current_time;
    }
  }

  /**
   *  Read the datestamp from the feature in the server.  If it is newer
   *  than the current datestamp (stored in the datestamp field of this
   *  Feature), then all cached information will be forgotten
   **/
  private void checkServerStamp () {
    final Datestamp server_datestamp = feature_handle.getDatestamp ();

    if (datestamp == null ||
        server_datestamp.value > (long)datestamp.getTime () / 1000) {

      datestamp = new Date ((long)server_datestamp.value * 1000);

      key = null;
      if (location != null) {
        old_location = location;
      }
      location = null;
      qualifiers = null;
    }
  }

  /**
   *  Make a Datestamp object from the given Date object.
   **/
  private static Datestamp makeServerStamp (final Date datestamp) {
    final Datestamp server_datestamp;

    if (datestamp == null) {
      server_datestamp = new Datestamp (0);
    } else {
      server_datestamp = new Datestamp ((int) (datestamp.getTime () / 1000));
    }

    return server_datestamp;
  }

  /**
   *  Return the EntryInformation object that was passed to the constructor.
   **/
  public EntryInformation getEntryInformation () {
    return entry_information; 
  }

  /**
   *  Returns true if and only if this Feature can't be changed and can't be
   *  removed from it's entry.
   **/
  public boolean isReadOnly () {
     return false;
  }

  /**
   *  The EntryInformation object that was passed to the constructor.
   **/
  final private EntryInformation entry_information;

  /**
   *  The NucFeature object that was passed to the constructor.
   **/
  private NucFeatureWriter feature_handle;

  /**
   *  The Key of this Feature (set by getKey () and checkStamp ())
   **/
  private Key key;

  /**
   *  The Location of this Feature (set by getLocation () and checkStamp ())
   **/
  private Location location;

  /**
   *  The Location of this Feature (set by getLocation () and checkStamp ())
   **/
  private QualifierVector qualifiers;

  /**
   *  The CorbaEntry object that contains this Feature as passed to the
   *  constructor.
   **/
  private RWCorbaEntry entry;

  /**
   *  A unique identifier for this feature.
   **/
  private final long id = id_counter++;

  /**
   *  This is incremented each time a constructor is called.
   **/
  private static long id_counter = 0;

  /**
   *  This is the "date stamp" of the feature. ie. the time it changed (as far
   *  as we know).  It will be checked (by calling
   *  feature_handle.getDatestamp()) in the get methods and will be passed to
   *  the set methods.
   **/
  private Date datestamp = null;

  /**
   *  The Location of this Feature before the last call to checkServerStamp ().
   **/
  private Location old_location = null;
}
