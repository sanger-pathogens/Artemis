/* CorbaFeature.java
 *
 * created: Wed Feb 10 1999
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1998,1999,2000,2001  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/CorbaFeature.java,v 1.2 2005-11-28 16:46:38 tjc Exp $
 **/

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;

import nsdb.NucFeature;
import type.NoResult;

/**
 *  This is a subclass of Feature that can read itself from the EMBL CORBA
 *  server.
 *
 *  @author Kim Rutherford
 *  @version $Id: CorbaFeature.java,v 1.2 2005-11-28 16:46:38 tjc Exp $
 **/

public class CorbaFeature extends EMBLObject
    implements ComparableFeature {
  /**
   *  Create a new CorbaFeature object.
   *  @param feature_handle The corba handle of the feature
   *  @exception InvalidRelationException Thrown if this Feature cannot contain
   *    one of the qualifiers that was read from the server.
   **/
  public CorbaFeature (final NucFeature feature_handle)
      throws LocationParseException, InvalidKeyException, NoResult,
             InvalidRelationException {
    this.feature_handle = feature_handle;

    final String location_string =
      feature_handle.getLocation ().getLocationString ();

    this.location = new Location (location_string);
    this.key = new Key (feature_handle.getKey ());
    grabQualifiers ();
  }

  /**
   *  Set the value of this object.
   *  @param key The new feature key
   *  @param location The Location object for the new feature
   *  @param qualifiers The qualifiers for the new feature
   *  @exception InvalidRelationException Thrown if this Feature cannot contain
   *    any of the given Qualifier objects.
   **/
  public void set (Key key,
                   Location location,
                   QualifierVector qualifiers)
      throws InvalidRelationException, ReadOnlyException {
    throw new ReadOnlyException ();
  }

  /**
   *  Set the owning Entry of this Feature.  Other objects should call
   *  setCorbaEntry () to change the owner.
   *  @param entry The Entry that now owns this Feature.
   **/
  private void setEntry (CorbaEntry entry) {
    this.entry = entry;
  }

  /**
   *  Set the owning CorbaEntry of this Feature.
   *  @param entry The CorbaEntry that now owns this Feature.
   **/
  void setCorbaEntry (final CorbaEntry entry) {
    setEntry (entry);
  }

  /**
   *  Set the key field of this object.
   *  @param key The new feature key
   **/
  public void setKey (Key key)
      throws ReadOnlyException {
    throw new ReadOnlyException ();
  }

  public void setLocation (Location location, Entry entry)
      throws ReadOnlyException {
    throw new ReadOnlyException ();
  }

  /**
   *  Set the location of this object.
   *  @param location The Location object for the new feature
   **/
  public void setLocation (Location location)
      throws ReadOnlyException {
    throw new ReadOnlyException ();
  }

  /**
   *  Set the qualifiers of this object discarding all the current qualifiers.
   *  @param qualifiers The qualifiers for the new feature
   *  @exception InvalidRelationException Throw if this Feature cannot contain
   *    any of the given Qualifier objects.
   **/
  public void setQualifiers (QualifierVector qualifiers)
      throws InvalidRelationException, ReadOnlyException {
    throw new ReadOnlyException ();
  }

  /**
   *  Add the given Qualifier to this Feature.  If this Feature contains a
   *  Qualifier with the same name as the new Qualifier it will be replaced.
   *  @param qualifier The new qualifier to add.
   *  @exception InvalidRelationException Throw if this Feature cannot contain
   *    the given Qualifier.
   **/
  public void setQualifier (Qualifier qualifier)
      throws InvalidRelationException, ReadOnlyException {
    throw new ReadOnlyException ();
  }

  /**
   *  Remove the Qualifier with the given name.  If there is no Qualifier with
   *  that name then return immediately.
   *  @param name The qualifier name to look for.
   **/
  public void removeQualifierByName (final String name)
      throws ReadOnlyException {
    throw new ReadOnlyException ();
  }

  /**
   *  Add the values from the given qualifier to the Qualifier object with the
   *  same name in this Feature or otherwise add a copy of the argument.
   *  @param qualifier This object contians name and values to add.
   *  @param values The values to add to the Qualifier.
   *  @return The Qualifier that was changed or created.
   *  @exception InvalidRelationException Thrown if this Feature cannot contain
   *    the given Qualifier.
   **/
  public Qualifier addQualifierValues (Qualifier qualifier)
      throws InvalidRelationException, ReadOnlyException {
    throw new ReadOnlyException ();
  }

  /**
   *  Return the unique identifier of this StreamFeature.
   **/
  public long getNumericID () {
    return id;
  }

  /**
   *  Return the key of this feature, as passed to the constructor.
   **/
  public Key getKey () {
    return key;
  }

  /**
   *  Return the Location of this Feature, as passed to the constructor.
   **/
  public Location getLocation () {
    return location;
  }

  /**
   *  Return a QualifierVector object containing the qualifiers for this
   *  feature.  This method does not return a copy of the qualifier vector so
   *  changing the vector will change the feature.
   **/
  public QualifierVector getQualifiers () {
    return qualifiers;
  }

  /**
   *  Return the Qualifier in this Feature with the given name or null if
   *  there no such Qualifier.
   **/
  public Qualifier getQualifierByName (final String name) {
    return getQualifiers ().getQualifierByName (name);
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
   *  will be an EmblStreamFeature and will not be in any Entry, so it should
   *  be explicitly added with a call to Entry.add ().
   **/
  public Feature copy () {
    return new EmblStreamFeature (this);
  }

  /**
   *  Returns true if and only if this Feature can't be changed and can't be
   *  removed from it's entry.
   **/
  public boolean isReadOnly () {
     return true;
  }

  /**
   *  Read the qualifiers from the given NucFeature, create an embl.Qualifier
   *  for each one and put them in the qualifiers vector.
   **/
  private void grabQualifiers ()
      throws InvalidRelationException {
    qualifiers = new QualifierVector ();

    try {
      final nsdb.NucFeaturePackage.Qualifier [] corba_qualifier_handles =
        feature_handle.getQualifiers ();

      for (int i = 0 ; i < corba_qualifier_handles.length ; ++i) {
        final nsdb.NucFeaturePackage.Qualifier this_corba_qualifier =
          corba_qualifier_handles[i];

        final StringVector values =
          getQualifierStringValues (this_corba_qualifier);

        EntryInformation entry_info = getEntryInformation ();

        final String qualifier_name = this_corba_qualifier.name;

        if (!entry_info.isValidQualifier (getKey (),
                                          qualifier_name)) {
          final String message =
            getKey () + " cannot have " + qualifier_name +
            " as a qualifier";
          throw new InvalidRelationException (message, getKey (),
                                              new Qualifier (qualifier_name,
                                                             values));
        }

        if (values.size () == 0) {
          qualifiers.setQualifier (new Qualifier (this_corba_qualifier.name));
        } else {
          qualifiers.setQualifier (new Qualifier (this_corba_qualifier.name,
                                                  values));
        }
      }
    } catch (NoResult e) {
      System.out.println ("unexpected NoResult exception while reading " +
                          "from corba " + e);
    }
  }

  /**
   *  Return (from corba) the values of the given qualifier as String objects.
   **/
  private StringVector
    getQualifierStringValues (nsdb.NucFeaturePackage.Qualifier qualifier) {

    final nsdb.NucFeaturePackage.QualifierValue_u [] corba_values =
      qualifier.values;

    final StringVector return_array = new StringVector ();

    for (int value_index = 0 ;
         value_index < corba_values.length ;
         ++value_index) {

      final nsdb.NucFeaturePackage.QualifierValue_u corba_value =
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
   *  Return the EntryInformation object from the Entry of this feature (if
   *  there is an Entry) or SimpleEntryInformation.default_entry_information
   *  otherwise.
   **/
  private EntryInformation getEntryInformation () {
    if (getEntry () == null) {
      return SimpleEntryInformation.getDefaultEntryInformation ();
    } else {
      return getEntry ().getEntryInformation ();
    }
  }

  /**
   *  The NucFeature object that was passed to the constructor.
   **/
  private NucFeature feature_handle;

  /**
   *  The Location of this Feature, set by the constructor.
   **/
  private Location location;

  /**
   *  The Key of this Feature, set by the constructor.
   **/
  private Key key;

  /**
   *  The qualifiers of this Feature, set by the first call to the
   *  getQualifiers () method.
   **/
  private QualifierVector qualifiers = null;

  /**
   *  The CorbaEntry object that contains this Feature as passed to the
   *  constructor.
   **/
  private CorbaEntry entry;

  /**
   *  A unique identifier for this feature.
   **/
  private final long id = id_counter++;

  /**
   *  This is incremented each time a constructor is called.
   **/
  private static long id_counter = 0;
}
