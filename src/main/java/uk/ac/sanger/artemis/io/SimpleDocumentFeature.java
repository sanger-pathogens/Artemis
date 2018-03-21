/* SimpleDocumentFeature.java
 *
 * created: Thu Feb 17 2000
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2000  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/SimpleDocumentFeature.java,v 1.4 2009-06-10 09:47:40 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;

/**
 *  SimpleDocumentFeature class
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: SimpleDocumentFeature.java,v 1.4 2009-06-10 09:47:40 tjc Exp $
 **/

abstract class SimpleDocumentFeature extends LineGroup
    implements Feature {
  /**
   *  Create a new DocumentFeature with the given DocumentEntry as it's
   *  owner.  The Entry can be null if this Feature doesn't have an owner
   *  (yet).
   **/
  public SimpleDocumentFeature (final DocumentEntry entry) {
    this.entry = entry;
  }

  /**
   *  Return the Entry object that contains this Feature.
   **/
  public Entry getEntry () {
    return entry;
  }

  /**
   *  If the entry is not null set its dirty flag.
   **/
  protected void setDirtyFlag () {
    if (getEntry () != null) {
      getDocumentEntry ().setDirtyFlag ();
    }
  }


  /**
   *  Return the DocumentEntry that contains this DocumentFeature.  A
   *  DocumentFeature will always be owned by a DocumentEntry.
   **/
  public DocumentEntry getDocumentEntry () {
    return entry;
  }

  /**
   *  Returns true if and only if this Feature can't be changed or can't be
   *  removed from it's entry.
   **/
  public boolean isReadOnly () {
    if (getEntry () != null && getEntry ().isReadOnly ()) {
      return true;
    } else {
      return false;
    }
  }

  /**
   *  Set the owning Entry of this Feature.  Other objects should call
   *  setDocumentEntry () to change the owner.
   *  @param entry The Entry that now owns this Feature.
   **/
  private void setEntry (final DocumentEntry entry) {
    this.entry = entry;
  }

  /**
   *  Set the owning DocumentEntry of this Feature.
   *  @param entry The Entry that now owns this Feature.
   **/
  protected void setDocumentEntry (final DocumentEntry entry)
      throws ReadOnlyException {
    setEntry (entry);
  }

  /**
   *  This is incremented each time a constructor is called.
   **/
  private static long id_counter = 0;

  /**
   *  Set the value of this object.
   *  @param key The new feature key - can be null if location or qualifiers
   *    is not null.
   *  @param location The Location object for the new feature - can be null if
   *    key or qualifiers is not null.
   *  @param qualifiers The qualifiers for the new feature - can be null if
   *    key or location is not null.
   *  @exception EntryInformationException Thrown if this Feature cannot
   *    contain the given Key, Qualifier or Key/Qualifier combination.
   **/
  public void set (final Key key,
                   final Location location,
                   final QualifierVector qualifiers)
      throws EntryInformationException, OutOfRangeException,
             ReadOnlyException {

    if (getEntry () != null && getEntry ().isReadOnly ()) {
      throw new ReadOnlyException ();
    }

    final Key key_to_test;

    if (key == null) {
      key_to_test = this.key;
    } else {
      key_to_test = key;
    }

    final QualifierVector qualifiers_to_test;

    if (qualifiers == null) {
      qualifiers_to_test = this.qualifiers;
    } else {
      qualifiers_to_test = qualifiers;
    }
    
    // save the Entry because the call to remove () will set it to null
    final Entry saved_entry = getEntry ();
      
    try {
      if (saved_entry != null) {
        // remove and then add the Feature because changing the Key may
        // change the position of the Feature in the FeatureTable eg. changing
        // CDS to CDS_motif can move the Feature because CDS is always sorted
        // before CDS_motif if they have the same start base
        saved_entry.remove (this);
      }
      
      if (key != null) {
        this.key = key;
      }
      if (qualifiers != null) {
        this.qualifiers = qualifiers;
      }
      if (location != null) {
        setLocation (location);
      }

      setDirtyFlag ();

    } finally {
      if (saved_entry != null) {
        try {
          saved_entry.add (this);
        } catch (EntryInformationException e) {
          saved_entry.forcedAdd (this);
          
          throw new Error ("internal error (feature key change may have " +
                           "been lost) - unexpected exception: " + e);
        }
      }
    }
  }

  /**
   *  Set the value of this object.
   *  @param entry The Entry that now owns this Feature.
   *  @param key The new feature key
   *  @param location The Location object for the new feature
   *  @param qualifiers The qualifiers for the new feature
   *  @exception EntryInformationException Thrown if this Feature cannot
   *    contain the given Key, Qualifier or Key/Qualifier combination.
   **/
  private void set (final DocumentEntry entry,
                    final Key key,
                    final Location location,
                    final QualifierVector qualifiers)
      throws EntryInformationException, ReadOnlyException {
    setDocumentEntry (entry);
    setKey (key);

    // throw away the cache
    first_base = -1;
    last_base = -1;
    this.location = location;

    setQualifiers (qualifiers);

    setDirtyFlag ();
  }

  /**
   *  Set the key field of this object.
   *  @param key The new feature key
   *  @exception EntryInformationException Thrown if this Feature cannot
   *    contain the given Key, Qualifier or Key/Qualifier combination.
   **/
  public void setKey (final Key key)
      throws EntryInformationException, ReadOnlyException {
    if (getEntry () != null && getEntry ().isReadOnly ()) {
      throw new ReadOnlyException ();
    }

    for (int i = 0 ; i < qualifiers.size () ; ++i) {
      final String this_name = ((Qualifier)qualifiers.elementAt(i)).getName();

      if (!getEntryInformation ().isValidQualifier (key, this_name)) {
        final String message =
          key + " cannot have /" + this_name +
          " as a qualifier";
        throw new InvalidRelationException(message, key,
                                           (Qualifier)qualifiers.elementAt(i));
      }
    }

    // save the Entry because the call to remove () will set it to null
    final Entry saved_entry = getEntry ();
      
    try {
      if (saved_entry != null) {
        // remove and then add the Feature because changing the Key may
        // change the position of the Feature in the FeatureTable eg. changing
        // CDS to CDS_motif can move the Feature because CDS is always sorted
        // before CDS_motif if they have the same start base
        saved_entry.remove (this);
      }
      
      this.key = key;
      
      setDirtyFlag ();
      
    } finally {
      if (saved_entry != null) {
        try {
          saved_entry.add (this);
        } catch (EntryInformationException e) {
          saved_entry.forcedAdd (this);
          
          throw new Error ("internal error (feature key change may have " +
                           "been lost) - unexpected exception: " + e);
        }
      }
    }
  }

  public void setLocation (final Location location)
      throws ReadOnlyException, OutOfRangeException {
    setLocation(location, null);
  }

  /**
   *  Set the location of this object.
   *  @param location The Location object for the new feature
   **/
  public void setLocation (final Location location, Entry saved_entry)
      throws ReadOnlyException, OutOfRangeException {

    if (getEntry () != null && getEntry ().isReadOnly ()) {
      throw new ReadOnlyException ();
    }

    // save the Entry because the call to remove () will set it to null
    if(saved_entry == null)
      saved_entry = getEntry ();
      
    try {
      if (saved_entry != null) {
        // remove and then add the Feature because changing the Location may
        // change the position of the Feature in the FeatureTable
        saved_entry.remove (this);
      }
      
      this.location = location;
      
      setDirtyFlag ();
      
      // throw away the cache
      first_base = -1;
      last_base = -1;

    } finally {

      if (saved_entry != null) {
        try {
          saved_entry.add (this);
        } catch (EntryInformationException e) {
          saved_entry.forcedAdd (this);
          
          first_base = -1;
          last_base = -1;

          throw new Error ("internal error (some feature qualifiers may " +
                           "have been lost) - unexpected exception: " + e);
        }
      }
    }
  }
  

  /**
   *  Set the qualifiers of this object.
   *  @param qualifiers The new qualifiers for the feature.  null means remove
   *    all qualifiers.
   *  @exception EntryInformationException Thrown if this Feature cannot
   *    contain the given Key, Qualifier or Key/Qualifier combination.
   **/
  public void setQualifiers (final QualifierVector qualifiers)
      throws EntryInformationException, ReadOnlyException {

    if (getEntry () != null && getEntry ().isReadOnly ()) {
      throw new ReadOnlyException ();
    }

    if (qualifiers == null) {
      this.qualifiers = new QualifierVector ();
    } else {
      for (int i = 0 ; i < qualifiers.size () ; ++i) {

        final String this_name = ((Qualifier)qualifiers.elementAt(i)).getName();

        if (!getEntryInformation ().isValidQualifier (getKey (), this_name)) {
          final String message =
            getKey () + " cannot have /" + this_name + " as a qualifier";
          throw new InvalidRelationException(message, getKey(),
                                             (Qualifier)qualifiers.elementAt(i));
        }
      }
      this.qualifiers = qualifiers.copy ();
    }

    setDirtyFlag ();
  }

  /**
   *  Add the given Qualifier to this Feature.  If this Feature contains a
   *  Qualifier with the same name as the new Qualifier it will be replaced.
   *  @param qualifier The new qualifier to add.
   *  @exception EntryInformationException Thrown if this Feature cannot
   *    contain the given Key, Qualifier or Key/Qualifier combination.
   **/
  public void setQualifier (final Qualifier qualifier)
      throws EntryInformationException, ReadOnlyException {
    if (getEntry () != null && getEntry ().isReadOnly ()) {
      throw new ReadOnlyException ();
    }

    getQualifiers ().setQualifier (qualifier);

    setDirtyFlag ();
  }

  /**
   *  Remove the Qualifier with the given name.  If there is no Qualifier with
   *  that name then return immediately.
   *  @exception EntryInformationException Thrown if a required qualifier is
   *    removed.
   *  @param name The qualifier name to look for.
   **/
  public void removeQualifierByName (final String name)
      throws EntryInformationException, ReadOnlyException {
    if (getEntry () != null && getEntry ().isReadOnly ()) {
      throw new ReadOnlyException ();
    }

    getQualifiers ().removeQualifierByName (name);

    setDirtyFlag ();
  }

  /**
   *  Add the values from the given qualifier to the Qualifier object with the
   *  same name in this Feature or otherwise add a copy of the argument.
   *  @param qualifier This object contains name and values to add.
   *  @return The Qualifier that was changed or created.
   *  @exception EntryInformationException Thrown if this Feature cannot
   *    contain the given Key, Qualifier or Key/Qualifier combination.
   **/
  public Qualifier addQualifierValues (final Qualifier qualifier)
      throws EntryInformationException, ReadOnlyException {
    if (getEntry () != null && getEntry ().isReadOnly ()) {
      throw new ReadOnlyException ();
    }

    setDirtyFlag ();

    return getQualifiers ().addQualifierValues (qualifier);
  }

  /**
   *  Return the unique identifier of this PublicDBStreamFeature.
   **/
  public long getNumericID () {
    return id;
  }

  /**
   *  Return the Key of this Feature, as passed to the constructor.
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
   *  Return the first base of this feature.
   **/
  public int getFirstBase () {
    if (first_base == -1) {
      first_base = getLocation ().getFirstBase ();
    }
    return first_base;
  }

  /**
   *  Return the last base of this feature.
   **/
  public int getLastBase () {
    if (last_base == -1) {
      last_base = getLocation ().getLastBase ();
    }
    return last_base;
  }

  /**
   *  Return a vector containing the qualifiers for this feature.  This method
   *  does not return a copy of the qualifier vector so changing the vector
   *  will change the feature.
   **/
  public QualifierVector getQualifiers () {
    return qualifiers;
  }

  /**
   *  Return the Qualifier in this Feature with the given name or null if
   *  there no such Qualifier.
   **/
  public Qualifier getQualifierByName (String name) {
    return getQualifiers ().getQualifierByName (name);
  }


  /**
   *  Return the EntryInformation object from the Entry of this feature (if
   *  there is an Entry) or SimpleEntryInformation.default_entry_information
   *  otherwise.
   **/
  public EntryInformation getEntryInformation () {
    if (getEntry () == null) {
      return SimpleEntryInformation.getDefaultEntryInformation ();
    } else {
      return getEntry ().getEntryInformation ();
    }
  }

  /**
   *  Return a String containing this Feature written in its native format.
   **/
//   public String toString () {
//     final StringWriter string_writer = new StringWriter ();

//     try {
//       writeToStream (string_writer);
//     } catch (IOException e) {
//       throw new Error ("internal error - unexpected exception: " + e);
//     }

//     return string_writer.toString ();
//   }

  /**
   *  Return the reference of a new copy of this Feature.
   **/
  public abstract Feature copy ();

  /**
   *  The DocumentEntry object that contains this Feature as passed to the
   *  constructor.
   **/
  private DocumentEntry entry;

  /**
   *  The is the key that was passed to the constructor
   **/
  private Key key;

  /**
   *  The is the Location object that was passed to the constructor
   **/
  private Location location;

  /**
   *  The is the vector of qualifiers that was passed to the constructor
   **/
  private QualifierVector qualifiers = new QualifierVector ();

  /**
   *  The first base of the location of this feature.  This is a cache.
   **/
  private int first_base = -1;

  /**
   *  The last base of the location of this feature.  This is a cache.
   **/
  private int last_base = -1;

  /**
   *  A unique identifier for this feature.
   **/
  private final long id = id_counter++;
}
