
package uk.ac.sanger.artemis.io;

import java.util.*;

import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import org.biojava.utils.ChangeVetoException;
import org.biojava.bio.Annotation;
import org.biojava.bio.symbol.LocationTools;

public class BioJavaFeature extends EMBLObject implements ComparableFeature
{
  org.biojava.bio.seq.Feature bioJavaFeature;

  private BioJavaEntry           nativeEntry;
  private Key             nativeKey;
  private Location        nativeLocation;
  private QualifierVector nativeQualifiers;

  private static long     id_counter = 0;
  private final long      id = id_counter++;

  public BioJavaFeature(org.biojava.bio.seq.Feature feature,
                        BioJavaEntry nativeEntry)
  {
    this.bioJavaFeature = feature;
    this.nativeEntry = nativeEntry;
  }

//   /**
//    *  Create a new BioJavaFeature object with the same key, location and
//    *  qualifiers as the given feature.  The feature will be added to
//    *  nativeEntry 
//    *  @param feature The feature to copy.
//    **/
//   BioJavaFeature(Feature feature,
//                  BioJavaEntry nativeEntry)
//   {
//     this.bioJavaFeature = new org.biojava.bio.seq.Feature;
//     this.nativeEntry = nativeEntry;
//   }

  public void set (Key key, Location loc, QualifierVector qualifiers)
      throws EntryInformationException, ReadOnlyException, OutOfRangeException
  {
    setKey (key);
    setLocation (loc);
    setQualifiers (qualifiers);
  }

  public Key getKey ()
  {
    return new Key (bioJavaFeature.getType());
  }

  public long getNumericID()
  {
    return id;
  }

  public Feature copy()
  {
    return new BioJavaFeature(bioJavaFeature, nativeEntry);
  }

  public Location getLocation()
  {
    return makeNativeLocation ();
  }

  public Entry getEntry()
  {
    return nativeEntry;
  }

  
  /**
   *  
   **/
  org.biojava.bio.seq.Feature getBioJavaFeature () {
    return bioJavaFeature;
  }

  void setBioJavaEntry(final BioJavaEntry entry)
  {
    this.nativeEntry = entry;
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


  public Qualifier getQualifierByName(String name)
      throws InvalidRelationException
  {
    org.biojava.bio.Annotation annotation = bioJavaFeature.getAnnotation();

    if (!annotation.containsProperty (name)) {
      return null;
    }

    Object value = annotation.getProperty (name);

    StringVector sv = new StringVector();

    if (Collection.class.isInstance(value)) {
      for (Iterator vi = ((Collection) value).iterator(); vi.hasNext();) {
        sv.add(vi.next().toString());
      }
    } else {
      sv.add(value.toString());
    }

    return new Qualifier (name, sv);
  }

  public QualifierVector getQualifiers()
  {
    org.biojava.bio.Annotation annotation = bioJavaFeature.getAnnotation();

    QualifierVector qualifiers = new QualifierVector();

    List keys = new ArrayList();
    keys.addAll(annotation.keys());
    Collections.sort(keys);

    for (Iterator ki = keys.iterator() ; ki.hasNext() ; ) {
      Object   key = ki.next();
      Object value = annotation.getProperty(key);

      if (key.equals (org.biojava.bio.seq.Feature.PROPERTY_DATA_KEY)) {
        continue;
      }

      if (value instanceof Collection) {
        StringVector sv = new StringVector();

        for (Iterator vi = ((Collection) value).iterator(); vi.hasNext();) {
          sv.add (vi.next().toString());
        }

        qualifiers.addQualifierValues (new Qualifier (key.toString(), sv));
      } else {
        if (value instanceof Boolean) {
          qualifiers.setQualifier (new Qualifier (key.toString ()));
        } else {
          if (value instanceof String) {
            qualifiers.addQualifierValues (new Qualifier (key.toString(), value.toString ()));
          }
        }
      }
    }
    return qualifiers;
  }

  /**
   *  
   **/
  public void removeQualifierByName(String name)
      throws EntryInformationException, ReadOnlyException
  {
    if (getEntry () == null) {
      return;
    }

    // save the Entry because the call to remove () will set it to null
    final BioJavaEntry saved_entry = nativeEntry;
      
    // remove and then add the Feature because changing the Key may
    // change the position of the Feature in the FeatureTable eg. changing
    // CDS to CDS_motif can move the Feature because CDS is always sorted
    // before CDS_motif if they have the same start base
    saved_entry.removeFromTable (this);

    try {
      bioJavaFeature.getAnnotation ().removeProperty (name);
    } catch (IllegalArgumentException e) {
      throw new EntryInformationException ("cannot remove qualifier: " + name);
    } catch (ChangeVetoException e) {
      throw new ReadOnlyException ("cannot remove qualifier: " + name);
    } finally {
      saved_entry.setDirtyFlag ();

      saved_entry.addToTable (this);
    }
  }

  public void setKey(Key key)
      throws EntryInformationException, ReadOnlyException
  {
    if (key == null) {
      return;
    }

    if (getEntry () == null) {
      return;
    }

    // save the Entry because the call to remove () will set it to null
    final BioJavaEntry saved_entry = nativeEntry;
    
    // remove and then add the Feature because changing the Key may
    // change the position of the Feature in the FeatureTable eg. changing
    // CDS to CDS_motif can move the Feature because CDS is always sorted
    // before CDS_motif if they have the same start base
    saved_entry.removeFromTable (this);
    

    try {
      bioJavaFeature.setType (key.toString ());
    } catch (ChangeVetoException e) {
      throw new ReadOnlyException ("cannot set key");
    } finally {
      saved_entry.setDirtyFlag ();

      saved_entry.addToTable (this);
    }
  }

  public void setLocation(Location location)
      throws OutOfRangeException, ReadOnlyException
  {
    if (location == null || getEntry () == null) {
      return;
    }

    // save the Entry because the call to remove () will set it to null
    final BioJavaEntry saved_entry = nativeEntry;
      
    // remove and then add the Feature because changing the Key may
    // change the position of the Feature in the FeatureTable eg. changing
    // CDS to CDS_motif can move the Feature because CDS is always sorted
    // before CDS_motif if they have the same start base
    saved_entry.removeFromTable (this);

//    final boolean complementFlag = loc.isComplement ();

    final org.biojava.bio.symbol.Location bioJavaLocation =
      makeBioJavaLocation (location);

    try {
      bioJavaFeature.setLocation (bioJavaLocation);
    } catch (ChangeVetoException e) {
      throw new ReadOnlyException ("cannot set location");
    } finally {
      saved_entry.setDirtyFlag ();

      saved_entry.addToTable (this);
    }
  }

  static org.biojava.bio.symbol.Location
    makeBioJavaLocation (final Location artemisLocation)
  {
    final ArrayList al = new ArrayList ();

    final RangeVector ranges = artemisLocation.getRanges ();

    for (int i = 0 ; i < ranges.size () ; ++i) {
      final Range range = ranges.elementAt (i);
      al.add (LocationTools.makeLocation (range.getStart (), range.getEnd ()));
    }
      
    return LocationTools.union (al);
  }


  public void setQualifier (final Qualifier qualifier)
      throws EntryInformationException, ReadOnlyException
  {
    // save the Entry because the call to remove () will set it to null
    final BioJavaEntry saved_entry = nativeEntry;
      
    if (getEntry () == null) {
      return;
    }

    // remove and then add the Feature because changing the Key may
    // change the position of the Feature in the FeatureTable eg. changing
    // CDS to CDS_motif can move the Feature because CDS is always sorted
    // before CDS_motif if they have the same start base
    saved_entry.removeFromTable (this);

    try {
      setQualifierInternal (qualifier);
    } catch (IllegalArgumentException e) {
      throw new EntryInformationException ("cannot remove qualifier: " +
                                           qualifier.getName ());
    } catch (ChangeVetoException e) {
      throw new ReadOnlyException ("cannot remove qualifier: " +
                                   qualifier.getName ());
    } finally {
      saved_entry.setDirtyFlag ();

      saved_entry.addToTable (this);
    }
  }

  private void setQualifierInternal (final Qualifier qualifier)
      throws IllegalArgumentException, ChangeVetoException {
    final Annotation annotation = bioJavaFeature.getAnnotation ();

    if (qualifier.getValues () == null) {
      annotation.setProperty (qualifier.getName (), new Boolean (true));
    } else {
      if (qualifier.getValues ().size () == 1) {
        annotation.setProperty (qualifier.getName (),
                                qualifier.getValues ().elementAt (0));
      } else {
        annotation.setProperty (qualifier.getName (),
                                qualifier.getValues ().asCollection ());
      }
    }
  }

  /**
   *
   **/
  private void clearAnnotation ()
      throws ChangeVetoException {
    org.biojava.bio.Annotation annotation = bioJavaFeature.getAnnotation();

    List keys = new ArrayList();

    keys.addAll(annotation.keys());

    for (Iterator ki = keys.iterator() ; ki.hasNext() ; ) {
      Object key = ki.next();

      if (key.equals (org.biojava.bio.seq.Feature.PROPERTY_DATA_KEY)) {
        continue;
      }

      annotation.removeProperty (key);
    }
  }

  /**
   *
   **/
  public void setQualifiers(QualifierVector qualifiers)
      throws EntryInformationException, ReadOnlyException
  {
    if (qualifiers == null) {
      return;
    }

    if (getEntry () == null) {
      return;
    }

    // save the Entry because the call to remove () will set it to null
    final BioJavaEntry saved_entry = nativeEntry;
      
    // remove and then add the Feature because changing the Key may
    // change the position of the Feature in the FeatureTable eg. changing
    // CDS to CDS_motif can move the Feature because CDS is always sorted
    // before CDS_motif if they have the same start base
    saved_entry.removeFromTable (this);

    
    try {
      clearAnnotation ();
      for (int i = 0 ; i < qualifiers.size () ; ++i) {
        setQualifierInternal (qualifiers.elementAt (i));
      }
    } catch (ChangeVetoException e) {
      throw new ReadOnlyException ("cannot set qualifiers");
    } finally {
      saved_entry.setDirtyFlag ();

      saved_entry.addToTable (this);
    }
  }

  public boolean isReadOnly () {
    return false;
  }

  private Location makeNativeLocation()
  {
    org.biojava.bio.symbol.Location locs = bioJavaFeature.getLocation();

    boolean complement = false;

    if (org.biojava.bio.seq.StrandedFeature.class.isInstance(bioJavaFeature))
      if (((org.biojava.bio.seq.StrandedFeature) bioJavaFeature).getStrand().getValue() == -1)
        complement = true;

    RangeVector ranges = new RangeVector();
    for (Iterator li = locs.blockIterator(); li.hasNext();)
      {
        org.biojava.bio.symbol.Location thisLoc =
          (org.biojava.bio.symbol.Location) li.next();

        int min = (thisLoc.getMin() == Integer.MIN_VALUE) ?
          1                                             :
          thisLoc.getMin();

        int max = (thisLoc.getMax() == Integer.MAX_VALUE) ?
          getEntry().getSequence().length()             :
          thisLoc.getMax();

        try
          {
            if (min != max)
              {
                if (complement)
                  {
                    ranges.insertElementAt(new Range(min, max), 0);
                  }
                else
                  {
                    ranges.add(new Range(min, max));
                  }
              }
            else
              {
                if (complement)
                  {
                    ranges.insertElementAt(new Range(min), 0);
                  }
                else
                  {
                    ranges.add(new Range(min));
                  }
              }
          }
        catch (OutOfRangeException ore)
          {
            System.err.println("Error converting BioJava range: "
                               + ore.getMessage());
          }
      }
    return new Location(ranges, complement);
  }
}
