
package uk.ac.sanger.artemis.io;

import java.io.*;
import java.util.Iterator;

import org.biojava.bio.BioException;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.seq.io.EmblProcessor;
import org.biojava.bio.seq.io.SequenceFormat;
import org.biojava.bio.seq.io.EmblLikeFormat;
// import org.biojava.bio.seq.io.GAMEFormat;
// import org.biojava.bio.seq.io.BSMLFormat;
import org.biojava.bio.seq.io.SequenceBuilderFactory;
import org.biojava.bio.seq.io.SmartSequenceBuilder;
import org.biojava.bio.seq.io.StreamReader;

import uk.ac.sanger.artemis.util.*;

public class BioJavaEntry implements Entry
{

  // tjc uk.ac.sanger.artemis.components/BioJavaEntrySource.java barfed 
  //     so I added this constructor
  public BioJavaEntry(final Document document,
                      final SequenceFormat sequenceFormat)
      throws IOException
  {
    this(null,document,sequenceFormat);
  }

  public BioJavaEntry(final EntryInformation entryInformation,
                      final Document document,
                      final SequenceFormat sequenceFormat)
      throws IOException
  {
    this.entryInformation = entryInformation;
    this.sequenceFormat = sequenceFormat;
    this.document = document;
    featTable = new StreamFeatureTable ();

    BufferedReader reader =
      new BufferedReader (document.getLinePushBackReader ());

    try {
      SequenceBuilderFactory sFact =
        new EmblProcessor.Factory(SmartSequenceBuilder.FACTORY);
      Alphabet alpha = DNATools.getDNA();
      SymbolTokenization rParser = alpha.getTokenization("token");
      SequenceFormat eFormat = new EmblLikeFormat();
      SequenceIterator seqIterator =
        new StreamReader(reader, eFormat, rParser, sFact);

      bioJavaSequence = seqIterator.nextSequence();
    } catch (BioException be) {
      be.printStackTrace ();
      throw new IOException("Error reading BioJava sequence: " + be);
    }

    setArtemisFeatures(bioJavaSequence);
    setArtemisSequence(bioJavaSequence);
  }

  public BioJavaEntry (final org.biojava.bio.seq.Sequence sequence) {
    setArtemisFeatures(bioJavaSequence);
    setArtemisSequence(bioJavaSequence);
  }

  public BioJavaEntry (final Entry old_entry) {
    final Sequence old_sequence = old_entry.getSequence ();
    final String old_sequence_str =
      old_sequence.getSubSequence (1, old_sequence.length ());
     
    try {
      bioJavaSequence =
        org.biojava.bio.seq.DNATools.createDNASequence (old_sequence_str, "dna");

      setArtemisSequence (bioJavaSequence);
    } catch (org.biojava.bio.symbol.IllegalSymbolException e) {
      // XXX
      throw new Error ("internal error - unexpected exception: " + e);
    }

    final FeatureEnumeration feature_enum = old_entry.features ();

    while (feature_enum.hasMoreFeatures ()) {
      final Feature old_feature = feature_enum.nextFeature ();

      try {
        createFeature (old_feature.getKey (),
                       old_feature.getLocation (),
                       old_feature.getQualifiers ());
      } catch (ReadOnlyException e) {
        throw new Error ("internal error - unexpected exception: " + e);
      } catch (OutOfRangeException e) {
        throw new Error ("internal error - unexpected exception: " + e);
      }
    }
  }


  public String getHeaderText()
  {
    return null;
  }

  public EntryInformation getEntryInformation()
  {
    return entryInformation;
  }

  /**
   *  Return the File reference that was passed to the constructor or null if
   *  none was passed.
   **/
  public Document getDocument () {
    return document;
  }

  /**
   *  Set the document to use when saving this DocumentEntry.
   **/
  public void setDocument (final Document document) {
    this.document = document;
  }

  public void save () throws IOException {
    save (getDocument ());
  }

  public void save (Document document) throws IOException {
    final PrintStream printStream =
      new PrintStream (document.getOutputStream ());

    sequenceFormat.writeSequence (bioJavaSequence, printStream);
  }

  public boolean hasUnsavedChanges () {
    return lastChangeTime != null;
  }

  public boolean isReadOnly () {
    return false;
  }

  /**
   *  Set the header of this Entry to be the given text.
   *  @return true if and only if the header was successfully set.  Not all
   *    Entry objects can change their header, so it is up to the calling
   *    function to check the return value.
   *  @exception IOException thrown if there is a problem reading the header
   *    from the String - most likely ReadFormatException.
   **/
  public boolean setHeaderText (final String newHeader)
      throws IOException {
    return false;
  }

  /**
   *  Set the name of this Entry - if possible (the return value will let the
   *  caller know).
   *  @return true if and only if the name was successfully set.  Not all
   *    Entry objects can change their name, so it is up to the calling
   *    function to check the return value.
   **/
  public boolean setName (final String name) {
    return false;
  }

  /**
   *  Create a new Feature object of an appropriate type in this Entry.
   *  @param key The new feature key
   *  @param location The Location object for the new feature
   *  @param qualifiers The qualifiers for the new feature (can be null if
   *    there are no qualifiers).
   *  @exception InvalidRelationException Thrown if this Feature cannot contain
   *    the given Qualifier.
   *  @exception EntryInformationException Thrown if a Feature in this Entry
   *    cannot contain the given Key, Qualifier or Key/Qualifier combination.
   **/
  public Feature createFeature (Key key,
                                Location location,
                                QualifierVector qualifiers)
      throws ReadOnlyException, OutOfRangeException {

    final org.biojava.bio.seq.StrandedFeature.Template template =
      new org.biojava.bio.seq.StrandedFeature.Template ();

    if (location.isComplement ()) {
      template.strand = org.biojava.bio.seq.StrandedFeature.NEGATIVE;
    } else {
      template.strand = org.biojava.bio.seq.StrandedFeature.POSITIVE;
    }

    template.annotation = new org.biojava.bio.SimpleAnnotation ();

    template.type = key.toString ();

    template.location = BioJavaFeature.makeBioJavaLocation (location);

    try {
      final org.biojava.bio.seq.Feature bioJavaFeature =
        bioJavaSequence.createFeature (template);

      final BioJavaFeature newFeature =
        new BioJavaFeature (bioJavaFeature, this);

      featTable.add (newFeature);

      setDirtyFlag ();

      return newFeature;
    } catch (org.biojava.utils.ChangeVetoException e) {
      throw new ReadOnlyException ("feature cannot be created");
    } catch (org.biojava.bio.BioException e) {
      // XXX - createFeature () should throw BioException
      throw new ReadOnlyException ("BioJava error: " + e);
    }
  }

  /**
   *  Add the given Feature to this Entry.  If the Feature is already in an
   *  Entry then Entry.remove () should be called on that Entry before calling
   *  Entry.add ().  An Error will be thrown otherwise.
   *  @exception ReadOnlyException If this entry is read only.
   *  @exception EntryInformationException Thrown if this Entry
   *    cannot contain the Key, Qualifier or Key/Qualifier combination of the
   *    given Feature.
   *  @return A reference that was passed to add (), if that Feature can be
   *    stored directly in this Entry, otherwise returns a reference to a new
   *    Feature, that is a copy of the argument.  The argument reference
   *    should not be used after the call to add (), unless the return
   *    reference happens to be the same as the argument.
   **/
  public Feature add (Feature feature)
      throws EntryInformationException, ReadOnlyException {
    return forcedAdd (feature);
  }

  /**
   *  Add the given Feature to this Entry.  If the Feature is already in an
   *  Entry then Entry.remove () should be called on that Entry before calling
   *  Entry.add ().  An Error will be thrown otherwise.  Invalid qualifiers
   *  will be quietly thrown away.  Features with invalid keys will not be
   *  added (and null will be returned).  "Invalid" means that the
   *  key/qualifier is non allowed to occur in an Entry of this type (probably
   *  determined by the EntryInformation object of this Entry).
   *  @exception ReadOnlyException If this entry is read only.
   *  @exception EntryInformationException Thrown if this Entry
   *    cannot contain the Key, Qualifier or Key/Qualifier combination of the
   *    given Feature
   *  @return A reference that was passed to add (), if that Feature can be
   *    stored directly in this Entry, otherwise returns a reference to a new
   *    Feature, that is a copy of the argument.  The argument reference
   *    should not be used after the call to add (), unless the return
   *    reference happens to be the same as the argument.  Returns null if and
   *    only if the new Feature has a key that is invalid for this Entry.
   **/
  public Feature forcedAdd (Feature feature)
      throws ReadOnlyException {
    if (feature.getEntry () != null) {
      throw new Error ("internal error - a feature must have one owner");
    }

    try {
      return createFeature (feature.getKey (),
                            feature.getLocation (),
                            feature.getQualifiers ());
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Remove the given Feature from this Entry.
   *  @return true if and only if the Feature was in this Entry.
   *  @exception ReadOnlyException If this entry is read only.
   **/
  public boolean remove (Feature feature)
      throws ReadOnlyException {
    final boolean reallyRemoved;

    if (featTable.contains (feature)) {
      reallyRemoved = true;
    } else {
      reallyRemoved = false;
    }

    featTable.remove (feature);

    if (reallyRemoved) {
      ((BioJavaFeature)feature).setBioJavaEntry (null);

      try
      {
        bioJavaSequence.removeFeature ((org.biojava.bio.seq.Feature) ((BioJavaFeature)feature).getBioJavaFeature ());
      }
      catch(org.biojava.utils.ChangeVetoException e) 
      {
        throw new ReadOnlyException("read only - feature cannot be removed");
      }
      catch(org.biojava.bio.BioException e)
      {
        //A nestable biological exception.
      }

      setDirtyFlag ();

      return true;
    } else {
      return false;
    }
  }

  /**
   *  Set lastChangeTime so that hasUnsavedChanges () will return true.
   *  lastChangeTime can be read with getLastChangeTime ().
   **/
  public void setDirtyFlag () {
    if (in_constructor) {
      // features are being added, but the entry hasn't really changed
    } else {
      final java.util.Calendar calendar = java.util.Calendar.getInstance ();

      lastChangeTime = calendar.getTime ();
    }
  }

  /**
   *  Return the Date when this Entry last changed or null if this Entry
   *  hasn't changed since the last save.
   **/
  public java.util.Date getLastChangeTime () {
    return lastChangeTime;
  }

  void removeFromTable (final BioJavaFeature feature) {
    featTable.remove (feature);
  }

  void addToTable (final BioJavaFeature feature) {
    featTable.add (feature);
  }

  public String getName()
  {
    return bioJavaSequence.getName();
  }

  public int getFeatureCount()
  {
    return featTable.getFeatureCount();
  }

  public Feature getFeatureAtIndex(final int index)
  {
    return featTable.getFeatureAtIndex (index);
  }

  public int indexOf(final Feature feature)
  {
    return featTable.indexOf (feature);
  }

  public boolean contains(final Feature feature)
  {
    return featTable.contains(feature);
  }

  public FeatureEnumeration features()
  {
    return featTable.features();
  }

  public FeatureVector getFeaturesInRange(final Range range)
      throws OutOfRangeException
  {
    return featTable.getFeaturesInRange(range);
  }

  public FeatureVector getAllFeatures()
  {
    return featTable.getAllFeatures ();
  }

  public Sequence getSequence()
  {
    return artemisSequence;
  }

  private void setArtemisFeatures(final FeatureHolder holder)
  {
    for (Iterator fi = holder.features(); fi.hasNext();) {
      org.biojava.bio.seq.Feature f =
        (org.biojava.bio.seq.Feature) fi.next();
      final BioJavaFeature biojavaFeature = new BioJavaFeature(f, this);

      featTable.add(biojavaFeature);
    }
  }


  private void setArtemisSequence(final SymbolList symbols) {
    artemisSequence = new BioJavaSequence(symbols);
  }


  /**
   *  Set to true in the constructor while features are added.  setDirtyFlag ()
   *  will do nothing while this is true.
   **/
  private boolean in_constructor = false;


  /**
   *  The Date when this Entry last changed or null if this Entry
   *  hasn't changed since the last save.  Set to null by save ().
   **/
  private java.util.Date lastChangeTime = null;

  /**
   *  The Document object that was passed to the constructor.
   **/
  private Document document = null;

  private FeatureTable featTable;
  private EntryInformation entryInformation;
  private Sequence artemisSequence;
  private org.biojava.bio.seq.Sequence bioJavaSequence;
  private SequenceFormat sequenceFormat;
  public void dispose()
  {
    // TODO Auto-generated method stub
    
  }
}
