/* Options.java
 *
 * created: Thu Dec 10 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/Options.java,v 1.15 2009-06-19 15:18:04 tjc Exp $
 **/

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.components.Splash;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.KeyVector;
import uk.ac.sanger.artemis.io.QualifierInfo;
import uk.ac.sanger.artemis.io.QualifierInfoVector;
import uk.ac.sanger.artemis.io.QualifierInfoException;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.SimpleEntryInformation;


import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Properties;
import java.util.Vector;


/**
 *  An object of this class is used to read, write and store all the
 *  configurable options in Diana.  The method Properties.load () is used to
 *  read in new options.
 *
 *  @author Kim Rutherford
 *  @version $Id: Options.java,v 1.15 2009-06-19 15:18:04 tjc Exp $
 **/

public class Options extends Properties 
{
  private static final long serialVersionUID = 1L;

  /**
   *  This is the object that will be returned by getOptions().
   **/
  static private Options options = new Options();

  /** The name of the "direct edit" option. */
  static private final String direct_edit_string = "direct_edit";

  /** The name of the "eukaryotic mode" option. */
  static private final String eukaryotic_mode_string = "organism_type";

  /** The name of the "highlight active entry" option. */
  static private final String highlight_active_entry_string =
    "highlight_active_entry";

  /** Cache for getDisplayQualifierNames(). */
  private static StringVector display_gene_qualifier_names = null;

  /** Cache for getSystematicQualifierNames(). */
  private static StringVector systematic_gene_qualifier_names = null;

  /** Cache for getAllGeneNames(). */
  private static StringVector all_gene_qualifier_names = null;

  /** Used as cache by readWritePossible(). */
  private static Boolean read_succeeded = null;

  /** The EntryInformation object to use for EMBL and GENBANK entries. */
  private static EntryInformation extended_entry_information;

  /** The EntryInformation object to use for EMBL and GENBANK entries. */
  private static EntryInformation db_entry_information;

  /** The default font that should be used for all windows. */
  private Font font = null;

  /**
   *  A map of colour numbers to Color object.  This is initialised by
   *  setDefaultColourMap().
   **/
  private Vector<Color> colour_map = null;

  /** A vector containing ExternalProgram objects that we can run. */
  private ExternalProgramVector external_programs = null;
  
  /** A vector containing ExternalProgram objects that we can run. */
  private ExternalProgramVector ncbi_programs = null;

  /** A vector of those objects listening for option change events. */
  private Hashtable option_listener_hash = new Hashtable();

  /** Set by getInvisibleQualifiers() and reset by resetCachedValues() */
  private StringVector invisible_qualifiers = null;
  
  public static String CACHE_PATH = 
      System.getProperty("user.home") + File.separatorChar +
      ".artemis" + File.separatorChar + "cache" + File.separatorChar;

  /** 
   *  Create a new Options object with default settings for the options.
   **/
  public Options() 
  {
    try 
    {
      put("font_size", "14");
      reset();
    } 
    catch(Throwable e) 
    {
      e.printStackTrace();
    }
  }

  /**
   *  Returns true if and only if the given String is not "no", "false", "n"
   *  or "f".
   **/
  private boolean getPropertyTruthValueInternal(final String value) 
  {
    if(value == null) 
      return false;

    final String lowercase_value = value.toLowerCase();

    if(lowercase_value.equals("false") ||
        lowercase_value.equals("f") || lowercase_value.equals("no") ||
        lowercase_value.equals("n"))  
      return false;
    else 
      return true;
  }

  /**
   *  Returns true if and only if the given property is set and is not "no",
   *  "false", "n" or "f".
   **/
  public boolean getPropertyTruthValue(final String name) 
  {
    return getPropertyTruthValueInternal(getProperty(name));
  }

  /**
   *  Names of qualifiers to search when attempting to find the primary or
   *  display name of a gene.
   **/
  public StringVector getDisplayQualifierNames() 
  {
    if(display_gene_qualifier_names == null) 
    {
      final String name_string = getProperty("display_name_qualifiers");
      
      if(name_string == null) 
      {
        display_gene_qualifier_names = new StringVector();
        display_gene_qualifier_names.add("gene");
      }
      else
        display_gene_qualifier_names = StringVector.getStrings(name_string);
    }
    return display_gene_qualifier_names;
  }

  /**
   *  Names of qualifiers to search when attempting to find the systematic
   *  name of a gene.
   **/
  public StringVector getSystematicQualifierNames() 
  {
    if(systematic_gene_qualifier_names == null) 
    {
      final String name_string = getProperty("systematic_name_qualifiers");
      
      if(name_string == null) 
      {
        systematic_gene_qualifier_names = new StringVector();
        systematic_gene_qualifier_names.add("gene");
      }
      else 
        systematic_gene_qualifier_names =
          StringVector.getStrings(name_string);
    }
    return systematic_gene_qualifier_names;
  }

  /**
   *  Return all possible gene names by combining the return values of
   *  getSystematicQualifierNames() and getDisplayQualifierNames()
   **/
  public StringVector getAllGeneNames() 
  {
    if(all_gene_qualifier_names == null) 
    {
      all_gene_qualifier_names = getSystematicQualifierNames().copy();
      all_gene_qualifier_names.add(getDisplayQualifierNames());
    }

    return all_gene_qualifier_names;
  }


  /**
   *  Read the options from the options file and uk.ac.sanger.artemis.ini
   **/
  private void readOptions() 
  {
    try 
    {
      final InputStream options_input_stream =
        Options.class.getResourceAsStream("/etc/options");

      if(options_input_stream == null) 
        return;

      load(options_input_stream);

      final boolean run_quietly =
        getPropertyTruthValueInternal(System.getProperty("run_quietly"));

      if(readWritePossible()) 
      {
        final String user_home = System.getProperty("user.home");

        // properties are read in order from these files.
        // the null is replaced by the file name given by the extra_options
        // system property (if it is set)
        final String [] standard_option_file_names = 
        {
          "Diana.ini",
          "options",
          "options.txt",
          "options.text",
          user_home + File.separator + ".artemis_options"
        };

        String [] option_file_names = standard_option_file_names;

        final Properties system_properties = System.getProperties();

        if(system_properties != null) 
        {
          final String extra_options_prop =
            system_properties.getProperty("extra_options");

          if(extra_options_prop !=null)
          {
            final StringVector extra_option_file_names =
              StringVector.getStrings(extra_options_prop, ":");

            option_file_names =
              new String [option_file_names.length +
                          extra_option_file_names.size()];

            for(int i = 0 ; i < extra_option_file_names.size() ; ++i) 
            {
              final String extra_option_file_name =
                (String)extra_option_file_names.elementAt(i);

              if(new File(extra_option_file_name).canRead()) 
                option_file_names[i] = extra_option_file_name;
              else 
              {
                System.err.println("warning: could not read options from \"" +
                                    extra_option_file_name + "\"");
                option_file_names[i] = null;
              }
            }

            for(int i = 0 ; i < standard_option_file_names.length ; ++i) 
              option_file_names[i + extra_option_file_names.size()] =
                standard_option_file_names[i];
          }
        }

        for(int i = 0 ; i < option_file_names.length ; ++i) 
        {
          final String this_options_file = option_file_names[i];

          if(this_options_file == null) 
            continue;

          final Document options_document =
            new FileDocument(new File(this_options_file));

          // read the option files if they exist
          if(options_document.readable()) 
          {

            final InputStream options_document_stream =
              options_document.getInputStream();

            /*if(!run_quietly) 
              System.err.println("reading options from \"" +
                                  this_options_file + "\"");*/

            Splash.appendToLog(this_options_file+" options read");
            load(options_document_stream);
          }
          else
            Splash.appendToLog(this_options_file+" not found");
        }
      }
    } 
    catch(IOException e) 
    {
      System.err.println("could not read an options file : " + e);
    }

    for(Enumeration enumProp = propertyNames(); enumProp.hasMoreElements();) 
    {
      final String property_name = (String)enumProp.nextElement();
      fireChangeEvent(property_name);
    }
  }

  /**
   *  Reset this Options object by setting all the options to their default
   *  values and then call readOptions().
   **/
  public void reset()
  {
    clear();

    setDefaultFeatureColours();
    setDefaultColourMap();

    readOptions();
    setChadoOptions();
    
    readSystemOptions();
    resetCachedValues();

/////
//  final Enumeration keys = propertyNames();
//  while (keys.hasMoreElements())
//  {
//    final String key   = (String)keys.nextElement();
//    final String value = getProperty(key);
//    if(key.startsWith("colo"))
//      System.out.println(key+"   \t"+value);   
//  }
/////

  }
  
  private void setChadoOptions()
  {
    final StringVector exon_models = getOptionValues("chado_exon_model");
    if(exon_models != null)
    {
      DatabaseDocument.EXONMODEL = (String)exon_models.elementAt(0);
      if(getProperty("colour_of_" + DatabaseDocument.EXONMODEL) == null)
        put("colour_of_" + DatabaseDocument.EXONMODEL, "7");
    }
    
    final StringVector chado_transcript = getOptionValues("chado_transcript");
    if(chado_transcript != null)
    {
      DatabaseDocument.TRANSCRIPT = (String)chado_transcript.elementAt(0);
      if(getProperty("colour_of_" + DatabaseDocument.TRANSCRIPT) == null)
        put("colour_of_" + DatabaseDocument.TRANSCRIPT, "1");
    }
    
    DatabaseDocument.CHADO_INFER_CDS = getPropertyTruthValue("chado_infer_CDS_UTR");
    if(DatabaseDocument.CHADO_INFER_CDS)
      DatabaseDocument.EXONMODEL = "exon";
  }

  /**
   *  Call load() from the super class and then reset the cached values for
   *  Font and Color objects.
   **/
  public void load(InputStream in_stream)
      throws IOException 
  {
    super.load(in_stream);

    resetCachedValues();
  }

  /**
   *  Return the reference of the global Options object.
   **/
  public static Options getOptions() 
  {
    return options;
  }

  /**
   *  Return true if and only if we are running on a Unix machine.
   **/
  public static boolean isUnixHost() 
  {
    if(System.getProperty("artemis.environment") != null &&
        System.getProperty("artemis.environment").equals("UNIX")) 
      return true;
    else 
      return false;
  }

  /**
   *  Return true if and only if Artemis is in "Noddy mode", meaning there
   *  will be lots of requesters asking for confirmation when delete/editing.
   **/
  public static boolean isNoddyMode() 
  {
    return !isBlackBeltMode();
  }

  /**
   *  Return true if and only if isNoddyMode() will return false.
   **/
  public static boolean isBlackBeltMode() 
  {
    return getOptions().getPropertyTruthValue("black_belt_mode");
  }

  /**
   *  Return true if and only if we are running in an applet.
   **/
  public static boolean readWritePossible() 
  {
    if(read_succeeded == null) 
    {
      try 
      {
        final File temp_file = File.createTempFile("dummy", "dummy");
        read_succeeded = new Boolean(true);
      }
      catch(Throwable _) 
      {
        read_succeeded = new Boolean(false);
      }
    }

    return read_succeeded.booleanValue();
  }

  /**
   *  Return a StringVector containing the values of the given option or null
   *  if there is no such option.  The values are separated by spaces in the
   *  properties file.
   *
   *  @see #getPropertyValues()
   **/
  public StringVector getOptionValues(final String option) 
  {
    return getPropertyValues(this, option);
  }

  /**
   *  Return a StringVector containing the values of the given option in the
   *  given Properties object or null if there is no such option.  The values
   *  are separated by spaces in the properties file.  For example, asking for
   *  the property foo in a Properties object where foo has the value: "thing1
   *  thing2 thing3", will return a vector containing three String objects:
   *  "thing1", "thing2", "thing3"
   **/
  public static StringVector getPropertyValues(final Properties properties,
                                               final String option) 
  {
    final String option_value = properties.getProperty(option);

    if(option_value == null) 
      return null;

    return StringVector.getStrings(option_value);
  }

  /**
   *  Return a QualifierInfoVector object containing one object for each
   *  qualifier given in extra_qualifiers option.
   **/
  public QualifierInfoVector getExtraQualifiers() 
  {
    return getExtraQualifiers("extra_qualifiers");
  }
  
  /**
   *  Return a QualifierInfoVector object containing one object for each
   *  qualifier given in extra_qualifiers option.
   **/
  public QualifierInfoVector getExtraGffQualifiers() 
  {
    return getExtraQualifiers("extra_qualifiers_gff");
  }
  
  /**
   *  Return a QualifierInfoVector object containing one object for each
   *  qualifier given in extra_qualifiers option.
   **/
  private QualifierInfoVector getExtraQualifiers(final String qualifier_options_flag) 
  {
    final QualifierInfoVector return_vector =
      new QualifierInfoVector();

    final StringVector extra_qualifiers_strings =
      getOptionValues(qualifier_options_flag);

    for(int i = 0 ; i < extra_qualifiers_strings.size() / 2 ; ++i) 
    {
      final String name = (String)extra_qualifiers_strings.elementAt(i * 2);
      final String type_string =
        (String)extra_qualifiers_strings.elementAt(i * 2 + 1);
      final int type = QualifierInfo.getQualifierTypeID(type_string);

      return_vector.add(new QualifierInfo(name, type, null, null, false));
    }

    return return_vector;
  }

  /**
   *  Return a Vector of the names of those qualifiers that shouldn't be shown
   *  by the QualifierChoice popup.
   **/
  public StringVector getInvisibleQualifiers(boolean isGFF) 
  {
    if(invisible_qualifiers == null) 
    { 
      invisible_qualifiers = getOptionValues("invisible_qualifiers");

      if(isGFF)
        invisible_qualifiers.add( getOptionValues("invisible_qualifiers_gff") );
      
      if(invisible_qualifiers == null) 
        invisible_qualifiers = new StringVector();
    }

    return invisible_qualifiers;
  }

  /**
   *  Return the path of the default sequence file or null if there is no
   *  default.  This will usually have been set in the uk.ac.sanger.artemis.ini file,
   *  eg. SEQUENCE_FILE='c5H10.seq'
   **/
  public String getDefaultSequenceFileName() 
  {
    final String property_string = getProperty("SEQUENCE_FILE");

    if(property_string == null) 
      return null;
    else 
    {
      if(property_string.length() == 0) 
        return null;

      // trim the quotes from the name
      if(property_string.startsWith("'") &&
          property_string.endsWith("'") &&
          property_string.length() > 2) 
        return property_string.substring(1, property_string.length() - 1);
      else 
        return property_string;
    }
  }

  /**
   *  Return the path of the default feature file or null if there is no
   *  default.  This will usually have been set in the uk.ac.sanger.artemis.ini file,
   *  eg. FEATURE_FILE='c5H10_embl.tab'
   **/
  public String getDefaultFeatureFileName() 
  {
    final String property_string = getProperty("FEATURE_FILE");

    if(property_string == null) 
      return null;
    else 
    {
      if(property_string.length() == 0) 
        return null;

      // trim the quotes from the name
      if(property_string.startsWith("'") &&
          property_string.endsWith("'") &&
          property_string.length() > 2) 
        return property_string.substring(1, property_string.length() - 1);
      else 
        return property_string;
    }
  }

  /**
   *  Return the minimum open reading frame size used when automatically
   *  creating ORFs with the "Mark Open Reading Frames" function.
   **/
  public int getMinimumORFSize() 
  {
    final Integer minimum_orf_size = getIntegerProperty("minimum_orf_size");

    if(minimum_orf_size == null) // default value
      return 100;
    else
      return minimum_orf_size.intValue();
  }

  /**
   *  Get the default colour for the given feature key, as specified in the
   *  options file.  If no colour is specified in the file then null is
   *  returned.
   *  @param key The feature key to get the colour of.
   **/
  public Color getDefaultFeatureColour(Key key) 
  {
    final Integer colour_integer = getIntegerProperty("colour_of_" + key);

    if(colour_integer == null) 
      return null;
    else 
      return getColorFromColourNumber(colour_integer.intValue());
  }

  /**
   *  Given a colour number (perhaps from a /colour qualifier) return an
   *  appropriate Color object.
   **/
  public Color getColorFromColourNumber(int colour_number) 
  {
    // first look up colour_map for speed.  if that fails we then try to look
    // for an appropriate property and turn into a Color object.
    if(colour_number < 0) 
      return null;

    if(colour_number >= colour_map.size() ||
        colour_map.elementAt(colour_number) == null) 
    {
      String col = getProperty("colour_" + colour_number);
      if(col != null)
        return parseColour(col);

      // there is no colour for this colour_number
      return null;
    }
    else 
      return colour_map.elementAt(colour_number);
  }

  /**
   *  Given a String containing three space separated integers (0-255), return
   *  a Color object.  The three integers represent red, green and blue
   *  respectively.  If a Color can't be parsed null is returned.
   **/
  private Color parseColour(String colour_string)
  {
    try 
    {
      // first get three integers from the String

      // trim any whitespace from the ends
      final StringVector value_strings =
        StringVector.getStrings(colour_string.trim());

      final int first_int =
        Integer.valueOf((String)value_strings.elementAt(0)).intValue();
      final int second_int =
        Integer.valueOf((String)value_strings.elementAt(1)).intValue();
      final int third_int =
        Integer.valueOf((String)value_strings.elementAt(2)).intValue();

      return new Color(first_int, second_int, third_int);

    }
    catch(NumberFormatException e) 
    {
      return null;
    }
  }

  /**
   *  Get the default Color for a colour
   **/
  public void setDefaultColourMap() 
  {
    put("colour_0", "255 255 255");  // white
    put("colour_1", "100 100 100");  // dark grey
    put("colour_2", "255 0 0");      // red
    put("colour_3", "0 255 0");      // green
    put("colour_4", "0 0 255");      // blue
    put("colour_5", "0 255 255");    // cyan
    put("colour_6", "255 0 255");    // magenta
    put("colour_7", "255 255 0");    // yellow
    put("colour_8", "152 251 152");  // pale green
    put("colour_9", "135 206 250");  // light sky blue
    put("colour_10", "255 165 0");   // orange
    put("colour_11", "200 150 100"); // brown
    put("colour_12", "255 200 200"); // pink
  }

  /**
   *  Return a vector containing the ExternalProgram objects of all the
   *  external programs that we can use.
   **/
  public ExternalProgramVector getExternalPrograms() 
  {
    if(external_programs == null) 
    {
      external_programs = new ExternalProgramVector();

      final StringVector protein_value_strings =
            getOptionValues("feature_protein_programs");

      if(protein_value_strings != null) 
      {
        for(int i = 0; i < protein_value_strings.size() / 2; ++i) 
        {
          final String program_name =
            (String)protein_value_strings.elementAt(i * 2);
          final String program_options =
            (String)protein_value_strings.elementAt(i * 2 + 1);

          final ExternalProgram program =
            new ExternalProgram(program_name, program_options,
                                ExternalProgram.AA_PROGRAM);

          external_programs.add(program);
        }
      }

      final StringVector dna_value_strings =
        getOptionValues("feature_dna_programs");

      if(dna_value_strings != null) 
      {
        for(int i = 0; i < dna_value_strings.size() / 2; ++i) 
        {
          final String program_name =
            (String)dna_value_strings.elementAt(i * 2);
          final String program_options =
            (String)dna_value_strings.elementAt(i * 2 + 1);

          final ExternalProgram program =
            new ExternalProgram(program_name, program_options,
                                       ExternalProgram.DNA_PROGRAM);

          external_programs.add(program);
        }
      }

      final StringVector application_value_strings =
        getOptionValues("application_programs");

      if(application_value_strings != null) 
      {
        for(int i = 0; i < application_value_strings.size(); ++i) 
        {
          final String program_name = (String)application_value_strings.elementAt(i);

          final ExternalProgram program =
            new ExternalProgram(program_name, null,
                                 ExternalProgram.APPLICATION);

          external_programs.add(program);
        }
      }
    }
    return external_programs;
  }

  
  
  /**
   *  Return a vector containing the ncbi ExternalProgram objects of all the
   *  external programs that we can use.
   **/
  public ExternalProgramVector getNCBIPrograms() 
  {
    if(ncbi_programs == null) 
    {
      ncbi_programs = new ExternalProgramVector();

      final StringVector protein_value_strings =
            getOptionValues("ncbi_protein_search");

      if(protein_value_strings != null) 
      {
        for(int i = 0; i < protein_value_strings.size() / 2; ++i) 
        {
          final String program_name =
            (String)protein_value_strings.elementAt(i * 2);
          final String program_options =
            (String)protein_value_strings.elementAt(i * 2 + 1);

          final ExternalProgram program =
            new ExternalProgram(program_name, program_options,
                                ExternalProgram.AA_PROGRAM);

          ncbi_programs.add(program);
        }
      }

      final StringVector dna_value_strings =
        getOptionValues("ncbi_dna_search");

      if(dna_value_strings != null) 
      {
        for(int i = 0; i < dna_value_strings.size() / 2; ++i) 
        {
          final String program_name =
            (String)dna_value_strings.elementAt(i * 2);
          final String program_options =
            (String)dna_value_strings.elementAt(i * 2 + 1);

          final ExternalProgram program =
            new ExternalProgram(program_name, program_options,
                                       ExternalProgram.DNA_PROGRAM);

          ncbi_programs.add(program);
        }
      }
    }
    return ncbi_programs;
  }

  /**
   *  Return a StringVector containing the bases of the possible start codons
   *  for prokaryotes.  This is stored in the prokaryotic_start_codons option
   *  in the options file.
   **/
//public StringVector getProkaryoticStartCodons() 
//{
//  final StringVector option_values =
//    getOptionValues("prokaryotic_start_codons");

//  for(int i = 0; i<option_values.size() ; ++i) 
//  {
//    final String new_value = option_values.elementAt(i).toLowerCase();
//    option_values.setElementAt(new_value, i);
//  }
//  return option_values;
//}

  /**
   *  Return a StringVector containing the bases of the possible eukaryotic
   *  start codons.  This is stored in the eukaryotic_start_codons option in
   *  the options file.
   **/
//public StringVector getEukaryoticStartCodons() 
//{
//  final StringVector option_values =
//    getOptionValues("eukaryotic_start_codons");

//  for(int i = 0; i<option_values.size() ; ++i) 
//  {
//    final String new_value = option_values.elementAt(i).toLowerCase();
//    option_values.setElementAt(new_value, i);
//  }

//  return option_values;
//}

  public StringVector getStartCodons()
  {
    final StringVector option_values;

    if(getProperty("start_codons") == null)
    {
      if(isEukaryoticMode())
        option_values = getOptionValues("eukaryotic_start_codons");
      else
        option_values = getOptionValues("prokaryotic_start_codons");
    }
    else
      option_values = getOptionValues("start_codons");

    for(int i = 0; i<option_values.size() ; ++i)
    {
      final String new_value = ((String)option_values.elementAt(i)).toLowerCase();
      option_values.setElementAt(new_value, i);
    }

    return option_values;
  }

  /**
   *  Return the default font that should be used for all windows.
   **/
  public Font getFont() 
  {
    return font;
  }

  /**
   *  Return the UIResource for the default font that should be used for all
   *  windows. 
   **/
  public javax.swing.plaf.FontUIResource getFontUIResource() 
  {
    return new javax.swing.plaf.FontUIResource(getFont());
  }

  /**
   *  Adds the specified event listener to receive option change events from
   *  this object.
   *  @param l the event change listener.
   **/
  public void addOptionChangeListener(OptionChangeListener l) 
  {
    option_listener_hash.put(l, this);
  }

  /**
   *  Removes the specified event listener so that it no longer receives
   *  option change events from this object.
   *  @param l the event change listener.
   **/
  public void removeOptionChangeListener(OptionChangeListener l) 
  {
    option_listener_hash.remove(l);
  }

  /**
   *  Read all the system properties and overwrite this Options object with
   *  those values.
   **/
  private void readSystemOptions() 
  {
    if(readWritePossible()) 
    {
      final Properties system_properties = System.getProperties();
      if(system_properties != null) 
      {
        final Enumeration enumeration = system_properties.keys();
        while(enumeration.hasMoreElements()) 
        {
          final String key   = (String)enumeration.nextElement();
          final String value = system_properties.getProperty(key);
          put(key, value);
        }
      }
    }
  }

  /**
   *  Set the default colour number for each a feature key.
   **/
  private void setDefaultFeatureColours() 
  { 
    final Object[] key_colour_map = 
    {
      "CDS", "5",
      "cds?", "7",
      "BLASTCDS", "2",
      "BLASTN_HIT", "6",
      "source", "0",
      "prim_tran", "0",
      "stem_loop", "2",
      "misc_feature", "3",
      "delta", "3",
      "repeat_region", "9",
      "repeat_unit", "9",
      "terminator", "3",
      "promoter", "3",
      "intron", "1",
      "exon", "7",
      DatabaseDocument.EXONMODEL, "7",
      "mRNA", "1",
      "tRNA", "8",
      "TATA", "3",
      "bldA", "2"
    };

    for(int i = 0 ; i < key_colour_map.length / 2 ; ++i) 
      put("colour_of_" + key_colour_map[i*2], key_colour_map[i*2+1]);
  }

  /**
   *  Clear all cached values (such as the font) and then recalculate.
   **/
  private void resetCachedValues() 
  {
    /*final*/ int font_size;

    try 
    {
      final Integer font_size_integer = getIntegerProperty("font_size");

      if(font_size_integer == null) 
        font_size = 12;
      else 
        font_size = font_size_integer.intValue();
    } 
    catch(NumberFormatException e) 
    {
      System.err.println("error in options file - " +
                          "font_size should be an integer");
      // a default value
      font_size = 14;
      put("font_size", String.valueOf(font_size));
    }

    if(getProperty("font_name") == null) 
      put("font_name", "Monospaced");

    font = new Font(getProperty("font_name"), Font.PLAIN, font_size);
    
    colour_map = new Vector<Color>(25);

    int colour_number = 0;

    while(true) 
    {
      final String colour_value_string =
        getProperty("colour_" + colour_number);

      if(colour_value_string == null) 
        // we know nothing about this colour number so we assume that there
        // aren't any more
        break;
      else 
      {
        final Color new_colour = parseColour(colour_value_string);

        if(new_colour == null) 
        {
          // we couldn't parse the colour
          System.err.println("error in options file could not understand " +
                              "this number: " + colour_value_string);
        } 
        else 
        {
          if(colour_number >= colour_map.size()) 
            colour_map.setSize(colour_number + 50);
          colour_map.setElementAt(new_colour, colour_number);
        }
      }
      ++colour_number;
    }

    invisible_qualifiers = null;

    display_gene_qualifier_names = null;
    systematic_gene_qualifier_names = null;
    all_gene_qualifier_names = null;
  }

  /**
   *  Return the value of a property/option as an Integer.  If the property
   *  does not exist or cannot be parsed as an Integer, the method returns
   *  null.
   **/
  public Integer getIntegerProperty(String property_name) 
  {
    final String property_string = getProperty(property_name);

    if(property_string == null) 
      return null;
    else 
      return Integer.valueOf(property_string);
  }

  /**
   *  Return the number of levels of undo to save or 0 if undo is disabled.
   **/
  public int getUndoLevels() 
  {
    final Integer undo_integer = getIntegerProperty("undo_levels");
    if(undo_integer == null) 
      return 0;
    else 
    {
      if(undo_integer.intValue() < 0) 
        return 0;
      else 
        return undo_integer.intValue();
    }
  }

  /**
   *  Enable direct editing if and only if the argument is true.  "Direct
   *  Editing" allows feature locations to be changed by dragging the ends of
   *  the features around with the mouse.
   **/
  public void setDirectEdit(final boolean direct_edit) 
  {
    if(direct_edit) 
    {
      if(!canDirectEdit())
      {
        put(direct_edit_string, "true");
        fireChangeEvent(direct_edit_string);
      }
    } 
    else 
    {
      if(canDirectEdit()) 
      {
        put(direct_edit_string, "false");
        fireChangeEvent(direct_edit_string);
      }
    }
  }

  /**
   *  Returns true if and only if direct editing is enabled.
   **/
  public boolean canDirectEdit() 
  {
    if(getPropertyTruthValue(direct_edit_string)) 
      return true;
    else 
      return false;
  }

  public void setGeneticCode(String table)
  {
    put("translation_table",table);
    fireChangeEvent("translation_table");
  }

  public void setDisplayNameQualifiers(final String display_name_qualifiers)
  {
    display_gene_qualifier_names = null;
    put("display_name_qualifiers", display_name_qualifiers);
    fireChangeEvent("display_name_qualifiers");
  }
  
  /**
   *  Names of qualifiers to search when attempting to find the systematic
   *  name of a gene.
   **/
  public void setSystematicQualifierNames(final String systematic_name_qualifiers) 
  {
    systematic_gene_qualifier_names = null;
    put("systematic_name_qualifiers", systematic_name_qualifiers);
    fireChangeEvent("systematic_name_qualifiers");
  }
  
  /**
   *  Set the organism type to eukaryotic if and only if the argument is
   *  true.  The other alternative is prokaryotic.
   **/
  public void setEukaryoticMode(final boolean eukaryotic_mode)
  {
    if(eukaryotic_mode) 
    {
      if(!isEukaryoticMode()) 
      {
        put(eukaryotic_mode_string, "eukaryotic");
        fireChangeEvent(eukaryotic_mode_string);
      }
    } 
    else 
    {
      if(isEukaryoticMode()) 
      {
        put(eukaryotic_mode_string, "prokaryotic");
        fireChangeEvent(eukaryotic_mode_string);
      }
    }
  }

  /**
   *  Returns true if and only if we should be using the eukaryotic settings.
   *  This is only the default.
   **/
  public boolean isEukaryoticMode() 
  {
    final String organism_type_prop = getProperty("organism_type");

    if(organism_type_prop == null ||
        organism_type_prop.equals("eukaryotic") ||
        organism_type_prop.equals("eukaryote") ||
        organism_type_prop.equals("euk")) 
      return true;
    else 
      return false;
  }

  /**
   *  Set whether or not to highlight the active entry.
   **/
  public void setHighlightActiveEntryFlag(final boolean highlight_active) 
  {
    if(highlight_active) 
    {
      if(!highlightActiveEntryFlag()) 
      {
        put(highlight_active_entry_string, "true");
        fireChangeEvent(highlight_active_entry_string);
      }
    } 
    else
    {
      if(highlightActiveEntryFlag()) 
      {
        put(highlight_active_entry_string, "false");
        fireChangeEvent(highlight_active_entry_string);
      }
    }
  }

  /**
   *  Returns true if and only if we should highlight the active entry in the
   *  display.
   **/
  public boolean highlightActiveEntryFlag() 
  {
    if(getPropertyTruthValue(highlight_active_entry_string)) 
      return true;
    else 
      return false;
  }

  /**
   *  Returns true if this version of Java is running on a GNU/Linux machine
   *  and is version 1.2.0 or later.
   **/
  public boolean isBuggyLinuxVM() 
  {
    if(!readWritePossible()) 
      // Java in a browser has problems, but not this problem
      return false;

    final String os_name = (String) System.getProperties().get("os.name");
    
    if(os_name.equals("Linux"))
    {
      final String java_version =
        (String) System.getProperties().get("java.version");

      if(java_version.startsWith("1.1.")) 
        return false;
      else 
        return true;
    }
    else 
      return false;
  }


  /**
   *  Return the EntryInformation object to use for EMBL and GENBANK entries.
   **/
  public static EntryInformation getDBEntryInformation() 
  {
    return db_entry_information;
  }

  /**
   *  The EntryInformation object that has all the information that
   *  getDBEntryInformation () would return, but which has the non-standard
   *  qualifiers and keys added.
   **/
  public static EntryInformation getArtemisEntryInformation() 
  {
    return extended_entry_information;
  }

  /**
   *  Return an EntryInformation object that is suitable for EMBL and GENBANK
   *  entries.
   **/
  private static EntryInformation makeEntryInformation()
      throws IOException, QualifierInfoException 
  {
    final InputStream feature_keys_stream =
      Options.class.getResourceAsStream("/etc/feature_keys");

    final InputStream qualifier_types_stream =
      Options.class.getResourceAsStream("/etc/qualifier_types");

    QualifierInfoVector qualifier_info_vector =
      readQualifierInfo(qualifier_types_stream, feature_keys_stream);

    final EntryInformation entry_information = new SimpleEntryInformation();

    for(int i = 0 ; i < qualifier_info_vector.size() ; ++i) 
    {
      final QualifierInfo qualifier_info =
        qualifier_info_vector.elementAt(i);

      entry_information.addQualifierInfo(qualifier_info);
    }

    entry_information.setEMBLFormat(true);

    return entry_information;
  }

  /**
   *  Return an EntryInformation object that is suitable for EMBL and GENBANK
   *  entries, and has some useful non-standard additions (specified by the
   *  options file).
   **/
  private static EntryInformation
    makeArtemisEntryInformation(final EntryInformation standard_entry_info)
      throws QualifierInfoException 
  {

    final StringVector extra_keys =
      getOptions().getOptionValues("extra_keys");

    final QualifierInfoVector extra_qualifiers =
      getOptions().getExtraQualifiers();

    final EntryInformation return_entry_information =
      new SimpleEntryInformation(standard_entry_info);

    for(int i = 0 ; i < extra_keys.size() ; ++i) 
    {
      final Key new_key = new Key((String)extra_keys.elementAt(i));
      return_entry_information.addKey(new_key);
    }

    for(int i = 0 ; i < extra_qualifiers.size() ; ++i) 
    {
      final QualifierInfo new_qualifier_info = extra_qualifiers.elementAt(i);
      return_entry_information.addQualifierInfo(new_qualifier_info);
    }

    // make a add qualifier info for each search program.  eg. add blastp_file
    // for blastp
    final ExternalProgramVector external_programs =
      getOptions().getExternalPrograms();

    for(int i = 0 ; i < external_programs.size() ; ++i) 
    {
      final ExternalProgram external_program = external_programs.elementAt(i);

      if(external_program.getType() == ExternalProgram.AA_PROGRAM ||
          external_program.getType() == ExternalProgram.DNA_PROGRAM) 
      {
        final QualifierInfo new_qualifier_info =
          new QualifierInfo(external_program.getName() + "_file",
                             QualifierInfo.QUOTED_TEXT,
                             null,
                             null,
                             true);

        return_entry_information.addQualifierInfo(new_qualifier_info);
      }
    }

    return_entry_information.setEMBLFormat(false);

    return return_entry_information;
  }

  /**
   *  Read the possible feature key and qualifier names and types from the two
   *  given streams (see etc/feature_keys and etc/qualifier_types for details
   *  on the formats).
   **/
  private static QualifierInfoVector
    readQualifierInfo(final InputStream qualifier_types_stream,
                       final InputStream feature_keys_stream)
      throws IOException 
  {

    final QualifierInfoVector return_vector = new QualifierInfoVector();

    Properties feature_properties = new Properties();
    final Properties qualifier_properties = new Properties();

    feature_properties.load(feature_keys_stream);
    qualifier_properties.load(qualifier_types_stream);

    // parse the feature_properties

    {
      final Properties new_feature_properties = new Properties();

      final Enumeration feature_enum = feature_properties.propertyNames();

      while(feature_enum.hasMoreElements()) {
        String current_feature_name = (String) feature_enum.nextElement();

        final StringVector property_values =
          Options.getPropertyValues(feature_properties, current_feature_name);

        new_feature_properties.put(current_feature_name, property_values);
      }

      feature_properties = new_feature_properties;
    }

    final Enumeration qualifier_enum = qualifier_properties.propertyNames();

    while(qualifier_enum.hasMoreElements()) 
    {
      String current_qualifier_name = (String) qualifier_enum.nextElement();

      final StringVector current_qualifier_values =
        Options.getPropertyValues(qualifier_properties,
                                   current_qualifier_name);

      final boolean once_only =
        current_qualifier_values.elementAt(0).equals("yes");

      final String type_string = (String)current_qualifier_values.elementAt(1);

      // find the keys for which this qualifier name is valid or required

      final KeyVector valid_keys = new KeyVector();
      final KeyVector required_keys = new KeyVector();

      final Enumeration features_enum = feature_properties.propertyNames();

      while(features_enum.hasMoreElements()) 
      {
        final String current_key_string = (String)features_enum.nextElement();

        final Key current_key = new Key(current_key_string);

        final StringVector current_feature_qualifiers =
          (StringVector) feature_properties.get(current_key_string);

        if(current_feature_qualifiers.contains(current_qualifier_name)) 
          valid_keys.add(current_key);
        else
          if(current_feature_qualifiers.contains("@" +
                                                   current_qualifier_name)) 
          {
            valid_keys.add(current_key);
            required_keys.add(current_key);
          }
      }

      final int type = QualifierInfo.getQualifierTypeID(type_string);

      final QualifierInfo qualifier_info =
        new QualifierInfo(current_qualifier_name, type, valid_keys,
                           required_keys, once_only);

      return_vector.add(qualifier_info);
    }

    return return_vector;
  }

  /**
   *  Send a change event to all the listeners.
   *  @param option_name The name of the option that changed.
   **/
  private void fireChangeEvent(final String option_name)
  {
    for(final Enumeration e = option_listener_hash.keys() ;
         e.hasMoreElements() ;) 
    {
      final OptionChangeEvent event =
        new OptionChangeEvent(this, option_name);

      final OptionChangeListener target =
        (OptionChangeListener) e.nextElement();

      target.optionChanged(event);
    }
  }

  static 
  {
    try 
    {
      db_entry_information = makeEntryInformation();
      extended_entry_information =
        makeArtemisEntryInformation(db_entry_information);
    }
    catch(QualifierInfoException e) 
    {
      System.err.println("could not initialise the embl package: " +
                          e.getMessage());
      System.exit(1);
    } 
    catch(IOException e) 
    {
      System.err.println("could not initialise the embl package: " +
                          e.getMessage());
      System.exit(1);
    }
  }

}
