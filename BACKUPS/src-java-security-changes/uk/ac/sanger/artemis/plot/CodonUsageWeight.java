/* CodonUsageWeight.java
 *
 * created: Tue Apr 13 1999
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1999  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/plot/CodonUsageWeight.java,v 1.1 2004-06-09 09:51:22 tjc Exp $
 */

package uk.ac.sanger.artemis.plot;


import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;

import uk.ac.sanger.artemis.util.*;

import java.io.*;
import java.util.*;

/**
 *  Objects of this class are the source of the codon usage weights that are
 *  used by the methods in the CodonWindowAlgorithm class.
 *  See Gribskov et al. (Nucl. Acids Res. 12(1); 539-549 (1984)).
 *  @author Kim Rutherford
 *  @version $Id: CodonUsageWeight.java,v 1.1 2004-06-09 09:51:22 tjc Exp $
 **/

public class CodonUsageWeight extends CodonWeight {
  /**
   *  Create a new CodonUsageWeight object.
   *  @param usage_file The File to read the codon usage stats from.
   *  @param strand The Strand that this object will be generating values from.
   *    This is needed so that we can get the frequencies of t, c, a and g.
   **/
  public CodonUsageWeight (final File usage_file, final Strand strand)
      throws IOException {
    this.usage_file = usage_file;

    makeSequenceData (strand);

    readFromFile ();
  }

  /**
   *  Returns the name of the file that the usage information was read from.
   **/
  public String getName () {
    return usage_file.getName ();
  }

  /**
   *  Return the codon score for the given codon.  The returned value will be
   *  between 0 and 1 (higher scores means a better value).
   *  @param codon XXX A lowercase string containing the bases of the codon to
   *    look up.
   **/
  public float getCodonValue (final char base1, final char base2,
                              final char base3) {
    final int index_of_base1 = Bases.getIndexOfBase (base1);
    final int index_of_base2 = Bases.getIndexOfBase (base2);
    final int index_of_base3 = Bases.getIndexOfBase (base3);

    if (index_of_base1 > 3 || index_of_base2 > 3 || index_of_base3 > 3) {
      // there is a non-base character in the sequence
      return 1.0F;
    }

    final int index =
      index_of_base1 * 16 + index_of_base2 * 4 + index_of_base3;

    final char translation_character =
      AminoAcidSequence.getCodonTranslation (base1, base2, base3);

    final int synonymous_codon_index =
      AminoAcidSequence.getSymbolIndex (translation_character);

//      System.out.println ("-----> " + data[index] + "   " +
//                          residue_data[synonymous_codon_index] + "  " +
//                          strand_data[index] + "  " +
//                          strand_residue_data[synonymous_codon_index]);

    return
      (data[index]/residue_data[synonymous_codon_index])/
      (strand_data[index]/strand_residue_data[synonymous_codon_index]);
  }

  /**
   *  Parse one line of the form:
   *  <pre>
   *   UUU 32.2( 48423)  UCU 30.5( 45913)  UAU 21.8( 32829)  UGU  8.9( 13371)
   *  </pre>
   *  Returns the four frequency count per 1000 from the line - in the
   *  example above it returns 32.2, 30.5, 21.8 and 8.9
   **/
  private float [] parseLine (String line)
      throws IOException {
    final float [] return_array = new float [4];

    line = line.trim ();

    int return_array_index = 0;

    final StringTokenizer tokeniser = new StringTokenizer (line, " ()");

//      System.out.println (line + ": " + tokeniser.countTokens ());

    if (tokeniser.countTokens () != 12) {
      final String message =
        "garbage codon usage data file at this point: " + line;
      throw new CodonUsageFormatException (message);
    }


    for (int i = 0 ; i < 4 ; ++i) {

      // ignore first token of the block - which is a codon string like "GCA"
      tokeniser.nextToken ();

      final String codon_value = tokeniser.nextToken ();

      return_array[i] = Float.valueOf (codon_value).floatValue ();

//        System.out.println (">" + return_array[i]);

      // ignore the occurrence count
      tokeniser.nextToken ();
    }

    return return_array;
  }

  /**
   *  Read the codon usage information from the given stream.  Each line of
   *  the input should be in this form:
   *  <pre>
   *   UUU 32.2( 48423)  UCU 30.5( 45913)  UAU 21.8( 32829)  UGU  8.9( 13371)
   *  </pre>
   *  The usage information is read into the data and residue_data arrays.
   **/
  private void readDataFromStream (final Document document)
      throws IOException {
    final Reader document_reader = document.getReader ();

    final BufferedReader buffered_reader =
      new BufferedReader (document_reader);

    int current_index = 0;

    while (true) {
      final String this_line = buffered_reader.readLine ();

      if (this_line == null) {
        break;
      }

      if (this_line.trim ().length () == 0) {
        // ignore blank lines
        continue;
      }

      if (current_index >= 64) {
        // we have all the values we need but there are more non blank lines
        final String message =
          "too many lines in the codon usage data file at this point: " +
          this_line;
        throw new CodonUsageFormatException (message);
      }

      // this will return exactly 4 values or will throw a
      // CodonUsageFormatException object
      final float [] line_data = parseLine (this_line);

      for (int i = 0 ; i < 4 ; ++i) {
        final int upper_index = (current_index & 0x0c) * 4;

        final int lower_index = current_index & 0x03;

        final int real_index = upper_index + lower_index + i*4;

        if (line_data[i] < 0.01F) {
          data[real_index] = 0.01F;
//           System.out.println ("zzzzz: " + line_data[i] + "  " + data[real_index]);
        } else {
          data[real_index] = line_data[i];
//           System.out.println ("foo: " + line_data[i] + "  " + data[real_index]);
        }

        final char translation_character =
          codon_translation_array[real_index];

        final int symbol_index =
          AminoAcidSequence.getSymbolIndex (translation_character);

        residue_data[symbol_index] += data[real_index];

//            System.out.println ("--> " + translation_character + "  " +
//                                symbol_index + "  " + line_data[i] + "  " +
//                                residue_data[symbol_index] + "  " + i +
//                                "  " + current_index + "  " +
//                                real_index);

      }

      ++current_index;
    }

//        for (int i = 0 ; i < 64 ; ++i) {
//          System.out.println (">>>> " + data[i]);
//        }

  }

  /**
   *  Read the codon usage information from the file_name that was passed to
   *  the constructor.
   **/
  private void readFromFile ()
      throws IOException {

    if (Options.readWritePossible ()) {
      // try the current directory first
      final Document current_dir_document = new FileDocument (usage_file);

      readDataFromStream (current_dir_document);
    }

    // try to read from the installation directory
//    readDataFromStream (Diana.getCodeDirectory ().append (file_name));
  }

  /**
   *  Fill in the strand_data and strand_residue_data arrays by looking at the
   *  given Strand.
   *  @param strand The Strand that this object will be generating values from.
   *    This is needed so that we can get the frequencies of t, c, a and g.
   **/
  private void makeSequenceData (final Strand strand) {

    for (int first_index = 0 ; first_index < 4 ; ++first_index) {
      final int first_base_count = getStrandBaseCount (first_index, strand);

      for (int second_index = 0 ; second_index < 4 ; ++second_index) {
        final int second_base_count =
          getStrandBaseCount (second_index, strand);

        for (int third_index = 0 ; third_index < 4 ; ++third_index) {
          final int third_base_count =
            getStrandBaseCount (third_index, strand);

          final int strand_data_index =
            first_index * 16 + second_index * 4 + third_index;

          strand_data [strand_data_index] =
            1.0F * first_base_count * second_base_count * third_base_count /
            strand.getSequenceLength ();

//            System.out.println ("---> " + strand_data [strand_data_index] +
//                                "  " + first_base_count);

          final char translation_character =
            AminoAcidSequence.codon_translation_array[strand_data_index];

          final int symbol_index =
            AminoAcidSequence.getSymbolIndex (translation_character);

          strand_residue_data [symbol_index] += strand_data[strand_data_index];

//            System.out.println ("foo: " + strand_residue_data [symbol_index]);
        }
      }
    }
  }

  /**
   *  Return a base count on the given strand.
   *  @param base_index T = 0, C = 1, A = 2 and G = 3.
   **/
  private int getStrandBaseCount (final int base_index, final Strand strand) {
    switch (base_index) {
    case 0:
      return strand.getTCount ();
    case 1:
      return strand.getCCount ();
    case 2:
      return strand.getACount ();
    case 3:
      return strand.getGCount ();
    }
    return 0;
  }

  /**
   *  This table is used for fast lookup of codon translations by
   *  readDataFromStream ().  There is one entry for each codon and the
   *  entries are in this order: TTT, TTC, TTA, TTG, TCT, TCC, ...
   **/
  private static char [] codon_translation_array = {
    'f', 's', 'y', 'c',
    'f', 's', 'y', 'c',
    'l', 's', '#', '*',
    'l', 's', '+', 'w',

    'l', 'p', 'h', 'r',
    'l', 'p', 'h', 'r',
    'l', 'p', 'q', 'r',
    'l', 'p', 'q', 'r',

    'i', 't', 'n', 's',
    'i', 't', 'n', 's',
    'i', 't', 'k', 'r',
    'm', 't', 'k', 'r',

    'v', 'a', 'd', 'g',
    'v', 'a', 'd', 'g',
    'v', 'a', 'e', 'g',
    'v', 'a', 'e', 'g',
  };

  /**
   * The 64 weight values should appear in this order: ttt, ttc, tta, ttg,
   * tct, ..., ggg.
   **/
  final private float [] data = new float [64];

  /**
   *  These are totals for each residue, derived from the data array above.
   **/
  final private float [] residue_data =
    new float [AminoAcidSequence.symbol_count];

  /**
   *  The 64 codon frequencies (per 1000) of the sequence that was passed to
   *  the constructor.
   **/
  final private float [] strand_data = new float [64];

  /**
   *  These are totals for each residue, derived from the strand_data array
   *  above.
   **/
  final private float [] strand_residue_data =
    new float [AminoAcidSequence.symbol_count];

  /**
   *  The File that was passed to the constructor.
   **/
  final private File usage_file;
}

