/* AminoAcidSequence.java
 *
 * created: Sat Dec 19 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/sequence/AminoAcidSequence.java,v 1.11 2007-07-09 12:38:51 tjc Exp $
 */

package uk.ac.sanger.artemis.sequence;

import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.StringVector;

/**
 *  Objects of this class represent a string of amino acids.
 *
 *  @author Kim Rutherford
 *  @version $Id: AminoAcidSequence.java,v 1.11 2007-07-09 12:38:51 tjc Exp $
 **/

public class AminoAcidSequence 
{
  /**
   *  Create a new AminoAcidSequence object from a string containing single
   *  character amino acids symbols.
   **/
  public AminoAcidSequence(String amino_acid_string) 
  {
    this.amino_acid_string = amino_acid_string; 
  }

  /**
   *  Translate a sequence of bases into the corresponding single letter amino
   *  acid codes.
   *  @param bases The bases to translated.  If the string length is not a
   *    multiple of three the last codon is incomplete and will not be
   *    translated.
   *  @param unknown_is_x If this parameter is true codons that contain
   *    ambiguous bases will be translated as 'x', if false they will be
   *    translated as '.'
   *  @return The translated sequence in one letter abbreviated form.
   **/
  public static AminoAcidSequence getTranslation(final String bases,
                                                 final boolean unknown_is_x) 
  {
//  this is set in Splash.java
//  setGeneCode();    
//
    final StringBuffer aa_buffer = new StringBuffer();
    final int number_of_codons = bases.length() / 3;

    for(int i = 0 ; i < number_of_codons * 3 ; i += 3) 
    {
      final char aa = getCodonTranslation(bases.charAt(i),
                                          bases.charAt(i+1),
                                          bases.charAt(i+2));
      if(aa == '.' && unknown_is_x) 
        aa_buffer.append ('x');
      else 
        aa_buffer.append (aa);
    }
    return new AminoAcidSequence(aa_buffer.toString());
  }

  /**
   *  Translate a sequence of bases into the corresponding single letter amino
   *  acid codes.
   *  @param bases The bases to translated.  If the string length is not a
   *    multiple of three the last codon is incomplete and will not be
   *    translated.
   *  @param unknown_is_x If this parameter is true codons that contain
   *    ambiguous bases will be translated as 'x', if false they will be
   *    translated as '.'
   *  @return The translated sequence in one letter abbreviated form.
   **/
  public static AminoAcidSequence getTranslation(final char[] bases,
                                                 final boolean unknown_is_x)
  {
//  this is set in Splash.java
//  setGeneCode();
//
    final StringBuffer aa_buffer = new StringBuffer();
    final int number_of_codons = bases.length / 3;

    for(int i = 0 ; i < number_of_codons * 3 ; i += 3)
    {
      final char aa = getCodonTranslation(bases[i],
                                          bases[i+1],
                                          bases[i+2]);
      if(aa == '.' && unknown_is_x)
        aa_buffer.append ('x');
      else
        aa_buffer.append (aa);
    }
    return new AminoAcidSequence(aa_buffer.toString());
  }

  /**
   *  Translate a sequence of bases into the corresponding single letter amino
   *  acid codes and appending 2 spaces after each amino acid character.
   *  @param bases The bases to translated.  If the string length is not a
   *    multiple of three the last codon is incomplete and will not be
   *    translated.
   *  @param unknown_is_x If this parameter is true codons that contain
   *    ambiguous bases will be translated as 'x', if false they will be
   *    translated as '.'
   *  @return The translated sequence in one letter abbreviated form.
   **/
  public static AminoAcidSequence getSpacedTranslation(final String bases,
                                               final boolean unknown_is_x)
  {
    final StringBuffer aa_buffer = new StringBuffer();
    final int number_of_codons = bases.length() / 3;
    for(int i = 0 ; i < number_of_codons * 3 ; i += 3)
    {
      final char aa = getCodonTranslation(bases.charAt(i),
                                          bases.charAt(i+1),
                                          bases.charAt(i+2));
      if(aa == '.' && unknown_is_x)
        aa_buffer.append('x');
      else
        aa_buffer.append(aa);
      aa_buffer.append("  ");
    }
    return new AminoAcidSequence(aa_buffer.toString());
  }

  /**
   *  Translate a sequence of bases into the corresponding single letter amino
   *  acid codes and appending 2 spaces after each amino acid character.
   *  @param bases The bases to translated.  If the string length is not a
   *    multiple of three the last codon is incomplete and will not be
   *    translated.
   *  @param unknown_is_x If this parameter is true codons that contain
   *    ambiguous bases will be translated as 'x', if false they will be
   *    translated as '.'
   *  @return The translated sequence in one letter abbreviated form.
   **/
  public static AminoAcidSequence getSpacedTranslation(final char bases[],
                                               final boolean unknown_is_x)
  { 
    final StringBuffer aa_buffer = new StringBuffer();
    final int number_of_codons = bases.length / 3;
    for(int i = 0 ; i < number_of_codons * 3 ; i += 3)
    {
      final char aa = getCodonTranslation(bases[i],
                                          bases[i+1],
                                          bases[i+2]);
      if(aa == '.' && unknown_is_x)
        aa_buffer.append('x');
      else
        aa_buffer.append(aa);
      aa_buffer.append("  ");
    }
    return new AminoAcidSequence(aa_buffer.toString());
  }

  /**
   *  Translate a single codon into the corresponding single letter amino acid
   *  code.
   *  @param codon_string A three character lowercase String containing the
   *    bases to translate.
   *  @return The translated sequence in one letter abbreviated form.  The
   *    return value will be '.' if the codon_string isn't a codon (eg. it
   *    contains more or less than three letters or the letters aren't from
   *    "CTAG")
   **/
  public static char getCodonTranslation(String codon_string) 
  {
    if(codon_string.length() < 3) 
      return '.';
    
    return getCodonTranslation(codon_string.charAt(0),
                               codon_string.charAt(1),
                               codon_string.charAt(2));
  }

  /**
   *  Translate a single codon into the corresponding single letter amino acid
   *  code.
   *  @param first_letter The first base of the codon.
   *  @param second_letter The second base of the codon.
   *  @param third_letter The third base of the codon.
   *  @return The translated sequence in one letter abbreviated form.  The
   *    return value will be '.' if the letters do not form a codon.
   **/
  public final static char getCodonTranslation(char first_letter,
                                               char second_letter,
                                               char third_letter)
  {
    final int first_index = getIndexOfBase(first_letter);
    if(first_index >= 4)
      return '.';

    final int second_index = getIndexOfBase(second_letter);
    if(second_index >= 4) 
      return '.';

    final int third_index = getIndexOfBase(third_letter);
    if(third_index >= 4) 
      return '.';

    final int codon_index = first_index * 16 + second_index * 4 + third_index;

    return codon_translation_array[codon_index];
  }

  /**
   *  Given a base letter return its index where t = 0, c = 1, a = 2, g = 3, 4
   *  otherwise.
   *  See letter_index.
   **/
  private final static int getIndexOfBase(final char base)    
  { 
    switch(base) 
    {
      case 'c':
        return 1;
      case 'a':
        return 2;
      case 'g':
        return 3;
      case 't':
      case 'u': 
        return 0;
    }
  
    return 4;
  }


  /**
   *  Return the number of units in this amino acid sequence.
   **/
  public int length() 
  {
    return amino_acid_string.length();
  }

  /**
   *  Return the one letter codon code of the codon at the given index
   *  (counting from zero).
   **/
  public char elementAt(final int index) 
  {
    return amino_acid_string.charAt(index);
  }

  /**
   *  Return the total molecular weight of the amino acids in this
   *  AminoAcidSequence..
   **/
  public float getMolecularWeight() 
  {
    float return_weight = 0;

    for(int i = 0 ; i < amino_acid_string.length() ; ++i) 
    {
      final char this_char = amino_acid_string.charAt(i);

      return_weight +=
        molecular_weights[getSymbolIndex(this_char)];
    }

    if(amino_acid_string.length() > 1) 
    {
      // need to take off the weight of a water molecule for each peptide bond
      return return_weight -
        molecular_weight_of_water * (amino_acid_string.length () - 1);
    }
    else
      return return_weight;
  }

  /**
   *  Return a string representation of this object.  This string will contain
   *  the one character amino acid codes for each acid in sequence.
   **/
  public String toString()
  {
    return amino_acid_string;
  }

  /**
   *  Search the subject_sequence for this AminoAcidSequence as a substring.
   *  'X' AAs are treated as wildcards in both sequences.
   **/
  public boolean checkForMatch(final AminoAcidSequence subject_sequence) 
  {
    final String subject_sequence_string = subject_sequence.toString();

    for(int subject_index = 0;
        subject_index < subject_sequence_string.length() -
          toString ().length () + 1;
         ++subject_index)
    {
      int query_index = 0;
      //boolean is_matching = true;
      for(; query_index < toString().length(); 
           ++query_index)
      {
        final char this_query_char =
                     toString().charAt(query_index);
        final char this_subject_char =
                     subject_sequence_string.charAt(subject_index + query_index);
        if(!aminoAcidMatches(this_subject_char,
                             this_query_char)) 
          break;
      }

      if(query_index == toString().length())
        return true;
    }

    return false;
  }

  /**
   *  Return true if and only if the two argument are the same AA or if one is
   *  an X.
   **/
  private static boolean aminoAcidMatches(final char aa_char1,
                                          final char aa_char2) 
  {
    if (aa_char1 == aa_char2) 
      return true;
    else 
    {
      if(aa_char1 == 'x' || aa_char2 == 'x')
        return true;
      else
        return false;
    }
  }

  /**
   *  Find the next occurrence of this seqeuence on either Strand of the given
   *  Bases object.  This method searches both strands simultaneously by
   *  searching the underlying Bases object directly.
   *  @param bases This holds the Strand objects to search.
   *  @param search_start_marker The match that will be returned will be
   *    after this base Marker position.  Position 1 on the reverse strand is
   *    considered to be after position 1 on the forward strand, forward
   *    position 2 is after reverse position 1, reverse position 2 is after
   *    forward position 2, etc.  This scheme allows the caller to iterate
   *    through all matches.
   *  @param search_end_position The search will not extend past this
   *    position.
   *  @param search_backwards If true the search will move from last base to
   *    first base, otherwise first to last.
   *  @return A MarkerRange covering the matching bases or null if there is no
   *    match in the given range.
   **/
  public MarkerRange findMatch(final Bases bases,
                               final Marker search_start_marker,
                               final boolean search_backwards,
                               final boolean search_fwd_strand,
                               final boolean search_bwd_strand) 
  {
    final String bases_string = bases.toString();

    // search the bases_string forward for the pattern_string and its
    // complement

    // the String index position in bases_string at which to start the search
    // for this pattern
    final int forward_search_start_index;

    // the String index position in bases_string at which to start the search
    // for the reverse complement of this position
    final int complement_search_start_index;

    if(search_backwards) 
    {
      if(search_start_marker == null) 
      {
        forward_search_start_index = bases.getLength () - 1;
        complement_search_start_index = bases.getLength () - 1;
      } 
      else 
      {
        complement_search_start_index =
                     search_start_marker.getRawPosition() - 2;
        if(search_start_marker.getStrand().isForwardStrand()) 
        {
          forward_search_start_index =
            search_start_marker.getRawPosition() - 2;
        } 
        else 
        {
          forward_search_start_index =
            search_start_marker.getRawPosition() - 1;
        }
      }
    } 
    else
    {
      if(search_start_marker == null) 
      {
        forward_search_start_index = 0;
        complement_search_start_index = 0;
      }
      else
      {
        forward_search_start_index = search_start_marker.getRawPosition();
        if(search_start_marker.getStrand().isForwardStrand()) 
        {
          complement_search_start_index =
            search_start_marker.getRawPosition() - 1;
        } 
        else
        {
          complement_search_start_index =
            search_start_marker.getRawPosition();
        }
      }
    }

    final int forward_search_result;
    if(search_fwd_strand)
      forward_search_result = searchFor(bases_string,
                                        forward_search_start_index,
                                        search_backwards);
    else
      forward_search_result = -1;

    final int complement_search_result;
    if(search_bwd_strand)
      complement_search_result = reverseComplementSearchFor(bases_string,
                                       complement_search_start_index,
                                       search_backwards);
    else
      complement_search_result = -1;

    final int match_first_base;
    final int match_last_base;

    final Strand match_strand;

    if(forward_search_result == -1)
    {
      // no match
      if(complement_search_result == -1) 
        return null;
    }

    if(search_backwards) 
    {
      // take the match that is closest to the end, or the complement match if
      // there is a tie
      if(complement_search_result != -1 &&
         (forward_search_result == -1 ||
          complement_search_result >= forward_search_result)) 
      {
        match_first_base =
          bases.getComplementPosition (complement_search_result + 1);
        match_last_base = match_first_base - (length () * 3 - 1);
        match_strand = bases.getReverseStrand ();
      }
      else
      {
        match_first_base = forward_search_result + 1;
        match_last_base = match_first_base + length () * 3 - 1;
        match_strand = bases.getForwardStrand ();
      }
    }
    else
    {
      // take the match that is closest to base 1, or the forward match if
      // there is a tie
      if(forward_search_result != -1 &&
         (complement_search_result == -1 ||
          forward_search_result <= complement_search_result)) 
      {
        match_first_base = forward_search_result + 1;
        match_last_base = match_first_base + length () * 3 - 1;
        match_strand = bases.getForwardStrand ();
      }
      else
      {
        match_first_base =
          bases.getComplementPosition (complement_search_result + 1);
        match_last_base = match_first_base - (length () * 3 - 1);
        match_strand = bases.getReverseStrand ();
      }
    }

    try 
    {
      return new MarkerRange(match_strand,
                             match_first_base,
                             match_last_base);
    } 
    catch (OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Search for this AminoAcidSequence in the given String of bases.  The
   *  String is treated as a sequence of bases and this AminoAcidSequence is
   *  searched for in each of the three reading frames.
   *  @param bases_string Search this String for the amino acid sequence.
   *  @param start_index This is the index in bases_string where the search
   *    should start.
   *  @param search_backwards If true the search will move from last base to
   *    first base, otherwise first to last.
   *  @return The index of the match or -1 if there is no match.
   **/
  private int searchFor(final String bases_string,
                       final int start_index,
                       final boolean search_backwards) 
  {
    if(search_backwards) 
      return searchBackwardFor(bases_string, start_index);
    else 
      return searchForwardFor(bases_string, start_index);
  }

  /**
   *  Search forward for this AminoAcidSequence in the given String of bases.
   *  The String is treated as a sequence of bases and this AminoAcidSequence
   *  is searched for in each of the three reading frames.
   *  @param bases_string Search this String for the amino acid sequence.
   *  @param start_index This is the index in bases_string where the search
   *    should start.
   *  @return The index of the match or -1 if there is no match.
   **/
  private int searchForwardFor(final String bases_string,
                              final int start_index) 
  {
    final int pattern_base_length = length() * 3;

    for(int base_index = start_index;
        base_index <= bases_string.length() - pattern_base_length ;
        ++base_index) 
    {
      boolean matched = true;

      for(int offset = 0 ; offset < length(); ++offset) 
      {
        final char search_aa = amino_acid_string.charAt(offset);

        // X matches any AA
        if(search_aa == 'x' || search_aa == 'X') 
          continue;

        final char base1 = bases_string.charAt(base_index + offset * 3 + 0);
        final char base2 = bases_string.charAt(base_index + offset * 3 + 1);
        final char base3 = bases_string.charAt(base_index + offset * 3 + 2);

        if(getCodonTranslation(base1, base2, base3) != search_aa)
        {
          matched = false;
          break;
        }
      }

      if(matched)
        return base_index;
    }

    return -1;
  }

  /**
   *  Search backward for this AminoAcidSequence in the given String of bases.
   *  The String is treated as a sequence of bases and this AminoAcidSequence
   *  is searched for in each of the three reading frames.
   *  @param bases_string Search this String for the amino acid sequence.
   *  @param start_index This is the index in bases_string where the search
   *    should start.
   *  @return The index of the match or -1 if there is no match.
   **/
  private int searchBackwardFor(final String bases_string,
                               int start_index) 
  {
    if(bases_string.length() - start_index < length() * 3) 
      start_index = bases_string.length() - length() * 3;

    for(int base_index = start_index; base_index >= 0;
        --base_index) 
    {
      boolean matched = true;

      for(int offset = 0 ; offset < length() ; ++offset) 
      {
        final char search_aa = amino_acid_string.charAt(offset);

        // X matches any AA
        if(search_aa == 'x' || search_aa == 'X') 
          continue;

        final char base1 = bases_string.charAt(base_index + offset * 3 + 0);
        final char base2 = bases_string.charAt(base_index + offset * 3 + 1);
        final char base3 = bases_string.charAt(base_index + offset * 3 + 2);

        if(getCodonTranslation(base1, base2, base3) != search_aa) 
        {
          matched = false;
          break;
        }
      }

      if(matched)
        return base_index;
    }

    return -1;
  }

  /**
   *  Search for this AminoAcidSequence in the reverse complement of the given
   *  String of bases.  The String is treated as a sequence of bases and this
   *  AminoAcidSequence is searched for in each of the three reading frames.
   *  @param bases_string Search this String for the amino acid sequence.
   *  @param start_index This is the index in bases_string where the search
   *    should start.
   *  @param search_backwards If true the search will move from last base to
   *    first base, otherwise first to last.
   *  @return The index of the match or -1 if there is no match.
   **/
  private int reverseComplementSearchFor(final String bases_string,
                                        final int start_index,
                                        final boolean search_backwards) 
  {
    if(search_backwards) 
      return reverseComplementSearchBackwardFor(bases_string, start_index);
    else 
      return reverseComplementSearchForwardFor(bases_string, start_index);
  }

  /**
   *  Search forward for this AminoAcidSequence in the reverse complement of
   *  the given String of bases.  The String is treated as a sequence of bases
   *  and this AminoAcidSequence is searched for in each of the three reading
   *  frames.
   *  @param bases_string Search this String for the amino acid sequence.
   *  @param start_index This is the index in bases_string where the search
   *    should start.
   *  @return The index of the match or -1 if there is no match.
   **/
  private int reverseComplementSearchForwardFor(final String bases_string,
                                               final int start_index) 
  {
    final int pattern_base_length = length() * 3;

    for(int base_index = start_index ;
        base_index <= bases_string.length() - pattern_base_length ;
        ++base_index) 
    {
      boolean matched = true;

      for(int offset = 0; offset < length (); ++offset)
      {
        final char base1 =
          Bases.complement(bases_string.charAt(base_index + offset * 3 + 0));
        final char base2 =
          Bases.complement(bases_string.charAt(base_index + offset * 3 + 1));
        final char base3 =
          Bases.complement(bases_string.charAt(base_index + offset * 3 + 2));

        final char amino_acid_char =
          amino_acid_string.charAt(amino_acid_string.length() - offset - 1);

        // X matches any AA
        if(amino_acid_char == 'x' || amino_acid_char == 'X')
          // X matches any AA
          continue;

        if(getCodonTranslation(base3, base2, base1) != amino_acid_char) 
        {
          matched = false;
          break;
        }
      }

      if(matched)
        return base_index;
    }

    return -1;
  }

  /**
   *  Search backward for this AminoAcidSequence in the reverse complement of
   *  the given String of bases.  The String is treated as a sequence of bases
   *  and this AminoAcidSequence is searched for in each of the three reading
   *  frames.
   *  @param bases_string Search this String for the amino acid sequence.
   *  @param start_index This is the index in bases_string where the search
   *    should start.
   *  @return The index of the match or -1 if there is no match.
   **/
  private int reverseComplementSearchBackwardFor(final String bases_string,
                                                int start_index) 
  {
    if(bases_string.length() - start_index < length() * 3) 
      start_index = bases_string.length() - length() * 3;

    for(int base_index = start_index; base_index >= 0;
         --base_index)
    {
      boolean matched = true;

      for(int offset = 0 ; offset < length() ; ++offset)
      {
        final char base1 =
          Bases.complement(bases_string.charAt(base_index + offset * 3 + 0));
        final char base2 =
          Bases.complement(bases_string.charAt(base_index + offset * 3 + 1));
        final char base3 =
          Bases.complement(bases_string.charAt(base_index + offset * 3 + 2));

        final char amino_acid_char =
          amino_acid_string.charAt(amino_acid_string.length() - offset - 1);

        // X matches any AA
        if(amino_acid_char == 'x' || amino_acid_char == 'X') 
          continue;

        if(getCodonTranslation(base3, base2, base1) != amino_acid_char)
        {
          matched = false;
          break;
        }
      }

      if(matched) 
        return base_index;
    }

    return -1;
  }

  /**
   *  Return true if and only if this sequence contains a stop codon.
   **/
  public boolean containsStopCodon()
  {
    for(int i = 0 ; i < amino_acid_string.length() ; ++i) 
    {
      final char this_char = amino_acid_string.charAt(i);

      if(isStopCodon (this_char))
        return true;
    }

    return false;
  }

  /**
   *  Return true if and only if the given amino acid symbol is the
   *  translation of a stop codon. ie #, * or +.
   **/
  public static boolean isStopCodon(final char amino_acid_char) 
  {
    if(amino_acid_char == '#' ||
       amino_acid_char == '*' ||
       amino_acid_char == '+') 
      return true;
    else 
      return false;
  }

  /**
   *  Return true if and only if the given one letter code symbol is the a
   *  legal amino acid or stop symbol.
   **/
  protected static boolean isLegalCodon(char one_letter_code)
  {
    one_letter_code = Character.toLowerCase(one_letter_code);
    switch(one_letter_code) 
    {
      case 'a': case 'r': case 'n': case 'd': case 'c': case 'q': case 'e':
      case 'g': case 'h': case 'i': case 'l': case 'k': case 'm': case 'f':
      case 'p': case 's': case 't': case 'w': case 'y': case 'v': case '*':
      case '#': case '+':
        return true;
      default:
        return false;
    }
  }

  /**
   *  This table is used for fast lookup of codon translations by
   *  getCodonTranslation().  There is one entry for each codon and the
   *  entries are in this order: TTT, TTC, TTA, TTG, TCT, TCC, ...
   **/
  final public static char [] codon_translation_array = {
    'f', 'f', 'l', 'l',
    's', 's', 's', 's',
    'y', 'y', '#', '+',
    'c', 'c', '*', 'w',

    'l', 'l', 'l', 'l',
    'p', 'p', 'p', 'p',
    'h', 'h', 'q', 'q',
    'r', 'r', 'r', 'r',

    'i', 'i', 'i', 'm',
    't', 't', 't', 't',
    'n', 'n', 'k', 'k',
    's', 's', 'r', 'r',

    'v', 'v', 'v', 'v',
    'a', 'a', 'a', 'a',
    'd', 'd', 'e', 'e',
    'g', 'g', 'g', 'g'
  };

  /**
   *  Used by getAminoAcidType().
   **/
  public final static int POLAR_UNCHARGED_AA    = 0;
  /**
   *  Used by getAminoAcidType().
   **/
  public final static int POSITIVELY_CHARGED_AA = 1;
  /**
   *  Used by getAminoAcidType().
   **/
  public final static int NEGATIVELY_CHARGED_AA = 2;
  /**
   *  Used by getAminoAcidType().
   **/
  public final static int HYDROPHOBIC_AA        = 3;
  /**
   *  Used by getAminoAcidType().
   **/
  public final static int SPECIAL_AA            = 4;
  /**
   *  Used by getAminoAcidType().
   **/
  public final static int STOP_AA               = 5;
  /**
   *  Used by getAminoAcidType().
   **/
  public final static int UNKNOWN_AA            = 6;
  /**
   *  Used by getAminoAcidType().
   **/
  public final static int ILLEGAL_AA            = 7;
  
  /**
   *  Returns one of POLAR_UNCHARGED_AA, POSITIVELY_CHARGED_AA,
   *  NEGATIVELY_CHARGED_AA, HYDROPHOBIC_AA, SPECIAL_AA or STOP_AA depending
   *  on the aa_char argument.
   **/
  public static int getAminoAcidType(final char aa_char) 
  {
    switch (aa_char) 
    {
      case 'S': case 'T': case 'N': case 'Q':
        return POLAR_UNCHARGED_AA;

      case 'K': case 'R': case 'H':
        return POSITIVELY_CHARGED_AA;

      case 'E': case 'D':
        return NEGATIVELY_CHARGED_AA;

      case 'A': case 'I': case 'L': case 'M': case 'F': case 'W': case 'V':
      case 'Y':
        return HYDROPHOBIC_AA;

      case 'C': case 'G': case 'P':
        return SPECIAL_AA;

      case '#': case '*': case '+':
        return STOP_AA;

      default:
        return ILLEGAL_AA;
    }
  }


  /**
   *  Return the one letter abbreviation for the given three letter amino acid
   *  name.
   *  @return A one letter code or -1 if three_letter_code can't be understood.
   **/
  public static char getOneLetterCode(final String three_letter_code) 
  {
    final String real_code =
      three_letter_code.substring(0, 1).toUpperCase() +
      three_letter_code.substring(1).toLowerCase();
    
    for(int i = 0 ; i < amino_acid_one_letter_names.length ; ++i) 
    {
      if(real_code.equals(amino_acid_abbreviated_names[i]))
        return amino_acid_one_letter_names[i];
    }

    return (char) -1;
  }

  /**
   *  Return the three letter abbreviation for the given one letter amino acid
   *  code.
   **/
  public static String getThreeLetterAbbreviation(char one_letter_code) 
  {
    one_letter_code = Character.toLowerCase(one_letter_code);
    for(int i = 0 ; i < amino_acid_one_letter_names.length ; ++i) 
    {
      if(one_letter_code == amino_acid_one_letter_names[i]) 
        return amino_acid_abbreviated_names[i];
    }

    throw new Error("internal error - illegal one letter amino acid code");
  }

  /**
   *  Return the three letter abbreviation for the given index code.
   **/
  public static String getThreeLetterAbbreviation(final int index) 
  {
    return amino_acid_abbreviated_names[index];
  }

  /**
   *  Return an integer from 0 to 22 representing the index of a codon
   *  symbol.
   **/
  public static int getSymbolIndex(char one_letter_code) 
  {
    one_letter_code = Character.toLowerCase(one_letter_code);
    switch(one_letter_code) 
    {
      case 'a': return 0;
      case 'r': return 1;
      case 'n': return 2;
      case 'd': return 3;
      case 'c': return 4;
      case 'q': return 5;
      case 'e': return 6;
      case 'g': return 7;
      case 'h': return 8;
      case 'i': return 9;
      case 'l': return 10;
      case 'k': return 11;
      case 'm': return 12;
      case 'f': return 13;
      case 'p': return 14;
      case 's': return 15;
      case 't': return 16;
      case 'w': return 17;
      case 'y': return 18;
      case 'v': return 19;
      case '*': return 20;
      case '#': return 21;
      case '+': return 22;
      case '.': return 23;
      case 'x': return 23;
      case 'u': return 24;
      default:
        throw new Error("Internal error - illegal one letter codon symbol: " +
                        one_letter_code);
    }
  }

  /**
   *  Given an index return a one letter codon symbol.  This method is the
   *  inverse of getSymbolIndex().
   **/
  public static char getSymbolFromIndex(final int index) 
  {
    return amino_acid_one_letter_names[index];
  }

  /**
   *  A String containg the amino acids symbols of this object.
   **/
  private String amino_acid_string = null;

  /**
   *  The three letter abbreviated names for the amino acids and stop codons.
   *  The names here correspond to the letter codes at the same indices in
   *  amino_acid_one_letter_names.
   **/
  private final static String [] amino_acid_abbreviated_names = 
  {
    "Ala", "Arg", "Asn", "Asp", "Cys",
    "Gln", "Glu", "Gly", "His", "Ile",
    "Leu", "Lys", "Met", "Phe", "Pro",
    "Ser", "Thr", "Trp", "Tyr", "Val",
    "Opl", "Ocr", "Amb", "---", "Sel"
  };


  /**
   *  The one letter abbreviated names for the amino acids and stop codons.
   *  The names here correspond to the three letter codes at the same indices
   *  in amino_acid_abbreviated_names.
   **/
  private final static char [] amino_acid_one_letter_names = 
  {
    'a', 'r', 'n', 'd', 'c',
    'q', 'e', 'g', 'h', 'i',
    'l', 'k', 'm', 'f', 'p',
    's', 't', 'w', 'y', 'v',
    '*', '#', '+', '.', 'u'
  };

  /**
   *  The molecular weights of the amino acids.  The values correspond to the
   *  three letter codes at the same indices in amino_acid_abbreviated_names.
   *  For example "Met" corresponds to a weight of 149.22
   **/
  private final static float [] molecular_weights = 
  {
    89.09F,  174.21F, 132.12F, 133.10F, 121.15F,
    146.15F, 147.13F, 75.07F,  155.16F, 131.18F,
    131.18F, 146.19F, 149.22F, 165.19F, 115.13F,
    105.09F, 119.12F, 204.22F, 181.19F, 117.15F,
    0.0F,    0.0F,    0.0F,    0.0F,    334.1F
  };

  /**
   *  The average molecular weight of all amino acids.
   **/
  //private final static float average_molecular_weight = 136.90F;

  /**
   *  The molecular weight of water.
   **/
  private final static float molecular_weight_of_water = 18.015F;

  /**
   *  The number of amino acid symbols, including stop codons and a symbol for
   *  unknown (---): 23.
   **/
  public final static int symbol_count = amino_acid_abbreviated_names.length;

  /**
   *  The number of amino acid symbols, not including stop codons: 20.
   **/
  //private final static int amino_acid_symbol_count = 20;


  public static void setGeneCode()
  {
    // if translation_table is in the options file use it to set
    // codon_translation_array
    final StringVector options_file_table =
             Options.getOptions().getOptionValues("translation_table");

    if(options_file_table != null) 
    {
      if(options_file_table.size () == 64) 
      {
        for(int i = 0 ; i < 64 ; ++i) 
        {
          final char new_table_char = Character.toUpperCase(
                      ((String)(options_file_table.elementAt(i))).charAt(0));

          if(isLegalCodon (new_table_char)) 
            codon_translation_array[i] = new_table_char;
          else 
            codon_translation_array[i] = '.';
        }
      }
    }
  }
}
