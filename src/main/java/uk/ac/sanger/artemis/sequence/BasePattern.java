/* BasePattern.java
 *
 * created: Sun Jan 10 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/sequence/BasePattern.java,v 1.1 2004-06-09 09:52:12 tjc Exp $
 */

package uk.ac.sanger.artemis.sequence;

import uk.ac.sanger.artemis.util.*;

/**
 *  A BasePattern is a String that describes a sequence of bases.
 *
 *  @author Kim Rutherford
 *  @version $Id: BasePattern.java,v 1.1 2004-06-09 09:52:12 tjc Exp $
 **/

public class BasePattern {
  /**
   *  Create a new BasePattern object.
   *  @param pattern_string The is a String representation of the pattern.
   **/
  public BasePattern (final String pattern_string)
      throws BasePatternFormatException {
    if (pattern_string.length () < 1) {
      throw new BasePatternFormatException ("pattern too short");
    }
    this.pattern_string = pattern_string.toLowerCase ();

    this.pattern_type = patternType (this.pattern_string);

    if (pattern_type == ILLEGAL_PATTERN) {
      throw new BasePatternFormatException ("illegal characters in pattern");
    }
  }

  /**
   *  A return value of patternType () - illegal characters in the pattern.
   **/
  public static int ILLEGAL_PATTERN = -1;

  /**
   *  A return value of patternType () - the pattern contains only the
   *  characters 'a', 't', 'g' and 'c'.
   **/
  private static int SIMPLE_PATTERN = 0;

  /**
   *  A return value of patternType () - the pattern contains only the
   *  characters a,t,g,c,r,y,k,m,s,w,n,b,d,h and v, ie IUC base codes.
   **/
  private static int IUC_PATTERN = 1;

  /**
   *  Return a String representation of this BasePattern.
   **/
  public String toString () {
    return pattern_string;
  }

  /**
   *  Returns true if and only if this pattern matches the given String.  The
   *  match must be exact so the pattern and the string must be the same
   *  length.
   **/
  public boolean matches (final String match_string) {
    if (match_string.length () == pattern_string.length () &&
        searchFor (match_string, 0) == 0) {
      return true;
    } else {
      return false;
    }
  }

  /**
   *  Find the next match of this pattern in either Strand of the given
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
  public MarkerRange findMatch (final Bases bases,
                                final Marker search_start_marker,
                                final int search_end_position,
                                final boolean search_backwards,
                                final boolean search_fwd_strand,
                                final boolean search_bwd_strand) {
    final String bases_string = bases.toString ();

    // search the bases_string forward for the pattern_string and its
    // complement

    // the String index position in bases_string at which to start the search
    // for this pattern
    final int forward_search_start_index;

    // the String index position in bases_string at which to start the search
    // for the reverse complement of this position
    final int complement_search_start_index;

    if (search_backwards) {
      if (search_start_marker == null) {
        forward_search_start_index = bases.getLength () - 1;
        complement_search_start_index = bases.getLength () - 1;
      } else {
        complement_search_start_index =
          search_start_marker.getRawPosition () - 2;
        if (search_start_marker.getStrand ().isForwardStrand ()) {
          forward_search_start_index =
            search_start_marker.getRawPosition () - 2;
        } else {
          forward_search_start_index =
            search_start_marker.getRawPosition () - 1;
        }
      }
    } else {
      if (search_start_marker == null) {
        forward_search_start_index = 0;
        complement_search_start_index = 0;
      } else {
        forward_search_start_index = search_start_marker.getRawPosition ();
        if (search_start_marker.getStrand ().isForwardStrand ()) {
          complement_search_start_index =
            search_start_marker.getRawPosition () - 1;
        } else {
          complement_search_start_index =
            search_start_marker.getRawPosition ();
        }
      }
    }
    
    final int forward_search_result;
    if(search_fwd_strand)
      forward_search_result = searchFor (bases_string, pattern_string,
                                         forward_search_start_index, search_backwards);
    else
      forward_search_result = -1;

    final int complement_search_result;
    if(search_bwd_strand)
      complement_search_result = searchFor (bases_string, Bases.reverseComplement (pattern_string),
                                            complement_search_start_index, search_backwards);
    else
      complement_search_result = -1;

    final int match_first_base;
    final int match_last_base;

    final Strand match_strand;

    if (forward_search_result == -1) {
      if (complement_search_result == -1) {
        // no match
        return null;
      }
    }

    if (search_backwards) {
      // take the match that is closest to the end, or the complement match if
      // there is a tie
      if (complement_search_result != -1 &&
          (forward_search_result == -1 ||
           forward_search_result != -1 &&
           complement_search_result >= forward_search_result)) {
        match_first_base =
          bases.getComplementPosition (complement_search_result + 1);
        match_last_base = match_first_base - (pattern_string.length () - 1);
        match_strand = bases.getReverseStrand ();
      } else {
        match_first_base = forward_search_result + 1;
        match_last_base = match_first_base + pattern_string.length () - 1;
        match_strand = bases.getForwardStrand ();
      }
    } else {
      // take the match that is closest to base 1, or the forward match if
      // there is a tie
      if (forward_search_result != -1 &&
          (complement_search_result == -1 ||
           complement_search_result != -1 &&
           forward_search_result <= complement_search_result)) {
        match_first_base = forward_search_result + 1;
        match_last_base = match_first_base + pattern_string.length () - 1;
        match_strand = bases.getForwardStrand ();
      } else {
        match_first_base =
          bases.getComplementPosition (complement_search_result + 1);
        match_last_base = match_first_base - (pattern_string.length () - 1);
        match_strand = bases.getReverseStrand ();
      }
    }

    if (match_last_base > search_end_position) {
      // there is no match within the range
      return null;
    }

    try {
      return new MarkerRange (match_strand,
                              match_first_base,
                              match_last_base);
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Find the all matches of this pattern in either Strand of the given
   *  Bases object.  This method searches both strands simultaneously by
   *  searching the underlying Bases object directly.
   *  @param bases This holds the Strand objects to search.
   *  @param search_start_marker The search will start at this Marker
   *    position.  See the comments on the search_start_marker argument to
   *    findMatch ().
   *  @param search_end_position The search will not extend past this base
   *    position.
   *  @return A MarkerRangeVector holding all the matches or null if there are
   *    no matches in the given range.
   **/
  public MarkerRangeVector findMatches (final Bases bases,
                                        final Marker search_start_marker,
                                        final int search_end_position) {
    final MarkerRangeVector return_vector = new MarkerRangeVector ();

    Marker current_position_marker = search_start_marker;

    while (true) {
      final MarkerRange new_match_position =
        findMatch (bases, current_position_marker, search_end_position, false, true, true);

      if (new_match_position == null) {
        break;
      } else {
        current_position_marker = new_match_position.getRawStart ();
        return_vector.add (new_match_position);
      }
    }

    return return_vector;
  }

  /**
   *  Search for this BasePattern in the given string.
   *  @param bases_string Search this String for the pattern.
   *  @param start_index This is the index in bases_string where the search
   *    should start.
   *  @return The index of the match or -1 if there is no match.
   **/
  private int searchFor (final String bases_string,
                        final int start_index) {
    return searchFor (bases_string, pattern_string, start_index, false);
  }

  /**
   *  Search for a pattern in a string.
   *  @param bases_string Search this String for the pattern.
   *  @param pattern_string This contains the pattern to search for and may
   *    include IUB base codes.
   *  @param start_index This is the index in bases_string where the search
   *    should start.
   *  @param search_backwards If true the search will move from last base to
   *    first base, otherwise first to last.
   *  @return The index of the match or -1 if there is no match.
   **/
  private int searchFor (final String bases_string,
                         final String pattern_string,
                         int start_index,
                         final boolean search_backwards) {

    if (search_backwards) {
      if (pattern_type == SIMPLE_PATTERN) {
        // indexOf () is much faster in this case
        return bases_string.lastIndexOf (pattern_string, start_index);
      }

      if (bases_string.length () - start_index < toString ().length ()) {
        start_index = bases_string.length () - toString ().length ();
      }

      for (int i = start_index ; i > 0 ; --i) {
        boolean match_failed = false;
        for (int pattern_index = 0 ;
             pattern_index < pattern_string.length () ;
             ++pattern_index) {
          if (charMatch (bases_string.charAt (i + pattern_index),
                         pattern_string.charAt (pattern_index))) {
            // OK, so continue with the inner loop
          } else {
            match_failed = true;
            break;
          }
        }

        if (match_failed) {
          // go around main loop again
        } else {
          // found a match
          return i;
        }
      }
    } else {
      if (pattern_type == SIMPLE_PATTERN) {
        // indexOf () is much faster in this case
        return bases_string.indexOf (pattern_string, start_index);
      }

      for (int i = start_index ;
           i < bases_string.length () - pattern_string.length () + 1;
           ++i) {
        boolean match_failed = false;
        for (int pattern_index = 0 ;
             pattern_index < pattern_string.length () ;
             ++pattern_index) {
          
          if (charMatch (bases_string.charAt (i + pattern_index),
                         pattern_string.charAt (pattern_index))) {
            // OK, so continue with the inner loop
          } else {
            match_failed = true;
            break;
          }
        }

        if (match_failed) {
          // go around main loop again
        } else {
          // found a match
          return i;
        }
      }
    }

    return -1;
  }

  /**
   *  Check a base to see if it matches an IUC base code.
   *  @param base_char The base character to be checked.
   *  @param pattern_char The single letter IUC base code to match the
   *    character against.
   **/
  private boolean charMatch (final char base_char,
                             final char pattern_char) {
    switch (base_char) {
    case 'c':
      switch (pattern_char) {
      case 'c':
      case 'y': case 'm': case 's': case 'n': case 'b': case 'h': case 'v':
        return true;
      }
      break;
    case 't':
      switch (pattern_char) {
      case 't':
      case 'y': case 'k': case 'w': case 'n': case 'b': case 'd': case 'h':
        return true;
      }
      break;
    case 'a':
      switch (pattern_char) {
      case 'a':
      case 'r': case 'm': case 'w': case 'n': case 'd': case 'h': case 'v':
        return true;
      }
      break;
    case 'g':
      switch (pattern_char) {
      case 'g':
      case 'r': case 'k': case 's': case 'n': case 'b': case 'd': case 'v':
        return true;
      }
      break;
    default:
      if (pattern_char == 'n') {
        return true;
      }
    }

    return false;
  }

  /**
   *  Return the pattern type of the argument.  Returns ILLEGAL_PATTERN if
   *  there are illegal characters in the pattern.  SIMPLE_PATTERN if the
   *  pattern contains only the characters 'a', 't', 'g' and 'c'.  IUC_PATTERN
   *  if the pattern contains only the characters a,t,g,c,r,y,k,m,s,w,n,b,d,h
   *  and v, ie IUC base codes.
   **/
  public static int patternType (String pattern_string) {
    boolean seen_iuc = false;

    for (int i = 0 ; i < pattern_string.length () ; ++i) {
      switch (pattern_string.charAt (i)) {
      case 'r': case 'y': case 'k': case 'm': case 's':
      case 'w': case 'n': case 'b': case 'd': case 'h': case 'v':
        seen_iuc = true;
        break;
      case 'a': case 't': case 'g': case 'c':
        // no problem
        break;
      default:
        // anything else is illegal
        return ILLEGAL_PATTERN;
      }
    }

    if (seen_iuc) {
      return IUC_PATTERN;
    } else {
      return SIMPLE_PATTERN;
    }
  }

  /**
   *  The pattern that was passed to the constructor.
   **/
  final String pattern_string;

  /**
   *  The type of this pattern, SIMPLE_PATTERN, IUC_PATTERN etc.
   **/
  final int pattern_type;
}


