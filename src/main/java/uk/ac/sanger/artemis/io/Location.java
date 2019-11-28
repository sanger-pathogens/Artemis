/* Location.java
 *
 * created: Mon Oct  5 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/Location.java,v 1.8 2009-02-02 16:05:32 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.io.LocationLexer.TokenEnumeration;

/**
 *  This class encapsulates the location field of a EMBL entry feature.
 *  As well as allowing access to the location string itself, this class has
 *  functions for parsing and manipulating the location.
 *
 *  @author Kim Rutherford
 *  @version $Id: Location.java,v 1.8 2009-02-02 16:05:32 tjc Exp $
 *
 */
public class Location 
{

  /**
   *  The canonical parsed version of the location string that was passed to
   *  the constructor.
   **/
  private LocationParseNode parse_tree = null;

  /**
   *  A cache of ranges, set and returned by getRanges ().
   **/
  private RangeVector ranges = null;

  /**
   *  A cache of the total range, set and returned by getTotalRange ().
   **/
  private Range total_range = null;

  /**
   *  Constructs a new Location object with that represents the value
   *  represented by the string.
   *  @exception LocationParseException Thrown if the location String does not
   *    represent a valid Location.
   **/
  public Location(final String location_string)
      throws LocationParseException 
  {
    parse_tree = getParseTree(location_string).getCanonical();
    
    if(parse_tree == null)
      throw new LocationParseException("invalid location", location_string);
  }

  /**
   *  Constructs a new Location object with that represents the given Range.
   *  @exception OutOfRangeException Thrown if the Range starts below one.
   **/
  public Location(final Range location_range)
      throws OutOfRangeException 
  {
    if(location_range.getStart() < 1) 
      throw new OutOfRangeException("location out of range: " +
                                    location_range.toString ());
    else
      parse_tree = new LocationParseNode(location_range);

    ranges = new RangeVector(location_range);
  }

  /**
   *  Constructs a new Location object with that represents the given Ranges
   *  @param ranges a non-empty vector containing the ranges for the new
   *    Location.
   *  @param complement true if and only if all the ranges should be
   *    complemented in the new location.
   **/
  public Location(final RangeVector ranges, final boolean complement) 
  {
    final LocationParseNodeVector vector = new LocationParseNodeVector();

    if(ranges.size() == 0) 
      throw new Error ("internal error - ranges.size () == 0");

    for(int i = 0 ; i<ranges.size() ; ++i) 
    {
      final LocationParseNode range_node = new LocationParseNode((Range)ranges.elementAt(i));

      if(complement)
      {
        final LocationParseNode complement_node =
          new LocationParseNode(LocationParseNode.COMPLEMENT, range_node);
        vector.addElement(complement_node);
      }
      else 
        vector.addElement (range_node);
    }

    if(vector.size() == 1) 
      parse_tree = vector.elementAt (0);
    else 
      parse_tree = new LocationParseNode(LocationParseNode.JOIN, vector);
  }

  /**
   *  This method translates the start and end of each Range in this Location
   *  into another coordinate system.  The Ranges will be truncated if
   *  necessary.
   *  @param constraint This contains the start and end base of the new
   *    coordinate system.  The position given by constraint.getStart () will
   *    be at postion/base 1 in the new coordinate system.
   *  @return the Location translated into the new coordinate system.  The
   *    Ranges in the Location will be truncated if thwy overlap either end of
   *    the constraint. null is returned if and only if all the Ranges lie
   *    outside of the new start and end.
   **/
  public Location truncate(final Range constraint) 
  {
    final boolean is_complement = isComplement();
    final RangeVector ranges    = getRanges();

    final RangeVector new_ranges = new RangeVector();

    for(int i = 0 ; i<ranges.size () ; ++i) 
    {
      final Range truncated_range =
        ((Range)ranges.elementAt(i)).truncate(constraint);

      if(truncated_range != null) 
        new_ranges.add(truncated_range);
    }

    if(new_ranges.size() > 0) 
      return new Location (new_ranges, is_complement);
    else 
      return null;
  }

  /**
   *  Returns the value of this Location as a string.  The value returned will
   *  be the same as string passed to the constructor (but will always look
   *  like join(complement(... if there is a complement and a join).
   **/
  public String toString() 
  {
    return getParsedLocation().toString();
  }

  /**
   *  Return true if and only if this Location is the same as the given
   *  Location.
   **/
  public boolean equals(final Location other_location) 
  {
    return toString().equals(other_location.toString());
  }


  /**
  *
  * Test whether this is a trans-spliced feature location.
  *
  */
  private boolean isTransSpliced()
  {
    final LocationParseNode top = getParsedLocation();
    if(top.getType () != LocationParseNode.JOIN)
      return false;

    final LocationParseNodeVector children = 
               top.getJoinChildren();
    
    int num_complement = 0;

    for(int i = children.size() - 1; i >= 0; --i)
    {
      final LocationParseNode child_node = children.elementAt(i);
      if(child_node.getType() == LocationParseNode.COMPLEMENT)
        num_complement++;
    }
   
    if(num_complement == 0 || num_complement == children.size())
      return false;

    return true;
  }

  /**
   * If lookAt5prime is set to true then only return true if the 5' end is 
   * partial otherwise only return true if the 3' end is partial.
   * @param lookAt5prime
   * @return
   */
  public boolean isPartial(final boolean lookAt5prime)
  {
    try
    {
      LocationParseNode top = getParsedLocation();
      if (top.getType() == LocationParseNode.COMPLEMENT)
          top = top.getComplementChild();

      if (top.getType() == LocationParseNode.JOIN || 
          top.getType() == LocationParseNode.ORDER)
      {
        LocationParseNodeVector nodes = top.getChildren();
        LocationParseNode node = nodes.elementAt( 
            (lookAt5prime ? 0 : nodes.size()-1));
        if (node.getType() == LocationParseNode.COMPLEMENT)
        {
          node = node.getComplementChild();
          if (node.getType() == LocationParseNode.RANGE
              && node.getRange() instanceof FuzzyRange)
            top = node;
        }
        else
          top = nodes.elementAt(0);
      }

      if (top.getType() == LocationParseNode.RANGE && 
          top.getRange() instanceof FuzzyRange)
      {
        final FuzzyRange fuzz = (FuzzyRange) top.getRange();
        if (fuzz.getStartObject() instanceof LowerInteger && (isComplement() ^ lookAt5prime))
          return true;
        if (fuzz.getEndObject() instanceof UpperInteger && !(isComplement() ^ lookAt5prime))
          return true;
      }
    }
    catch (Exception e){}
    return false;
  }

  
  /**
   *  Returns the value of this Location as a string with the complement (if
   *  any) as the top node and the join or order (if any) below it.
   *  eg. "complement(join(..." not "join(complement(..."
   **/
  public String toStringShort() 
  {
    final LocationParseNode top = getParsedLocation();

    // cope with trans-spliced locations
    if(isTransSpliced())
      return toString();

    if((top.getType () == LocationParseNode.JOIN &&
        top.getJoinChildren ().elementAt (0).getType () ==
        LocationParseNode.COMPLEMENT) ||
       (top.getType () == LocationParseNode.ORDER &&
        top.getOrderChildren ().elementAt (0).getType () ==
        LocationParseNode.COMPLEMENT)) 
    {
      final StringBuffer return_buffer;
      final LocationParseNodeVector children;

      if (top.getType () == LocationParseNode.JOIN) 
      {
        return_buffer = new StringBuffer("complement(join(");
        children = top.getJoinChildren ();
      } 
      else 
      {
        return_buffer = new StringBuffer("complement(order(");
        children = top.getOrderChildren();
      }


      // reverse the nodes because in a join(complement( location the
      // exons a in reverse order
      for(int i = children.size() - 1; i >= 0; --i) 
      {
        final LocationParseNode child_node = children.elementAt(i);

        if(child_node.getType() == LocationParseNode.COMPLEMENT)
        {
          final LocationParseNode complement_child =
            child_node.getComplementChild ();

          return_buffer.append (complement_child.toString ());
        }
        else   // probably an ENTRY_RANGE
          return_buffer.append (child_node.toString ());

        if(i != 0) 
          return_buffer.append(',');
      }

      return_buffer.append("))");

      return return_buffer.toString();
    } 
    else 
    {
      // does not contain a join and a complement, so there is no need to swap
      // them
      return toString ();
    }
  }

  /**
   *  This method returns true if and only if this Location is a complement.
   **/
  public boolean isComplement() 
  {
	int type = getParsedLocation ().getType ();
	
    if(type == LocationParseNode.COMPLEMENT) 
    {
      // COMPLEMENT(RANGE)
      return true;
    }
    else 
    {
      if(type == LocationParseNode.JOIN || type == LocationParseNode.ORDER) 
      {
        // a join/order must have at least one child node
        if(getParsedLocation().getChildren().elementAt(0).getType() ==
            LocationParseNode.COMPLEMENT) 
        {
          // JOIN(COMPLEMENT(RANGE)...) or ORDER(COMPLEMENT(RANGE)...)
          return true;
        }
        else 
        {
          // JOIN(RANGE,...) or ORDER(RANGE,...)
          return false;
        }
      }
      else 
      {
        // the head node is not a JOIN, COMPLEMENT or ORDER so it must be a
        // RANGE
        return false;
      }
    }
  }

  /**
   *  Return the reference of a new copy of this Location.
   **/
  public Location copy()
  {
    return new Location(getParsedLocation().copyClean());
  }

  /**
   *  Change the given range in this Location.  If there are two or more
   *  identical ranges in the location the first one will be changed.
   *  @param old_range The Range to change
   *  @param new_range The Range to replace old_range with.
   *  @return a new Location that will be a copy of this Location with the
   *    given Range changed.
   **/
  public Location changeRange(final Range old_range, final Range new_range) 
  {
    final Location new_location = new Location(getParsedLocation().copy());
    new_location.getParsedLocation().changeRange(old_range, new_range);

    return new_location;
  }

  /**
   *  Add the given Range in this Location.
   *  @return a new Location that will be a copy of this Location with the
   *    given Range added.
   **/
  public Location addRange(final Range new_range) 
  {
    final Location new_location = new Location(getParsedLocation().copy());
    final LocationParseNode new_node = new LocationParseNode(new_range);

    new_location.parse_tree =
      new_location.getParsedLocation().addRangeNode(new_node);

    return new_location;
  }

  /**
   *  Remove the given Range in this Location.  If there are two or more
   *  Ranges in this Location that are identical to remove_range the first
   *  will be removed.
   *  @return a new Location that will be a copy of this Location with the
   *    given LocationParseNode removed.
   **/
  public Location removeRange(final Range remove_range) 
  {
    final Location new_location = new Location(getParsedLocation().copy());

    new_location.parse_tree =
      new_location.getParsedLocation().removeRange(remove_range);

    return new_location;
  }

  /**
   *  Return a reversed and complemented copy of this Location.
   *  @param sequence_length The length of the sequence that this Location is
   *    associated with.
   *  @return a new Location that will be a reversed and complemented copy of
   *    this Location.  The user_data fields in the returned Location will be
   *    the same as the original.
   **/
  public Location reverseComplement(final int sequence_length) 
  {
    return reverseComplement(sequence_length, 1);
  }

   /**
   *  Return a reversed and complemented copy of this Location.
   *  @param sequence_length The length of the sequence that this Location is
   *    associated with.
   *  @return a new Location that will be a reversed and complemented copy of
   *    this Location.  The user_data fields in the returned Location will be
   *    the same as the original.
   **/
  public Location reverseComplement(final int sequence_length, 
                                    final int offset)
  {
    final Location new_location = new Location(getParsedLocation().copy());

    new_location.parse_tree =
      new_location.getParsedLocation().reverseComplement(sequence_length, offset);

    new_location.parse_tree = new_location.parse_tree.getCanonical();
    return new_location;
  }


  /**
   *  Return a Range that spans the whole Location.
   **/
  public Range getTotalRange() 
  {
    try
    {
      if(total_range == null) 
      {
        final RangeVector ranges = getRanges();

        int lowest_so_far  = -1;
        int highest_so_far = -1;
        
        for(int i = 0 ; i<ranges.size() ; ++i) 
        {
          final int this_start = ((Range)ranges.elementAt(i)).getStart();
          final int this_end   = ((Range)ranges.elementAt(i)).getEnd();
          
          if(lowest_so_far == -1 || this_start < lowest_so_far) 
            lowest_so_far = this_start;
          
          if(this_end > highest_so_far) 
            highest_so_far = this_end;
        }

        total_range = new Range(lowest_so_far, highest_so_far);
      }

      return total_range;
    } 
    catch (OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Return a Vector containing all the Ranges of this Location (which could
   *  be none).
   **/
  public RangeVector getRanges()
  {
    if(ranges != null)
      return ranges;

    ranges = new RangeVector();

    final LocationParseNode parsed_location;

    parsed_location = getParsedLocation();

    // this will contain nodes that will either be COMPLEMENT(RANGE) or
    // RANGE
    LocationParseNodeVector range_nodes;

    // strip out the JOIN or ORDER node if any - must be the top node because
    // we called getCanonical ()
    if(parsed_location.getType() == LocationParseNode.JOIN) 
      range_nodes = parsed_location.getJoinChildren ();
    else 
    {
      if(parsed_location.getType() == LocationParseNode.ORDER) 
        range_nodes = parsed_location.getOrderChildren();
      else 
      {
        range_nodes = new LocationParseNodeVector();
        range_nodes.addElement(parsed_location);
      }
    }

    // loop over complement_nodes extracting the Range object the nodes.
    for(int i = 0 ; i<range_nodes.size() ; ++i) 
    {
      final LocationParseNode current_node;
      current_node = range_nodes.elementAt(i);

      final Range current_node_range;

      if(current_node.getType() == LocationParseNode.COMPLEMENT) 
      {
        final LocationParseNode child = current_node.getComplementChild();
        if(child.getType() == LocationParseNode.RANGE) 
          current_node_range = child.getRange ();
       else 
       {
          // probably an ENTRY_RANGE node
          continue;
        }
      }
      else 
      {
        if(current_node.getType() == LocationParseNode.RANGE) 
          current_node_range = current_node.getRange();
        else 
        {
          // probably an ENTRY_RANGE node
          continue;
        }
      }

      ranges.add(current_node_range);
    }

    return ranges;
  }

  /**
   *  Return the first base that this Location covers or -1 if the location
   *  does not refer to any bases in this sequence.
   **/
  public int getFirstBase()
  {
    return getTotalRange().getStart();
  }

  /**
   *  Return the last base that this Location covers or -1 if the location
   *  does not refer to any bases in this sequence.
   **/
  public int getLastBase()
  {
    return getTotalRange().getEnd();
  }

  /**
   *  Return a complemented copy of this Location.  If this Location is
   *  already a complement then a non complemented copy is returned.
   **/
  public Location getComplement() 
  {
    return new Location(getParsedLocation().getNodeComplement());
  }

  /**
   *  Returns a parse tree for this Location.  Changes to the parse tree will
   *  change the original Location.
   *  @exception LocationParseException Thrown if a parse error occurs.
   **/
  public LocationParseNode getParsedLocation()
  {
    return parse_tree;
  }

  public void setParsedLocation(LocationParseNode parse_tree)
  {
    this.parse_tree = parse_tree;
  }


  /**
   *  Create a new Location object from the given LocationParseNode tree.
   **/
  private Location(LocationParseNode parse_tree) 
  {
    this.parse_tree = parse_tree;
  }

  /**
   *  Returns a parse tree for this Location.  Changes to the parse tree will
   *  not effect the original Location.
   *  @exception LocationParseException Thrown if a parse error occurs.
   */
  private LocationParseNode getParseTree(String location_string)
      throws LocationParseException 
  {
    final LocationLexer lexer     = new LocationLexer(location_string);
    final TokenEnumeration enumTk = lexer.getTokens();
    final LocationParseNode join  = parseJoinContents(enumTk);

    if(enumTk.peekElement() != null) 
      throw new LocationParseException("garbage at the end of the " +
                                       "location string", enumTk);
    return join;
  }

  /**
   *  Parse and return a location in the form of a LocationParseNode.  A
   *  location is of the form: complement(...), join(...), order(...) or
   *  <range> (where <range> is described in parseRange()).  If the next
   *  tokens in the TokenEnumeration are not a location then an exception is
   *  thrown.  */
  private static LocationParseNode parseLocation(TokenEnumeration enumTk)
      throws LocationParseException 
  {
    // failed to parse a range - try parsing a functional(join, etc.)

    final LocationParseNode parsed_functional = parseFunctional(enumTk);

    if(parsed_functional != null)
       return parsed_functional;

    final LocationParseNode entry_range = parseEntryRange(enumTk);

    if(entry_range != null)
      return entry_range;

    final LocationParseNode parsed_range = parseRange(enumTk);

    if(parsed_range != null) 
       return parsed_range;
    else 
      throw new LocationParseException("expected a range or a functional",
                                       enumTk);
  }

  /**
   *  Attempt to parse and return a Range with an entry label in front (like
   *  J00194:(100..202) or J00193:hladr).  If the next tokens do not look like
   *  an entry range then it immediately returns null, otherwise an attempt is
   *  made to parse the next tokens as an entry range.  If the parse fails
   *  part way through then a LocationParseException is thrown.
   **/
  private static LocationParseNode parseEntryRange(TokenEnumeration enumTk)
      throws LocationParseException 
  {
    if(!(enumTk.peekElement() instanceof String)) 
      return null;

    final String entry_string = (String)enumTk.nextElement();

    if(!enumTk.eatToken(':'))
      throw new LocationParseException("parse error after reading \"" +
                                        entry_string + "\"", enumTk);

    final LocationParseNode entry_location = parseLocation(enumTk);

     return new LocationParseNode(entry_string, entry_location);
  }

  /**
   *  Attempt to parse and return the contents of a JOIN node (or the top
   *  level node).
   **/
  private static LocationParseNode parseJoinContents(TokenEnumeration enumTk)
      throws LocationParseException 
  {
    final LocationParseNodeVector return_ranges =
      new LocationParseNodeVector();

    final LocationParseNode parsed_range = parseLocation(enumTk);

    if(parsed_range == null)
      return null;

    return_ranges.addElement(parsed_range);

    final Object peek_element = enumTk.peekElement();

    while(enumTk.eatToken(',')) 
    {
      LocationParseNode new_child = parseLocation(enumTk);
      return_ranges.addElement(new_child);
    }

    if(return_ranges.size() > 1) 
       return new LocationParseNode(LocationParseNode.JOIN, return_ranges);
    else 
       return return_ranges.elementAt(0);
  }

  /**
   *  Attempt to parse and return a Range.  If the next tokens do not look
   *  like a range then it immediately returns null, otherwise an attempt is
   *  made to read a range.  If the parse fails part way through what looks
   *  like a range then a LocationParseException is thrown.
   **/
  private static LocationParseNode parseRange(TokenEnumeration enumTk)
      throws LocationParseException 
  {
    boolean is_lower_unbounded_range = false;
    boolean is_upper_unbounded_range = false;

    // returns Integer, UpperInteger, LowerInteger, Range or null
    Object start_object = parseRangeBound(enumTk);

    if(start_object == null) 
      return null;

    final Object middle_token = enumTk.peekElement();

    if(!enumTk.eatToken("..") && !enumTk.eatToken('^')) 
       return new LocationParseNode(FuzzyRange.makeRange(start_object));

    if(enumTk.peekElement() == null) 
      throw new LocationParseException("location ends in the middle of " +
                                        "a range", enumTk);

    if(start_object instanceof UpperInteger) 
      throw new LocationParseException("range cannot start with: " +
                                       start_object, enumTk);

    Object end_object = parseRangeBound(enumTk);

    if(end_object == null) 
      throw new LocationParseException("unexpected characters in location",
                                       enumTk);

    if(end_object instanceof LowerInteger) 
      throw new LocationParseException("a range cannot end with: " +
                                       end_object, enumTk);

    if(middle_token instanceof Character &&
       ((Character)middle_token).charValue() == '^') 
    {
      if(start_object instanceof Integer && end_object instanceof Integer) 
      {
        try 
        {
          final Range new_range =
            new BetweenRange(((Integer) start_object).intValue(),
                             ((Integer) end_object).intValue());
           return new LocationParseNode(new_range);
        }
        catch(OutOfRangeException e) 
        {
          throw new LocationParseException("a range must start before it " +
                                           "ends", enumTk);
        }
      } 
      else 
        throw new LocationParseException("a range that contains a " +
                                          "'^' must start and end with a " +
                                          "plain number", enumTk);
    }

    try
    {
      return new LocationParseNode(FuzzyRange.makeRange(start_object,
                                                        end_object));
    }
    catch(OutOfRangeException e) 
    {
      throw new LocationParseException("a range must start before it ends",
                                        enumTk);
    }
  }

  /**
   *  Parse the start or end of a range.
   **/
  private static Object parseRangeBound(TokenEnumeration enumTk) throws
      LocationParseException 
  {
    final Object peek_element = enumTk.peekElement();
    
    if(peek_element instanceof Integer)
    {
      if(((Integer) peek_element).intValue() < 1) 
        throw new LocationParseException("range bounds must be " +
                                          "greater than 0", enumTk);
      return enumTk.nextElement();
    }

    if(peek_element instanceof LowerInteger) 
    {
      if(((LowerInteger) peek_element).getPosition() < 1) 
        throw new LocationParseException("range bounds must be " +
                                          "greater than 0", enumTk);

      return enumTk.nextElement();
    }

    if(peek_element instanceof UpperInteger)
    {
      if(((UpperInteger)peek_element).getPosition() < 1)
        throw new LocationParseException("range bounds must be " +
                                          "greater than 0", enumTk);
      return enumTk.nextElement();
    }

    // try to parse a range of the form (100.200)
    if(enumTk.eatToken ('(')) 
    {
      if(enumTk.peekElement() instanceof Integer) 
      {
        final Integer start = (Integer)enumTk.nextElement();

        if(enumTk.eatToken('.')) 
        {
          if(enumTk.peekElement() instanceof Integer) 
          {
            final Integer end = (Integer)enumTk.nextElement();

            int start_value = start.intValue();
            int end_value   = end.intValue();

            // make sure the ranges are the correct way around
            if(start_value > end_value) 
            {
              final int tmp = start_value;
              start_value = end_value;
              end_value = tmp;
            }

            final Object return_object;

            if(start_value < 1) 
              throw new LocationParseException("range bounds must be " +
                                               "greater than 0", enumTk);

            if(start_value == end_value) 
              return_object = new Integer(start_value);
            else 
            {
              try
              {
                return_object = new Range(start_value, end_value);
              } 
              catch(OutOfRangeException e) 
              {
                throw new Error("internal error - unexpected exception: " +e);
              }
            }

            if(enumTk.eatToken(')')) 
               return return_object;
            else
            {
              final String message =
                "expected a closing parenthesis after reading: (" +
                start.intValue() + "." + end.intValue();
              throw new LocationParseException(message, enumTk);
            }
          }
          else 
            throw new LocationParseException("expected an integer", enumTk);
        } 
        else 
          throw new LocationParseException("expected a '.'", enumTk);
      } 
      else 
        throw new LocationParseException("expected an integer", enumTk);
    }

    return null;
  }

  /**
   *  Attempt to parse and return a functional.  If the next tokens do not
   *  look like a functional then it immediately returns null, otherwise an
   *  attempt is made to read a functional.  If the parse fails part way
   *  through what looks like a functional then a LocationParseException is
   *  thrown.
   */
  private static LocationParseNode parseFunctional(TokenEnumeration enumTk)
      throws LocationParseException 
  {
    Object peeked_element = enumTk.peekElement();

    if(peeked_element instanceof String &&
       ("complement".equals((String) peeked_element) ||
        "join".equals((String) peeked_element) ||
        "order".equals((String) peeked_element))) 
    {
      // remove the first token from the enumTk
      String functional_name = (String)enumTk.nextElement();

      if(!enumTk.eatToken('(')) 
        throw new LocationParseException("expected(", enumTk);

      if("complement".equals(functional_name)) 
      {
        // special case for complement - it should have only one argument
        LocationParseNode child = parseJoinContents(enumTk);

        if(!enumTk.eatToken(')')) {
//          throw new LocationParseException ("expected )", enumTk);
        }

        // create a COMPLEMENT type node
         return new LocationParseNode(LocationParseNode.COMPLEMENT, child);
      } 
      else
      {
        // parseLocation() will parse range,range,range as a JOIN
        LocationParseNode join = parseJoinContents(enumTk);

        if(!enumTk.eatToken(')')) {
//          throw new LocationParseException("expected )", enumTk);
        }

        if(join.getType() == LocationParseNode.JOIN) 
        {
          // create a JOIN or ORDER type node
          if("join".equals(functional_name)) 
            return new LocationParseNode(LocationParseNode.JOIN,
                                         join.getChildren());
          else 
            return new LocationParseNode(LocationParseNode.ORDER,
                                          join.getChildren());
        }
        else 
          return join;
      }
    } 
    else  // not a functional
      return null;
  }

  /**
  *
  * Given a range check which strand it is on.
  *
  */
  public boolean isComplement(Range range)
  {
    final LocationParseNode parse_tree = getParsedLocation();
 
    if(parse_tree.getType() != LocationParseNode.JOIN)
      return isComplement();

    final LocationParseNodeVector children =
                                parse_tree.getJoinChildren();

    for(int i = children.size() - 1; i >= 0; --i)
    {
      LocationParseNode child_node = children.elementAt(i);

      Range child_range =  null;
      if(child_node.getType() == LocationParseNode.COMPLEMENT)
      {
      	if(child_node.getComplementChild().getType() == LocationParseNode.RANGE)
         	child_range = child_node.getComplementChild().getRange();
      }
      else if(child_node.getType() == LocationParseNode.RANGE)
        child_range = child_node.getRange();

      if((child_node.getType() == LocationParseNode.RANGE ||
          child_node.getType() == LocationParseNode.COMPLEMENT) &&
          (child_range != null && child_range.equals(range)))
      {
        if(child_node.getType() == LocationParseNode.COMPLEMENT)
          return true;
        else
          return false;
      }
    }

    return isComplement();
  }
}
