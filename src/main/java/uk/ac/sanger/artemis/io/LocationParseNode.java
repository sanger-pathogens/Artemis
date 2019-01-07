/* LocationParseNode.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/LocationParseNode.java,v 1.7 2005-11-15 14:02:37 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import javax.swing.JOptionPane;
import uk.ac.sanger.artemis.util.OutOfRangeException;

/**
 *  LocationParseNode class represents one node in the parse tree tree of a
 *  Location.  It is a utility class for EMBL.Location.
 *
 *  @author Kim Rutherford
 *  @version $Id: LocationParseNode.java,v 1.7 2005-11-15 14:02:37 tjc Exp $
 *
 **/

class LocationParseNode extends EMBLObject
{

  /** Node type for an unknown or error node. */
  final static public int UNKNOWN     = 0;

  /** Node type for complement (...)          */
  final static public int COMPLEMENT  = 1;

  /** Node type for join (...)                */
  final static public int JOIN        = 2;

  /** Node type for order (...)               */
  final static public int ORDER       = 3;

  /** Node type for number..number            */
  final static public int RANGE       = 4;

  /** Node type for a location reference to another entry */
  final static public int ENTRY_RANGE = 5;

  /**
   *  The node type can be COMPLEMENT, JOIN, ORDER, RANGE or ENTRY_RANGE
   **/
  private int node_type;

  /**
   *  If this LocationParseNode has a node type of JOIN or ORDER this member
   *  will point to an LocationParseNode object containing the children of
   *  this node.  If the node type is COMPLEMENT this member will point to an
   *  LocationParseNode object containing the child node of the complement.
   *  If the node type is RANGE this member will point to a Range object.
   *  If the node type is ENTRY_RANGE this will point to a LocationParseNode
   *  of type RANGE and the entry_name member this object will be set.
   **/
  private Object child;

  /**
   *  If this LocationParseNode has a node type of ENTRY_RANGE this member
   *  will give the name of the entry, otherwise it will be null.
   **/
  private String entry_name;

  /**
   *  Create a new LocationParseNode object with a type given by the first
   *  argument and sub nodes given by the second.  The type should be JOIN or
   *  ORDER.
   **/
  LocationParseNode(final int node_type,
                    final LocationParseNodeVector children) 
  {
    if(!(node_type == JOIN || node_type == ORDER)) 
      throw new Error("LocationParseNode constructor was called with the wrong type");

    this.node_type = node_type;
    this.child = children;

    if(children.size() < 1) 
      throw new Error("A functional must have at least one argument");
  }

  /**
   *  Create a new LocationParseNode object of type COMPLEMENT with child node
   *  given by argument.
   *  @param type must be COMPLEMENT
   **/
  LocationParseNode(final int type,
                    final LocationParseNode child) 
  {
    if(type != COMPLEMENT) 
      throw new Error("LocationParseNode constructor was called with the wrong type");

    this.node_type = COMPLEMENT;
    this.child = child;
  }

  /**
   *  Create a new LocationParseNode object of type ENTRY_RANGE with the given
   *  node as the child.
   **/
  LocationParseNode(final String entry_name, final LocationParseNode node) 
  {
    this.node_type  = ENTRY_RANGE;
    this.child      = node;
    this.entry_name = entry_name;
  }

  /**
   *  Create a new LocationParseNode object of type RANGE with range
   *  given by argument.
   **/
  public LocationParseNode(final Range range) 
  {
    this.node_type = RANGE;
    this.child = range;
  }

  /**
   *  Return the type of this node.
   **/
  public int getType() 
  {
    return node_type;
  }

  /**
   *  Returns the value of this LocationParseNode and it's children as a
   *  string.  The method will call toString() on the children of this node
   *  and return those strings as part of the result.
   **/
  public String toString() 
  {
    final String return_string;
    switch(getType())
    {
        case COMPLEMENT:
        return_string = "complement(" + toStringChild() + ")";
        break;
      case JOIN:
        return_string = "join(" + toStringChildren() + ")";
        break;
      case ORDER:
        return_string = "order(" + toStringChildren() + ")";
        break;
      case RANGE:
        return_string = child.toString();
        break;
      case ENTRY_RANGE:
        return_string = entry_name + ":" + toStringChild();
        break;

      default:
        return_string = null;
    }

    return return_string;
  }

  /**
   *  Find the node containing old_range and replace the Range with new_range.
   *  @return true if the node was successfully replaced.
   **/
  boolean changeRange(final Range old_range,
                      final Range new_range) 
  {
    switch(getType())
    {
      case ENTRY_RANGE:
      case COMPLEMENT:
        return ((LocationParseNode) child).changeRange(old_range, new_range);
      case JOIN:
      case ORDER:
        final LocationParseNodeVector children = getChildren();

        for(int i = 0; i < children.size(); ++i) 
        {
          final LocationParseNode this_child = children.elementAt(i);
  
          final boolean return_value =
            this_child.changeRange(old_range, new_range);
  
          if(return_value) 
            return true;
        }
        return false;
      case RANGE:
        if(getRange().equals(old_range))
        {
          setRange(new_range);
          return true;
        } 
        else 
          return false;
      default:
        return false;
    }
  }

  /**
   *  Make a JOIN node from the given RANGE node.
   **/
  private LocationParseNode makeJoinFromRange(final LocationParseNode node) 
  {
    final LocationParseNodeVector children = new LocationParseNodeVector();
    children.addElement(node);

    return new LocationParseNode(JOIN, children);
  }

  /**
   *  Add a RANGE node to the current JOIN or ORDER node.
   **/
  private void addRangeToJoin(final LocationParseNode new_node) 
  {
    final LocationParseNodeVector children = getChildren();
    final int new_range_start = new_node.getRange().getStart();

    if(children.elementAt(0).getType() == COMPLEMENT)
    {
      final LocationParseNode complement_node =
        new LocationParseNode(COMPLEMENT, new_node);

      for(int i = children.size() - 1 ; i >= 0 ; --i) 
      {
        final LocationParseNode this_child = children.elementAt(i);

        if(this_child.getType() != COMPLEMENT) 
          continue;

        final LocationParseNode complement_child =
          this_child.getComplementChild();

        if(complement_child.getType() == RANGE)
        {
          final Range child_range = complement_child.getRange();

          if(child_range.getStart() > new_range_start) 
          {
            children.insertElementAt(complement_node, i + 1);
            return;
          }
        }
      }
      children.insertElementAt(complement_node, 0);
    }
    else
    {
      for(int i = 0 ; i < children.size() ; ++i) 
      {
        final LocationParseNode this_child = children.elementAt(i);

        if(this_child.getType() == RANGE &&
           this_child.getRange().getStart() > new_range_start) 
        {
          children.insertElementAt(new_node, i);
          return;
        }

      }
      children.addElementAtEnd(new_node);
    }
  }

  /**
   *  Add new_node to this LocationParseNode (or its children).  new_node must
   *  be a RANGE node.
   *  @return a changed tree if the node was successfully added or the old
   *    tree otherwise.
   **/
  LocationParseNode addRangeNode(final LocationParseNode new_node) 
  {
    switch(getType())
    {
      case ENTRY_RANGE:
        return this;
      case RANGE:
      case COMPLEMENT:
        return makeJoinFromRange(this).addRangeNode(new_node);
      case JOIN:
      case ORDER:
        addRangeToJoin(new_node);
        return this;
      default:
        throw new Error("internal error - unknown location node type");
    }
  }

  /**
   *  Remove a RANGE node to the current JOIN or ORDER node.
   *  @return a new COMPLEMENT or RANGE node if this node currently has two
   *    children.  returns this node otherwise.
   **/
  private LocationParseNode removeRangeFromJoin(final Range remove_range) 
  {
    final LocationParseNodeVector children = getChildren();

    if(children.size() > 2) 
    {
      for(int i = 0 ; i < children.size() ; ++i) 
      {
        final LocationParseNode this_child = children.elementAt(i);

        if(this_child.getType() == COMPLEMENT) 
        {
          final LocationParseNode complement_child =
            this_child.getComplementChild();
          if(complement_child.getType() == RANGE &&
             complement_child.getRange().equals(remove_range)) 
          {
            children.removeElement(this_child);
            return this;
          }
        } 
        else 
        {
          if(this_child.getType() == RANGE &&
              this_child.getRange().equals(remove_range)) 
          {
            children.removeElement(this_child);
            return this;
          }
        }
      }
    } 
    else
    {
      for(int i = 0 ; i < children.size() ; ++i) 
      {
        final LocationParseNode this_child = children.elementAt(i);

        if(this_child.getType() == COMPLEMENT) 
        {
          final LocationParseNode complement_child =
            this_child.getComplementChild();
          if(complement_child.getType() == RANGE &&
              complement_child.getRange().equals(remove_range)) 
          {
            children.removeElement(this_child);
            return children.elementAt(0);
          }
        } 
        else
        {
          if(this_child.getType() == RANGE &&
              this_child.getRange().equals(remove_range)) 
          {
            children.removeElement(this_child);
            return children.elementAt(0);
          }
        }
      }
    }
    return this;
  }

  /**
   *  Remove the node containing remove_range from this LocationParseNode (or
   *  its children).
   *  @return a changed tree if the node was successfully removed or the old
   *    tree otherwise.
   **/
  LocationParseNode removeRange(final Range remove_range) 
  {
    switch(getType()) 
    {
      case ENTRY_RANGE:
        return this;
      case RANGE:
      case COMPLEMENT:
        throw new Error("internal error - inconsistent location");
      case JOIN:
      case ORDER:
        return removeRangeFromJoin(remove_range);
      default:
        throw new Error("internal error - unknown location node type");
    }
  }

  /**
   *  Return a reversed and complemented copy of this Location.
   *  @param sequence_length The length of the sequence that this Location is
   *    associated with.
   *  @param offset this is set to zero if the whole sequence is being
   *    operated on or to the start of the region being reverse complemented.
   *  @return a reversed and complemented tree.
   **/
  LocationParseNode reverseComplement(final int sequence_length,
                                      final int offset) 
  {
    try
    {
      switch(getType())
      {
        case ENTRY_RANGE:
          return this;
        case RANGE:
        {
          final Range range = getRange();
          final int start = sequence_length - (range.getEnd() - offset + 1) + offset;
          final int end   = sequence_length - (range.getStart() - offset + 1) + offset;

//        System.out.println("LocationParseNode.reverseComplement() HERE "+ 
//                  start+ ".."+ end +"   sequence_length="+sequence_length+
//                 "   range.getEnd()="+range.getEnd()+ "   range.getStart()="+range.getStart()+
//                 "   offset="+offset+ " "+
//                 start_old+ ".."+end_old);

          final Range new_range = new Range(start, end);
          final LocationParseNode new_range_node =
            new LocationParseNode(new_range);
          return new LocationParseNode(COMPLEMENT, new_range_node);
        }
        case COMPLEMENT:
        {
          final LocationParseNode child = getComplementChild();

          if(child.getType() == RANGE) 
          {
            final Range range = child.getRange();
            final int start = sequence_length - (range.getEnd() - offset + 1) + offset;
            final int end   = sequence_length - (range.getStart() - offset + 1) + offset;

            final Range new_range = new Range(start, end);
//            new Range(sequence_length - (range.getEnd() - offset + 1) + offset,
//                      sequence_length - (range.getStart() - offset + 1) + offset);
            return new LocationParseNode(new_range);
          }
          else
            return child;
        }
        case JOIN:
        case ORDER:
        {
          final LocationParseNodeVector children = getChildren();

          final LocationParseNodeVector new_children =
            new LocationParseNodeVector();

          for(int i = 0 ; i < children.size() ; ++i) 
          {
            final LocationParseNode this_child = children.elementAt(i);

            final LocationParseNode new_child =
              this_child.reverseComplement(sequence_length, offset);
            new_children.addElementAtEnd(new_child);
          }

          child = new_children;

          return this;
        }
        default:
          throw new Error("internal error - unknown location node type: " +
                         getType());
      }
    } 
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Returns a "canonical" version of the parse tree starting at the current
   *  node.  In this context a canonical location is one where any join
   *  or order functional is first in the location string.  In a canonical
   *  parse tree only the head node can be a JOIN or ORDER node.  Also in a
   *  "join(complement(..." location the exons must be listed in high to low
   *  order ie. (30..40,10..20) not (10..20,30..40).  Here is an example: <p>
   *  The canonical version of this string: <p>
   *  "complement(join(100..200,400..500))" is: <p>
   *  "join(complement(400..500),complement(100..200))" <p>
   *  The canonical version is less readable but is more standard.
   *  @return A canonical tree or null if the location is nonsense, for
   *    example "join(complement(join(..."
   **/
  LocationParseNode getCanonical() 
  {
    if(!isValid()) 
      return null;

    // this is a new parse tree (a copy of tree starting at this node) with
    // the JOIN or ORDER node (if any) at the top.
    final LocationParseNode new_tree = getCanonicalComplement();

    if(new_tree.getType() != JOIN && new_tree.getType() != ORDER) 
      return new_tree;

    // now put the children of the JOIN or ORDER node in the correct order
    final LocationParseNodeVector new_children =
      (LocationParseNodeVector)new_tree.child;

    if(new_children.size() >= 2)
    {
      // first move all ENTRY_RANGE nodes to the end

      final LocationParseNodeVector entry_range_nodes =
        new LocationParseNodeVector();

      for(int i = new_children.size() - 1 ; i >= 0 ; --i) 
      {
        final LocationParseNode node = new_children.elementAt(i);
        if(node.getType() == ENTRY_RANGE) 
        {
          new_children.removeElementAt(i);
          entry_range_nodes.addElementAtEnd(node);
        }
        else
        {
          if(node.getType() == COMPLEMENT &&
             node.getComplementChild().getType() == ENTRY_RANGE) 
          {
            new_children.removeElementAt(i);
            entry_range_nodes.addElementAtEnd(node);
          }
        }
      }

      if(new_children.size() > 0)
      {
      
        // the new head node is a JOIN or ORDER and the children are RANGE
        // nodes so put the children in order - ascending left to right
        // or
        // the children are COMPLEMENT nodes so put the children in order 
        // - ascending right to left

        // bubble sort in place putting the RANGE with the smallest start
        // base first

        boolean change_happened = true;
        while(change_happened) 
        {
          change_happened = false;
          int this_start;
          int next_start;

          for(int i = 0 ; i < new_children.size() - 1 ; ++i) 
          {
            final LocationParseNode this_node = new_children.elementAt(i);
            final LocationParseNode next_node = new_children.elementAt(i+1);

            if(this_node.getType() == RANGE)
              this_start = this_node.getRange().getStart();
            else
              this_start = this_node.getComplementChild().getRange().getEnd();

            if(next_node.getType() == RANGE)
              next_start = next_node.getRange().getStart();
            else
              next_start = next_node.getComplementChild().getRange().getEnd();

            if(this_node.getType() != RANGE && next_node.getType() != RANGE)
            {
              if(this_start < next_start)      // swap the two nodes
              {
                new_children.setElementAt(next_node, i);
                new_children.setElementAt(this_node, i + 1);
                change_happened = true;
              }
            }
            else if(this_start > next_start)   // swap the two nodes
            {
              new_children.setElementAt(next_node, i);
              new_children.setElementAt(this_node, i + 1);
              change_happened = true;
            }
          }
        }
      }

      for(int i = 0 ; i < entry_range_nodes.size() ; ++i) 
        new_children.addElementAtEnd(entry_range_nodes.elementAt(i));
    } 
    else 
      throw new Error("internal error - a JOIN should have > 1 child");

    return new_tree;
  }

  /**
   *  Return true if and only if the tree starting at this node is valid.
   *  "valid" means it has at most two levels of functionals and the lowest
   *  level consists only of RANGE nodes.  Also the children of a JOIN or
   *  ORDER node must be all RANGE nodes or all COMPLEMENT nodes - there can't
   *  be a mixture.
   **/
  public boolean isValid() 
  {
    // 2 means check that we have at most two levels.
    // UNKNOWN means the parent isn't a JOIN, ORDER or COMPLEMENT node (there
    // is no parent at all.)
    return isValid(2, UNKNOWN);
  }

  /**
   *  Return the child node of this node.  This node must be a COMPLEMENT type
   *  or an Error will be thrown.
   **/
  public LocationParseNode getComplementChild()
  {
    if(getType() != COMPLEMENT) 
      throw new Error("in LocationParseNode.getComplementChild() - node " +
                      "is not a COMPLEMENT");
    
    return (LocationParseNode)child;
  }

  /**
   *  Return the children of this node.  If this node is not a JOIN or ORDER
   *  type then an error will be thrown.
   **/
  public LocationParseNodeVector getChildren()
  {
    if(!(getType() == JOIN || getType() == ORDER)) 
      throw new Error("in LocationParseNode.getChildren() - node " +
                      "is not a JOIN or ORDER");
    
    return (LocationParseNodeVector)child;
  }

  /**
   *  Return the children of this node.  This node must be a JOIN type node or
   *  an Error will be thrown.
   **/
  public LocationParseNodeVector getJoinChildren()
  {
    if(getType() != JOIN) 
      throw new Error("in LocationParseNode.getJoinChildren() - node " +
                      "is not a JOIN");
    
    return (LocationParseNodeVector) child;
  }

  /**
   *  Return the children of this node.  This node must be a JOIN type node or
   *  an Error will be thrown.
   **/
  public LocationParseNodeVector getOrderChildren()
  {
    if(getType() != ORDER) 
      throw new Error("in LocationParseNode.getOrderChildren() - node " +
                       "is not an ORDER");
    
    return(LocationParseNodeVector) child;
  }

  /**
   *  Return the child of this ENTRY_RANGE node.  This node must be an
   *  ENTRY_RANGE node or an Error will be thrown.  The returned node will be
   *  of type RANGE.
   **/
  public LocationParseNode getEntryRangeChild()
  {
    if(getType() != ENTRY_RANGE) 
      throw new Error("in LocationParseNode.getEntryRangeChild() - node " +
                      "is not an ENTRY_RANGE");
    
    return (LocationParseNode)child;
  }

  /**
   *  Return the Range object for this node.  This node must be a RANGE type
   *  node or an Error will be thrown.
   **/
  public Range getRange() 
  {
    if(getType() != RANGE)
      throw new Error("in LocationParseNode.getRange() - node " +
                      "is not a RANGE " + getType());
    
    return (Range)child;
  }


  /**
   *  Internal method that implements copy() and copyClean().
   *  @param copy_user_data user_data (from EMBLObject) will be copied if and
   *    only if this parameter is true.
   **/
  private LocationParseNode copyInternal(final boolean copy_user_data) 
  {
    final LocationParseNode new_node;

    switch(getType())
    {
      case COMPLEMENT:
      {
        final LocationParseNode new_child =
          getComplementChild().copyInternal(copy_user_data);

        new_node = new LocationParseNode(getType(), new_child);
      }
      break;
      case JOIN:
      case ORDER:
        // copy the children
        final LocationParseNodeVector children_copy =
          new LocationParseNodeVector();
        final LocationParseNodeVector children = (LocationParseNodeVector)child;
        for(int i = 0; i < children.size(); ++i) 
        {
          final LocationParseNode new_child =
            children.elementAt(i).copyInternal(copy_user_data);
          children_copy.addElement(new_child);
        }
        new_node = new LocationParseNode(getType(), children_copy);
        break;
      case RANGE:
        new_node = new LocationParseNode(getRange().copy());
        break;
      case ENTRY_RANGE:
      {
        final LocationParseNode new_child =
          ((LocationParseNode)child).copyInternal(copy_user_data);
        new_node = new LocationParseNode(entry_name, new_child);
      }
      break;
      default:
        throw new Error("internal error - unknown location node type");
    }

    if(copy_user_data)
      new_node.setUserData(getUserData());
    else 
      new_node.setUserData(null);
   
    return new_node;
  }

  /**
   *  Return a copy of this node.  All children of this node will be copied as
   *  well.  No user_data will be copied.
   **/
  LocationParseNode copyClean()
  {
    return copyInternal(false);
  }

  /**
   *  Return a copy of this node.  All children of this node and all user_data
   *  fields will be copied as well.
   **/
  LocationParseNode copy() 
  {
    return copyInternal(true);
  }

  /**
   *  Set the value of a RANGE node.  If this node is not a RANGE node an
   *  Error will be thrown.
   **/
  void setRange(final Range range) 
  {
    if(getType() != RANGE) 
      throw new Error("in LocationParseNode.getRange() - node " +
                      "is not a RANGE " + getType());
    this.child = range;
  }

  /**
   *  Returns the result of calling toString() on each of the children of
   *  this node separated by ",".
   **/
  private String toStringChildren() 
  {
    String return_string = "";

    final LocationParseNodeVector children = (LocationParseNodeVector) child;

    if(getType() == JOIN || getType() == ORDER) 
    {
      for(int i = 0; i<children.size(); ++i) 
      {
        if(0 != i) 
          return_string = return_string + ",";

        //        System.out.println("adding:" + children.elementAt(i).toString());
        return_string = return_string + children.elementAt(i).toString();
      }
    }
    else 
      throw new Error ("LocationParseNode.toStringChildren(): Illegal type");

    return return_string;
  }

  /**
   *  Turn an Object into a String.  The Object must be an instance of String,
   *  LocationParseNode or Range.
   **/
  private String toStringChild() 
  {
    if(child instanceof LocationParseNode) 
      return ((LocationParseNode)child).toString();
    else 
      throw new Error("LocationParseNode.toStringChild() was called on " +
                      "the wrong type of Object");
  }

  /**
   *  Return a complemented copy of a node.  ie. COMPLEMENT(THING) will become
   *  THING and THING will become COMPLEMENT(THING).
   **/
  LocationParseNode getNodeComplement() 
  {
    if(getType() == COMPLEMENT) 
      return getComplementChild().copy();
    else 
      return new LocationParseNode(COMPLEMENT, copy());
  }

  /**
   *  Return a copy of the tree starting at this node with any
   *  "COMPLEMENT(JOIN(..." nodes turned to "JOIN(COMPLEMENT(..."
   **/
  private LocationParseNode getCanonicalComplement() 
  {
    if(getType() == COMPLEMENT &&
       (getComplementChild().getType() == JOIN ||
        getComplementChild().getType() == ORDER)) 
    {
      final LocationParseNode complement_child = getComplementChild();

      // complement each child of the JOIN or ORDER
      final LocationParseNodeVector complemented_children =
        new LocationParseNodeVector();

      final LocationParseNodeVector complement_children =
        (LocationParseNodeVector)complement_child.child;

      for(int i = 0 ; i < complement_children.size() ; ++i) 
      {
        final LocationParseNode complemented_child =
          complement_children.elementAt(i).getNodeComplement();
        complemented_children.addElement(complemented_child);
      }

      return new LocationParseNode(complement_child.getType(),
                                    complemented_children);
    }
    else   // no changes needed so just make a copy
      return this.copy();
  }

  /**
   *  Return true if and only if the tree starting at this node is valid.
   *  @see #isValid()
   *  @param max_levels The maximum number of levels below this node - more
   *    than this number of levels means that the tree is invalid.
   **/
  private boolean isValid(final int max_levels, final int parent_type) 
  {
    if(getType() == RANGE) 
    {
      // a RANGE a valid possibility at any level
      return true;
    }

    if(getType() == ENTRY_RANGE) 
    {
      return getEntryRangeChild().isValid();
    }

    if(max_levels == 0) 
    {
      // we are at the lowest level and this isn't a RANGE so this tree is
      // nonsense
      return false;
    } 
    else 
    {
      switch(getType()) 
      {
        case JOIN:
        case ORDER:
          if(parent_type == UNKNOWN || parent_type == COMPLEMENT) 
          {
            // parent has the correct type

            // this is set to true if any of the children are RANGE nodes
            boolean found_range_child = false;

            // this is set to true if any of the children are COMPLEMENT nodes
            boolean found_complement_child = false;
  
            final LocationParseNodeVector children =
              (LocationParseNodeVector) child;

            for(int i = 0; i < children.size(); ++i) 
            {
              final LocationParseNode child_node = children.elementAt(i);
              final boolean return_val =
                child_node.isValid(max_levels - 1, getType());

              // one of the children is invalid so this node is invalid
              if(!return_val) 
                return false;

              if(child_node.getType() == RANGE)
                found_range_child = true;
  
              if(child_node.getType() == COMPLEMENT) 
                found_complement_child = true;
            }

            // all the children were OK

            if(found_range_child && found_complement_child) 
            {
              // a mixture of RANGE and COMPLEMENT nodes is not valid
              //return false;
              JOptionPane.showMessageDialog(null, "Contains trans-spliced feature!\n"+toString(),
                                            "Trans-spliced site found",
                                            JOptionPane.INFORMATION_MESSAGE);
              return true;
            } 
            else 
              return true;
          }
          else   // parent has the wrong type
            return false;

        case COMPLEMENT:
          if(parent_type == UNKNOWN ||
             parent_type == JOIN ||
             parent_type == ORDER) 
          {
            return getComplementChild().isValid(max_levels - 1, getType());
          } 
          else    // parent has the wrong type (ie. it is a COMPLEMENT)
            return false;

        default:
          // anything else is invalid
          return false;
      }
    }
  }

}
