/* DisplayAdjustmentEvent.java
 *
 * created: Tue Dec 15 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/DisplayAdjustmentEvent.java,v 1.3 2005-12-02 14:58:57 tjc Exp $
 **/

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.ChangeEvent;

/**
 *  This event is sent when a FeatureDisplay is scrolled.
 *
 *  @author Kim Rutherford
 *  @version $Id: DisplayAdjustmentEvent.java,v 1.3 2005-12-02 14:58:57 tjc Exp $
 **/

public class DisplayAdjustmentEvent extends ChangeEvent 
{

  /** The new start base, as passed to the constructor. */
  private int start_base;

  /** The new end base, as passed to the constructor. */
  private int end_base;

  /**
   *  The width in bases of the display, as passed to the constructor.
   **/
  private int width_in_bases;

  /**
   *  The width of each base on the display, as passed to the constructor.
   **/
  private float base_width;

  /** The scale factor, as passed to the constructor. */
  private int scale_factor;

  private int drop_position;

  /**
   *  True if and only if the FeatureDisplay is drawing in
   *  reverse complement mode.
   **/
  private boolean rev_comp_display;

  /**
   *  The type of event.  One of: SCALE_ADJUST_EVENT, SCROLL_ADJUST_EVENT or
   *  ALL_CHANGE_ADJUST_EVENT.
   **/
  private int type;

  /**
   *  The type of DisplayAdjustmentEvent where the scale (only) has changed
   **/
  final static public int SCALE_ADJUST_EVENT = 0;

  /**
   *  The type of DisplayAdjustmentEvent where the display has scrolled but
   *  the scale has stayed the same.
   **/
  final static public int SCROLL_ADJUST_EVENT = 1;

  /**
   *  The type of DisplayAdjustmentEvent where the display has been reverse
   *  complemented.
   **/
  final static public int REV_COMP_EVENT = 2;

  /**
   *  The type of DisplayAdjustmentEvent where the display has scrolled and
   *  the the scale has changed.
   **/
  final static public int ALL_CHANGE_ADJUST_EVENT = 3;
  
  /**
   *  The type of DisplayAdjustmentEvent where the display has been reverse
   *  complemented.
   **/
  final static public int CONTIG_REV_COMP_EVENT = 4;

  final static public int CONTIG_REORDER = 5;
  
  /** change tabix indexed sequence */
  final static public int IDX_SEQUENCE_CHANGE = 6;

  public DisplayAdjustmentEvent(Object source,
                                int start_base, int end_base,
                                int drop_position, int type)
  {
    super(source);
    this.start_base    = start_base;
    this.end_base      = end_base;
    this.drop_position = drop_position;
    this.type          = type;
  }

  /**
   *  Create a new DisplayAdjustmentEvent.
   *  @param source The Object that generated the event - probably a component.
   *  @param start_base The new start base on the display.
   *  @param end_base The new end base on the display.
   *  @param width_in_bases The width of the drawing area in bases.  We need
   *    this because if the user scrolls to the end of the sequence, the last
   *    visible base (end_base) may not be at the right of the screen.  This
   *    might be greater than end_base - start_base.
   *  @param base_width The width in pixels of a base on screen.
   *  @param scale_factor This is the scale factor use by the FeatureDisplay
   *    component.  A factor of zero means the full translation will be
   *    visible.  At higher scale factors only stop codons are visible, and
   *    a bigger number will mean more bases are visible.
   *  @param type the type of event: SCALE_ADJUST_EVENT, SCROLL_ADJUST_EVENT
   *    or ALL_CHANGE_ADJUST_EVENT.
   **/
  public DisplayAdjustmentEvent(Object source,
                                int start_base, int end_base,
                                int width_in_bases, float base_width,
                                int scale_factor,
                                boolean rev_comp_display, int type)
  {
    super(source);
    this.start_base     = start_base;
    this.end_base       = end_base;
    this.width_in_bases = width_in_bases;
    this.base_width     = base_width;
    this.scale_factor   = scale_factor;
    this.rev_comp_display = rev_comp_display;
    this.type           = type;
  }

  /**
   *  Return the new start base to display, as passed to the constructor.
   **/
  public int getStart() 
  {
    return start_base;
  }

  /**
   *  Return the new end base to display, as passed to the constructor.
   **/
  public int getEnd() 
  {
    return end_base;
  }

  public int getDropPosition()
  {
    return drop_position;
  }

  /**
   *  Return the width in bases of the display, as passed to the constructor.
   **/
  public int getWidthInBases() 
  {
    return width_in_bases;
  }

  /**
   *  Return the width of a base on the display, as passed to the constructor.
   **/
  public float getBaseWidth() 
  {
    return base_width;
  }

  /**
   *  Return the scale factor that was passed to the constructor.
   **/
  public int getScaleFactor() 
  {
    return scale_factor;
  }

  /**
   *  Return true if and only if the FeatureDisplay is drawing in reverse
   *  complement mode.
   **/
  public boolean isRevCompDisplay() 
  {
    return rev_comp_display;
  }

  /**
   *  Return the type that was passed to the constructor.
   **/
  public int getType() 
  {
    return type;
  }

}


