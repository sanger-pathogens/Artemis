/* PlotMouseListener.java
 *
 * created: Wed Sep 13 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/PlotMouseListener.java,v 1.1 2004-06-09 09:47:12 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

/**
 *  Implemented by those objects that need to mouse events from a Plot object.
 *  The coordinates of the mouse click are translated to base/aa positions for
 *  ease of use.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: PlotMouseListener.java,v 1.1 2004-06-09 09:47:12 tjc Exp $
 **/

public interface PlotMouseListener {
  /**
   *  Called when the user clicks somewhere on the plot canvas.
   *  @param position the base/amino acid position of the click.  This is
   *    -1 if and only if the click was outside the graph (eg. in the label at
   *    the top)
   **/
  void mouseClick (final int position);

  /**
   *  Called when the user drags the mouse over the plot.
   *  @param drag_start_position The base/amnino acid position where the drag
   *    started or -1 if the drag was started outside the graph.
   *  @param current_position the base/amino acid position of the click.
   *    This is -1 if and only if the user has dragged the mouse out of
   *    the graph (eg. in the label at the top)
   **/
  void mouseDrag (final int drag_start_position,
                  final int current_position);

  /**
   *  Called when the user double-clicks somewhere on the plot.
   *  @param position the base/amino acid position of the click.  This is
   *    -1 if and only if the click was outside the graph (eg. in the label at
   *    the top)
   **/
  void mouseDoubleClick (final int position);
}
