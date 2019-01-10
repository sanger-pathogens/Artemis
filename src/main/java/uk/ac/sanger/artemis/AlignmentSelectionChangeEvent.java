/* AlignmentSelectionChangeEvent.java
 *
 * created: Tue Feb 13 2001
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2001  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/AlignmentSelectionChangeEvent.java,v 1.1 2004-06-09 09:44:08 tjc Exp $
 */

package uk.ac.sanger.artemis;

/**
 *  Events of this type are send from AlignmentViewer objects whenever the list
 *  of selected/highlighted hits changes.
 *
 *  @author Kim Rutherford
 *  @version $Id: AlignmentSelectionChangeEvent.java,v 1.1 2004-06-09 09:44:08 tjc Exp $
 **/

public class AlignmentSelectionChangeEvent extends java.util.EventObject {
  /**
   *  Create a new SelectionChangeEvent object.
   *  @param source The source of the event.
   *  @param selection The selected AlignMatch objects after the change.
   **/
  public AlignmentSelectionChangeEvent (final Object source,
                                        final AlignMatchVector selection) {
    super (source);
    this.selected_matches = selection;
  }

  /**
   *  The selected matches.
   **/
  private AlignMatchVector selected_matches = null;
}
