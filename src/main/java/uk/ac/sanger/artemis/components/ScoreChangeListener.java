/* ScoreAdjustmentListener.java
 *
 * created: Thu Oct 21 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/ScoreChangeListener.java,v 1.1 2004-06-09 09:47:28 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

/**
 *  The listener interface for receiving score change events. 
 *
 *  @author Kim Rutherford
 *  @version $Id: ScoreChangeListener.java,v 1.1 2004-06-09 09:47:28 tjc Exp $
 **/

public interface ScoreChangeListener {
  /**
   *  Invoked when the value of the score has changed.
   **/
  void scoreChanged (final ScoreChangeEvent event);
}
