/* MarkerChangeListener.java
 *
 * created: Sun Mar 21 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/sequence/MarkerChangeListener.java,v 1.1 2004-06-09 09:52:19 tjc Exp $
 */

package uk.ac.sanger.artemis.sequence;

/**
 *  The MarkerChangeListener interface is implemented by those classes that
 *  need to know when a Marker changes it's position or Strand.
 *
 *  @author Kim Rutherford
 *  @version $Id: MarkerChangeListener.java,v 1.1 2004-06-09 09:52:19 tjc Exp $
 **/

public interface MarkerChangeListener extends uk.ac.sanger.artemis.ChangeListener {
  /**
   *  Invoked when a Marker changes it's position or Strand.
   **/
  void markerChanged (MarkerChangeEvent event);
}


