/* ProcessWatcherListener.java
 *
 * created: Tue Feb 29 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/ProcessWatcherListener.java,v 1.1 2004-06-09 09:47:15 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

/**
 *  This interface is implemented by those classes that need to be notified
 *  when a Process finishes.  The Process must be watched with a
 *  ProcessWatcher object.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: ProcessWatcherListener.java,v 1.1 2004-06-09 09:47:15 tjc Exp $
 **/

public interface ProcessWatcherListener {
  /**
   *  Invoked by a ProcessWatcher object when a Process finishes.
   **/
  void processFinished (final ProcessWatcherEvent event);
}
