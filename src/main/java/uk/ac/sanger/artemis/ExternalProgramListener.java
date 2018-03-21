/* ExternalProgramListener.java
 *
 * created: Wed Aug  6 2003
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2003  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/ExternalProgramListener.java,v 1.1 2004-06-09 09:44:31 tjc Exp $
 */

package uk.ac.sanger.artemis;

/**
 *  This interface is implemented by those classes that need to listen for
 *  changes to ExternalProgram objects.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: ExternalProgramListener.java,v 1.1 2004-06-09 09:44:31 tjc Exp $
 **/

public interface ExternalProgramListener {
  /**
   *  Called when the status of an ExternalProgram changes.
   **/
  void statusChanged (final ExternalProgramEvent event);
}
