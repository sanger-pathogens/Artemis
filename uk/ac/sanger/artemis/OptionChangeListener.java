/* OptionChangeListener.java
 *
 * created: Tue Jun 12 2001
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/OptionChangeListener.java,v 1.1 2004-06-09 09:44:59 tjc Exp $
 */

package uk.ac.sanger.artemis;

/**
 *  The OptionChangeListener interface is implemented by those classes that
 *  need to listen for changes to Option objects.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: OptionChangeListener.java,v 1.1 2004-06-09 09:44:59 tjc Exp $
 **/

public interface OptionChangeListener extends ChangeListener {
  /**
   *  Invoked when an Option is changed.
   **/
  void optionChanged (OptionChangeEvent event);
}
