/* DisplayComponent.java (formally SelectionDisplayer.java)
 *
 * created: Fri Nov 13 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/DisplayComponent.java,v 1.2 2008-10-30 15:25:24 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.GotoEventSource;

import javax.swing.JFrame;

/**
 *  Interface discribing those methods common to all the classes in
 *  uk.ac.sanger.artemis.components that can display EntryGroup objects.
 *
 *  @author Kim Rutherford
 *  @version $Id: DisplayComponent.java,v 1.2 2008-10-30 15:25:24 tjc Exp $
 **/

public interface DisplayComponent
{
  /**
   *  Return an object that implements the GotoEventSource interface.
   **/
  GotoEventSource getGotoEventSource ();

  /**
   *  Return the reference of the JFrame that owns this component.
   **/
  JFrame getParentFrame ();
}
