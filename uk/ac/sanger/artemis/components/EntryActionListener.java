/* EntryActionListener.java
 *
 * created: Mon Jan 31 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/EntryActionListener.java,v 1.1 2004-06-09 09:46:21 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;

import java.awt.event.*;

/**
 *  This class is an implementation of the ActionListener interface that can
 *  remember an entry.
 *
 *  @author Kim Rutherford
 *  @version $Id: EntryActionListener.java,v 1.1 2004-06-09 09:46:21 tjc Exp $
 **/

abstract public class EntryActionListener implements ActionListener {
  /**
   *  Make a new EntryActionListener from the given Entry.
   **/
  EntryActionListener (final EntryEdit entry_edit,
                       final Entry entry) {
    this.entry = entry;
    this.entry_edit = entry_edit;
  }

  abstract public void actionPerformed (final ActionEvent event);

  public Entry getEntry () {
    return entry;
  }

  public EntryEdit getEntryEdit () {
    return entry_edit;
  }

  final private Entry entry;
  final private EntryEdit entry_edit;
}

