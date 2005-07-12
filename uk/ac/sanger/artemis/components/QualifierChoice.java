/* QualifierChoice.java
 *
 * created: Tue Sep  7 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/QualifierChoice.java,v 1.2 2005-07-12 14:44:31 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;

import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.InvalidKeyException;
import uk.ac.sanger.artemis.io.EntryInformation;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;
import javax.swing.event.*;

/**
 *  This is a Choice component that shows only the qualifier names for a given
 *  key.
 *
 *  @author Kim Rutherford
 *  @version $Id: QualifierChoice.java,v 1.2 2005-07-12 14:44:31 tjc Exp $
 **/

public class QualifierChoice extends JComboBox {
  /**
   *  Create a new QualifierChoice component for the given Key with the given
   *  qualifier as the default.
   *  @param entry_information The object to get the list of possible
   *    qualifiers from.
   *  @param default_qualifier The name of the Qualifier that should be shown
   *    initially.  If null the first (alphabetically) is selected.
   **/
  public QualifierChoice (final EntryInformation entry_information,
                          final Key key, final String default_qualifier) {
    this.entry_information = entry_information;
    this.key = key;

    if (default_qualifier != null &&
        entry_information.isValidQualifier (key, default_qualifier)) {
      this.default_qualifier = default_qualifier;
    } else {
      this.default_qualifier = null;
    }

    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    int screen_height = screen.height;

    final int MAX_VISIBLE_ROWS;
    if(screen_height < 1024)
      MAX_VISIBLE_ROWS = 20;
    else
      MAX_VISIBLE_ROWS = 30;
    
    setMaximumRowCount (MAX_VISIBLE_ROWS);

    setEditable(true);

    update ();
  }

  /**
   *  Change the qualifiers shown in this component to be those of the given
   *  Key.
   **/
  public void setKey (final Key key) {
    if (this.key != key) {
      this.key = key;

      update ();
    }
  }

  /**
   *  Select the given qualifier_name.
   **/
  private void setSelectedQualifierByName (final String qualifier_name) {
    final int index = indexOf (qualifier_name);

    if (index == -1) {
      // add the key
      addItem (qualifier_name);
      setSelectedItem (qualifier_name);
    } else {
      setSelectedIndex (index);
    }
  }

  /**
   *  Return the index in the Choice component of the given qualifier_name.
   **/
  private int indexOf (final String qualifier_name) {
    for (int i = 0 ; i < getItemCount () ; ++i) {
      if (getItemAt (i).equals (qualifier_name)) {
        return i;
      }
    }
    return -1;
  }

  /**
   *  Update the Choice to refect the current Key.
   **/
  private void update () {
    removeAllItems ();

//    final StringVector common_qualifiers = new StringVector ();
//    final StringVector uncommon_qualifiers = new StringVector ();
    
    StringVector qualifier_names =
      entry_information.getValidQualifierNames (key);

    if (qualifier_names == null) {
      qualifier_names = new StringVector ("note");
    }

    if (default_qualifier != null &&
        !qualifier_names.contains (default_qualifier)) {
      qualifier_names.add (default_qualifier);
    }

    qualifier_names.sort ();

//    final StringVector invisible_qualifiers =
//      Options.getOptions ().getInvisibleQualifiers ();

    for (int i = 0 ; i < qualifier_names.size () ; ++i) {
      final String qualifier_name = qualifier_names.elementAt (i);
//      if (!invisible_qualifiers.contains (qualifier_name)) {
        addItem (qualifier_name);
//      }
    }

    if (default_qualifier == null) {
      if (indexOf ("note") != -1) {
        setSelectedQualifierByName ("note");
      } else {
        if (indexOf ("locus_tag") != -1) {
          setSelectedQualifierByName ("locus_tag");
        } else {
          setSelectedIndex (0);
        }
      }
    } else {
      setSelectedQualifierByName (default_qualifier);
    }
  }

  /**
   *  The Key that was passed to the constructor.
   **/
  private Key key = null;

  /**
   *  The qualifier name that was passed to the constructor.
   **/
  private String default_qualifier = null;

  /**
   *  The EntryInformation object that was passed to the constructor.
   **/
  private EntryInformation entry_information;
}

