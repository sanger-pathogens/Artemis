/* QualifierTextArea.java
 *
 * created: Tue Oct 23 2001
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/QualifierTextArea.java,v 1.1 2004-06-09 09:47:19 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.Options;

import uk.ac.sanger.artemis.io.QualifierParseException;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.EmblStreamFeature;

import java.awt.*;
import java.io.*;
import javax.swing.*;

/**
 *  This component is a TextArea that understands qualifiers.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: QualifierTextArea.java,v 1.1 2004-06-09 09:47:19 tjc Exp $
 **/

public class QualifierTextArea extends JTextArea {
  /**
   *  Create a new QualifierTextArea containing no text.
   **/
  public QualifierTextArea () {
    super ((Options.getOptions ().getPropertyTruthValue ("alicat_mode") ||
            Options.getOptions ().getPropertyTruthValue ("val_mode") ?
            40 :
            18),
           81);
    setLineWrap (true);
    setBackground (Color.white);

    setDragEnabled(true);
  }

  /**
   *  Parse and return the qualifiers in this TextArea in a QualifierVector.
   **/
  public QualifierVector
    getParsedQualifiers (final EntryInformation entry_information)
      throws QualifierParseException {
    final String qualifier_string = getText ();
    return getQualifiersFromString (qualifier_string,
                                    entry_information);
  }

  /**
   *  Return a QualifierVector containing the qualifiers from a String.
   *  @param qual_string contains the qualifiers to parse
   */
  private static QualifierVector
    getQualifiersFromString (final String qual_string,
                             final EntryInformation entry_information)
      throws QualifierParseException {

    try {
      final StringReader string_reader = new StringReader (qual_string);

      final QualifierVector embl_qualifiers =
        EmblStreamFeature.readQualifiers (string_reader,
                                          entry_information);

      string_reader.close ();

      return embl_qualifiers;
    } catch (IOException exception) {
      throw (new QualifierParseException (exception.getMessage ()));
    }
  }
}
