/* ComparisonDataFactory.java
 *
 * created: Thu Jul 15 1999
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 1999-2002  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/ComparisonDataFactory.java,v 1.1 2004-06-09 09:44:16 tjc Exp $
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.util.LinePushBackReader;

import java.io.*;

/**
 *  This class contains the method readComparisonData (), which returns an
 *  appropriate ComparisonData object for a given Document.
 *
 *  @author Kim Rutherford
 *  @version $Id: ComparisonDataFactory.java,v 1.1 2004-06-09 09:44:16 tjc Exp $
 **/

public class ComparisonDataFactory {
  /**
   *  This method creates an appropriate ComparisonData object from a Document.
   **/
  static public ComparisonData readComparisonData (Document data_document)
      throws IOException {
    
    final Reader in_file = data_document.getReader ();

    final LinePushBackReader pushback_reader =
      new LinePushBackReader (in_file);

    final String line = pushback_reader.readLine ();

    if (line == null) {
      throw new IOException ("End of file while reading from: " +
                             data_document);
    }

    pushback_reader.pushBack (line);

    if (MSPcrunchComparisonData.formatCorrect (line)) {
      return new MSPcrunchComparisonData (pushback_reader);
    } else {
      if (SSAHAComparisonData.formatCorrect (line)) {
        return new SSAHAComparisonData (pushback_reader);
      } else {
        if (BlastM8ComparisonData.formatCorrect (line)) {
          return new BlastM8ComparisonData (pushback_reader);
        } else {
          if (MegaBlastComparisonData.formatCorrect (line)) {
            return new MegaBlastComparisonData (pushback_reader);
          } else {
//      if (tokenizer.countTokens () < 8) {
//        return new MUMmerComparisonData (pushback_reader);
//      } else {
            throw new IOException ("cannot understand the comparison file format");
//      }
          }
        }
      }
    }
  }
}

