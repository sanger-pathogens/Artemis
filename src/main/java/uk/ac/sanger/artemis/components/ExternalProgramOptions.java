/* ExternalProgramOptions.java
 *
 * created: Mon Oct  4 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/ExternalProgramOptions.java,v 1.3 2009-05-13 15:45:04 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;

/**
 *  This component allows the user to set the options of an ExternalProgram.
 *
 *  @author Kim Rutherford
 *  @version $Id: ExternalProgramOptions.java,v 1.3 2009-05-13 15:45:04 tjc Exp $
 **/

public class ExternalProgramOptions {
  /**
   *  Create a new ExternalProgramOptions object for the given
   *  ExternalProgram.
   **/
  public ExternalProgramOptions (final ExternalProgram external_program) {
    if (external_program.getType () == ExternalProgram.AA_PROGRAM ||
        external_program.getType () == ExternalProgram.DNA_PROGRAM) {

      final TextRequester requester =
        new TextRequester ("Options for " + external_program.getName () + ":",
                           18, external_program.getProgramOptions ());

      requester.addTextRequesterListener (new TextRequesterListener () {
        public void actionPerformed (final TextRequesterEvent event) {
          final String requester_text = event.getRequesterText ().trim ();
          if (requester_text.length () > 0) {
            external_program.setProgramOptions (requester_text);
          }
        }
      });

      requester.setVisible(true);
    } else {
      throw new Error ("internal error - please hit the programmer");
    }
  }
}
