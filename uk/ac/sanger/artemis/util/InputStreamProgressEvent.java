/* InputStreamProgressEvent.java
 *
 * created: Thu Jun  8 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/util/InputStreamProgressEvent.java,v 1.1 2004-06-09 09:53:02 tjc Exp $
 */

package uk.ac.sanger.artemis.util;

/**
 *  This event is sent when progress is made while reading from an
 *  InputStream.  "progress" means some more bytes have been read. 
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: InputStreamProgressEvent.java,v 1.1 2004-06-09 09:53:02 tjc Exp $
 **/

public class InputStreamProgressEvent {
  /**
   *  Call the constructor with the char_count parameter set to this if the
   *  end of stream has been reached.
   **/
  public final static int EOF = -1;

  /**
   *  Create a new InputStreamProgressEvent.
   *  @param char_count The total number of chars that have been read so far.
   *    Call the constructor with this parameter set to EOF at end of stream.
   **/
  public InputStreamProgressEvent (final int char_count) {
    this.char_count = char_count;
  }

  /**
   *  Return the char_count that was passed to the constructor.  A return
   *  value of EOF means that end of file has been hit.
   **/
  public int getCharCount () {
    return char_count;
  }

  /**
   *  The count that was passed to the constructor.
   **/
  private int char_count;
}
