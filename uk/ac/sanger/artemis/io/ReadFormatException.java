/* ReadFormatException.java
 *
 * created: Tue Oct 13 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/ReadFormatException.java,v 1.1 2004-06-09 09:50:22 tjc Exp $
 */

package uk.ac.sanger.artemis.io;
import java.io.IOException;

/**
 *  An object of this type is thrown while reading if the format of the input
 *  stream is incorrect.
 *
 *  @author Kim Rutherford
 *  @version $Id: ReadFormatException.java,v 1.1 2004-06-09 09:50:22 tjc Exp $
 *
 **/
public class ReadFormatException extends IOException {
  /**
   *  This constructor creates a ReadFormatException with no detail message.
   **/
  private ReadFormatException () {
    super ();
  }

  /**
   *  This constructor creates a ReadFormatException with the given String
   *  as the message.
   *  @param message the detail message
   **/
  public ReadFormatException (String message) {
    this (message, 0);
  }

  /**
   *  This constructor creates a ReadFormatException with the given String
   *  as the message and a line number.
   *  @param message The detail message.
   *  @param line_number The line number at which the error occured.
   **/
  public ReadFormatException (String message, int line_number) {
    super (message);
    this.line_number = line_number;
  }

  /**
   *  Return the line number of the exception.
   **/
  public int getLineNumber () {
    return line_number;
  }

  /**
   *  The line number of the error.
   **/
  private int line_number = 0;
}
