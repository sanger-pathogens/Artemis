/* LocationParseException.java
 *
 * created: Wed Oct  7 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/LocationParseException.java,v 1.2 2004-10-29 09:36:24 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.io.LocationLexer.*;

/**
 *  LocationParseException objects are thrown if a parse error occurs during
 *  the parsing of a location string.
 *
 *  @author Kim Rutherford
 *  @version $Id: LocationParseException.java,v 1.2 2004-10-29 09:36:24 tjc Exp $
 *
 */

public class LocationParseException extends ReadFormatException 
{
  /**
   *  Create a new LocationParseException with the given String as the
   *  message.
   *  @param message The detail message.
   *  @param enumTk An enumeration containing the next tokens to be read.  This
   *    is used to give the user an indication of where in the location string
   *    the parsed error happens.
   **/
  public LocationParseException(String message, TokenEnumeration enumTk)
  {
    super("Parse error at this point: " + enumTk.toString() +
          ": " + message);
  }

  /**
   *  Create a new LocationParseException with the given String as the
   *  message.
   *  @param message The detail message.
   *  @param location_string The String in which the error occured.
   **/
  public LocationParseException(String message, String location_string) 
  {
    super("Parse error in this location: " + location_string +
          ": " + message);
  }
}


