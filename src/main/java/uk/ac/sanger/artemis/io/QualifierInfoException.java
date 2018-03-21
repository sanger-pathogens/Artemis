/* QualifierInfoException.java
 *
 * created: Sun Feb 21 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/QualifierInfoException.java,v 1.1 2004-06-09 09:50:09 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

/**
 *  Objects of this class are thrown when as incorrect type string is passed
 *  to the QualifierInfo constructor or when a QualifierInfo object is
 *  re-added to an EntryInformation object.
 *
 *  @author Kim Rutherford
 *  @version $Id: QualifierInfoException.java,v 1.1 2004-06-09 09:50:09 tjc Exp $
 **/

public class QualifierInfoException extends Exception {
  /**
   *  Create a new QualifierInfoException object with the given String as the
   *  message.
   **/
  public QualifierInfoException (final String message) {
    super (message);
  }
}


