/* CorbaEntrySource.java
 *
 * created: Wed Jun  7 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/CorbaEntrySource.java,v 1.1 2004-06-09 09:46:11 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import java.io.*;
import java.net.*;
import java.awt.*;

import javax.swing.*;

/**
 *  This class contains the methods common to all EntrySource implementations
 *  that read from a CORBA server.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: CorbaEntrySource.java,v 1.1 2004-06-09 09:46:11 tjc Exp $
 **/

abstract public class CorbaEntrySource {
  /**
   *  Create a new CorbaEntrySource using the given URL as the location of the
   *  server IOR.
   *  @param frame The component that created this EntrySource.  (Used for
   *    requesters.)
   *  @param ior_url_string A String containing the URL of the IOR for the
   *    server.
   **/
  public CorbaEntrySource (final JFrame frame,
                           final String ior_url_string)
      throws MalformedURLException //, IOException
  {
    this.frame = frame;
    this.ior_url = new URL (ior_url_string);
  }

  /**
   *  Finds and returns a stringified object reference (IOR) for a CORBA
   *  server.  The IOR is read from the URL that was passed to the
   *  constructor.
   *  @return The stringified IOR.
   *  @exception IOException Thrown if an error occurs while reading.
   **/
  public String getIOR ()
      throws IOException {
    //  get URL of IOR
    final BufferedReader in =
      new BufferedReader (new InputStreamReader (ior_url.openStream()));
    return in.readLine();
  }

  /**
   *  Return the JFrame that was passed to the constructor.
   **/
  public JFrame getFrame () {
    return frame;
  }

  /**
   *  The URL containing the IOR of the CORBA server.
   **/
  private URL ior_url = null;

  /**
   *  The JFrame that was passed to the constructor.
   **/
  private JFrame frame = null;
}
