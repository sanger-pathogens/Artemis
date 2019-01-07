/* DocumentFactory.java
 *
 * created: Fri Jul 11 2003
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2003  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/util/DocumentFactory.java,v 1.3 2006-03-02 19:46:41 tjc Exp $
 */

package uk.ac.sanger.artemis.util;

import java.net.*;
import java.io.*;

import uk.ac.sanger.artemis.components.filetree.RemoteFileNode;

/**
 *  A Factory for Document objects.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: DocumentFactory.java,v 1.3 2006-03-02 19:46:41 tjc Exp $
 **/

public class DocumentFactory {
  /**
   *  Create a new Document of the appropriate type from the given String.
   *  eg. FileDocument for files or a URLDocument for a URL.
   **/
  public static Document makeDocument (final String source_string) {
    if (source_string.indexOf ("://") != -1) {
      try {
        return new URLDocument (new URL (source_string));
      } catch (MalformedURLException e) {
        return new FileDocument (new File (source_string));
      }
    } 
    else 
    {
      File file = new File (source_string);
      if(file.exists())                    // assume a local file
        return new FileDocument(file);
      else                                 // assume a remote file
      {
        int index;
        String parent = source_string;

        if( (index = parent.lastIndexOf("/")) > -1)
          parent = parent.substring(0,index);
        else
          parent = file.getParent();

        RemoteFileNode node = new RemoteFileNode("", file.getName(), null,
                                                 parent, false);
        return new RemoteFileDocument(node);
      }  
    }
  }
}
