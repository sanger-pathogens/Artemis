/* StreamFeature.java (formally ReaderFeature.java)
 *
 * created: Wed Dec 30 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/StreamFeature.java,v 1.1 2004-06-09 09:50:33 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import java.io.IOException;
import java.io.Writer;

/**
 *  This is an implementation of Feature that can read and write itself to a
 *  stream.
 *  @author Kim Rutherford
 *  @version $Id: StreamFeature.java,v 1.1 2004-06-09 09:50:33 tjc Exp $
 **/

public interface StreamFeature extends Feature {
  /**
   *  Write this Feature to the given stream.
   *  @param writer The stream to write to.
   *  @exception IOException thrown if there is an io problem while writing
   *    the Feature.
   **/
  void writeToStream (final Writer writer)
      throws IOException;
}
