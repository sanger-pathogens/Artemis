/* InputStreamProgressListener.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/util/InputStreamProgressListener.java,v 1.2 2005-04-01 16:08:23 tjc Exp $
 */

package uk.ac.sanger.artemis.util;

/**
 *  This interface is implemented by classes which need to know when progress
 *  is made while reading from an InputStream.  "progress" means some more
 *  bytes have been read.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: InputStreamProgressListener.java,v 1.2 2005-04-01 16:08:23 tjc Exp $
 **/

public interface InputStreamProgressListener 
{
  void progressMade(final InputStreamProgressEvent event);
  void progressMade(final String progress);
}
