/* AbstractReader
 *
 * created: July 2011
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2010  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or(at your option) any later version.
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
 */

package uk.ac.sanger.artemis.components.variant;

import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

public abstract class AbstractVCFReader
{
  protected abstract String[] getSeqNames();
  protected abstract String getFileName();
  
  protected static int readInt(final InputStream is) throws IOException {
    byte[] buf = new byte[4];
    is.read(buf);
    return ByteBuffer.wrap(buf).order(ByteOrder.LITTLE_ENDIAN).getInt();
  }
  
  
  protected static float readFloat(final InputStream is) throws IOException {
    byte[] buf = new byte[4];
    is.read(buf);
    return ByteBuffer.wrap(buf).order(ByteOrder.LITTLE_ENDIAN).getFloat();
  }

  protected static long readLong(final InputStream is) throws IOException {
    byte[] buf = new byte[8];
    is.read(buf);
    return ByteBuffer.wrap(buf).order(ByteOrder.LITTLE_ENDIAN).getLong();
  }
}