/* ByteBuffer.java
 *
 * created: 2005
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2005  Genome Research Limited
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
 */

package uk.ac.sanger.artemis.util;

/**
*
* Buffer for appending to a byte array.
*
*/
public class ByteBuffer
{

  private byte buff[];
  private static int byteCapacity = 128;
  private int count = 0;

  /** 
  * Creates new ByteBuffer
  */
  public ByteBuffer() 
  {
  }
 
  public int size()
  {
    return count;
  }
   
  /**
  * Appends the subarray of the <CODE>byte</CODE> array. 
  * @param b the array to be appended
  * @param off the offset to the start of the array
  * @param len the length of bytes to append
  * @return a reference to this <CODE>ByteBuffer</CODE> object
  */
  public void append(byte b[], int off, int len) 
  {
    int newcount = count + len;

    if(buff == null || newcount > buff.length) 
    {
      if(buff == null)
      {
        byteCapacity = newcount;
        buff = new byte[byteCapacity];
      }
      else
      {
        byteCapacity = Math.max(buff.length << 1, newcount);
        byte newbuff[] = new byte[byteCapacity];
        System.arraycopy(buff, 0, newbuff, 0, count);
        buff = newbuff;
      }
    }
    System.arraycopy(b, off, buff, count, len);

    count = newcount;
  }

  /**
  * Appends the byte array to the end of the buffer.
  */
  public void append(byte b[])
  {
    append(b, 0, b.length);
  }

  /**
  * Convert the string to a byte array and appends 
  * to the end of the buffer.
  */
  public void append(String s)
  {
    byte b[] = s.getBytes();
    append(b);
  }

  /**
  * Appends the ByteBuffer to the end of the buffer.
  */
  public void append(ByteBuffer newbuff)
  {
    append(newbuff.getBytes());
  }

  /**
  * Get the byte[] that has been filled.
  * @return the byte array
  */
  public byte[] getBytes()
  { 
    byte[] newbuff = new byte[count];
    System.arraycopy(buff, 0, newbuff, 0, count);
    return newbuff;
  }
}

