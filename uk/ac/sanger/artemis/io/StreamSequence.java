/* StreamSequence.java (formally ReaderSequence.java)
 *
 * created: Wed Dec 30 1998
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1998-2005  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/StreamSequence.java,v 1.15 2008-12-11 15:43:31 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import java.io.IOException;
import java.io.Writer;

/**
 *  This is an implementation of Sequence that can read and write itself to a
 *  stream.
 *
 *  Sequence stored in 4 bit chunks.
 *
 *  @author Kim Rutherford
 *  @version $Id: StreamSequence.java,v 1.15 2008-12-11 15:43:31 tjc Exp $
 **/

public abstract class StreamSequence
    extends LineGroup implements Sequence 
{
  
  /**
   *  Return a new StreamSequence object that is a copy of this one.
   **/
  abstract public StreamSequence copy();

  /**
   *  Return the sequence type (one of EMBL_FORMAT, RAW_FORMAT, FASTA_FORMAT,
   *  etc.)
   **/
  abstract public int getFormatType();

  /**
   *  Contains the sequence data for this object.  It will contain the bases
   *  of the sequence with no spaces after the Feature constructor finishes.
   **/
  private byte[] sequencePacked;

  /** Count of the a bases in the sequence. */
  private int a_count = 0;

  /** Count of the c bases in the sequence. */
  private int c_count = 0;

  /** Count of the g bases in the sequence. */
  private int g_count = 0;

  /** Count of the t bases in the sequence. */
  private int t_count = 0;

  /** char array for returning sequence chunks */
  private char[] dst = null;

  private int sequence_length;
  private char bases[];

  /**
   *  Return a the given range of bases as a String.  Returns an empty
   *  sequence if the end position is less than the start position.
   *  @param start The start base of the range.
   *  @param end The end base of the range.
   **/
  public String getSubSequence(int start, int end) 
  {
    if(end < start)    // empty range
      return "";
    else 
    {
      final char[] c = getCharSubSequence(start,end);
      if(end > length())
        end = length();

      return new String(c,0,end-start+1);
    }
  }

  public void forceReset()
  {
    dst = null;
  }

  public char[] getCharSubSequence(int start, int end) 
  {
    char[] this_dst = null;
    if(end-start > 1000)
      this_dst = dst;

    int dst_length = 0;
    if(this_dst != null)
      dst_length = this_dst.length;

    if(this_dst == null || dst_length < end-start+1 || end >= length()) 
    {
      if(end-start > 1000)
      {
        dst = new char[end-start+1];
        this_dst = dst;
      }
      else
        this_dst = new char[end-start+1];
      
      dst_length = this_dst.length;
//    dst = new char[end-start+3];
//    System.out.println("REALLOCATE "+ this_dst.length);
    }

//  int packStart = Math.round( (float)start/2.f ) - 1;
    int packStart = (start - 1) >> 1;
    int packEnd   = (int) Math.round( (double)end/2.d ); // end/2;
    
    int count = 0;
    byte currStorageUnit;
    int index1;
    int index2;

    // skip first four bits
    if(start % 2 == 0)
    {
      currStorageUnit = sequencePacked[packStart];
      index1 = (int)(currStorageUnit & 0x000F);
      this_dst[count] = bases[index1];

      packStart++;
      packEnd++;
      count++;
    }

    for(int i=packStart; i <= packEnd && count < dst_length && i < sequencePacked.length; i++) 
    {
      currStorageUnit = sequencePacked[i];
      index1 = (int)(currStorageUnit & 0x000F);
      index2 = (int)( (currStorageUnit >> 4) & 0x000F);

      try
      {
        this_dst[count] = bases[index2];
      }
      catch(ArrayIndexOutOfBoundsException oob)
      {
//      System.out.println("start "+start+" end "+end+"  packStart "+packStart+" packEnd "+packEnd);
//      System.out.println("count "+count+"  this_dst.length "+this_dst.length+"  index2 "+index2);
//      oob.printStackTrace();
      }
      count++;

      if(count < dst_length)
      {
        this_dst[count] = bases[index1]; 
      }

      count++;
    }

    return this_dst;
  }

  public char[] getCharSequence()
  {
    char dst[] = new char[length()];
    int packEnd = (int)Math.round( (double)length()/2.d );
    int count = 0;
    byte currStorageUnit = 0;

    for(int i=0; i < packEnd; i++)
    {
      currStorageUnit = sequencePacked[i];
      int index1 = (int)(currStorageUnit & 0x000F);
      currStorageUnit = (byte) (currStorageUnit>>4);
      int index2 = (int)(currStorageUnit & 0x000F);

      dst[count]   = bases[index2];
      count++;

      if(count < dst.length)
        dst[count] = bases[index1];

//    System.out.print(Packing.unpack(index2));
//    if(count < dst.length)
//      System.out.print(Packing.unpack(index1));

      count++;
    }
//  System.out.println("\n"+currStorageUnit);

    return dst;  
  }

  public char charAt(final int i)
  {
    final int packStart = (i-1) >> 1;
    final byte currStorageUnit = sequencePacked[packStart];
    final int index;

    if(i % 2 == 0)
      index = (int)(currStorageUnit & 0x000F);
    else
      index = (int)( (currStorageUnit >> 4) & 0x000F);
  
    return bases[index];
  }

  public void setFromChar(final char dna[])
  { 
    sequence_length = dna.length;
    int numBytes = (int)Math.round( (double)sequence_length/2.d );

    sequencePacked = new byte[numBytes];
    setFromChar(dna, 0, 0, sequence_length);
    setCounts(dna);
  }

  /** 
   *
   *  Set this sequence to hold the bases in the given byte array.
   *
   **/
  private void setFromChar(final char dna[], int offset, 
                          final int bit_shift,
                          final int new_sequence_length)
  {
    int offsetSize    = offset >> 1;
    int bytePointer   = offset >> 1;
    int symbolPointer = 0;
    int numBytes      = (int)Math.round( (double)dna.length/2.d );

    // filled last unit if = 0
    int filledLastUnit = dna.length & 0x0001;
    byte currByte;
    if(bit_shift != 0)
    {
//    System.out.println("Bit Shifty..."+offset);
      currByte = sequencePacked[sequence_length>>1];
      currByte = (byte)(currByte | Packing.pack(dna[symbolPointer]));
      sequencePacked[bytePointer] = currByte;
      bytePointer++;
      symbolPointer++;
      if(filledLastUnit == 0)
        offsetSize += 1;
    }
    else 
      currByte = 0; 

    for(int i=bytePointer; i<numBytes+offsetSize; i++) 
    {
      // each byte consists of 4 bit nibbles
      // process each separately
      for(int j=0; j < 2; j++) 
      {
        currByte = (byte)(currByte<<4 | Packing.pack(dna[symbolPointer]));

        symbolPointer++;

        if(j == 0 && symbolPointer == dna.length) 
        {
          currByte = (byte)(currByte<<4);
          break;
        }
      }
      sequencePacked[bytePointer] = currByte;
//    System.out.println(" "+symbolPointer+" bytePointer " + bytePointer + " to " +
//                          Integer.toString((int) currByte, 2) + " decimal " +
//                          Integer.toString((int) currByte, 10) + " hex " +
//                          Integer.toString((int) currByte, 16));
      bytePointer++;
      currByte = 0;
    }

    if(bases == null)
      bases = Packing.bases;

  }

  protected void appendChar(final char dna[])
  {
    int newlength = sequence_length + dna.length;
    int numBytes  = (int)Math.round( (double)newlength/2.d );
    if(numBytes > capacity())
      expandCapacity(numBytes);

//  System.out.println("capacity "+capacity()+" numBytes "+numBytes);
    // filled last unit if = 0
    int filledLastUnit = sequence_length & 0x0001;

    setFromChar(dna, sequence_length, filledLastUnit, newlength);
    sequence_length = newlength;
  }

  protected void setSequencePackingCapacity(final int n)
  {
    int numBytes  = Math.round( n/2.f );
    sequencePacked = new byte[numBytes];
  }


  /**
  * This implements the expansion semantics of ensureCapacity but is
  * unsynchronized for use internally by methods which are already
  * synchronized.
  *
  * @see java.lang.StringBuffer#ensureCapacity(int)
  */
  private void expandCapacity(int minimumCapacity) 
  {
    int newCapacity = (sequencePacked.length + 1) * 2;
    if(newCapacity < 0) 
      newCapacity = Integer.MAX_VALUE;
    else if(minimumCapacity > newCapacity) 
      newCapacity = minimumCapacity;
	
//  System.out.println("EXPANDING.... "+newCapacity);
    byte newValue[] = new byte[newCapacity];
    System.arraycopy(sequencePacked, 0, newValue, 0, sequencePacked.length);
    sequencePacked = newValue;
  }

  /**
  * Returns the current capacity of the String buffer. The capacity
  * is the amount of storage available for newly inserted
  * characters; beyond which an allocation will occur.
  *
  * @return  the current capacity of this string buffer.
  */
  private synchronized int capacity() 
  {
    return sequencePacked.length;
  }


  /**
   *  Write this Sequence to the given stream.
   *  @param writer The stream to write to.
   **/
  public abstract void writeToStream(final Writer writer)
      throws IOException;

  /**
   *  Returns the length of the sequence in bases.
   **/
  public int length() 
  {
    return sequence_length;
  }

  /**
   *  Return the count of c bases in the whole sequence.
   **/
  public int getCCount()  
  {
    return c_count;
  }

  /**
   *  Return the count of g bases in the whole sequence.
   **/
  public int getGCount()
  {
    return g_count;
  }

  /**
   *  Return the count of a bases in the whole sequence.
   **/
  public int getACount() 
  {
    return a_count;
  }

  /**
   *  Return the count of t bases in the whole sequence.
   **/
  public int getTCount() 
  {
    return t_count;
  }

  /**
   *  Return the count of non-g,c,t,a bases in the whole sequence.
   **/
  public int getOtherCount() 
  {
    return
      length() - (getCCount() + getACount() + getTCount() + getGCount());
  }


  /**
   *  Set the a_count, c_count, t_count and g_count variables.
   **/
  protected void setCounts()
  {
    final int len = length();
    final int packEnd = (int)Math.round( (double)len/2.d );
    int count = 0;
    byte currStorageUnit;
    int index1;
    int index2;

    for(int i=0; i < packEnd; i++)
    {
      currStorageUnit = sequencePacked[i];
      index1 = currStorageUnit & 0x000F;
      index2 = (currStorageUnit>>4) & 0x000F;

      counter(bases[index2]);
      count++;

      if(count < len)
        counter(bases[index1]);

      count++;
    }
  }


  /**
   *  Set the a_count, c_count, t_count and g_count variables.
   **/
  private void setCounts(char[] sequence_chars) 
  {
    a_count = c_count = t_count = g_count = 0;

    for(int i = 0 ; i < length(); ++i) 
      counter(sequence_chars[i]);
  }

  private void counter(char c)
  {
    switch(c)
    {
      case 'a':
        ++a_count;
        break;
      case 'c':
        ++c_count;
        break;
      case 'g':
        ++g_count;
        break;
      case 't':
        ++t_count;
        break;
      default:
        break;
    }
  }
}
