/* EmptySequence.java
 * This file is part of Artemis
 *
 * Copyright (C) 2007  Genome Research Limited
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

package uk.ac.sanger.artemis.io;

import org.biojava.bio.symbol.IllegalSymbolException;
import uk.ac.sanger.artemis.util.ReadOnlyException;

/**
 * Used when retrieving just part of a sequence with associated features.
 **/

public class PartialSequence implements Sequence 
{
  private char[] sequence;
  /** Count of the a bases in the sequence. */
  private int a_count = 0;

  /** Count of the c bases in the sequence. */
  private int c_count = 0;

  /** Count of the g bases in the sequence. */
  private int g_count = 0;

  /** Count of the t bases in the sequence. */
  private int t_count = 0;

  public PartialSequence(char[] sequence)
  {
    try
    {
      setFromChar(sequence);
    }
    catch(IllegalSymbolException e)
    {
      e.printStackTrace();
    }
    catch(ReadOnlyException e)
    {
      e.printStackTrace();
    }
  }
  
  public char charAt(int i)
  {
    return sequence[i];
  }

  public int getACount()
  {
    return a_count;
  }

  public int getCCount()
  {
    return c_count;
  }

  public char[] getCharSubSequence(int start, int end)
  {
    int subSeqLength = end-start+1;
    char subSequence[] = new char[subSeqLength];
    int count = 0;
    for(int i=start;i<=end;i++)
      subSequence[count] = sequence[i];
    
    return subSequence;
  }

  public int getGCount()
  {
    return g_count;
  }

  public int getOtherCount()
  {
    return
      length() - (getCCount() + getACount() + getTCount() + getGCount());
  }

  public String getSubSequence(int start, int end)
  {
    return new String(getCharSubSequence(start, end));
  }
  
  public char[] getSequence()
  {
    return sequence;
  }

  public int getTCount()
  {
    return t_count;
  }

  public int length()
  {
    return sequence.length;
  }

  public void setFromChar(char[] sequence) throws ReadOnlyException, IllegalSymbolException
  { 
    this.sequence = sequence;
    setCounts();
  }
  
  /**
   *  Set the a_count, c_count, t_count and g_count variables.
   **/
  private void setCounts() 
  {
    a_count = c_count = t_count = g_count = 0;

    for(int i = 0 ; i < length(); ++i) 
      counter(sequence[i]);
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
