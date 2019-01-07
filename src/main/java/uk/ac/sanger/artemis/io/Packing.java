/* Packing.java
 *
 * created: Jun 2005
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

package uk.ac.sanger.artemis.io;

class Packing
{
  static char bases[] = { 'u',
                          'a',
                          'g',
                          'm',
                          'c',
                          'r',
                          'w',
                          's',
                          't',
                          'y',
                          'k',
                          'b',
                          'd',
                          'h',
                          'v',
                          'n' };

  public Packing(char dna[])
  {
  }

  protected static byte pack(char c) 
  {
    if(c == 'a')
      return 1;
    else if(c == 'g')
      return 2;
    else if(c == 'c')
      return 4;
    else if(c == 't')
      return 8;
    else if(c == 'n')
      return 15;

    for(int i = 0; i < bases.length; i++)
    {
      if(c == bases[i])
        return (byte)i;
    }

    return 15;
  }

  protected static char unpack(int i)
  {
    return bases[i];
  }

  public byte wordSize() 
  {
    return 4;
  }


}
