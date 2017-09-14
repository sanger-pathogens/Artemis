/* 
 * This file is part of Artemis
 *
 * Copyright (C) 2014  Genome Research Limited
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
package uk.ac.sanger.artemis.components.alignment;

class ReadCount
{
  protected float senseCnt = 0;
  protected float antiCnt  = 0;
  ReadCount(float[] cnt, boolean isFwd)
  {
    if(isFwd)
    {
      senseCnt = cnt[0];
      antiCnt  = cnt[1];
    }
    else
    {
      senseCnt = cnt[1];
      antiCnt  = cnt[0];
    }
  }
}