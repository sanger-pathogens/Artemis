/*
 * created: 2010
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

class VCFRecord
{
  protected String ID;
  protected String alt;
  protected String ref;
  protected String filter;
  protected String seqID;
  protected int pos;
  protected float quality;
  protected String info;
  protected String format;
  protected String data[][];
  
  public String toString()
  {
    return seqID+"\t"+pos+"\t"+ID+"\t"+ref+"\t"+alt+"\t"+quality+
           "\t"+filter+"\t"+info+"\t"+format+"\t"+getSampleDataString();
  }
  
  /**
   * Return the sample data as a tab-delimited string
   * @return
   */
  private String getSampleDataString()
  {
    StringBuffer buff = new StringBuffer();
    for(int i=0; i<data.length; i++)       // loop over samples
    {
      for(int j=0; j<data[i].length; j++)  // loop over values
      {
        buff.append(data[i][j]);
        if(j<data[i].length-1)
          buff.append(":");
      }
      if(i<data.length-1)
        buff.append("\t");
    }
    return buff.toString();
  }
}