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
  protected String chrom;
  protected int pos;
  protected float quality;
  protected String info;
  protected String format;
  protected String data[][];
  
  public String toString()
  {
    return chrom+"\t"+pos+"\t"+ID+"\t"+ref+"\t"+alt+"\t"+quality+
           "\t"+filter+"\t"+info+"\t"+format+"\t"+getSampleDataString();
  }
  
  protected int getNumAlleles()
  {
    if (alt.equals(".")) 
      return 1;

    return alt.split(",").length+1;
  }
  
  /**
   * Parse a VCF line and return a VCFRecord
   * @param line
   * @return
   */
  protected static VCFRecord parse(String line)
  {
    VCFRecord rec = new VCFRecord();
    String parts[] = line.split("\\t");

    rec.chrom = parts[0];
    rec.pos   = Integer.parseInt(parts[1]);
    rec.ID    = parts[2];
    rec.ref   = parts[3];
    rec.alt   = parts[4];
    rec.quality = Float.parseFloat(parts[5]);
    rec.filter  = parts[6];
    rec.info    = parts[7];
    
    if(parts.length > 9)
    {
      rec.format  = parts[8].trim();
      int nsamples = parts.length-9;
      int nfmt = rec.format.split(":").length;
      
      rec.data = new String[nsamples][nfmt];
      for(int i=0; i<nsamples; i++)
      {
        String data[] = parts[9+i].split(":");
        rec.data[i] = data;
      }
    }
    return rec;
  }
  
  /**
   * Return the sample data as a tab-delimited string
   * @return
   */
  private String getSampleDataString()
  {
    if(data == null)
      return "";
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