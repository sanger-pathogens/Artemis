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

import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.RangeVector;
import uk.ac.sanger.artemis.sequence.AminoAcidSequence;
import uk.ac.sanger.artemis.sequence.Bases;

class VCFRecord
{
  private String chrom;
  private int pos;
  private String ID;
  private String ref;
  private String alt;
  private float quality;
  private String filter;
  private String info;
  private String format;
  private String data[][];
  private short synFlag = -1;
  

  /**
   * Return the string representation of the VCF record as a
   * tab-delimited string.
   */
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
    
    try
    {
      rec.quality = Float.parseFloat(parts[5]);
    }
    catch(NumberFormatException e)
    {
      rec.quality = 0.f;
    }
    
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
   * For example DP or MQ
   * @param key
   * @return
   */
  protected String getInfoValue(String key)
  {
    String parts[] = info.split(";");
    for(int i=0; i<parts.length; i++)
      if(parts[i].startsWith(key+"="))
        return parts[i].substring(key.length()+1);
    return null;
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

  /**
   * @return the chrom
   */
  protected String getChrom()
  {
    return chrom;
  }

  /**
   * @param chrom the chrom to set
   */
  protected void setChrom(String chrom)
  {
    this.chrom = chrom;
  }

  /**
   * @return the pos
   */
  protected int getPos()
  {
    return pos;
  }

  /**
   * @param pos the pos to set
   */
  protected void setPos(int pos)
  {
    this.pos = pos;
  }

  /**
   * @return the iD
   */
  protected String getID()
  {
    return ID;
  }

  /**
   * @param iD the iD to set
   */
  protected void setID(String iD)
  {
    ID = iD;
  }

  /**
   * @return the ref
   */
  protected String getRef()
  {
    return ref;
  }

  /**
   * @param ref the ref to set
   */
  protected void setRef(String ref)
  {
    this.ref = ref;
  }

  /**
   * @return the alt
   */
  protected String getAlt()
  {
    return alt;
  }

  /**
   * @param alt the alt to set
   */
  protected void setAlt(String alt)
  {
    this.alt = alt;
  }

  /**
   * @return the quality
   */
  protected float getQuality()
  {
    return quality;
  }

  /**
   * @param quality the quality to set
   */
  protected void setQuality(float quality)
  {
    this.quality = quality;
  }

  /**
   * @return the filter
   */
  protected String getFilter()
  {
    return filter;
  }

  /**
   * @param filter the filter to set
   */
  protected void setFilter(String filter)
  {
    this.filter = filter;
  }

  /**
   * @return the info
   */
  protected String getInfo()
  {
    return info;
  }

  /**
   * @param info the info to set
   */
  protected void setInfo(String info)
  {
    this.info = info;
  }

  /**
   * @return the format
   */
  protected String getFormat()
  {
    return format;
  }

  /**
   * @param format the format to set
   */
  protected void setFormat(String format)
  {
    this.format = format;
  }

  /**
   * @return the data
   */
  protected String[][] getData()
  {
    return data;
  }

  /**
   * @param data the data to set
   */
  protected void setData(String[][] data)
  {
    this.data = data;
  }


  /**
   * @param features
   * @param basePosition
   * 0 if non-synonymous;
   * 1 if synonymous;  
   * 2 if non-synonymous and creates a stop codon
   */
  protected short getSynFlag(FeatureVector features, int basePosition)
  {
    if(synFlag == -1)
      this.synFlag = isSynonymous(features, basePosition);
    return synFlag;
  }

  /**
   * @param features
   * @param basePosition
   * @return
   * 0 if non-synonymous;
   * 1 if synonymous;  
   * 2 if non-synonymous and creates a stop codon
   * 3 not within a gene
   */
  private short isSynonymous(FeatureVector features, int basePosition)
  {
    char variant = getAlt().toLowerCase().charAt(0);
    int intronlength = 0;
    Range lastRange = null;
    
    for(int i = 0; i<features.size(); i++)
    {
      Feature feature = features.elementAt(i);
      if(feature.getRawFirstBase() < basePosition && feature.getRawLastBase() > basePosition)
      {
        RangeVector ranges = feature.getLocation().getRanges();
        intronlength = 0;

        for(int j=0; j< ranges.size(); j++)
        {
          Range range = (Range) ranges.get(j);
          
          if(j > 0)
          {
            if(feature.isForwardFeature())
              intronlength+=range.getStart()-lastRange.getEnd()-1;
            else
              intronlength+=lastRange.getStart()-range.getEnd()-1;
            
            if(intronlength < 0)
              intronlength = 0;
          }
          
          if(range.getStart() < basePosition && range.getEnd() > basePosition)
          {
            int mod;
            int codonStart;
            
            if(feature.isForwardFeature())
            {
              mod = (basePosition-feature.getRawFirstBase())%3;
              codonStart = basePosition-feature.getRawFirstBase()-mod;
            }
            else
            {
              mod = (feature.getRawLastBase()-basePosition)%3;
              codonStart = feature.getRawLastBase()-basePosition-mod;
            }

            codonStart-=intronlength;
            
            try
            {
              if(codonStart+3 > feature.getBases().length())
                return 0;
              char codon[] = feature.getBases().substring(codonStart,
                  codonStart + 3).toLowerCase().toCharArray();

              // String oldBase = new String(codon);
              char aaRef = AminoAcidSequence.getCodonTranslation(codon[0],
                  codon[1], codon[2]);

              if(!feature.isForwardFeature())
                variant = Bases.complement(variant);
              codon[mod] = variant;
              char aaNew = AminoAcidSequence.getCodonTranslation(codon[0],
                  codon[1], codon[2]);

              if (aaNew == aaRef) 
                return 1;
              else if(AminoAcidSequence.isStopCodon(aaNew))
                return 2;
              else
                return 0;
            }
            catch(Exception e)
            {
              for(int k=0; k<ranges.size(); k++)
                System.out.println(k+" "+ ((Range)ranges.get(k)).getStart() );
              
              System.out.println(feature.getIDString()+"  "+codonStart+" "+intronlength+" basePosition="+basePosition+" segment="+range.getStart()+".."+range.getEnd()+" mod="+mod);
              throw new RuntimeException(e);
            }
          }

          lastRange = range;
        }
      }
    }
    
    return 3;
  }
  
}