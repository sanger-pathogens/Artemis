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

import java.util.List;
import java.util.regex.Pattern;

import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.RangeVector;
import uk.ac.sanger.artemis.sequence.AminoAcidSequence;
import uk.ac.sanger.artemis.sequence.Bases;

public class VCFRecord
{
  //private static Logger logger = Logger.getLogger(VCFRecord.class);
  private String chrom;
  private int pos;
  private String ID;
  private String ref;
  private VariantBase var;
  private float quality;
  private String filter;
  private String info;
  private String infos[];
  private String format;
  private String genotypeData[][];
  private short synFlag = -1;
  private boolean markAsNewStop = false;

  protected static Pattern MULTI_ALLELE_PATTERN = Pattern.compile(
      "^[AGCTNMRWSYKBDHVagctnmrwsykbdhv]+,[AGCTNMRWSYKBDHVagctnmrwsykbdhv,]+$");
  protected static Pattern COLON_PATTERN = Pattern.compile(":");
  protected static Pattern SEMICOLON_PATTERN = Pattern.compile(";");
  protected static Pattern TAB_PATTERN = Pattern.compile("\\t");

  /**
   * Return the string representation of the VCF record as a
   * tab-delimited string.
   */
  public String toString()
  {
    return chrom+"\t"+pos+"\t"+ID+"\t"+ref+"\t"+var.toString()+"\t"+quality+
           "\t"+filter+"\t"+info+"\t"+format+"\t"+getSampleDataString();
  }
 
  
  /**
   * Parse a VCF line and return a VCFRecord
   * @param line
   * @return
   */
  protected static VCFRecord parse(final String line, int nsamples)
  {
    final VCFRecord rec = new VCFRecord();
    final String parts[] = split(line, "\t", 9+nsamples);
    //final String parts[] = TAB_PATTERN.split(line);
    
    rec.chrom = parts[0];
    rec.pos   = Integer.parseInt(parts[1]);
    rec.ID    = parts[2];
    rec.ref   = parts[3];
    rec.var   = new VariantBase(rec, parts[4]);
    
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
      rec.format  = (parts[8]).trim();
      final int nfmt = countOccurrences(rec.format, ':')+1; //rec.format.split(":").length;
      nsamples = parts.length-9;
          
      rec.genotypeData = new String[nsamples][nfmt];
      for(int i=0; i<nsamples; i++)
      {
        //rec.genotypeData[i] = COLON_PATTERN.split(parts[9+i]);
        rec.genotypeData[i] = split(parts[9+i], ":", nfmt);
      }
    }
    return rec;
  }
  
  protected static int countOccurrences(final String str, final char search)
  {
    int count = 0;
    for(int i=0; i < str.length(); i++)
    {
      if(str.charAt(i) == search)
        count++;
    }
    return count;
  }

  /**
   * Split a string into an array
   * @param arg
   * @param splitChar
   * @param nsize
   * @return
   */
  protected static String[] split(final String argStr, final String splitChar, final int nsize)
  {
    final String str[] = new String[nsize];
    String value;

    int ind1 = 0;
    int ind2;
    int count = 0;
    int argLen = argStr.length();

    while(ind1 < argLen)
    {
      ind2 = argStr.indexOf(splitChar,ind1);
      if(ind2 == ind1)
      {
        ind1++;
        continue;
      }

      if(ind2 < 0)
        ind2 = argLen;
 
      value = argStr.substring(ind1,ind2);
      ind1 = ind2+1;

      str[count] = value;
      count++;
    }
    
    // shrink array if there are fewer elements
    if(count < nsize)
    {
      String tmp[] = new String[count];
      System.arraycopy( str, 0, tmp, 0, count );
      return tmp;
    }
   
    return str;
  }

  /**
   * For example DP or MQ
   * @param key
   * @return
   */
  protected String getInfoValue(String key)
  {
    if(infos == null)
      infos = SEMICOLON_PATTERN.split(info);
    for(int i=0; i<infos.length; i++)
      if(infos[i].startsWith(key+"="))
        return infos[i].substring(key.length()+1);
    return null;
  }
  
  /**
   * Test if a INFO flag key is present
   * @param key
   * @return
   */
  protected boolean containsInfoFlag(String key)
  {
    if(infos == null)
      infos = SEMICOLON_PATTERN.split(info);
    for(int i=0; i<infos.length; i++)
      if(infos[i].equals(key))
        return true;
    return false;
  }
  
  /**
   * Get genotype values for a given sample.
   * @param sampleIndex
   * @return
   */
  protected String getFormatValueForSample(int sampleIndex)
  {
    if(getFormat() == null)
      return null;
    final StringBuffer buff = new StringBuffer();
    for(int i=0; i<genotypeData[sampleIndex].length; i++)  // loop over values
    {
      buff.append(genotypeData[sampleIndex][i]);
      if(i<genotypeData[sampleIndex].length-1)
        buff.append(":");
    }
    return buff.toString();
  }
  
  /**
   * Get genotype values for a given key within a given sample.
   * @param key
   * @param sampleIndex
   * @return
   */
  protected String getFormatValueForSample(String key, int sampleIndex)
  {
    final String fmtStr[] = getFormatValues(key);
    if(fmtStr == null)
      return null;
    return fmtStr[sampleIndex];
  }

  /**
   * Get genotype values for a given key
   * @param key
   * @return
   */
  protected String[] getFormatValues(final String key)
  {
    if(getFormat() == null)
      return null;
    final String fmts[] = COLON_PATTERN.split(getFormat());

    for(int i=0; i<fmts.length; i++)
    {
      if(fmts[i].equals(key))
      {
        final String keyData[] = new String[genotypeData.length];
        for(int j=0; j<genotypeData.length; j++)
        {
          if(genotypeData[j].length == fmts.length)
            keyData[j] = genotypeData[j][i];
        }
        return keyData;
      }
    }
    return null;
  }
  
  /**
   * Return the sample data as a tab-delimited string
   * @return
   */
  protected String getSampleDataString()
  {
    if(genotypeData == null)
      return "";
    StringBuffer buff = new StringBuffer();
    for(int i=0; i<genotypeData.length; i++)       // loop over samples
    {
      for(int j=0; j<genotypeData[i].length; j++)  // loop over values
      {
        buff.append(genotypeData[i][j]);
        if(j<genotypeData[i].length-1)
          buff.append(":");
      }
      if(i<genotypeData.length-1)
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
  protected VariantBase getAlt()
  {
    return var;
  }

  /**
   * @param alt the alt to set
   */
  protected void setAlt(String alt)
  {
    this.var = new VariantBase(this, alt);
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
    if(filter.equals(":"))
      this.filter = ".";
    else
      this.filter = filter;
  }
  
  protected void appendFilter(String filter)
  {
    if(  getFilter().length() == 0 || 
        (getFilter().length() == 1 && getFilter().equals(".")) ||
        (getFilter().length() == 3 && getFilter().equals("PASS")))
      this.filter = filter;
    else
      this.filter += ";" + filter;
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
  protected String[][] getGenoTypeData()
  {
    return genotypeData;
  }

  /**
   * @param data the data to set
   */
  protected void setGenoTypeData(String[][] data)
  {
    this.genotypeData = data;
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
  
  protected short getSynFlag(List<CDSFeature> features, int basePosition)
  {
    //logger.info("getSynFlag(List<CDSFeature>) current syn : " + synFlag + " size? " + features.size());
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
      
    char variant = getAlt().toString().toLowerCase().charAt(0);
    for(int i = 0; i<features.size(); i++)
    {
      Feature feature = features.elementAt(i);
      short isSyn = checkSyn(new CDSFeature(feature), basePosition, variant);
      if(isSyn > - 1)
    	return isSyn;
    }

    return 3;
  }
  
  private short isSynonymous(List<CDSFeature> features, int basePosition) 
  {
      char variant = getAlt().toString().toLowerCase().charAt(0);
      for(CDSFeature feature : features)
      {
        short isSyn = checkSyn(feature, basePosition, variant);
        if(isSyn > - 1)
          return isSyn;
      }

      return 3;
  }
  
  protected static short checkSyn(CDSFeature gfeat, int basePosition, char variant)
  {
    //logger.info("CDSFEATURE\t"+gfeat);
    //logger.info("BASEANDVARIANT\t"+basePosition + "\t" + variant);
    if(gfeat.firstBase < basePosition && gfeat.lastBase > basePosition)
    {
      RangeVector ranges = gfeat.ranges;
      for(int j=0; j< ranges.size(); j++)
      {
        Range range = (Range) ranges.get(j);
        if(j > 0)
        {
          if(gfeat.isFwd)
            gfeat.intronlength+=range.getStart()-gfeat.lastRange.getEnd()-1;
          else
            gfeat.intronlength+=gfeat.lastRange.getStart()-range.getEnd()-1;
            
          if(gfeat.intronlength < 0)
            gfeat.intronlength = 0;
        }
          
        if(range.getStart() < basePosition && range.getEnd() > basePosition)
        {
          int mod;
          int codonStart;
            
          if(gfeat.isFwd)
          {
            mod = (basePosition-gfeat.firstBase-gfeat.intronlength)%3;
            codonStart = basePosition-gfeat.firstBase-gfeat.intronlength-mod;
          }
          else
          {
            mod = (gfeat.lastBase-basePosition-gfeat.intronlength)%3;
            codonStart = gfeat.lastBase-basePosition-gfeat.intronlength-mod;
          }
            
          try
          {
            if(codonStart+3 > gfeat.bases.length())
              return 0;
            char codon[] = gfeat.bases.substring(codonStart,
                codonStart + 3).toLowerCase().toCharArray();

            char aaRef = AminoAcidSequence.getCodonTranslation(codon[0],
                  codon[1], codon[2]);
            //logger.info(String.format("%d %d %s%s%s", mod, codonStart, codon[0],codon[1],codon[2]));
            if(!gfeat.isFwd)
              variant = Bases.complement(variant);
            codon[mod] = variant;
            //logger.info(String.format("%d %d %s%s%s", mod, codonStart, codon[0],codon[1],codon[2]));
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
              
            System.out.println(gfeat.feature.getIDString()+"  "+codonStart+" "+gfeat.intronlength+" basePosition="+basePosition+" segment="+range.getStart()+".."+range.getEnd()+" mod="+mod);
            throw new RuntimeException(e);
          }
        }

        gfeat.lastRange = range;
      }
    }
    return -1;
  }
  
  protected boolean isMarkAsNewStop() 
  {
      return markAsNewStop;
  }

  protected void setMarkAsNewStop(boolean markAsNewStop) 
  {
      this.markAsNewStop = markAsNewStop;
  }
}