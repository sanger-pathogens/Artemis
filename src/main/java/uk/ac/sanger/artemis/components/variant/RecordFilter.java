/*
 * created: 2011
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2011  Genome Research Limited
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

  public class RecordFilter
  {
    /** associated header line */
    private HeaderLine hLine;
    /** number of possible values */
    private int NUMBER;
    /** min/max filter values (for integer) */
    protected int minIVal[];
    protected int maxIVal[];
    
    /** min/max filter values (for float) */
    protected float minFVal[];
    protected float maxFVal[];
    
    RecordFilter(HeaderLine hLine, int NUMBER)
    {
      this.hLine = hLine;
      this.NUMBER = NUMBER;
      
      if (hLine.getType().equals("Integer"))
      {
        minIVal = new int[NUMBER];
        maxIVal = new int[NUMBER];
        
        for(int i=0; i<NUMBER; i++)
        {
          minIVal[i] = Integer.MIN_VALUE;
          maxIVal[i] = Integer.MAX_VALUE;
        }
      }
      else if(hLine.getType().equals("Float"))
      {
        minFVal = new float[NUMBER];
        maxFVal = new float[NUMBER];
        
        for(int i=0; i<NUMBER; i++)
        {
          minFVal[i] = Float.MIN_VALUE;
          maxFVal[i] = Float.MAX_VALUE;
        }
      }
    }
    
    protected boolean pass(final VCFRecord record, final String valStr[], final AbstractVCFReader vcfReader)
    {
      String numStr = hLine.getNumberString();
      if(numStr == null)
        return pass(valStr);

      // numStr - if this is not a number it can be:
      // '.' - number of possible values varies, is unknown, or is unbounded
      // 'A' - one value per alternate allele
      // 'G' - one value for each possible genotype

/*      int nvals = 0;
      if (numStr.equals("A"))
        nvals = record.getAlt().getNumAlleles();
      else if(numStr.equals("G"))
        nvals = vcfReader.getNumberOfSamples();
      else
        nvals = valStr.length;*/
      
      for (int i = 0; i < valStr.length; i++)
      {
        try
        {
          if (hLine.getType().equals("Integer"))
          {
            int val = Integer.parseInt(valStr[i]);
            if (val < minIVal[0] || val > maxIVal[0])
              return false;
          }
          else if (hLine.getType().equals("Float"))
          {
            float val = Float.parseFloat(valStr[i]);
            if (val < minFVal[0] || val > maxFVal[0])
              return false;
          }
        }
        catch(NumberFormatException nfe)
        {
          return false;
        }
      }
      
      return true;
    }
    
   
    /**
     * For a fixed number of values check the min and max
     * values.
     * @param valStr
     * @return
     */
    private boolean pass(final String valStr[])
    {
      for (int i = 0; i < NUMBER; i++)
      {
        if (hLine.getType().equals("Integer"))
        {
          int val = Integer.parseInt(valStr[i]);
          if (val < minIVal[i] || val > maxIVal[i])
            return false;
        }
        else if (hLine.getType().equals("Float"))
        {
          float val = Float.parseFloat(valStr[i]);
          if (val < minFVal[i] || val > maxFVal[i])
            return false;
        }
      }
      return true;
    }
    
    protected HeaderLine getHeaderLine()
    {
      return hLine;
    }
    
    public String toString()
    { 
      if(hLine.isFlag())
        return hLine.getDescription();
      
      StringBuffer buff = new StringBuffer(hLine.getHeaderTypeStr());
      buff.append(": ");
      buff.append(hLine.getDescription());
      buff.append(", ");
      
      final String id = (hLine.getHeaderType() == HeaderLine.FORMAT_LINE ? "sample" : "") +
          hLine.getID();
      
      for(int i=0; i<NUMBER; i++)
      {
        if(i > 0)
          buff.append(" : ");
        if (hLine.getType().equals("Integer"))
        {
          if(minIVal[i] > Integer.MIN_VALUE) 
          {  
            buff.append(id+" < "+minIVal[i]);
            if(maxIVal[i] < Integer.MAX_VALUE)
              buff.append(" ");
          }
          if(maxIVal[i] < Integer.MAX_VALUE)
            buff.append(id+" > "+maxIVal[i]);
        }
        else if (hLine.getType().equals("Float"))
        {
          if(minFVal[i] > Float.MIN_VALUE)
          {
            buff.append(id+" < "+minFVal[i]);
            if(maxFVal[i] < Float.MAX_VALUE)
              buff.append(" ");
          }
          if(maxFVal[i] < Float.MAX_VALUE)
            buff.append(id+" > "+maxFVal[i]);;
        }
      }
      return buff.toString();
    }
  }