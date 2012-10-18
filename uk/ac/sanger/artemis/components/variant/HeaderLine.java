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

import java.util.Hashtable;


public class HeaderLine
{
  private int headerType;
  private String ID;
  private String type;
  private String description;
  private int number = 1;
  private String numberString = null;
  private String headerTypeStr;
  private String origLine;
  private boolean isFlag = false;
  
  protected static final String[] filterFlagStr = new String[]{
    "INFO", "FORMAT", "FILTER", "FILTER_QUAL", "FILTER_NV_FLAG",
    "FILTER_OVERLAP_FLAG", "FILTER_MULTALL_FLAG", 
    "FILTER_INS", "FILTER_DEL",
    "FILTER_NONSYN", "FILTER_SYN", "FILTER_MANUAL", "FILTER_HOMOZYG"};
  protected static final int INFO_LINE   = 0;
  protected static final int FORMAT_LINE = 1;
  protected static final int FILTER_LINE = 2;
  protected static final int FILTER_QUAL = 3;
  protected static final int FILTER_NV_FLAG = 4;
  protected static final int FILTER_OVERLAP_FLAG = 5;
  protected static final int FILTER_MULTALL_FLAG = 6;
  protected static final int FILTER_INS = 7;
  protected static final int FILTER_DEL = 8;
  protected static final int FILTER_NONSYN = 9;
  protected static final int FILTER_SYN = 10;
  protected static final int FILTER_MANUAL = 11;
  protected static final int FILTER_HOMOZYG = 12;
  
  public HeaderLine(final String origLine, String headerTypeStr, Hashtable<String, String> lineHash)
  {
    this.origLine = origLine;
    this.ID = lineHash.get("ID");
    this.type = lineHash.get("Type");
    this.description = lineHash.get("Description");
    this.headerTypeStr = headerTypeStr;
    
    for(int i=0; i<filterFlagStr.length; i++)
    {
      if(headerTypeStr.equals(filterFlagStr[i]))
        headerType = i;
    }

    if(type !=null && type.equals("Flag"))
      isFlag = true;
    
    String numStr = lineHash.get("Number");
    if (numStr != null)
    {
      try
      {
        number = Integer.parseInt(numStr);
      }
      catch(NumberFormatException e)
      {
        // '.' - number of possible values varies, is unknown, or is unbounded
        // 'A' - one value per alternate allele
        // 'G' - one value for each possible genotype (more relevant to the FORMAT tags)
        numberString = numStr;
        number = 1;
      }
    }
  }
  
  /**
   * Type of line : FORMAT, INFO, FILTER
   * @return
   */
  protected int getHeaderType()
  {
    return headerType;
  }
  
  protected String getHeaderTypeStr()
  {
    return headerTypeStr;
  }
  
  protected String getID()
  {
    return ID;
  }
  
  protected String getType()
  {
    return type;
  }
  
  protected boolean isFlag()
  {
    return isFlag;
  }
  
  protected int getNumber()
  {
    return number;
  }
  
  protected String getDescription()
  {
    return description;
  }
  
  /**
   * 
   * @return
   */
  protected String getNumberString()
  {
    return numberString;
  }
  
  public String toString()
  {
    return origLine;
  }
}