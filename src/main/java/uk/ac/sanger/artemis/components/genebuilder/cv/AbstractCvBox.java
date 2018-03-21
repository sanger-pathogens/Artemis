/* AbstractCvBox.java
 *
 * created: 2007
 *
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
 **/

package uk.ac.sanger.artemis.components.genebuilder.cv;

import java.util.Calendar;
import java.util.Date;

import uk.ac.sanger.artemis.io.QualifierVector;

/**
 * Component for controlled vocabulary qualifier, e.g
 * GO, controlled_curation, product
 */
abstract class AbstractCvBox
{

  /**
   * Strip out the value of a field of interest from a qualifier string
   * 
   * @param fieldName
   * @param qualifierString
   * @return
   */
  protected static String getField(final String fieldName, final String qualifierString)
  {
    String field = "";
    
    int ind1 = qualifierString.toLowerCase().indexOf(fieldName.toLowerCase());
    int ind2 = qualifierString.indexOf(";", ind1);
    
    int len = fieldName.length();

    if(ind2 > ind1 && ind1 > -1)
      field = qualifierString.substring(ind1+len,ind2);
    else if(ind1 > -1)
      field = qualifierString.substring(ind1+len);
    
    return field;
  }
  
  /**
   * Strip out the value of a field of interest from a qualifier string
   * by splitting on '='.
   * @param fieldName
   * @param qualifierString
   * @return
   */
  public static String getFieldIgnoreSeparator(final String fieldName, final String qualifierString)
  {
    String[] completeValues = qualifierString.split(";[\\S&&[^,;=]]+=");   
    StringBuffer buff = null;
    for(int i=0; i<completeValues.length; i++)
    {
      if(qualifierString.indexOf(fieldName+"="+completeValues[i])>-1)
      {
        if(buff == null)
          buff = new StringBuffer();
        else
          buff.append("; ");
        buff.append(completeValues[i]);
      }
    }
    if(buff == null)
      return getField(fieldName, qualifierString);
    return buff.toString();
  }
  
  /**
   * Strip out the value of a field of interest from a qualifier string
   * 
   * @param fieldName
   * @param qualifierString
   * @return
   */
  protected String changeField(final String fieldName, 
                               final String newFieldStr,
                               String newQualifierString)
  {
    int ind1 = newQualifierString.toLowerCase().indexOf(fieldName.toLowerCase());
    int ind2 = newQualifierString.indexOf(";", ind1);
    
    int len = fieldName.length();

    if(ind2 > ind1 && ind1 > -1)
    {
      if(newFieldStr.equals(""))
        newQualifierString =
          newQualifierString.substring(0, ind1) +
          newQualifierString.substring(ind2);
      else
        newQualifierString =
          newQualifierString.substring(0, ind1+len) +
          newFieldStr +
          newQualifierString.substring(ind2);
    }
    else if(ind1 > -1)
    {
      if(newFieldStr.equals(""))
        newQualifierString =
          newQualifierString.substring(0, ind1);
      else
        newQualifierString =
          newQualifierString.substring(0, ind1+len) +
          newFieldStr;
    }
    else
    {
      if(!newFieldStr.equals(""))
        newQualifierString = newQualifierString + ";" + 
                             fieldName + newFieldStr;
    }
    
    return newQualifierString;
  }
  
/*  protected String replace(String str, String pattern, String replace) 
  { 
    int s = 0; 
    int e = 0; 
    StringBuffer result = new StringBuffer(); 
    while ((e = str.indexOf(pattern, s)) >= 0) 
    { 
      result.append(str.substring(s, e)); 
      result.append(replace); 
      s = e+pattern.length(); 
    } 
    result.append(str.substring(s)); 
    return result.toString(); 
  } */
  
  /**
   * Get a Date object from a date string in the format 20061129
   * @param dateStr
   * @return
   */
  protected Date getDate(String dateStr)
  {
    Calendar cal = Calendar.getInstance();
    cal.clear();
    if(dateStr == null ||
       dateStr.equals("") ||
       dateStr.length() != 8)
      return null;
    
    int year  = Integer.parseInt(dateStr.substring(0,4));
    int month = Integer.parseInt(dateStr.substring(4, 6))-1;
    int day   = Integer.parseInt(dateStr.substring(6,8));
    
    cal.set(year,month,day);
    return cal.getTime();
  }
  
  protected abstract boolean isQualifierChanged();
  protected abstract void updateQualifier(final QualifierVector qv);
}