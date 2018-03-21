/* MultipleAlleleVariant
 *
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

public class MultipleAlleleVariant
{
  private static String IUB_2[][] = 
  { 
      {"m", "a", "c"},
      {"r", "a", "g"},
      {"w", "a", "t"},
      {"w", "a", "u"},
      {"s", "c", "g"},
      {"y", "c", "t"},
      {"y", "c", "u"},
      {"k", "g", "t"},
      {"k", "g", "u"}
  };

  private static String IUB_3[][] = 
  { 
      {"b", "c", "g", "t"},
      {"b", "c", "g", "u"},
      {"d", "a", "g", "t"},
      {"d", "a", "g", "u"},
      {"h", "a", "c", "t"},
      {"h", "a", "c", "u"},
      {"v", "a", "c", "g"}
  };
  
  /**
   * M (A or C)             B (C, G or T/U)
   * R (A or G)             D (A, G or T/U)
   * W (A or T/U)           H (A, C or T/U)
   * S (C or G)             V (A, C or G)
   * Y (C or T/U)
   * K (G or T/U)
   * @param base
   * @return
   */
  protected static String getIUBCode(VCFRecord record)
  {
    String alt = record.getAlt().toString();
    String alleles[] = alt.toLowerCase().split(",");
    String pl;
    if ((pl = record.getFormatValueForSample("PL", 0)) != null && pl.split(",").length == 3 &&
         pl.split(",")[1].equals("0")) 
    {
      // include ref
      String[] temp = new String[alleles.length+1];
      System.arraycopy(alleles, 0, temp, 0, alleles.length);
      alleles = temp;
      alleles[alleles.length-1] = record.getRef().toLowerCase();
    }

    boolean isSNP = true;
    for(int i=0; i<alleles.length; i++)
    {
      alleles[i] = alleles[i];
      if(alleles[i].length() > 1)
        isSNP = false;
    }

    if(isSNP && alleles.length == 2)
    {
      for(int i=0; i<IUB_2.length; i++)
      {
        if( (IUB_2[i][1].equals(alleles[0]) && IUB_2[i][2].equals(alleles[1])) ||
            (IUB_2[i][1].equals(alleles[1]) && IUB_2[i][2].equals(alleles[0])) )
          return IUB_2[i][0];
      }
    }

    if(isSNP && alleles.length == 3)
    {
      //System.out.println(record.toString());
      for(int i=0; i<IUB_3.length; i++)
      {
        if( (IUB_3[i][1].equals(alleles[0]) && IUB_3[i][2].equals(alleles[1]) && IUB_3[i][3].equals(alleles[2])) ||
            (IUB_3[i][1].equals(alleles[1]) && IUB_3[i][2].equals(alleles[0]) && IUB_3[i][3].equals(alleles[2])) ||
            (IUB_3[i][1].equals(alleles[2]) && IUB_3[i][2].equals(alleles[0]) && IUB_3[i][3].equals(alleles[1])) ||
            (IUB_3[i][1].equals(alleles[2]) && IUB_3[i][2].equals(alleles[1]) && IUB_3[i][3].equals(alleles[0])) ||
            (IUB_3[i][1].equals(alleles[0]) && IUB_3[i][2].equals(alleles[2]) && IUB_3[i][3].equals(alleles[1])) ||
            (IUB_3[i][1].equals(alleles[1]) && IUB_3[i][2].equals(alleles[2]) && IUB_3[i][3].equals(alleles[0])))
          return IUB_3[i][0];
      }
    }
    
    if(isSNP)
      return "n";

    return null;
  }

}