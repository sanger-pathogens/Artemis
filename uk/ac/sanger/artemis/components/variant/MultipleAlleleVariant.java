
package uk.ac.sanger.artemis.components.variant;

public class MultipleAlleleVariant
{

  private static String IUB[][] = 
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
  
  /**
   * M (A or C)
   * R (A or G)
   * W (A or T/U)
   * S (C or G)
   * Y (C or T/U)
   * K (G or T/U)
   * @param base
   * @return
   */
  protected static String getIUBCode(VCFRecord record)
  {
    String alt = record.getAlt().toString();
    String alleles[] = alt.split(",");
    String pl;
    if ((pl = record.getFormatValue("PL")) != null && pl.split(",").length == 3 &&
         pl.split(",")[1].equals("0")) 
    {
      // include ref
      String[] temp = new String[alleles.length+1];
      System.arraycopy(alleles, 0, temp, 0, alleles.length);
      alleles = temp;
      alleles[alleles.length-1] = record.getRef();
    }
    
    boolean isSNP = true;
    for(int i=0; i<alleles.length; i++)
    {
      alleles[i] = alleles[i].toLowerCase();
      if(alleles[i].length() > 1)
        isSNP = false;
    }
    
    if(isSNP && alleles.length == 2)
    {
      for(int i=0; i<IUB.length; i++)
      {
        if(IUB[i][1].equals(alleles[0]) && IUB[i][2].equals(alleles[1]) ||
           IUB[i][1].equals(alleles[1]) && IUB[i][2].equals(alleles[0]))
          return IUB[i][0];
      }
    }
    return null;
  }

}