
package uk.ac.sanger.artemis.components.variant;

import java.util.regex.Pattern;

public class VariantBase
{
  private VCFRecord record;
  private String alt;
  protected static Pattern COMMA_PATTERN = Pattern.compile(",");

  public VariantBase(VCFRecord record, String alt)
  {
    this.record = record;
    this.alt = alt;
  }

  public String toString()
  {
    return alt;
  }

  protected int length()
  {
    return alt.length();
  }

  /**
   * Is this a deletion type.
   * @param variant
   * @return
   */
  protected boolean isDeletion(boolean vcf_v4)
  {
    if (vcf_v4)
    {
      if (alt.length() < record.getRef().length() && !(alt.indexOf(",") > -1))
        return true;
    }
    else if (alt.indexOf("D") > -1)
      return true;
    return false;
  }

  /**
   * Is this an insertion type.
   * 
   * @param variant
   * @return
   */
  protected boolean isInsertion(boolean vcf_v4)
  {
    if (vcf_v4)
    {
      if (alt.length() > record.getRef().length() && !(alt.indexOf(",") > -1))
        return true;
    }
    else if (alt.indexOf("I") > -1)
      return true;
    return false;
  }

  protected boolean isMultiAllele(int sampleIndex)
  {
    if (VCFRecord.MULTI_ALLELE_PATTERN.matcher(alt).matches())
      return true;

    // look at probability of each genotype (PL) information as well
    if(sampleIndex < 0)
      sampleIndex = 0;
    final String pl = record.getFormatValueForSample("PL", sampleIndex);
    if(pl != null)
    {
      final String pls[] = COMMA_PATTERN.split(pl);
      if(pls.length == 3 && pls[1].equals("0")) // middle value is zero, e.g.
        return true;
    }

    return false;
  }

  protected int getNumAlleles()
  {
    return VCFRecord.countOccurrences(alt, ',')+2;
    //return COMMA_PATTERN.split(alt).length + 1;
    //return alt.split(",").length + 1;
  }

  protected int getNumberOfIndels(boolean vcf_v4)
  {
    if (vcf_v4)
    {
      if (alt.equals("."))
        return record.getRef().length();
      return Math.abs(record.getRef().length() - alt.length());
    }

    int index = alt.indexOf("D");
    if(index < 0)
      index = alt.indexOf("I");
    int ndel = 0;
    try
    {
      ndel = Integer.parseInt(alt.substring(index + 1));
    }
    catch (NumberFormatException e)
    {
      System.out.println(alt);
      e.printStackTrace();
    }
    return ndel;
  }

  protected boolean isNonVariant()
  {
    if ((alt.equals(".") ||  alt.equals("X")) && record.getRef().length() == 1)
      return true;
    return false;
  }
}