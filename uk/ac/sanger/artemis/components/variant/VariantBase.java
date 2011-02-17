
package uk.ac.sanger.artemis.components.variant;

public class VariantBase
{
  private VCFRecord record;
  private String alt;

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

  protected boolean isMultiAllele()
  {
    if (VCFRecord.MULTI_ALLELE_PATTERN.matcher(alt).matches())
      return true;

    // look at probability of each genotype (PL) information as well
    String pl;
    if ((pl = record.getFormatValue("PL")) != null && pl.split(",").length == 3 &&
         pl.split(",")[1].equals("0")) // middle value is zero, e.g.
      return true;

    return false;
  }

  protected int getNumAlleles()
  {
    if (alt.equals("."))
      return 1;

    return alt.split(",").length + 1;
  }

  protected int getNumberOfDeletions(boolean vcf_v4)
  {
    if (vcf_v4)
    {
      if (alt.equals("."))
        return record.getRef().length();
      return record.getRef().length() - alt.length();
    }

    int index = alt.indexOf("D");
    int ndel = 0;
    try
    {
      ndel = Integer.parseInt(alt.substring(index + 1));
    }
    catch (NumberFormatException e)
    {
      e.printStackTrace();
    }
    return ndel;
  }

  protected boolean isNonVariant()
  {
    if (alt.equals(".") && record.getRef().length() == 1)
      return true;
    return false;
  }
}