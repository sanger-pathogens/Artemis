package uk.ac.sanger.artemis.chado;

import java.io.Serializable;


public class FeatureCvTermProp implements Serializable 
{

  private int featureCvTermPropId;
  private CvTerm cvTerm;
  private FeatureCvTerm featureCvTerm;
  private String value;
  private int rank;

  // Property accessors
  public int getFeatureCvTermPropId()
  {
    return this.featureCvTermPropId;
  }

  public void setFeatureCvTermPropId(int featureCvTermPropId)
  {
    this.featureCvTermPropId = featureCvTermPropId;
  }

  public CvTerm getCvterm()
  {
    return this.cvTerm;
  }

  public void setCvterm(CvTerm cvTerm)
  {
    this.cvTerm = cvTerm;
  }

  public FeatureCvTerm getFeatureCvTerm()
  {
    return this.featureCvTerm;
  }

  public void setFeatureCvterm(FeatureCvTerm featureCvTerm)
  {
    this.featureCvTerm = featureCvTerm;
  }

  public String getValue()
  {
    return this.value;
  }

  public void setValue(String value)
  {
    this.value = value;
  }

  public int getRank()
  {
    return this.rank;
  }

  public void setRank(int rank)
  {
    this.rank = rank;
  }

}


