package uk.ac.sanger.artemis.chado;

import java.io.Serializable;


public class FeatureCvTermDbXRef implements Serializable
{ 
  private int featureCvTermDbXRefId;
  private DbXRef dbXRef;
  private FeatureCvTerm featureCvTerm;

  public int getFeatureCvTermDbXRefId()
  {
    return this.featureCvTermDbXRefId;
  }

  public void setFeatureCvTermDbXRefId(int featureCvTermDbXRefId)
  {
    this.featureCvTermDbXRefId = featureCvTermDbXRefId;
  }

  public DbXRef getDbXRef()
  {
    return this.dbXRef;
  }

  public void setDbXRef(DbXRef dbXRef)
  {
    this.dbXRef = dbXRef;
  }

  public FeatureCvTerm getFeatureCvTerm()
  {
    return this.featureCvTerm;
  }

  public void setFeatureCvTerm(FeatureCvTerm featureCvTerm)
  {
    this.featureCvTerm = featureCvTerm;
  }

}


