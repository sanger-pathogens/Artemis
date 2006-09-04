package uk.ac.sanger.artemis.chado;


import java.io.Serializable;
import java.util.HashSet;
import java.util.Set;


public class FeatureCvTerm implements Serializable 
{

    // Fields 
  private int featureCvTermId;
     
  private CvTerm cvTerm;
  private Feature feature;
  private Pub pub;
  private boolean not;
  private Set featureCvTermProps = new HashSet(0);
  private Set featureCvTermPubs  = new HashSet(0);
  private Set featureCvTermDbXRefs = new HashSet(0);

   
  // Property accessors
  public int getFeatureCvTermId()
  {
    return this.featureCvTermId;
  }

  public void setFeatureCvTermId(int featureCvTermId)
  {
    this.featureCvTermId = featureCvTermId;
  }

  public CvTerm getCvTerm()
  {
    return this.cvTerm;
  }

  public void setCvTerm(CvTerm cvTerm)
  {
    this.cvTerm = cvTerm;
  }

  public Feature getFeature()
  {
    return this.feature;
  }

  public void setFeature(Feature feature)
  {
    this.feature = feature;
  }

  public Pub getPub()
  {
    return this.pub;
  }

  public void setPub(Pub pub)
  {
    this.pub = pub;
  }

  public boolean isNot()
  {
    return this.not;
  }

  public void setNot(boolean not)
  {
    this.not = not;
  }

  public Set getFeatureCvTermProps()
  {
    return this.featureCvTermProps;
  }

  public void setFeatureCvtermprops(Set featureCvTermProps)
  {
    this.featureCvTermProps = featureCvTermProps;
  }

  public Set getFeatureCvTermPubs()
  {
    return this.featureCvTermPubs;
  }

  public void setFeatureCvTermPubs(Set featureCvTermPubs)
  {
    this.featureCvTermPubs = featureCvTermPubs;
  }

  public Set getFeatureCvTermDbXRefs()
  {
    return this.featureCvTermDbXRefs;
  }

  public void setFeatureCvTermDbxrefs(Set featureCvTermDbXRefs)
  {
    this.featureCvTermDbXRefs = featureCvTermDbXRefs;
  }

}
