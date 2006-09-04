package uk.ac.sanger.artemis.chado;


import java.io.Serializable;
import java.util.HashSet;
import java.util.Set;


public class Pub implements Serializable 
{  
  private Set phylotreePubs = new HashSet(0);
  private Set phylonodePubs = new HashSet(0);
  private int pubId;
  private CvTerm cvTerm;
  private String title;
  private String volumeTitle;
  private String volume;
  private String seriesName;
  private String issue;
  private String pyear;
  private String pages;
  private String miniRef;
  private String uniqueName;
  private Boolean obsolete;
  private String publisher;
  private String pubPlace;
  private Set pubAuthors = new HashSet(0);
  private Set pubRelationshipsForObjectId = new HashSet(0);
  private Set pubDbXRefs = new HashSet(0);
  private Set featureCvTerms = new HashSet(0);
  private Set featureRelationshipPubs = new HashSet(0);
  private Set featurePubs = new HashSet(0);
  private Set featurePropPubs = new HashSet(0);
  private Set featureSynonyms = new HashSet(0);
  private Set featureCvTermPubs = new HashSet(0);
  private Set featureRelationshipPropPubs = new HashSet(0);
  private Set pubProps = new HashSet(0);
  private Set pubRelationshipsForSubjectId = new HashSet(0);
  private Set featureLocPubs = new HashSet(0);
  
  // Property accessors
  public int getPubId()
  {
    return this.pubId;
  }

  public void setPubId(int pubId)
  {
    this.pubId = pubId;
  }

  public CvTerm getCvTerm()
  {
    return this.cvTerm;
  }

  public void setCvTerm(CvTerm cvTerm)
  {
    this.cvTerm = cvTerm;
  }

  public String getTitle()
  {
    return this.title;
  }

  public void setTitle(String title)
  {
    this.title = title;
  }

  public String getVolumeTitle()
  {
    return this.volumeTitle;
  }

  public void setVolumeTitle(String volumeTitle)
  {
    this.volumeTitle = volumeTitle;
  }

  public String getVolume()
  {
    return this.volume;
  }

  public void setVolume(String volume)
  {
    this.volume = volume;
  }

  public String getSeriesName()
  {
    return this.seriesName;
  }

  public void setSeriesName(String seriesName)
  {
    this.seriesName = seriesName;
  }

  public String getIssue()
  {
    return this.issue;
  }

  public void setIssue(String issue)
  {
    this.issue = issue;
  }

  public String getPyear()
  {
    return this.pyear;
  }

  public void setPyear(String pyear)
  {
    this.pyear = pyear;
  }

  public String getPages()
  {
    return this.pages;
  }

  public void setPages(String pages)
  {
    this.pages = pages;
  }

  public String getMiniRef()
  {
    return this.miniRef;
  }

  public void setMiniRef(String miniRef)
  {
    this.miniRef = miniRef;
  }

  public String getUniqueName()
  {
    return this.uniqueName;
  }

  public void setUniqueName(String uniqueName)
  {
    this.uniqueName = uniqueName;
  }

  public Boolean getObsolete()
  {
    return this.obsolete;
  }

  public void setObsolete(Boolean obsolete)
  {
    this.obsolete = obsolete;
  }

  public String getPublisher()
  {
    return this.publisher;
  }

  public void setPublisher(String publisher)
  {
    this.publisher = publisher;
  }

  public String getPubPlace()
  {
    return this.pubPlace;
  }

  public void setPubPlace(String pubPlace)
  {
    this.pubPlace = pubPlace;
  }

  public Set getPubAuthors()
  {
    return this.pubAuthors;
  }

  public void setPubAuthors(Set pubAuthors)
  {
    this.pubAuthors = pubAuthors;
  }

  public Set getPubRelationshipsForObjectId()
  {
    return this.pubRelationshipsForObjectId;
  }

  public void setPubRelationshipsForObjectId(Set pubRelationshipsForObjectId)
  {
    this.pubRelationshipsForObjectId = pubRelationshipsForObjectId;
  }

  public Set getPubDbXRefs()
  {
    return this.pubDbXRefs;
  }

  public void setPubDbXRefs(Set pubDbXRefs)
  {
    this.pubDbXRefs = pubDbXRefs;
  }

  public Set getFeatureCvTerms()
  {
    return this.featureCvTerms;
  }

  public void setFeatureCvTerms(Set featureCvTerms)
  {
    this.featureCvTerms = featureCvTerms;
  }

  public Set getFeatureRelationshipPubs()
  {
    return this.featureRelationshipPubs;
  }

  public void setFeatureRelationshipPubs(Set featureRelationshipPubs)
  {
    this.featureRelationshipPubs = featureRelationshipPubs;
  }

  public Set getFeaturePubs()
  {
    return this.featurePubs;
  }

  public void setFeaturePubs(Set featurePubs)
  {
    this.featurePubs = featurePubs;
  }

  public Set getFeaturePropPubs()
  {
    return this.featurePropPubs;
  }

  public void setFeaturePropPubs(Set featurePropPubs)
  {
    this.featurePropPubs = featurePropPubs;
  }

  public Set getFeatureSynonyms()
  {
    return this.featureSynonyms;
  }

  public void setFeatureSynonyms(Set featureSynonyms)
  {
    this.featureSynonyms = featureSynonyms;
  }

  public Set getFeatureCvTermPubs()
  {
    return this.featureCvTermPubs;
  }

  public void setFeatureCvTermPubs(Set featureCvTermPubs)
  {
    this.featureCvTermPubs = featureCvTermPubs;
  }

  public Set getFeatureRelationshipPropPubs()
  {
    return this.featureRelationshipPropPubs;
  }

  public void setFeatureRelationshipPropPubs(Set featureRelationshipPropPubs)
  {
    this.featureRelationshipPropPubs = featureRelationshipPropPubs;
  }

  public Set getPubProps()
  {
    return this.pubProps;
  }

  public void setPubProps(Set pubProps)
  {
    this.pubProps = pubProps;
  }

  public Set getPubRelationshipsForSubjectId()
  {
    return this.pubRelationshipsForSubjectId;
  }

  public void setPubRelationshipsForSubjectId(Set pubRelationshipsForSubjectId)
  {
    this.pubRelationshipsForSubjectId = pubRelationshipsForSubjectId;
  }

  public Set getFeatureLocPubs()
  {
    return this.featureLocPubs;
  }

  public void setFeatureLocPubs(Set featureLocPubs)
  {
    this.featureLocPubs = featureLocPubs;
  }

  public Set getPhylotreePubs()
  {
    return this.phylotreePubs;
  }

  public void setPhylotreePubs(Set phylotreePubs)
  {
    this.phylotreePubs = phylotreePubs;
  }

  public Set getPhylonodePubs()
  {
    return this.phylonodePubs;
  }

  public void setPhylonodePubs(Set phylonodePubs)
  {
    this.phylonodePubs = phylonodePubs;
  }

}
