package uk.ac.sanger.artemis.io;

public interface LazyQualifierValue
{
  public abstract String getString();
  public abstract void setForceLoad(boolean lazyLoad);
  public abstract boolean isLazyLoaded();
}
