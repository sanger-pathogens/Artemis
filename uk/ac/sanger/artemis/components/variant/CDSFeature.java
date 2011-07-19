package uk.ac.sanger.artemis.components.variant;

import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.RangeVector;

class CDSFeature
{
  protected boolean isFwd;
  protected RangeVector ranges;
  protected int firstBase;
  protected int lastBase;
  protected int intronlength = 0;
  protected String bases;
  protected String id;
  protected Range lastRange = null;

  public CDSFeature(boolean isFwd,
		  	    	RangeVector ranges,
		  			int firstBase,
		  			int lastBase,
		  			String bases)
  {
    this.isFwd = isFwd;
    this.ranges = ranges;
    this.firstBase = firstBase;
    this.lastBase = lastBase;
    this.bases = bases;
 }
  
  public CDSFeature(Feature feature)
  {
    this.isFwd     = feature.isForwardFeature();
    this.ranges    = feature.getLocation().getRanges();
    this.firstBase = feature.getRawFirstBase();
    this.lastBase  = feature.getRawLastBase();
    this.bases     = feature.getBases();
    this.id        = feature.getIDString();
  } 
}