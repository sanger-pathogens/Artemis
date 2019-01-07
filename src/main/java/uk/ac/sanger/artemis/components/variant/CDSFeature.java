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
  protected Feature feature;
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
    this.feature   = feature;
  } 
  
  public String toString() {
      StringBuilder sb = new StringBuilder();
      final String sep = "\t";
      sb.append(isFwd);
      sb.append(sep);
      sb.append(firstBase);
      sb.append(sep);
      sb.append(lastBase);
      sb.append(sep);
      String b = bases;
      if (b.length() > 10) 
          b = b.substring(0, 10);
      sb.append(b);
      sb.append(sep);
      sb.append(ranges.size());
      sb.append(sep);
      for (int i = 0; i < ranges.size(); i++) {
          Range r = (Range) ranges.get(i);
          sb.append(r.getStart()+"-"+r.getEnd());
          sb.append(sep);
      }
      return sb.toString();
  }
}