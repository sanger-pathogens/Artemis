/*
 * Copyright (C) 2008  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 *  @author: Tim Carver
 */
 
package uk.ac.sanger.artemis.circular;

import java.awt.Color;
import java.io.IOException;
import java.io.Writer;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Vector;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureKeyPredicate;
import uk.ac.sanger.artemis.FeatureKeyQualifierPredicate;
import uk.ac.sanger.artemis.FeaturePredicate;
import uk.ac.sanger.artemis.io.Key;

public class Track
{
  private double position = .9d;
  private float size = 10.f;
  private FeaturePredicate featurePredicate;
  private boolean showForward = true;
  private boolean showReverse = true;
  private boolean notQualifier = false;
  private boolean any = false;
  private String keyStr;
  private String qualifier;
  private String qualifierValue;
  private Entry entry;
  private Color colour;
  
  public Track(double position, Entry entry)
  {
    this.position = position;
    this.entry = entry;
  }
  
  public Track(double position, 
               String keyStr,
               String qualifier,
               boolean notQualifier,
               boolean showForward,
               boolean showReverse,
               Entry entry)
  {
    this.position = position;
    this.showForward = showForward;
    this.showReverse = showReverse;
    this.keyStr = keyStr;
    this.qualifier = qualifier;
    this.notQualifier = notQualifier;
    this.entry = entry;
    
    if(keyStr != null)
    {
      final Key key = new Key(keyStr);
      if(qualifier != null)
        featurePredicate = 
          new FeatureKeyQualifierPredicate(key, qualifier, isNotQualifier());
      else
        featurePredicate = new FeatureKeyPredicate(key);
    }
    
    if(featurePredicate == null)
      any = true;
  }
  
  
  public Track(double position, 
               String keyStr,
               boolean showForward, 
               boolean showReverse,
               Entry entry)
  {
    this(position, keyStr, null, true, showForward, showReverse, entry);
  }
  
  
  public double getPosition()
  {
    return position;
  }

  public void setPosition(double position)
  {
    this.position = position;
  }

  /**
   * Test if this feature is drawn on this track
   * @param this_feature
   * @return
   */
  public boolean isOnTrack(Feature this_feature)
  {
    if(getEntry() != null && !getEntry().contains(this_feature))
      return false;
    
    if(isAny())
      return true;
    
    if(featurePredicate == null)
      return false;
    
    if(featurePredicate.testPredicate(this_feature))
    {
      if( (this_feature.isForwardFeature()  && isShowForward()) ||
          (!this_feature.isForwardFeature() && isShowReverse()) )
        return true;
    }
    return false;
  }

  public boolean isShowForward()
  {
    return showForward;
  }

  public void setShowForward(boolean showForward)
  {
    this.showForward = showForward;
  }

  public boolean isShowReverse()
  {
    return showReverse;
  }

  public void setShowReverse(boolean showReverse)
  {
    this.showReverse = showReverse;
  }

  public boolean isAny()
  {
    return any;
  }

  public void setAny(boolean any)
  {
    this.any = any;
  }

  public FeaturePredicate getFeaturePredicate()
  {
    return featurePredicate;
  }

  public void setFeaturePredicate(FeaturePredicate featurePredicate)
  {
    this.featurePredicate = featurePredicate;
  }

  public String getKeyStr()
  {
    return keyStr;
  }

  public void setKeyStr(String keyStr)
  {
    this.keyStr = keyStr;
  }

  public String getQualifier()
  {
    return qualifier;
  }

  public void setQualifier(String qualifier)
  {
    this.qualifier = qualifier;
  }

  public boolean isNotQualifier()
  {
    return notQualifier;
  }

  public void setNotQualifier(boolean notQualifier)
  {
    this.notQualifier = notQualifier;
  }

  public String getQualifierValue()
  {
    return qualifierValue;
  }

  public void setQualifierValue(String qualifierValue)
  {
    this.qualifierValue = qualifierValue;
  }

  public float getSize()
  {
    return size;
  }

  public void setSize(float size)
  {
    this.size = size;
  }

  public Entry getEntry()
  {
    return entry;
  }

  public void setEntry(Entry entry)
  {
    this.entry = entry;
  }

  public Color getColour()
  {
    return colour;
  }

  public void setColour(Color colour)
  {
    this.colour = colour;
  }
  
  protected void setPropertiesFromTemplate(final String line)
  {
    final String properties[] = line.split("\t");

    setPosition(Double.parseDouble(properties[0]));
    setSize(Float.parseFloat(properties[1]));
    setShowForward(Boolean.parseBoolean(properties[2]));
    setShowReverse(Boolean.parseBoolean(properties[3]));
    setNotQualifier(Boolean.parseBoolean(properties[4]));
    setAny(Boolean.parseBoolean(properties[5]));
    
    if(properties[6].equals("null"))
      setKeyStr(null);
    else
      setKeyStr(properties[6]);
    if(properties[7].equals("null"))
      setQualifier(null);
    else
      setQualifier(properties[7]);
    if(properties[8].equals("null"))
      setQualifierValue(null);
    else
      setQualifierValue(properties[8]);
    if(properties[9].equals("null"))
      setColour(null);
    else
    {
      String colourRGB[] = properties[9].split(":");
      Color colour = new Color(Integer.parseInt(colourRGB[0]),
                               Integer.parseInt(colourRGB[1]),
                               Integer.parseInt(colourRGB[2]));
      setColour(colour);
    }
  }
  
  /**
   * Write the track properties out
   * @param writer
   * @throws IOException
   */
  protected void write(final Writer writer) throws IOException
  { 
    writer.write(position+"\t"+
                 size+"\t"+
                 showForward+"\t"+
                 showReverse+"\t"+
                 isNotQualifier()+"\t"+
                 any+"\t"+
                 keyStr+"\t"+
                 qualifier+"\t"+
                 qualifierValue+"\t"+
                 ( colour != null ? 
                     colour.getRed()+":"+colour.getGreen()+":"+colour.getBlue() : null)+"\t"+
                 ( entry != null ? entry.getName()+"\t"+entry.getRootDocument() : null )+
                 "\n");
  }
  
  /**
   * Write a header line for the track properties
   * @param writer
   * @throws IOException
   */
  protected static void writeHeader(final Writer writer, final DNADraw dna) throws IOException
  {
    DateFormat dateFormat = new SimpleDateFormat("dd/MM/yyyy HH:mm:ss");      
    writer.write("## DNA Plot :: track template (created: "+
        dateFormat.format(new Date())+")\n");
    writer.write("# line attributes: start="+dna.getStart()+
                                     " end="+dna.getEnd()+" " +
                                "line_size="+dna.getLineSize()+
                                " circular="+dna.isCircular());
                                
    if(!dna.isCircular())
      writer.write(" line_height="+dna.getLineHeight()+
                   " bases_per_line="+dna.getBasesPerLine());
    
    writer.write("\n# tick marks: major="+dna.getTickInterval()+
                              " minor="+dna.getMinorTickInterval()+"\n");
    
    // graphs
    Vector<Graph> userGraphs = dna.getUserGraphs();
    for(int i=0; i<userGraphs.size(); i++)
    {
      Graph userGraph = userGraphs.get(i);
      if(dna.containsGraph(userGraph))
      {
        String fileName = ((UserGraph)userGraph).getFileName();  
        writer.write("# User Graph: "+userGraph.getOptionsStr()+
                     " file_name="+fileName+"\n");
      }
    }
  
    if(dna.getGcGraph() != null && dna.containsGraph(dna.getGcGraph()))
      writer.write("# GC Graph: "+dna.getGcGraph().getOptionsStr()+"\n");
    if(dna.getGcSkewGraph() != null && dna.containsGraph(dna.getGcSkewGraph()))
      writer.write("# GC Skew Graph: "+dna.getGcSkewGraph().getOptionsStr()+"\n");
    
    writer.write(
        "# Columns are:\n"+
        "# POS  - track position\n"+
        "# SIZE - track size\n"+
        "# FWD  - show forward strand features\n"+
        "# REV  - show reverse strand features\n"+
        "# NOT  - use NOT\n"+
        "# ANY  - show any features\n"+
        "# KEY  - show features of this key\n"+
        "# QUAL - show features with this qualifier\n"+
        "# VAL  - show features with this qualifier value(s)\n"+
        "# COL  - colour for this track e.g. 255:0:0 (R:G:B) or NULL\n"+
        "# NAME - file entry name or null\n"+
        "# DIR  - root directory for this file\n#\n");
    writer.write(
        "#POS\t"+
        "SIZE\t"+
        "FWD\t"+
        "REV\t"+
        "NOT \t"+
        "ANY\t"+
        "KEY\t"+
        "QUAL\t"+
        "VAL\t"+
        "COL\t"+
        "NAME\t"+
        "DIR\n");
  }
}