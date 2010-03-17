/* UserDataAlgorithm.java
 *
 * created: Wed May 10 2000
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2000  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/plot/UserDataAlgorithm.java,v 1.14 2009-07-22 12:51:54 tjc Exp $
 */

package uk.ac.sanger.artemis.plot;

import uk.ac.sanger.artemis.sequence.*;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.ReadFormatException;

import java.awt.Color;
import java.awt.GridLayout;
import java.io.*;
import java.util.HashMap;
import java.util.regex.Pattern;

import javax.swing.ButtonGroup;
import javax.swing.JCheckBox;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

/**
 *  Objects of this class have one useful method - getValues (), which takes a
 *  range of bases and returns a single floating point number.  The number is
 *  calculated by averaging the values from a data file.  The Strand to use is
 *  set in the constructor.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: UserDataAlgorithm.java,v 1.14 2009-07-22 12:51:54 tjc Exp $
 **/

public class UserDataAlgorithm extends BaseAlgorithm
{
  /** A base per line file format */
  public static int BASE_PER_LINE_FORMAT  = 1;
  
  /** Base position is specified in the first column file format */
  public static int BASE_SPECIFIED_FORMAT = 2;
  
  /** Wiggle format */
  public static int WIGGLE_VARIABLE_STEP_FORMAT = 3;
  public static int WIGGLE_FIXED_STEP_FORMAT = 4;
  
  public static int BLAST_FORMAT = 5;
  public static int MSPCRUNCH_BLAST_FORMAT = 6;
  
  /** The data read by the constructor - for BASE_PER_LINE_FORMAT */
  private float data[][] = null;
  
  /** The data read by the constructor - for BASE_SPECIFIED_FORMAT 
                                       - and WIGGLE_VARIABLE_STEP_FORMAT
                                       - and WIGGLE_FIXED_STEP_FORMAT */
  private HashMap<Integer, Float[]> dataMap;
  
  /** The maximum value in the data array. */
  private float data_max = Float.MIN_VALUE;

  /** The minimum value in the data array. */
  private float data_min = Float.MAX_VALUE;
  
  /** Default window size */
  private int default_window_size = 3;

  /** The average calculated by readData (). */
  private float average_value = 0;

  /** The value returned by getValueCount (). */
  private int number_of_values;
  
  private boolean logTransform;
  
  /** Format type for this instance */
  public int FORMAT = BASE_PER_LINE_FORMAT;
  
  private LineAttributes lines[];
  
  public Wiggle wiggle[];
  
  /**
   *  Create a new UserDataAlgorithm object. This reads a file
   *  which can be one of two types of formats:
   *  a. one line of values per base.
   *  b. the first column specifies the base position with
   *     subsequent columns being values.
   *  @param strand The strand to do the calculation on.
   *  @param document The Document to read the data from.
   *  @param logTransform true if the log transformation is to be
   *  shown.
   **/
  public UserDataAlgorithm (final Strand strand, final Document document, 
                            final boolean logTransform)
      throws IOException 
  {
    super (strand, "User algorithm from " + document.getName (), "user");

    this.logTransform = logTransform;
    final Reader document_reader = document.getReader ();

    LinePushBackReader pushback_reader = new LinePushBackReader (document_reader);
    String first_line = pushback_reader.readLine (); 

    Pattern dataPattern = Pattern.compile("^\\s*([\\d\\.-]+\\s*)+$");
    Pattern blastPattern = Pattern.compile(
      "^(\\S+\\t+){2}[\\d\\.]+\\t+(\\d+\\t+){7}\\S+\\t+(\\s*\\d+)$");
    Pattern mspCrunchPattern = Pattern.compile(
        "^\\d+\\s[\\d\\.]+(\\s\\d+){2}\\s\\D\\S+(\\s\\d+){2}\\s\\D\\S+.*");

    if(dataPattern.matcher(first_line).matches())
      FORMAT = BASE_PER_LINE_FORMAT;
    else if(blastPattern.matcher(first_line).matches())
      FORMAT = BLAST_FORMAT;
    else if(mspCrunchPattern.matcher(first_line).matches())
      FORMAT = MSPCRUNCH_BLAST_FORMAT;
    else
    { 
      StringBuffer header = new StringBuffer(first_line+"\n");

      while(!dataPattern.matcher(first_line).matches())
      {
        first_line = pushback_reader.readLine ().trim();
        header.append(first_line+"\n");
      }
      
      FORMAT = parseHeader(header);
    }

    final Pattern patt = Pattern.compile("\\s+");
    String tokens[] = patt.split(first_line);
    
    if (tokens.length < 1) 
      throw new ReadFormatException ("unknown file type");

    this.number_of_values = tokens.length;
    pushback_reader.pushBack (first_line);
    
    if(FORMAT == BASE_PER_LINE_FORMAT)
      data = new float [strand.getSequenceLength ()][tokens.length];
    
    if(FORMAT == BASE_SPECIFIED_FORMAT ||
       FORMAT == BASE_PER_LINE_FORMAT)
      readData(pushback_reader);
    else if(FORMAT == BLAST_FORMAT ||
            FORMAT == MSPCRUNCH_BLAST_FORMAT)
      readBlast(pushback_reader);
    else
      readWiggle(pushback_reader);
    pushback_reader.close();
    document_reader.close();
  }

  /**
   *  Read all from buffered_reader into data.
   **/
  private void readData (final LinePushBackReader pushback_reader)
      throws IOException
  {
    String line = null;
    int count = 0;
    int countAll = 0;
    int estimate_window_size = Integer.MAX_VALUE;
    final int seqLength = getStrand ().getSequenceLength ();
    final Pattern patt = Pattern.compile("\\s+");
    
    while ((line = pushback_reader.readLine ()) != null)
    {
      if (count >= seqLength) 
        throw new ReadFormatException ("too many values in input file");

      String tokens[] = patt.split(line); 
      if (FORMAT == BASE_PER_LINE_FORMAT && tokens.length != data[0].length)
        throw new ReadFormatException ("line has the wrong number of fields:\n"+line);
      
      int base = 0;
      Float line_data[] = new Float[tokens.length-1];
      for (int i = 0 ; i < tokens.length ; ++i)
      {
        try 
        {
          if( FORMAT == BASE_SPECIFIED_FORMAT &&
              i == 0)
          {
            int last_base = base;
            base = (int) Float.parseFloat(tokens[i]);
            
            if(base > seqLength)
              throw new ReadFormatException (
                  "a base position is greater than the sequence length:\n"+line);

            if((base - last_base) < estimate_window_size &&
               (base - last_base) > 0)
              estimate_window_size = base - last_base;
            if(dataMap == null)
              dataMap = new HashMap<Integer, Float[]>();
            continue;
          }
            
          float value = Float.parseFloat(tokens[i]);
          if(logTransform)
            value = (float) Math.log(value+1);

          if (value > data_max) 
            data_max = value;
          if (value < data_min)
            data_min = value;
            
          if(FORMAT == BASE_PER_LINE_FORMAT)
            data[count][i] = value;
          else
          {
            line_data[i-1] = value;
            countAll++;
          }
          average_value += value;
        } 
        catch (NumberFormatException e) 
        {
          throw new ReadFormatException ("cannot understand this number: " +
                                         tokens[i] + " - " +e.getMessage ());
        }
      }
      ++count;
      if(FORMAT == BASE_SPECIFIED_FORMAT)
        dataMap.put(base, line_data);
    }

    if (FORMAT == BASE_PER_LINE_FORMAT)
      average_value /= data[0].length * seqLength;
    else
    {
      average_value = average_value/countAll;
      if(estimate_window_size != Integer.MAX_VALUE)
        default_window_size = estimate_window_size;
    }
  }
  
  /**
   *  Read all from buffered_reader into data.
   **/
  private void readWiggle (final LinePushBackReader pushback_reader)
      throws IOException
  {
    String line = null;
    int count = 0;
    int stepCount = 0;
    final int seqLength = getStrand ().getSequenceLength ();
    final Pattern patt = Pattern.compile("\\s+");

    dataMap = new HashMap<Integer, Float[]>();
    this.number_of_values = 1;
    
    while ((line = pushback_reader.readLine ()) != null)
    {
      if(line.startsWith("track"))
      {
        parseTrackLine(line);
        continue;
      }
      else if(line.startsWith("variableStep ") ||
              line.startsWith("fixedStep"))
      {
        FORMAT = parseWiggle(line);
        stepCount = 0;
        this.number_of_values++;
        continue;
      }
      else if(line.startsWith("#"))
        continue;
      
      String tokens[] = patt.split(line.trim()); 
      int base = 0;

      int valueIndex = 0;
      try 
      {
        if(FORMAT == WIGGLE_VARIABLE_STEP_FORMAT)
        {
          base = (int) Float.parseFloat(tokens[0]);
          valueIndex = 1;
        }
        else
        {
          base  = wiggle[wiggle.length-1].start + 
           (stepCount*wiggle[wiggle.length-1].step);
        }
        
        if(base > seqLength)
          throw new ReadFormatException (
                "the base position ("+base+") is greater than the sequence length:\n"+line);
    
        float value = Float.parseFloat(tokens[valueIndex]);
        if(logTransform)
          value = (float) Math.log(value+1);

        if (value > data_max) 
          data_max = value;
        if (value < data_min)
          data_min = value;
        
        final Float valueArray[] = new Float[number_of_values];;
        if(dataMap.containsKey(base))
        {
          Float oldValues[] = dataMap.get(base);
          for(int i=0; i<oldValues.length; i++)
            valueArray[i] = oldValues[i];
        }

        valueArray[number_of_values-1] = value;
        dataMap.put(base, valueArray);
  
        count++;
        stepCount++;
        average_value += value;
      } 
      catch (NumberFormatException e) 
      {
        throw new ReadFormatException ("cannot understand this number: " +
                                       tokens[valueIndex] + " - " +e.getMessage ());
      } 
    }

    average_value = average_value/count;
    default_window_size = 1;
  }

  
  /**
   *  Read all from buffered_reader into data.
   **/
  private void readBlast (final LinePushBackReader pushback_reader)
      throws IOException
  {
    String line = null;
    int count = 0;
    int lineNum = 0;
    final int seqLength = getStrand ().getSequenceLength ();
    final Pattern patt;
    
    if(FORMAT == BLAST_FORMAT)
      patt = Pattern.compile("\\t+");
    else
      patt = Pattern.compile("\\s");
    
    dataMap = new HashMap<Integer, Float[]>();
    this.number_of_values = 1;
    
    int coordIndexStart = 6;
    int coordIndexEnd   = 7;
    
    if(FORMAT == MSPCRUNCH_BLAST_FORMAT)
    {
      coordIndexStart = 2;
      coordIndexEnd   = 3;
    }
    
    
    while ((line = pushback_reader.readLine ()) != null)
    { 
      String tokens[] = patt.split(line.trim());
      
      if(lineNum == 0)
      {
        final JPanel message = new JPanel(new GridLayout(2,1));

        String queryStr = (FORMAT == BLAST_FORMAT) ? tokens[0] : tokens[4];
        String subjStr  = (FORMAT == BLAST_FORMAT) ? tokens[1] : tokens[7];
        
        JCheckBox query = new JCheckBox(queryStr, true);
        message.add(query);
        JCheckBox subj  = new JCheckBox(subjStr, false);
        message.add(subj);
        ButtonGroup group = new ButtonGroup();
        group.add(query);
        group.add(subj);
        
        JOptionPane.showConfirmDialog(null, message, 
            "Use Coordinates From", JOptionPane.OK_OPTION);
        if(subj.isSelected())
        {
          if(FORMAT == BLAST_FORMAT)
          {
            coordIndexStart = 8;
            coordIndexEnd = 9;
          }
          else
          {
            coordIndexStart = 5;
            coordIndexEnd = 6;
          }
        }
      }
      lineNum++;
      
      int startBase = Integer.parseInt(tokens[coordIndexStart]);
      int endBase   = Integer.parseInt(tokens[coordIndexEnd]);
      float value;
      if(FORMAT == BLAST_FORMAT)
        value = Float.parseFloat(tokens[11]);
      else
        value = Float.parseFloat(tokens[0]);
      
      int valueIndex = 0;
      try 
      {
        if(startBase > seqLength || endBase > seqLength)
          throw new ReadFormatException (
                "the location ("+startBase+".."+endBase+
                ") is outside than the sequence length:\n"+line);

        if(logTransform)
          value = (float) Math.log(value+1);

        if (value > data_max) 
          data_max = value;
        if (value < data_min)
          data_min = value;
        
        if(startBase > endBase)
        {
          int tmpStart = startBase;
          startBase = endBase;
          endBase   = tmpStart;
        }
            
        final Float valueArray[] = new Float[number_of_values];
        for (int base = startBase; base <= endBase; base++)
        {
          if (dataMap.containsKey(base))
          {
            Float oldValues[] = dataMap.get(base);
            if(oldValues[0] < value)
            {
              valueArray[number_of_values - 1] = value;
              dataMap.put(base, valueArray);
            }
          }
          else
          {
            valueArray[number_of_values - 1] = value;
            dataMap.put(base, valueArray);
          }
          count++;
        }

        average_value += value;
      } 
      catch (NumberFormatException e) 
      {
        throw new ReadFormatException ("cannot understand this number: " +
                                       tokens[valueIndex] + " - " +e.getMessage ());
      } 
    }

    average_value = average_value/count;
    default_window_size = 1;
    
    FORMAT = BLAST_FORMAT;
  }

  /**
   *  Return the value of the function between a pair of bases.
   *  @param start The start base (included in the range).
   *  @param end The end base (included in the range).
   *  @param values The one return value for this algorithm is returned in
   *    this array.
   **/
  public void getValues (int start, int end, final float [] values) 
  {
    final int value_count = getValueCount ();

    if(getStrand ().getDirection() == Bases.REVERSE)
    {
      int tstart = start;
      int tend   = end;
      end   = getStrand().getBases().getComplementPosition(tstart);
      start = getStrand().getBases().getComplementPosition(tend);
    }
    
    if(FORMAT == BASE_SPECIFIED_FORMAT ||
       FORMAT == WIGGLE_VARIABLE_STEP_FORMAT ||
       FORMAT == WIGGLE_FIXED_STEP_FORMAT ||
       FORMAT == BLAST_FORMAT)
    {
      for (int i = 0 ; i < value_count ; ++i) 
      {
        values[i] = 0;
        int count = 0;
        for (int base = start ; base <= end ; ++base) 
        {
          if(dataMap.containsKey(base) && 
             ((Float[])dataMap.get(base)).length > i &&
             ((Float[])dataMap.get(base))[i] != null)
          {
            values[i] += ((Float[])dataMap.get(base))[i];
            count++;
          }
        }

        if(count > 1)
          values[i] = values[i]/count;
      }
    }
    else
    {
      for (int i = 0 ; i < value_count ; ++i) 
      {
        values [i] = 0;
        for (int base = start ; base <= end ; ++base) 
          values [i] += data[base - 1][i] / (end - start + 1);
      }
    }
  }

  /**
   * Determine the graph file format.
   * Read the line colour from the header, there should be
   * one per line and space separated.
   * @param line
   * @throws ReadFormatException 
   */
  private int parseHeader(StringBuffer headerText) throws ReadFormatException
  {
    FORMAT = BASE_SPECIFIED_FORMAT;
    
    BufferedReader reader = new BufferedReader(
        new StringReader(headerText.toString()));
      
    String line = null;
    try
    {
      while ((line = reader.readLine()) != null) 
      {
        if((line.indexOf("colour") == -1 && 
            line.indexOf("color")  == -1 ) ||
            line.startsWith("track"))
        {
          if(line.startsWith("track "))
            parseTrackLine(line);
          else  if(line.startsWith("variableStep ") ||
                   line.startsWith("fixedStep"))
            FORMAT = parseWiggle(line);

          continue;
        }

        int index = line.indexOf("colour");
        if (index == -1)
          index = line.indexOf("color");

        index = line.indexOf(" ", index + 1);
        line = line.substring(index).trim();
        String rgbValues[] = line.split(" ");

        lines = new LineAttributes[rgbValues.length];
        for (int j = 0; j < rgbValues.length; j++)
          lines[j] = new LineAttributes(LineAttributes.parse(rgbValues[j]));
      }
    }
    catch (NumberFormatException e) 
    {
      throw new ReadFormatException ("cannot understand this number: " +
                                     line + " - " +e.getMessage ());
    } 
    catch (Exception e)
    {
      e.printStackTrace();
    }
    return FORMAT;
  }
  
  /**
   * http://genome.ucsc.edu/goldenPath/help/hgWiggleTrackHelp.html
   * 
   * All options are placed in a single line separated by spaces:
   *
   *  track type=wiggle_0 name=track_label description=center_label
   *        visibility=display_mode color=r,g,b altColor=r,g,b
   *        priority=priority autoScale=on|off
   *        gridDefault=on|off maxHeightPixels=max:default:min
   *        graphType=bar|points viewLimits=lower:upper
   *        yLineMark=real-value yLineOnOff=on|off
   *        windowingFunction=maximum|mean|minimum smoothingWindow=off|2-16
   *
   * @param trackLine
   */
  private void parseTrackLine(String trackLine)
  {
    String colour = "0,0,0";
    int beginIndex = trackLine.indexOf(" color=");
    if(beginIndex > -1)
    {
      beginIndex+=7;
      int endIndex   = trackLine.indexOf(" ", beginIndex);
      colour = trackLine.substring(beginIndex, endIndex);
    }
    
    incrementLines(LineAttributes.parse(colour));
  }
  
  private void incrementLines(Color colour)
  {
    LineAttributes line = new LineAttributes(colour);
    
    if(lines == null)
      lines = new LineAttributes[1];
    else
    {
      LineAttributes linesTmp[] = new LineAttributes[lines.length];
      System.arraycopy(lines, 0, linesTmp, 0, lines.length);
      lines = new LineAttributes[linesTmp.length+1];
      System.arraycopy(linesTmp, 0, lines, 0, linesTmp.length);
    }
    lines[lines.length-1] = line;
  }
  
  /**
   * Wiggle formats  (default: span=1) :
   * variableStep  chrom=chrN  [span=windowSize]
   * fixedStep  chrom=chrN  start=position  step=stepInterval  [span=windowSize]
   * @param line
   * @return
   */
  private int parseWiggle(String line) throws NumberFormatException
  {
    if(line.startsWith("variableStep "))
      FORMAT = WIGGLE_VARIABLE_STEP_FORMAT;
    else if (line.startsWith("fixedStep "))
      FORMAT = WIGGLE_FIXED_STEP_FORMAT;
    
    if(wiggle == null)
    {
      wiggle = new Wiggle[lines.length];
      for(int i=0; i<wiggle.length; i++)
        wiggle[i] = new Wiggle();
    }
    else
    {
      Wiggle wiggleTmp[] = new Wiggle[wiggle.length];
      System.arraycopy(wiggle, 0, wiggleTmp, 0, wiggle.length);
      wiggle = new Wiggle[wiggleTmp.length+1];
      System.arraycopy(wiggleTmp, 0, wiggle, 0, wiggleTmp.length);
      wiggle[wiggle.length-1] = new Wiggle();
    }
    
    if(wiggle.length > lines.length)
      incrementLines(lines[lines.length-1].getLineColour());
    
    if(FORMAT == WIGGLE_FIXED_STEP_FORMAT)
    {
      wiggle[wiggle.length-1].start = 
        Integer.parseInt(getSubString(" start=" ,line));

      wiggle[wiggle.length-1].step = 
        Integer.parseInt(getSubString(" step=" ,line));
    }
    
    int beginIndex = line.indexOf(" span=");
    if(beginIndex > -1)
    {
      wiggle[wiggle.length-1].span = 
        Integer.parseInt(getSubString(" span=" ,line));
    }
    
    return FORMAT;
  }


  /**
   * Find the value of a key within a string.
   * @param key
   * @param line
   * @return
   */
  private String getSubString(String key, String line)
  {
    int beginIndex = line.indexOf(key)+key.length();
    int endIndex   = line.indexOf(" ", beginIndex);
    if(endIndex == -1)
      endIndex = line.length();
    return line.substring(beginIndex, endIndex);
  }


  /**
   * Return any LineAttributes read from the header (for
   * BASE_SPECIFIED_FORMAT).
   * @return
   */
  public LineAttributes[] getLineAttributes()
  {
    return lines;
  }
  
  /**
   *  Return the number of values a call to getValues () will return - one
   *  in this case.
   **/
  public int getValueCount () 
  {
    if(FORMAT == BASE_SPECIFIED_FORMAT)
      return number_of_values - 1;
    return number_of_values;
  }

  /**
   *  Return the default or optimal window size.
   *  @return null is returned if this algorithm doesn't have optimal window
   *    size.
   **/
  public Integer getDefaultWindowSize () 
  {
    return new Integer (default_window_size);
  }

  /**
   *  Return the default maximum window size for this algorithm.
   *  @return null is returned if this algorithm doesn't have maximum window
   *    size.
   **/
  public Integer getDefaultMaxWindowSize ()
  {
    return new Integer (100);
  }

  /**
   *  Return the default minimum window size for this algorithm.
   *  @return null is returned if this algorithm doesn't have minimum window
   *    size.
   **/
  public Integer getDefaultMinWindowSize () 
  {
    return new Integer (1);
  }

  /**
   *  Return the default or optimal step size.
   *  @return null is returned if this algorithm doesn't have optimal step
   *    size.
   **/
  public Integer getDefaultStepSize (int window_size)
  {
    if (window_size > 10) 
      return new Integer (window_size / 10);
    else 
      return null;
  }

  /**
   *  Return the maximum value of this algorithm.
   **/
  protected Float getMaximumInternal ()
  {
    return new Float (data_max);
  }

  /**
   *  Return the minimum value of this algorithm.
   **/
  protected Float getMinimumInternal () 
  {
    return new Float (data_min);
  }

  /**
   *  Return the average value of function over the whole strand.
   **/
  public Float getAverage () 
  {
    return new Float (average_value);
  }
  
  public int getWiggleStart(int index)
  {
    return wiggle[index].start;
  }
  
  public int getWiggleSpan(int index)
  {
    return wiggle[index].span;
  }
  
  public boolean isWiggleFormat()
  {
    if(FORMAT == WIGGLE_VARIABLE_STEP_FORMAT || 
       FORMAT == WIGGLE_FIXED_STEP_FORMAT)
      return true;
    return false;
  }
  
  class Wiggle
  {
    int start;
    int step;
    int span = 0;
  }
}
