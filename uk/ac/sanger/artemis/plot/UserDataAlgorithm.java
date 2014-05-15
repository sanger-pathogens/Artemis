/* UserDataAlgorithm.java
 * This file is part of Artemis
 *
 * Copyright (C) 2000-2012  Genome Research Limited
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
 */

package uk.ac.sanger.artemis.plot;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.Strand;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.FileDocument;
import uk.ac.sanger.artemis.util.LinePushBackReader;
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.components.variant.TabixReader;
import uk.ac.sanger.artemis.io.IndexFastaStream;
import uk.ac.sanger.artemis.io.ReadFormatException;
import java.awt.Color;
import java.awt.GridLayout;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Pattern;
import javax.swing.ButtonGroup;
import javax.swing.JCheckBox;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

import net.sf.samtools.util.BlockCompressedInputStream;

/**
 *  Objects of this class have one useful method - getValues (), which takes a
 *  range of bases and returns a single floating point number.  The number is
 *  calculated by averaging the values from a data file.  The Strand to use is
 *  set in the constructor.
 *  @author Kim Rutherford
 **/

public class UserDataAlgorithm extends BaseAlgorithm
{
  /** A base per line file format */
  public static int BASE_PER_LINE_FORMAT  = 1;
  /** Base position is specified in the first column file format */
  public static int BASE_SPECIFIED_FORMAT = 2;
  /** Wiggle format */
  public static int WIGGLE_VARIABLE_STEP_FORMAT = 3;
  public static int WIGGLE_FIXED_STEP_FORMAT    = 4;
  public static int BLAST_FORMAT                = 5;
  public static int MSPCRUNCH_BLAST_FORMAT      = 6;
  public static int TABIX_INDEXED_FORMAT        = 7;
  
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
  private Wiggle wiggle[];
  private TabixIdxGraph idxReader;
  private static org.apache.log4j.Logger logger4j = 
      org.apache.log4j.Logger.getLogger(UserDataAlgorithm.class);
  
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
  public UserDataAlgorithm (final Strand strand, final Document doc, 
                            final boolean logTransform)
      throws IOException 
  {
    super (strand, 
        (logTransform ? "User log plot " + doc.getName () : 
                        "User plot " + doc.getName ()), "user");
    this.logTransform = logTransform;

    if(isIndexed(doc) && doc.getInputStream() instanceof BlockCompressedInputStream)
    {
      setAlgorithmName(
          (logTransform ? "Indexed log plot " + doc.getName () :
                          "Indexed plot " + doc.getName ()));
      FORMAT = TABIX_INDEXED_FORMAT;
      idxReader = new TabixIdxGraph(
          ((FileDocument) doc).getFile().getAbsolutePath());
      number_of_values = idxReader.getNumberOfValues();
      return;
    }
    
    final Reader doc_reader = doc.getReader ();
    final LinePushBackReader pushback_reader = new LinePushBackReader (doc_reader);
    String firstLn = pushback_reader.readLine (); 

    final Pattern dataPattern = Pattern.compile("^\\s*([\\d\\.-]+\\s*)+$");
    final Pattern blastPattern = Pattern.compile(
      "^(\\S+\\t+){2}[\\d\\.]+\\t+(\\d+\\t+){7}\\S+\\t+(\\s*\\d+)$");
    final Pattern mspCrunchPattern = Pattern.compile(
        "^\\d+\\s[\\d\\.]+(\\s\\d+){2}\\s\\D\\S+(\\s\\d+){2}\\s\\D\\S+.*");

    if(dataPattern.matcher(firstLn).matches())
      FORMAT = BASE_PER_LINE_FORMAT;
    else if(blastPattern.matcher(firstLn).matches())
      FORMAT = BLAST_FORMAT;
    else if(mspCrunchPattern.matcher(firstLn).matches())
      FORMAT = MSPCRUNCH_BLAST_FORMAT;
    else
    { 
      StringBuffer header = new StringBuffer(firstLn+"\n");
      while(!dataPattern.matcher(firstLn).matches())
      {
        firstLn = pushback_reader.readLine ().trim();
        header.append(firstLn+"\n");
      }
      
      FORMAT = parseHeader(header);
    }

    final Pattern patt = Pattern.compile("\\s+");
    String tokens[] = patt.split(firstLn);
    
    if (tokens.length < 1) 
      throw new ReadFormatException ("unknown file type");

    this.number_of_values = tokens.length;
    pushback_reader.pushBack (firstLn);
    
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
    doc_reader.close();
  }
  
  /**
   * Test if the tabix (.tbi) index is present.
   * @param doc
   * @return
   */
  private static boolean isIndexed(Document doc)
  {
    if(doc instanceof FileDocument)
      return (new File(
          ((FileDocument)doc).getFile().getAbsolutePath() + ".tbi")).exists();
    return false;
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
    
    final boolean useEstimate = 
        Options.getOptions ().getIntegerProperty (getAlgorithmShortName () +
        "_default_window_size") == null;
    if(!useEstimate)
      default_window_size = Options.getOptions ().getIntegerProperty (getAlgorithmShortName () +
          "_default_window_size");
    
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
          throw new ReadFormatException ("at line number: "+(count+1)+" cannot understand this number: "+
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
      if(estimate_window_size != Integer.MAX_VALUE && useEstimate)
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
    else if(FORMAT == TABIX_INDEXED_FORMAT)
      idxReader.getValues(start, end, values);
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
            line.indexOf("color")  == -1 &&
            line.indexOf("label ") == -1) ||
            line.startsWith("track"))
        {
          if(line.startsWith("track "))
            parseTrackLine(line);
          else  if(line.startsWith("variableStep ") ||
                   line.startsWith("fixedStep"))
            FORMAT = parseWiggle(line);

          continue;
        }

        int idx = line.indexOf("colour");
        if (idx == -1)
          idx = line.indexOf("color");

        if (idx != -1)
        {
          idx = line.indexOf(" ", idx + 1);
          line = line.substring(idx).trim();
          String rgbValues[] = line.split(" ");

          lines = new LineAttributes[rgbValues.length];
          for (int j = 0; j < rgbValues.length; j++)
            lines[j] = new LineAttributes(LineAttributes.parse(rgbValues[j]));
        }
        else if(lines != null)
        {
          idx = line.indexOf("label ");
          if (idx != -1)
          {
            line = line.substring(idx+6).trim();
            String labels[] = line.split(" ");
            for (int j = 0; j < labels.length; j++)
              lines[j].setLabel(labels[j]);
          }
        }
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
  private void parseTrackLine(final String trackLine)
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
  
  private void incrementLines(final Color c)
  {
    LineAttributes line = new LineAttributes(c);
    
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
    if(Options.getOptions ().getIntegerProperty (getAlgorithmShortName () +
       "_default_max_window") == null)
      return new Integer (100);
    else
      return super.getDefaultMaxWindowSize ();
  }

  /**
   *  Return the default minimum window size for this algorithm.
   *  @return null is returned if this algorithm doesn't have minimum window
   *    size.
   **/
  public Integer getDefaultMinWindowSize () 
  {
    if(Options.getOptions ().getIntegerProperty (getAlgorithmShortName () +
       "_default_min_window") == null)
      return new Integer (1);
    else
      return super.getDefaultMinWindowSize();
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
  
  protected boolean isWiggleFormat()
  {
    if(FORMAT == WIGGLE_VARIABLE_STEP_FORMAT || 
       FORMAT == WIGGLE_FIXED_STEP_FORMAT)
      return true;
    return false;
  }
  
  public void readIndexValues(boolean recalculate_flag, Entry seqEntry, int start, int end)
  {
    if(start<1)
      start = 1;
    idxReader.readValuesForRange(recalculate_flag, seqEntry, start, end);
  }
  
  class Wiggle
  {
    int start;
    int step;
    int span = 0;
  }
  
  class TabixIdxGraph
  {
    private TabixReader reader;
    /** number of columns with values */
    private int nValues = 1;
    private float[][] rvalues;
    private boolean startColIsEndCol = true;
    private int sbeg, send;
    
    TabixIdxGraph(final String fName) throws IOException
    {
      this.reader = new TabixReader(fName);
      if(getStartColumn() != getEndColumn())
        startColIsEndCol = false;

      final StringBuffer headerBuffer = new StringBuffer();
      String metaChar = String.valueOf(reader.getCommentChar());
      String hdr;
      while( (hdr = reader.readLine() ) != null && hdr.startsWith(metaChar) )
        headerBuffer.append(hdr + "\n");
      //System.out.println("Header "+headerBuffer.toString());

      // assume base end column is last column before columns of values
      nValues = reader.readLine().split("\\t").length - getEndColumn();
    }
    
    /**
     * Return the values between a pair of bases
     * @param start
     * @param end
     * @param values
     */
    private void getValues(int start, int end, final float[] values)
    {
      for (int i = 0 ; i < getNumberOfValues() ; ++i) 
      {
        values [i] = 0;
        for (int base = start ; base <= end ; ++base) 
          values [i] += rvalues[base - this.sbeg][i] / (end - start + 1);
      }
    }
    
    /**
     * Get the sequence name
     * @param seqEntry - sequence entry
     * @return
     */
    private String getReferenceName(Entry seqEntry)
    {
      String refStr = null;
      if(seqEntry.getEMBLEntry().getSequence() instanceof IndexFastaStream)
        refStr = 
           ((IndexFastaStream)seqEntry.getEMBLEntry().getSequence()).getContig();
      else if(seqEntry.getHeaderText() != null)
      {
        final String hdr = seqEntry.getHeaderText();
        int idx = hdr.indexOf("ID   ");
        if (idx == -1)
        {
          idx = hdr.indexOf("LOCUS       "); // genbank
          if (idx > -1)
            refStr = hdr.substring(idx + 12).split("[;\\s]")[0];
        }
        else
          refStr = hdr.substring(idx + 5).split("[;\\s]")[0];
      }

      final String seqNames[] = reader.getSeqNames();
      if(refStr == null || !Arrays.asList(seqNames).contains(refStr))
      {
        logger4j.debug(refStr+" NOT FOUND IN "+reader.getFileName()+
            " SET TO DEFAULT "+seqNames[0]);
        refStr = reader.getSeqNames()[0];
      }
      return refStr;
    }

    /**
     * Read the values in a range
     * @param recalculate_flag
     * @param seqEntry - sequence entry
     * @param start - start base
     * @param end   - end base
     */
    private void readValuesForRange(boolean recalculate_flag, Entry seqEntry, int start, int end)
    {
      if(!recalculate_flag && (end <= start || (start == this.sbeg && end == this.send)) )
        return;

      this.sbeg = start;
      this.send = end;
      rvalues = new float[end-start+1][getNumberOfValues()];
      final String r = getReferenceName(seqEntry)+":"+start+"-"+end;

      try
      {
        final TabixReader.Iterator tbxIt = reader.query(r);
        String ln;
        while( tbxIt != null &&  (ln = tbxIt.next()) != null )
        {
          StringVector parts = StringVector.getStrings(ln, "\t", true);
          final int base;
          if(startColIsEndCol)
            base = Integer.parseInt((String)parts.get(getStartColumn()-1)) - start;
          else
          {
            int b = Integer.parseInt(parts.get(getStartColumn()-1));
            int e = Integer.parseInt(parts.get(getEndColumn()-1));
            base = b + ((b - e)/2) - start;
          }
          for(int i=0; i<rvalues[base].length; i++)
          {
            float val = Float.parseFloat( parts.get(i+getEndColumn()) );
            if(logTransform)
              val = (float) Math.log(val+1);

            if (val > data_max) 
              data_max = val;
            if (val < data_min)
              data_min = val;
            rvalues[base][i] = val;
          }
        }
      }
      catch (IOException e)
      {
        logger4j.debug("IOException READING RANGE "+r+" FROM "+reader.getFileName());
        e.printStackTrace();
      }
      catch (NumberFormatException e)
      {
        logger4j.debug("NumberFormatException READING RANGE "+r+" FROM "+reader.getFileName());
        e.printStackTrace();
      }
    }
    
    /**
     * Number of columns with values
     * @return
     */
    private int getNumberOfValues()
    {
      return nValues;
    }
    
    private int getStartColumn()
    {
      return reader.getStartColumn();
    }
    
    private int getEndColumn()
    {
      return reader.getEndColumn();
    }
  }
}
