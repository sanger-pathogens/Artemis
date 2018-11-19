/* BlastServerTableComparisonData.java
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2018  Genome Research Limited
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
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.util.LinePushBackReader;

import java.io.*;
import java.util.List;
import java.util.StringTokenizer;

/**
 *  This class implements the SimpleComparisonData interface 
 *  for blast web site hit table output.
 *
 *  @author kp11
 */
public class BlastWebSiteHitTableComparisonData extends SimpleComparisonData {
  
  /** Min number of fields in web site blastn/tblastx file. */
  public static final int MIN_NUM_FIELDS = 12;
  
  /** Max number of fields in web site blastn/tblastx file. */
  public static final int MAX_NUM_FIELDS = 14;
  
  /** Comparison file descriptive type. */
  public static final String TYPE = "Blast web site hit table comparison data";
  
  /**
   *  Create a new BlastWebSiteHitTableComparisonData by reading from the given
   *  LinePushBackReader.
   */
  public BlastWebSiteHitTableComparisonData (final LinePushBackReader stream)
      throws IOException {
    super (stream);
  }

  /**
   *  Create a new, empty instance of BlastServerTableComparisonData.
   */
  public BlastWebSiteHitTableComparisonData () {
    
  }

  /**
   *  Returns a new, empty instance of this type of object;
   */
  protected SimpleComparisonData getNewSimpleComparisonData () {
    return new BlastWebSiteHitTableComparisonData ();
  }


  /**
   *  Make an AlignMatch object from the given String.
   */
  private static AlignMatch makeMatchFromStringStatic (String line)
      throws IOException {
    
	if (line.trim ().length () == 0 || line.startsWith ("#")) {
      return null;
    }

    final StringTokenizer tokenizer = new StringTokenizer (line, "\t");

    int numTokens = tokenizer.countTokens ();
    if (
    		numTokens < BlastWebSiteHitTableComparisonData.MIN_NUM_FIELDS ||
    		numTokens > BlastWebSiteHitTableComparisonData.MAX_NUM_FIELDS) {
      
      final String message = "while reading " + 
    		  BlastWebSiteHitTableComparisonData.TYPE +
    		  					": unexpected number of fields for this line: " + line;
      throw new ComparisonDataParseException (message);
      
    }
    
    // Parse fields from line...
    
    // throw away the query name
    tokenizer.nextToken ();
    
    // throw away the subject name
    tokenizer.nextToken ();
    
    // % ident
    final String percentIdentToken = tokenizer.nextToken ();
    
    // throw away alignment length
    tokenizer.nextToken ();
    
    // throw away mismatches
    tokenizer.nextToken ();
    
    // throw away gap opens
    tokenizer.nextToken ();
    
    // Query start/end
    final String qStartToken = tokenizer.nextToken ();
    final String qEndToken = tokenizer.nextToken ();
    
    final String sStartToken = tokenizer.nextToken ();
    final String sEndToken   = tokenizer.nextToken ();
    
    // throw away evalue
    tokenizer.nextToken ();
    
    // % bit score
    final String scoreToken = tokenizer.nextToken ();
    
    // And ignore all the end fields for tblastx

    try {
    	
      final int score = (int)(Float.valueOf (scoreToken).floatValue ());
      final int percentIdent = (int)(Float.valueOf (percentIdentToken).floatValue ());
      final int qStart = Integer.valueOf (qStartToken).intValue ();
      final int qEnd   = Integer.valueOf (qEndToken).intValue ();
      final int sStart = Integer.valueOf (sStartToken).intValue ();
      final int sEnd   = Integer.valueOf (sEndToken).intValue ();

      return makeAlignMatch (sStart, sEnd, qStart, qEnd, score,
                             percentIdent);
      
    } catch (NumberFormatException e) {
      throw new IOException ("while reading " +
    		  BlastWebSiteHitTableComparisonData.TYPE +
                             ": failed to parse a number from this string: " +
                             e.getMessage ());
    }
  }

  /**
   *  Make an AlignMatch object from the given String.  The String must be in
   *  a format appropriate for this object.
   */
  @Override
  protected AlignMatch makeMatchFromString (final String line)
      throws IOException {
    return makeMatchFromStringStatic (line);
  }

  /**
   *  Returns true if and only if the given line is in the correct format for
   *  this type of ComparisonData.  This should be as strict as possible.
   */
  public static boolean formatCorrect (final List<String> headers) {
	
	 boolean result = false;
	 
	 if (headers.size() >= 4) {
		 
		 String blastTypeHeader = headers.get(0);
		 String iterationHeader = headers.get(1);
		 String queryHeader = headers.get(2);
		 String ridHeader = headers.get(3);
		 
		 if (
				 (blastTypeHeader.startsWith ("# tblastx") || blastTypeHeader.startsWith ("# blastn")) &&
			 	 iterationHeader.startsWith ("# Iteration:") &&
				 queryHeader.startsWith ("# Query:") &&
				 ridHeader.startsWith ("# RID:")) {
				 		
			result = true;	 
		 }
	 }
	 
	 return result;
  }
  
}
