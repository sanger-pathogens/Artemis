/* ComparisonDataFactory.java
 *
 * created: Thu Jul 15 1999
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 1999-2018  Genome Research Limited
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

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.util.LinePushBackReader;

import java.io.*;
import java.util.LinkedList;
import java.util.List;

import org.apache.log4j.Logger;

/**
 *  This class contains the method readComparisonData (), which returns an
 *  appropriate ComparisonData object for a given Document.
 *
 *  @author Kim Rutherford
 *  @version $Id: ComparisonDataFactory.java,v 1.1 2004-06-09 09:44:16 tjc Exp $
 **/

public class ComparisonDataFactory {
	
  /** Logging instance. */
  private static Logger logger4j = Logger.getLogger(ComparisonDataFactory.class);
  
  /**
   *  This method creates an appropriate ComparisonData object from a Document.
   */
  static public ComparisonData readComparisonData (Document data_document)
      throws IOException {
    
    final Reader in_file = data_document.getReader ();
    String fileName  = data_document.getName();

    final LinePushBackReader pushback_reader =
      new LinePushBackReader (in_file);

    
    String line = null;
    List<String> headers = null;
    try {
    	line = peekFirstLine(pushback_reader, fileName);
    	headers = readHeaders(pushback_reader, fileName);
    } catch (IOException e) {
    	
    	try {
    		// close the reader
    		pushback_reader.close();
    	} catch (IOException ioe) {
    		// Ignore
    	}
    	
    	throw e;
    }
    

    if (BlastWebSiteHitTableComparisonData.formatCorrect (headers)) {
        logger4j.info("Loading Blast web site hit table comparison file: " + fileName);
        return new BlastWebSiteHitTableComparisonData (pushback_reader);
    }
    
    if (MSPcrunchComparisonData.formatCorrect (line)) {
      logger4j.info("Loading crunch comparison file: " + fileName);
      return new MSPcrunchComparisonData (pushback_reader);
    } else {
      if (SSAHAComparisonData.formatCorrect (line)) {
    	logger4j.info("Loading SSAHA comparison file: " + fileName);
        return new SSAHAComparisonData (pushback_reader);
      } else {
        if (MegaBlastComparisonData.formatCorrect (line)) {
    	  logger4j.info("Loading mega blast comparison file: " + fileName);
          return new MegaBlastComparisonData (pushback_reader);
        } else {
          if (BlastM8ComparisonData.formatCorrect (line)) {
        	logger4j.info("Loading Blast m8 comparison file: " + fileName);
            return new BlastM8ComparisonData (pushback_reader);
          } else {
        	
        	try {
        		// close the reader
        		pushback_reader.close();
        	} catch (IOException ioe) {
        		// Ignore
        	}
        	
        	logger4j.info("Failed to load ACT comparison file: " + fileName);
        	throw new IOException ("cannot understand the comparison file format");
          }
        }
      }
    }
      
  }
  
  protected static List<String> readHeaders(LinePushBackReader reader, String fileName) throws IOException {
	  
	  List<String> headerList = new LinkedList<String>();
	  String line = null;
	  boolean finished = false;
	  
	  do {
		  line = reader.readLine();
		  
		  if (line == null) {
			  throw new IOException (
					  "End of file while reading from: " +
                      fileName);
		  }
		  
		  if (line.startsWith("#")) {
			  headerList.add(line);
		  } else {
			  finished = true;
		  }
		  
	  } while (!finished);
	  
	  reader.pushBack(line);
	  
	  return headerList;
  }
  
  protected static String peekFirstLine(LinePushBackReader reader, String fileName) throws IOException {
	  
	  final String line = reader.readLine ();

	  if (line == null) {
    	
		  throw new IOException (
				  "End of file while reading from: " +
                  fileName);
	  }
	  
	  reader.pushBack(line);
	  
	  return line;
	    
  }
}

