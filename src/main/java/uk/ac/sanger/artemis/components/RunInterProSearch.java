/* RunPfamSearch.java
 *
 * created: 2009
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2009  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or(at your option) any later version.
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
 **/

package uk.ac.sanger.artemis.components;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLEncoder;
import java.net.UnknownHostException;

import javax.swing.JOptionPane;

import org.apache.log4j.Logger;

import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.editor.BrowserControl;

/**
 * InterPro sequence search functionality.
 * 
 * @author kp11
 */
public class RunInterProSearch
{
  /** The InterPro web site URL */
  protected String searchURL = Options.getOptions().getProperty("interpro_search_url");
  
  /** Title of dialog to display when a search error is encountered. */
  protected static final String SEARCH_ERROR_DIALOG_TITLE = "InterPro search error";
  
  /** Sequence to search on. */
  private String sequence;
  
  /** Logging instance. */
  private static Logger logger = Logger.getLogger(RunInterProSearch.class);
  
  /**
   * Default constructor.
   */
  public RunInterProSearch()
  {
    // Do nothing
  }
  
  /**
   * Constructor.
   * @param sequence - protein sequence to search on.
   */
  public RunInterProSearch(final String sequence)
  {
    this.sequence = sequence;
  }
  
  /**
   * Constructor.
   * @param sequence - protein sequence to search on.
   * @paramsearchURL - the InterPro search URL.
   */
  public RunInterProSearch(final String sequence, final String searchURL)
  {
    this.sequence = sequence;
    this.searchURL = searchURL;
  }
  
  /**
   * Run a search.
   */
  public void start()
  { 
	OutputStreamWriter wr = null;
	BufferedReader rd = null;
	  
    try
    {
      // Construct POST data
      String data = URLEncoder.encode("queryString", "UTF-8") + "="
          + URLEncoder.encode(sequence, "UTF-8");

      // Send data to search URL

      HttpURLConnection conn = getConnection();
      
      conn.setDoOutput(true);
      conn.setInstanceFollowRedirects(false);
      
      wr = new OutputStreamWriter(conn.getOutputStream());
      wr.write(data);
      wr.flush();

      // Get the response
      rd = new BufferedReader(
          new InputStreamReader(conn.getInputStream()));

      // Check response status
      int status = conn.getResponseCode();
      if (status != 200 && status != 302) 
      {
    	  // [not Accepted]
    	  
    	  logger.error("ERROR: HTTP error received from EBI InterPro website: " + status);
    	  showError("HTTP error received from EBI InterPro website: " + status);
          return;
      }
      
      // Get the redirection URL
      String location = conn.getHeaderField( "Location" );
      
      // Display page in browser
      displayURL(location);
      
    }
    catch (UnknownHostException e)
    { 
    	logger.error("ERROR: Failed to send InterPro query due to an unknown host error: " + e.getMessage());
    	showError("Cannot contact address " + e.getMessage());
    }
    catch (Exception e)
    { 
    	logger.error("ERROR: Failed to send InterPro query: " + e.getMessage());
    	showError("Error : " + e.getMessage());
    }
    finally
    {
    	// Close streams
    	
    	if (wr != null)
    	{
    		try
    		{
    			wr.close();
    		}
    		catch (Exception e)
    		{
    			// Ignore
    		}
    	}
    	

    	if (rd != null)
    	{
    		try
    		{
    			rd.close();
    		}
    		catch (Exception e)
    		{
    			// Ignore
    		}
    	}
    }
  }
  
  /**
   * Create and return an HTTP connection object.
   * @return HttpURLConnection
   * @throws IOException
   */
  protected HttpURLConnection getConnection() throws IOException
  {
	  URL url = new URL(searchURL);

      return (HttpURLConnection)url.openConnection();
  }
  
  /**
   * Invoke a browser to display the specified URL.
   * @param url String
   */
  protected void displayURL(String url)
  {
	  BrowserControl.displayURL(url);
  }
  
  /**
   * Display an error message to the user.
   * @param message String
   */
  protected void showError(String message)
  {
	  JOptionPane.showMessageDialog(null,
			  message, SEARCH_ERROR_DIALOG_TITLE, JOptionPane.INFORMATION_MESSAGE);
  }

  /**
   * Get method for the InterPro search URL.
   * @return String - url
   */
  public String getSearchURL()
  {
	  return searchURL;
  }

  /**
   * The amino acid sequence to search on.
   * @return String - amino acid letters
   */
  public String getSequence()
  {
	  return sequence;
  }

  /**
   * Return the InterPro URL.
   * @param searchURL String
   */
  public void setSearchURL(String searchURL)
  {
	this.searchURL = searchURL;
  }
	
  /**
   * Set the amino acid sequence to send.
   * @param sequence String
   */
  public void setSequence(String sequence)
  {
	this.sequence = sequence;
  }
}