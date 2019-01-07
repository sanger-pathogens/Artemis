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
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;

import javax.swing.JOptionPane;

import uk.ac.sanger.artemis.editor.BrowserControl;

/**
 * Pfam and Rfam sequence search
 */
public class RunPfamSearchThread extends Thread
{
  protected static String pfamUrl = "http://pfam.xfam.org/search/sequence";
  protected static String rfamUrl = "http://rfam.xfam.org/search/sequence";
  private String searchURL = pfamUrl;
  
  private String residues;
  
  public RunPfamSearchThread(final String residues)
  {
    this.residues = residues;
  }
  
  public RunPfamSearchThread(final String residues, final String searchURL)
  {
    this.residues = residues;
    this.searchURL = searchURL;
  }
  
  public void run()
  {
    boolean isPfam = searchURL.equals(pfamUrl);
    
    ProgressBarFrame progress = null;

    try
    {
      // Construct data
      String data = URLEncoder.encode("seq", "UTF-8") + "="
          + URLEncoder.encode(residues, "UTF-8");
      data += "&" + URLEncoder.encode("output", "UTF-8") + "="
          + URLEncoder.encode("xml", "UTF-8");

      // Send data
      URL url = new URL(searchURL);

      URLConnection conn = url.openConnection();
      conn.setDoOutput(true);
      
      if(!isPfam)
        conn.addRequestProperty("Accept", "text/xml");
      
      OutputStreamWriter wr = new OutputStreamWriter(conn.getOutputStream());
      wr.write(data);
      wr.flush();

      // Get the response
      BufferedReader rd = new BufferedReader(
          new InputStreamReader(conn.getInputStream()));
      String urlResults = searchURL + (isPfam ? "/results?" : "/");

      String line;
      String eta = "5";
      while ((line = rd.readLine()) != null)
      {
        int index;
        if((index = line.indexOf("job_id=")) > -1)
        {
          line = line.substring(index+8, line.length()-2);
          if(isPfam)
            urlResults = urlResults.concat("jobId="+line);
          else
            urlResults = urlResults.concat(line);
        }
        else if((index = line.indexOf("<estimated_time>")) > -1)
        {
          eta = line.substring(index+16);
          index = eta.indexOf("<");
          if(index > -1)
            eta = eta.substring(0, index);
        }
        else if((index = line.indexOf("<h2>Error</h2>")) > -1)
        {
          // Got an error tag in the page so read the 
          // error message in the next line...
          String next  = rd.readLine();
          String error = "Pfam search error";
         
          if (next != null)
          {
        	  // Get rid of any double quotes and whitespace
        	  String errorLine = next.replace("\"", "").trim();
        	  if (errorLine.length() > 0)
        	  {
        		  error = errorLine;
        	  }
          }
          
          JOptionPane.showMessageDialog(null,
        		  error, "Pfam search error", JOptionPane.INFORMATION_MESSAGE);
          return;
        }
      }
      wr.close();
      rd.close();
      
      int waitTime = Integer.parseInt(eta);
      progress = new ProgressBarFrame(waitTime, (isPfam ? "Pfam" : "Rfam"));
      URL result = new URL(urlResults);
      Thread.sleep(waitTime*1000);
      
      int cnt = 0;
      while(((HttpURLConnection) result.openConnection()).getResponseCode() == 204 &&
            cnt < 500)
      {
        cnt++;
        Thread.sleep(500);
      }

      BrowserControl.displayURL(urlResults);
    }
    catch (Exception e)
    { 
    	e.printStackTrace(); 
    }
    finally
    {
    	if (progress != null)
    	{
    		progress.dispose();
    	}
    }
  }
}