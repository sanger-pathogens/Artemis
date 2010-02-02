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

public class RunPfamSearchThread extends Thread
{
  private static String pfamUrl = "http://pfam.sanger.ac.uk/search/sequence";

  private String residues;
  
  public RunPfamSearchThread(final String residues)
  {
    this.residues = residues;
  }
  
  public void run()
  {
    try
    {
      // Construct data
      String data = URLEncoder.encode("seq", "UTF-8") + "="
          + URLEncoder.encode(residues, "UTF-8");
      data += "&" + URLEncoder.encode("output", "UTF-8") + "="
          + URLEncoder.encode("xml", "UTF-8");

      // Send data
      URL url = new URL(pfamUrl);

      URLConnection conn = url.openConnection();
      conn.setDoOutput(true);
      OutputStreamWriter wr = new OutputStreamWriter(conn.getOutputStream());
      wr.write(data);
      wr.flush();

      // Get the response
      BufferedReader rd = new BufferedReader(
          new InputStreamReader(conn.getInputStream()));
      String urlResults = pfamUrl+"/results?";
      String line;
      String eta = "5";
      while ((line = rd.readLine()) != null)
      {
        int index;
        if((index = line.indexOf("job_id=")) > -1)
        {
          line = line.substring(index+8, line.length()-2);
          urlResults = urlResults.concat("jobId="+line);
        }
        else if((index = line.indexOf("<estimated_time>")) > -1)
        {
          eta = line.substring(index+16);
          index = eta.indexOf("<");
          if(index > -1)
            eta = eta.substring(0, index);
        }
        else if((index = line.indexOf("<error>")) > -1)
        {
          line = line.substring(index+7);
          index = line.indexOf("<");
          if(index > -1)
            line = line.substring(0, index);
          JOptionPane.showMessageDialog(null,
              line, "Pfam search error", JOptionPane.INFORMATION_MESSAGE);
          return;
        }
      }
      wr.close();
      rd.close();
      
      int waitTime = Integer.parseInt(eta);
      ProgressBarFrame progress = new ProgressBarFrame(waitTime, "Pfam");
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
      progress.dispose();
    }
    catch (Exception e){}
  }
}