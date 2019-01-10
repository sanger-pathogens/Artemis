/* 
 * This file is part of Artemis
 *
 * Copyright (C) 2017  Genome Research Limited
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
package uk.ac.sanger.artemis.components;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.net.URL;
import java.net.URLConnection;

import org.junit.Test;

import uk.ac.sanger.artemis.components.RunBlastAtNCBI;


/**
 * Run some Blast queries.
 * As this test class fires queries at the NCBI database it should only
 * really be run on an adhoc basis to test this functionality is
 * working as expected.
 * 
 * @author kp11
 *
 */
public class RunBlastAtNCBITest {

	private String sendRequestToBlast(String data) throws Exception
	{
		
		URL url = new URL("https://blast.ncbi.nlm.nih.gov/Blast.cgi");
		URLConnection conn = url.openConnection();
		conn.setDoOutput(true);
		
		System.out.println("RunBlastAtNCBITest: Established connection to NCBI");
		
		OutputStreamWriter wr = new OutputStreamWriter(conn.getOutputStream());
		wr.write(data);
		wr.flush();

		System.out.println("RunBlastAtNCBITest: Data request flushed to target");
		
		// Get the response
		BufferedReader rd = new BufferedReader(new InputStreamReader(conn.getInputStream()));
		String urlResults = url + "?RID=";
		String line;
		String eta = "10";
		while ((line = rd.readLine()) != null) {
			int index;
			if ((index = line.indexOf("RID =")) > -1) {
				line = line.substring(index + 5).trim();
				urlResults = urlResults.concat(line);
				if (line.equals("")) {
					org.junit.Assert.fail("Call to Blast URL failed - no result returned");
				}
			} else if ((index = line.indexOf("RTOE =")) > -1)
				eta = line.substring(index + 6).trim();
		}
		
		System.out.println("RunBlastAtNCBITest: About to close connection");
		
		wr.close();
		rd.close();
		
		System.out.println("RunBlastAtNCBITest: Connection closed");
		
		//int waitTime = Integer.parseInt(eta);
		//Thread.sleep(waitTime * 1000);
		
		return urlResults + "&CMD=Get&OLD_BLAST=false";
	}

	@Test
	public void testBlastNRequest() throws Exception 
	{

		System.out.println("RunBlastAtNCBITest: About to construct blastn request string...");
		
		String programName =  "blastn";
		String residues = "aggctgttttccacagatttcacagtattggttcaaatggtcaaaaattgttttaaccagt";
		
		String request = RunBlastAtNCBI.constructRequest(
				programName,	
				residues,
			  	"nr",										// Database
			  	"500", 										// Hits
			  	"None",										// Filter - "None", "Low Complexity", "Human Repeats", "Masked"
			  	"10.0",										// Expect
			  	"plain",										// Service - "plain" or "megablast"
			  	(programName.equals("blastn") ? "5": "11"), 	// Gap Open
			  	(programName.equals("blastn") ? "2": "1") 	// Gap Close
		);
		
		System.out.println("RunBlastAtNCBITest: Request string => " + request);
		
		String resultUrl = sendRequestToBlast(request);
		
		System.out.println("RunBlastAtNCBITest: Blastn request completed. Result URL => " + resultUrl);
		
		assertNotNull(resultUrl);
		assertTrue("Check results URL looks valid using a reg exp", resultUrl.matches("^https[:]//blast[\\.]ncbi[\\.]nlm[\\.]nih[\\.]gov/Blast[\\.]cgi[\\?].*"));
		
	}
	
	/*
	 * This test seems to take a few 10s of seconds to be able to view the actual results
	 * in a browser if you need to.
	 */
	@Test
	public void testBlastPRequest() throws Exception 
	{

		System.out.println("RunBlastAtNCBITest: About to construct blastp request string...");
		
		String programName =  "blastp";
		String residues = "ATHIEDLHNITSNQLYETYRTEKLSTSQLLLDSTVXTIDKNLSQHDQVLREDRLR";
		
		String request = RunBlastAtNCBI.constructRequest(
				programName,	
				residues,
			  	"nr",										// Database
			  	"500", 										// Hits
			  	"None",										// Filter - "None", "Low Complexity", "Human Repeats", "Masked"
			  	"10.0",										// Expect
			  	"plain",										// Service - "plain" or "megablast"
			  	(programName.equals("blastn") ? "5": "11"), 	// Gap Open
			  	(programName.equals("blastn") ? "2": "1") 	// Gap Close
		);
		
		System.out.println("RunBlastAtNCBITest: Request string => " + request);
		
		String resultUrl = sendRequestToBlast(request);
		
		System.out.println("RunBlastAtNCBITest: Blastp request completed. Result URL => " + resultUrl);
		
		assertNotNull(resultUrl);
		assertTrue("Check results URL looks valid using a reg exp", resultUrl.matches("^https[:]//blast[\\.]ncbi[\\.]nlm[\\.]nih[\\.]gov/Blast[\\.]cgi[\\?].*"));
		

	}
}
