/*
 * Copyright (C) 2009  Genome Research Limited
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
package uk.ac.sanger.artemis.circular.digest;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class EmbossTableParser
{
	List<CutSite> list = new ArrayList<CutSite>();
	private int length;

	public List<CutSite> parse(BufferedReader br) throws IOException
	{
		String line;
		while ((line = br.readLine()) != null)
		{
			String trim = line.trim();
			if (trim.startsWith("#"))
			{
				if (trim.indexOf("Sequence:") != -1)
				{
					int pos = trim.indexOf("to: ");
					String l = trim.substring(pos + 4);
					length = Integer.parseInt(l);
				}
				continue;
			}
			
			if (trim.length() == 0)
				continue;
			
			String[] parts = trim.split("\\s+");
			if (parts.length < 3)
				continue;
			
			if ("Start".equals(parts[0]))
				continue;
			
			CutSite cutSite;
			
			if(parts.length > 8) // new EMBOSS format
			  cutSite = new CutSite(
					parts[3], parts[5], parts[6], parts[7], parts[8], parts[2]);
			else
			  cutSite = new CutSite(
			      parts[2], parts[4], parts[5], parts[6], parts[7], "+");
			list.add(cutSite);
		}
		return list;
	}

	public int getLength()
	{
		return this.length;
	}

}
