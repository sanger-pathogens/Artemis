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

public class CutSite
{
	private String enzymeName;
	private int fivePrime;
	private int threePrime;
	private int fivePrimeRev;
	private int threePrimeRev;
	private boolean forward = false;
  private boolean highlighted = false;

	CutSite(final String enzymeName,
			    String fivePrimeStr, String threePrimeStr,
			    String fivePrimeRevStr, String threePrimeRevStr,
			    String strand)
	{
		this.enzymeName = enzymeName;
		this.fivePrime  = Integer.parseInt(fivePrimeStr);
		this.threePrime = Integer.parseInt(threePrimeStr);
		if(!fivePrimeRevStr.equals("."))
			fivePrimeRev = Integer.parseInt(fivePrimeRevStr);
		if(!threePrimeRevStr.equals("."))
			threePrimeRev = Integer.parseInt(threePrimeRevStr);
		
		if(strand.equals("+"))
			forward = true;
	}

	public int getThreePrime()
	{
		return threePrime;
	}

	public int getFivePrime()
	{
		return fivePrime;
	}

	public int getFivePrimeRev()
	{
		return fivePrimeRev;
	}

	public int getThreePrimeRev()
	{
		return threePrimeRev;
	}

	public boolean isForward()
	{
		return forward;
	}
	
	public boolean isHighlighted()
	{
		return highlighted;
	}

	public void setHighlighted(boolean highlighted)
	{
		this.highlighted = highlighted;
	}
	
	public String getEnzymeName()
	{
		return enzymeName;
	}
}
