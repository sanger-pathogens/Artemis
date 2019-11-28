/* 
 * This file is part of Artemis
 *
 * Copyright (C) 2019  Genome Research Limited
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
package uk.ac.sanger.artemis.io;

import static org.junit.Assert.*;

import org.junit.Test;


/**
 * JUnit test for the Location class.
 * 
 * @author kp11
 *
 */
public class LocationTest
{
	
	/**
	 * Test the isComplement() method.
	 * @throws Exception
	 */
	@Test
	public void testIsComplement1() throws Exception
	{
		// Given
		Location loc = new Location("complement(3..10)");
		
		// When
		boolean result = loc.isComplement();
		
		// Then
		assertTrue(result);
	}
	
	/**
	 * Test the isComplement() method.
	 * @throws Exception
	 */
	@Test
	public void testIsComplement2() throws Exception
	{
		// Given
		Location loc = new Location("join(184264..184375,187655..188211,196610..196906,199841..199927)");
		
		// When
		boolean result = loc.isComplement();
		
		// Then
		assertFalse(result);
	}
	
	/**
	 * Test the isComplement() method.
	 * @throws Exception
	 */
	@Test
	public void testIsComplement3() throws Exception
	{
		// Given
		Location loc = new Location("order(3..10)");
		
		// When
		boolean result = loc.isComplement();
		
		// Then
		assertFalse(result);
	}
	
	/**
	 * Test the isComplement() method.
	 * @throws Exception
	 */
	@Test
	public void testIsComplement4() throws Exception
	{
		// Given
		Location loc = new Location("join(complement(701966..702762))");
		
		// When
		boolean result = loc.isComplement();
		
		// Then
		assertTrue(result);
	}
	
	/**
	 * Test the isComplement() method.
	 * @throws Exception
	 */
	@Test
	public void testIsComplement5() throws Exception
	{
		// Given
		Location loc = new Location("join(complement(249124..251322),complement(242210..243004))");
		
		// When
		boolean result = loc.isComplement();
		
		// Then
		assertTrue(result);
	}
	
	/**
	 * Test the isComplement() method, for range
	 * @throws Exception
	 */
	@Test
	public void testIsComplement6() throws Exception
	{
		// Given
		Location loc = new Location("10..30");
		
		// When
		boolean result = loc.isComplement();
		
		// Then
		assertFalse(result);
	}
	
	/**
	 * Test the isComplement() method, for range
	 * @throws Exception
	 */
	@Test
	public void testIsComplement7() throws Exception
	{
		// Given
		Location loc = new Location("order(complement(701966..702762))");
		
		// When
		boolean result = loc.isComplement();
		
		// Then
		assertTrue(result);
	}
	
}
