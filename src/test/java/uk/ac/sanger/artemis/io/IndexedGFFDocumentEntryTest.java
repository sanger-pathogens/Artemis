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

import static org.junit.Assert.fail;

import java.io.IOException;
import org.junit.Test;

import uk.ac.sanger.artemis.components.variant.NoFeaturesException;
import uk.ac.sanger.artemis.components.variant.TabixReader;

/**
 * Unit test for the IndexedGFFDocumentEntry class.
 * It currently only covers new functionality/changes.
 * 
 * @author kp11
 *
 */
public class IndexedGFFDocumentEntryTest
{

	public static final String TEST_INDEXED_GFF_FILE = "/data/gff/indexed-gff.gff.gz";
	
	/**
	 * Test the query method.
	 * @throws IOException 
	 */
	@Test
	public void testQuery() throws IOException 
	{	
		TabixReader reader = new TabixReader(
				IndexedGFFDocumentEntryTest.class.getResource(TEST_INDEXED_GFF_FILE).getFile());
		
		try 
		{
			reader.query(-1, 1, 19);
			fail("Expected a NoFeaturesException to be thrown");
		} 
		catch (NoFeaturesException e)
		{
			// expected
		}
	}
}
