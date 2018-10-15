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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.List;

import org.junit.Test;

import uk.ac.sanger.artemis.circular.digest.CutSite;
import uk.ac.sanger.artemis.circular.digest.EmbossTableParser;

public class EmbossTableParserTest
{
  /**
   * Test that restrict output is parsed
   */
  @Test
  public void testParser() throws IOException
  {
    final InputStream inputStream = EmbossTableParserTest.class
        .getResourceAsStream("/data/foo.restrict");

    InputStreamReader reader = new InputStreamReader(inputStream);
    EmbossTableParser etp = new EmbossTableParser();

    List<CutSite> cutSites = etp.parse(new BufferedReader(reader));
    CutSite firstCutSite = cutSites.get(0);

    assertEquals("Number of cut sites", cutSites.size(), 4);
    assertEquals("Enzyme name", firstCutSite.getEnzymeName(), "HindIII");
    assertEquals("3prime", firstCutSite.getThreePrime(), 85);
    assertEquals("5prime", firstCutSite.getFivePrime(), 81);
    assertEquals("3prime-rev", firstCutSite.getThreePrimeRev(), 0);
    assertEquals("5prime-rev", firstCutSite.getFivePrimeRev(), 0);
    assertTrue("Cut site strand", firstCutSite.isForward());
  }
}
