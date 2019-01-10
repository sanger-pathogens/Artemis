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

import static org.junit.Assert.*;

import java.io.File;
import org.junit.Test;

public class CircularGenomeControllerTest
{
  /**
   * Test that EMBOSS_ROOT is set
   */
  @Test
  public void testEmbossRoot()
  {
    assertNotNull("EMBOSS_ROOT not set", System.getProperty("EMBOSS_ROOT"));
  }

  /**
   * Test that restrict can be found
   */
  @Test
  public void testEmbossInstalled()
  {
    String restrictPath = System.getProperty("EMBOSS_ROOT") + "/bin/restrict";
    File restrict = new File(restrictPath);
    assertTrue("restrict not found at: " + restrictPath, restrict.exists());
  }
}