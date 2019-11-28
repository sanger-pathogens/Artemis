/* 
 * This file is part of Artemis
 *
 * Copyright (C) 2013  Genome Research Limited
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

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;

import java.io.File;
import java.net.URL;

import junit.framework.Assert;

import org.junit.Test;

public class GffToEMBLTest
{
  @Test
  /**
   * Test writing EMBL format from GFF files.
   */
  public void testGffToEMBL()
  {
    testGffConversion("/data/gff/Pf3D7_01_02_v3.gff.gz");
    testGffConversion("/data/gff/test.gff.gz");
    testGffConversion("/data/gff/test_boundary.gff.gz");
  }
  
  private void testGffConversion(String gff)
  {
    URL gffURL = Utils.class.getResource(gff);
    File gffFile = new File(gffURL.getFile());
    final File outDir = gffFile.getParentFile();
    final boolean emblSubmission = false;
    final boolean flatten = true;
    final boolean gzip = true;
    
    final File result = new File(outDir.getAbsolutePath()+File.separator+gffFile.getName()+".embl.gz");
    assertFalse("EMBL file exists: "+result.getAbsolutePath(),  result.exists());
    if(result.exists())
      return;
    
    new GffToEMBL(gffURL.getFile(), outDir.getAbsolutePath(),
        emblSubmission, flatten, gzip);
    
    assertTrue("EMBL not created: "+result.getAbsolutePath(),  result.exists());
    if(result.exists() && result.isFile())
      result.delete();
  }
}
