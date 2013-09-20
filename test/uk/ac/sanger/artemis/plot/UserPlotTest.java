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
package uk.ac.sanger.artemis.plot;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.plot.LineAttributes;
import uk.ac.sanger.artemis.plot.UserDataAlgorithm;
import uk.ac.sanger.artemis.util.FileDocument;
import uk.ac.sanger.artemis.io.Utils;

import uk.ac.sanger.artemis.sequence.Strand;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.net.URL;

import junit.framework.Assert;
import org.junit.Test;

public class UserPlotTest
{
  @Test
  /**
   * Base position user plot.
   */
  public void basePosition()
  {
    URL gffURL = UserPlotTest.class.getResource("/data/plot/base_position.plot.gz");
    EntryGroup entryGrp = Utils.getEntryGroup("/data/MAL1.embl.gz");
    final Strand fwdStrand = entryGrp.getBases ().getForwardStrand ();
    final FileDocument doc = new FileDocument (new File(gffURL.getFile()));
    
    try
    {
      final UserDataAlgorithm alg = new UserDataAlgorithm (fwdStrand, doc, false);
      float [] values = new float [alg.getValueCount()];
      int start = 278;
      alg.getValues(start, start+10, values);
      assertEquals("Number of plots",values.length,6);
      assertTrue("Value of plot 1 at base position 278",values[0]==804.99f);
    }
    catch(IOException e)
    {
      Assert.fail(e.getMessage());
    }
  }
  
  @Test
  /**
   * Base position (with labels) user plot.
   */
  public void labelBasePosition()
  {
    URL gffURL = UserPlotTest.class.getResource("/data/plot/base_position_labels.plot.gz");
    EntryGroup entryGrp = Utils.getEntryGroup("/data/MAL1.embl.gz");
    final Strand fwdStrand = entryGrp.getBases ().getForwardStrand ();
    final FileDocument doc = new FileDocument (new File(gffURL.getFile()));
    
    try
    {
      final UserDataAlgorithm alg = new UserDataAlgorithm (fwdStrand, doc, false);
      final LineAttributes lines[] = alg.getLineAttributes();
      assertEquals("Number of lines",lines.length,6);
      assertEquals("Number of plots",lines[0].getLabel(),"lab1");

      float [] values = new float [alg.getValueCount()];
      int start = 278;
      alg.getValues(start, start+10, values);
      assertEquals("Number of plots",values.length,6);
      assertTrue("Value of plot 1 at base position 278",values[0]==804.99f);
    }
    catch(IOException e)
    {
      Assert.fail(e.getMessage());
    }
  }

}
