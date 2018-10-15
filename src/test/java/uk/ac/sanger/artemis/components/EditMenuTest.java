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
package uk.ac.sanger.artemis.components;

import static org.junit.Assert.assertTrue;

import java.awt.GraphicsEnvironment;
import javax.swing.JFrame;

import junit.framework.Assert;

import org.junit.Test;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.FeatureSegmentVector;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.io.Entry;
import uk.ac.sanger.artemis.io.Utils;
import uk.ac.sanger.artemis.io.ValidateFeatureTest;


public class EditMenuTest 
{
  @Test
  /**
   * Split a gene model and undo
   */
  public void gffUnMergeAndUndo()
  {
    // ignore if in headless mode with no x11
    if(GraphicsEnvironment.getLocalGraphicsEnvironment().isHeadless())
      return;
        
    final EntryGroup egrp = Utils.getEntryGroup("/data/Pf3D7_01_02_v3.gff.gz");

    FeatureVector features = egrp.getAllFeatures();
    final Selection selection = new Selection(null);
    Feature f = Utils.getCDSFeatureByIdPrefix("PF3D7_0103500.1", features);
    final FeatureSegmentVector segments = f.getSegments();
    selection.add(f);
    selection.add(segments.elementAt(0));
    selection.add(segments.elementAt(1));
      
    // validate GFF features before split and undo and after
    // to ensure the gene model structures remain consistent
    ValidateFeatureTest.testAll(egrp);
    EditMenu.unmergeFeature(null, selection, egrp);
    EditMenu.undo(null, selection, egrp);
    ValidateFeatureTest.testAll(egrp);
  }
}