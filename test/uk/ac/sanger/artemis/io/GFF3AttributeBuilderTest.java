/* 
 * This file is part of Artemis
 *
 * Copyright (C) 2014  Genome Research Limited
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
import static org.junit.Assert.assertEquals;
import junit.framework.Assert;

import org.junit.Test;
import org.junit.Before;

import uk.ac.sanger.artemis.io.GFF3AttributeAggregator;
import uk.ac.sanger.artemis.io.GFF3AttributeBuilder;
import uk.ac.sanger.artemis.io.GFF3Encoder;
import uk.ac.sanger.artemis.util.StringVector;

public class GFF3AttributeBuilderTest {
  private GFF3Encoder                    enc;
  private String[]                       invals  = { "foo", "bar" };
  private String[]                       invalsC = { "foo,bar", "baz,quux" };
  private static GFF3AttributeAggregator testProc;

  @Before
  public void setUp() {
    enc = new GFF3Encoder();
    testProc = new GFF3AttributeAggregator() {
      @Override
      public String process(StringVector values) {
        StringBuilder buffer = new StringBuilder();
        if (values != null && values.size() > 0) {
          for (int value_index = 0; value_index < values.size(); ++value_index) {
            buffer.append(">>"
                + GFF3Encoder.encode(values.elementAt(value_index)) + "<<");
            if (value_index < (values.size()) - 1)
              buffer.append("|");
          }
        }
        return buffer.toString();
      }
    };
  }

  @Test
  /**
   * Tests adding one attribute.
   */
  public void testAdd1() {
    GFF3AttributeBuilder ab = new GFF3AttributeBuilder();
    StringVector in = new StringVector(invals);

    ab.add("attr1", in);
    assertEquals(ab.toString(), "attr1=foo,bar");
  }

  @Test
  /**
   * Tests adding two different attributes.
   */
  public void testAdd2() {
    GFF3AttributeBuilder ab = new GFF3AttributeBuilder();
    StringVector in = new StringVector(invals);

    ab.add("attr1", in);
    ab.add("attr2", in);
    assertEquals(ab.toString(), "attr1=foo,bar;attr2=foo,bar");
  }

  @Test
  /**
   * Tests adding attributes with custom aggregators.
   */
  public void testAddWithAggs() {
    GFF3AttributeBuilder ab = new GFF3AttributeBuilder();
    StringVector in = new StringVector(invals);
    ab.setAggregator("attr1", testProc);
    ab.add("attr1", in);
    assertEquals(ab.toString(), "attr1=>>foo<<|>>bar<<");
    ab.add("attr1", in);
    assertEquals(ab.toString(), "attr1=>>foo<<|>>bar<< >>foo<<|>>bar<<");
  }

  @Test
  /**
   * Tests adding attributes (encoded values) with custom aggregators.
   */
  public void testAddWithAggsCommas() {
    GFF3AttributeBuilder ab = new GFF3AttributeBuilder();
    StringVector in = new StringVector(invalsC);
    ab.add("attr1", in);
    assertEquals(ab.toString(), "attr1=foo%2Cbar,baz%2Cquux");
    ab.add("attr1", in);
    assertEquals(ab.toString(),
        "attr1=foo%2Cbar,baz%2Cquux foo%2Cbar,baz%2Cquux");
    ab = new GFF3AttributeBuilder();
    ab.setAggregator("attr1", testProc);
    ab.add("attr1", in);
    assertEquals(ab.toString(), "attr1=>>foo%2Cbar<<|>>baz%2Cquux<<");
    ab.add("attr1", in);
    assertEquals(ab.toString(),
        "attr1=>>foo%2Cbar<<|>>baz%2Cquux<< >>foo%2Cbar<<|>>baz%2Cquux<<");
  }

  @Test
  /**
   * Tests the ignoring of attribute fields in the output.
   */
  public void testIgnore() {
    GFF3AttributeBuilder ab = new GFF3AttributeBuilder();
    StringVector in = new StringVector(invals);

    ab.ignore("attr1");
    ab.add("attr1", in);
    ab.add("attr2", in);
    assertEquals(ab.toString(), "attr2=foo,bar");
    ab.unignore("attr1");
    assertEquals(ab.toString(), "attr1=foo,bar;attr2=foo,bar");
  }

  @Test
  /**
   * Tests the handling of duplicate attributes.
   */
  public void testAddMultiAttr() {
    GFF3AttributeBuilder ab = new GFF3AttributeBuilder();
    StringVector in = new StringVector(invals);

    ab.add("attr1", in);
    ab.add("attr1", in);
    assertEquals(ab.toString(), "attr1=foo,bar foo,bar");
  }

  @Test
  /**
   * Tests the handling of duplicate attributes, with delimiter.
   */
  public void testAddMultiAttrWithGlue() {
    GFF3AttributeBuilder ab = new GFF3AttributeBuilder();
    StringVector in = new StringVector(invals);

    ab.setGlue("attr1", "X");
    ab.add("attr1", in);
    ab.add("attr1", in);
    ab.add("attr1", in);
    assertEquals(ab.toString(), "attr1=foo,barXfoo,barXfoo,bar");
  }

  @Test
  /**
   * Tests cloning of attributes to separate keys with default aggregator.
   */
  public void testClone() {
    GFF3AttributeBuilder ab = new GFF3AttributeBuilder();
    StringVector in = new StringVector(invals);

    ab.setClone("attr1", "brand_new");
    ab.add("attr1", in);
    assertEquals(ab.toString(), "attr1=foo,bar;brand_new=foo,bar");
    ab.add("brand_new", in);
    assertEquals(ab.toString(), "attr1=foo,bar;brand_new=foo,bar foo,bar");
  }

  @Test
  /**
   * Tests cloning of attributes to separate keys with different aggregators.
   */
  public void testCloneWithAggs() {
    GFF3AttributeBuilder ab = new GFF3AttributeBuilder();
    StringVector in = new StringVector(invals);
    ab.setClone("attr1", "brand_new");
    ab.setAggregator("brand_new", testProc);
    ab.add("attr1", in);
    assertEquals(ab.toString(), "attr1=foo,bar;brand_new=>>foo<<|>>bar<<");
    ab.add("attr1", in);
    assertEquals(ab.toString(),
        "attr1=foo,bar foo,bar;brand_new=>>foo<<|>>bar<< >>foo<<|>>bar<<");
  }

  @Test
  /**
   * Tests mapping/cloning of attributes to separate keys with different aggregators.
   */
  public void testMappingAndCloneWithAggs1() {
    GFF3AttributeBuilder ab = new GFF3AttributeBuilder();
    StringVector in = new StringVector(invals);

    ab.setMapping("attr1", "aaaa");
    ab.setClone("aaaa", "brand_new");
    ab.setAggregator("brand_new", testProc);
    ab.add("attr1", in);
    assertEquals(ab.toString(), "aaaa=foo,bar;brand_new=>>foo<<|>>bar<<");
    ab.add("attr1", in);
    assertEquals(ab.toString(),
        "aaaa=foo,bar foo,bar;brand_new=>>foo<<|>>bar<< >>foo<<|>>bar<<");
  }

  @Test
  /**
   * Tests mapping/cloning of attributes to separate keys with different aggregators.
   */
  public void testMappingAndCloneWithAggs2() {
    GFF3AttributeBuilder ab = new GFF3AttributeBuilder();
    StringVector in = new StringVector(invals);

    ab.setMapping("attr1", "aaaa");
    ab.setClone("attr1", "brand_new");
    ab.setAggregator("brand_new", testProc);
    ab.add("attr1", in);
    assertEquals(ab.toString(), "aaaa=foo,bar;brand_new=>>foo<<|>>bar<<");
    ab.add("attr1", in);
    assertEquals(ab.toString(),
        "aaaa=foo,bar foo,bar;brand_new=>>foo<<|>>bar<< >>foo<<|>>bar<<");
  }

  @Test
  /**
   * Tests mapping one attribute to a new name.
   */
  public void testMapping() {
    GFF3AttributeBuilder ab = new GFF3AttributeBuilder();
    StringVector in = new StringVector(invals);

    ab.setMapping("attr1", "brand_new");
    ab.add("attr1", in);
    ab.add("attr2", in);
    assertEquals(ab.toString(), "attr2=foo,bar;brand_new=foo,bar");
  }

  @Test
  /**
   * Tests mapping one attribute to a new name with custom target aggregator.
   */
  public void testMappingWithAggs() {
    GFF3AttributeBuilder ab = new GFF3AttributeBuilder();
    StringVector in = new StringVector(invals);
    ab.setMapping("attr1", "brand_new");
    ab.setAggregator("brand_new", testProc);
    ab.add("attr1", in);
    ab.add("attr2", in);
    assertEquals(ab.toString(), "attr2=foo,bar;brand_new=>>foo<<|>>bar<<");
  }

  @Test
  /**
   * Tests mapping one attribute to a new name with custom target aggregator.
   */
  public void testMappingCollision() {
    GFF3AttributeBuilder ab = new GFF3AttributeBuilder();
    StringVector in = new StringVector(invals);

    ab.setMapping("attr1", "attr2");
    ab.add("attr1", in);
    ab.add("attr2", in);
    assertEquals(ab.toString(), "attr2=foo,bar foo,bar");
  }

  @Test
  /**
   * Tests mapping one attribute to a new name with custom target aggregator.
   */
  public void testMappingCollisionWithAggs() {
    GFF3AttributeBuilder ab = new GFF3AttributeBuilder();
    StringVector in = new StringVector(invals);

    ab.setMapping("attr1", "attr2");
    ab.setAggregator("attr2", testProc);
    ab.add("attr1", in);
    ab.add("attr2", in);
    assertEquals(ab.toString(), "attr2=>>foo<<|>>bar<< >>foo<<|>>bar<<");
    ab = new GFF3AttributeBuilder();
    ab.setMapping("attr1", "attr2");
    ab.setAggregator("attr1", testProc);
    ab.add("attr1", in);
    ab.add("attr2", in);
    assertEquals(ab.toString(), "attr2=>>foo<<|>>bar<< foo,bar");
  }
}
