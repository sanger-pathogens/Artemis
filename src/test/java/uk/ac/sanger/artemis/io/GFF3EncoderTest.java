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

import uk.ac.sanger.artemis.io.GFF3Encoder;

public class GFF3EncoderTest {
  private GFF3Encoder enc;

  @Before
  public void setUp() {
    enc = new GFF3Encoder();
  }

  @Test
  /**
   * Tests whether there are missing map characters.
   */
  public void testMapEqualLength() {
    assertEquals(GFF3Encoder.mapChars.length, GFF3Encoder.mappedEscapes.length);
  }

  @Test
  /**
   * Tests whether the mapping is correct according to GFF3 spec.
   */
  public void testMapChars() {
    assertEquals(enc.encode("%"), "%25");
    assertEquals(enc.decode("%25"), "%");
    assertEquals(enc.encode("&"), "%26");
    assertEquals(enc.decode("%26"), "&");
    assertEquals(enc.encode(","), "%2C");
    assertEquals(enc.decode("%2C"), ",");
    assertEquals(enc.encode(";"), "%3B");
    assertEquals(enc.decode("%3B"), ";");
    assertEquals(enc.encode("="), "%3D");
    assertEquals(enc.decode("%3D"), "=");
    assertEquals(enc.encode("\t"), "%09");
    assertEquals(enc.decode("%09"), "\t");
    assertEquals(enc.encode("\n"), "%0A");
    assertEquals(enc.decode("%0A"), "\n");
    assertEquals(enc.encode("\r"), "%0D");
    assertEquals(enc.decode("%0D"), "\r");
  }

  @Test
  /**
   * Tests the decoding of escaped characters in GFF files.
   */
  public void testDecode() {
    for (int i=0; i < GFF3Encoder.mappedEscapes.length; i++) {
      assertEquals(enc.decode("test"+GFF3Encoder.mappedEscapes[i]+"foo"),
                   "test"+GFF3Encoder.mapChars[i]+"foo");
      assertEquals(enc.decode("test%"+GFF3Encoder.mappedEscapes[i]+"foo"),
                   "test%"+GFF3Encoder.mapChars[i]+"foo");
      assertEquals(enc.decode("test1"+GFF3Encoder.mappedEscapes[i]+"foo," +
                              "test2"+GFF3Encoder.mappedEscapes[i]),
                   "test1"+GFF3Encoder.mapChars[i]+"foo," + 
                   "test2"+GFF3Encoder.mapChars[i]);
    }
  }

  @Test
  /**
   * Tests the encoding of escaped characters in GFF files.
   */
  public void testEncode() {
    for (int i=0; i < GFF3Encoder.mappedEscapes.length; i++) {
      assertEquals(enc.encode("test"+GFF3Encoder.mapChars[i]+"foo"),
                   "test"+GFF3Encoder.mappedEscapes[i]+"foo");
      assertEquals(enc.encode("test%"+GFF3Encoder.mapChars[i]+"foo"),
                   "test%25"+GFF3Encoder.mappedEscapes[i]+"foo");
    }
  }

}
